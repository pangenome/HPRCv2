#!/usr/bin/env python3
"""
Dense visualization of pangenome variation and non-reference coverage w.r.t CHM13.

HPRCv2 vs CHM13 coverage + variant density tracks.

Tracks (densities per 100kb) plotted below the baseline
  - SVs (large): INS, DEL, INV  (HGSVC3 VCFs, CHM13 coords)
  - MEI: mobile element insertions (CSV, CHM13 coords)
  - SNV: single-nucleotide variants (VCF, CHM13 coords)
  - INDEL: small insertions+deletions from 'indel_alt' VCF (aggregated INS+DEL)

Coverage curves (log10 depth) plotted above the baseline
  - Total, Non-reference, Nested

Usage: e.g. python annotated_karyogram.py --input bin_coverage_spread.tsv --show-ins 
--show-del --show-inv --show-indel --show-mei --show-snv

Currently pulling from: https://www.internationalgenome.org/data-portal/data-collection/hgsvc3

"""

import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from textwrap import dedent
from pathlib import Path

# ---------------- I/O helpers ----------------

def read_vcf_positions(vcf_path, chroms_keep):
    """
    Fast VCF loader: returns DataFrame ['chromosome','pos','svtype'].
    Parses POS and SVTYPE from INFO. Filters to chroms_keep.
    """
    if not vcf_path:
        return pd.DataFrame(columns=["chromosome", "pos", "svtype"])
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        print(f"[warn] VCF not found: {vcf_path}", file=sys.stderr)
        return pd.DataFrame(columns=["chromosome", "pos", "svtype"])

    # Find the '#CHROM' header
    with vcf_path.open() as f:
        header = None
        for line in f:
            if line.startswith("#CHROM"):
                header = line.lstrip("#").strip().split("\t")
                break
    if header is None:
        raise RuntimeError(f"Could not find #CHROM header in {vcf_path}")

    usecols = ["#CHROM", "POS", "INFO"] if "#CHROM" in header else ["CHROM", "POS", "INFO"]
    df = pd.read_csv(
        vcf_path,
        sep="\t",
        comment="#",
        names=header,
        usecols=[c for c in usecols if c in header],
        dtype={"POS": np.int64},
        engine="c",
    )
    chrom_col = "#CHROM" if "#CHROM" in df.columns else "CHROM"
    df = df.rename(columns={chrom_col: "chromosome", "POS": "pos"})

    df = df[df["chromosome"].isin(chroms_keep)].copy()
    if df.empty:
        return pd.DataFrame(columns=["chromosome", "pos", "svtype"])

    svtype = df["INFO"].str.extract(r"(?:^|;)SVTYPE=([^;]+)")[0]
    df["svtype"] = svtype.fillna("NA")
    return df[["chromosome", "pos", "svtype"]]


def read_mei_positions(csv_path, chroms_keep):
    """
    Read MEI callset CSV: needs CHROM, POS (1-based).
    Returns DataFrame ['chromosome','pos','svtype'] with svtype='MEI'.
    """
    if not csv_path:
        return pd.DataFrame(columns=["chromosome", "pos", "svtype"])
    csv_path = Path(csv_path)
    if not csv_path.exists():
        print(f"[warn] MEI CSV not found: {csv_path}", file=sys.stderr)
        return pd.DataFrame(columns=["chromosome", "pos", "svtype"])

    usecols = ["CHROM", "POS"]
    df = pd.read_csv(csv_path, usecols=usecols, dtype={"CHROM": "string", "POS": np.int64})
    df = df.rename(columns={"CHROM": "chromosome", "POS": "pos"})
    df["chromosome"] = df["chromosome"].astype(str)
    df = df[df["chromosome"].isin(chroms_keep)].copy()
    if df.empty:
        return pd.DataFrame(columns=["chromosome", "pos", "svtype"])

    df["svtype"] = "MEI"
    return df[["chromosome", "pos", "svtype"]]


def sv_density_by_bin(df_sv, bin_size, sv_classes=None, aggregate_label=None):
    """
    Bin positions to windows and count per class.
    Returns tidy DF: ['chromosome','bin_start','svclass','count','log10_density'].

    If sv_classes is None => use the svtype column as-is.
    If aggregate_label is provided => aggregate all selected rows into a single class name.
    """
    if df_sv.empty:
        return pd.DataFrame(columns=["chromosome","bin_start","svclass","count","log10_density"])

    work = df_sv.copy()
    if sv_classes is not None:
        work = work[work["svtype"].isin(sv_classes)].copy()
        if work.empty:
            return pd.DataFrame(columns=["chromosome","bin_start","svclass","count","log10_density"])

    work["bin_start"] = ((work["pos"] - 1) // bin_size) * bin_size

    if aggregate_label:
        grouped = (
            work.groupby(["chromosome", "bin_start"], observed=True)
                .size()
                .reset_index(name="count")
        )
        grouped["svclass"] = aggregate_label
    else:
        grouped = (
            work.groupby(["chromosome", "bin_start", "svtype"], observed=True)
                .size()
                .reset_index(name="count")
                .rename(columns={"svtype":"svclass"})
        )

    grouped["log10_density"] = np.log10(grouped["count"].astype(float) + 0.01)
    return grouped


def _parse_track_list(s):
    """Parse comma/space separated list to canonical track names."""
    if not s:
        return set()
    tokens = [t.strip().lower() for t in s.replace(",", " ").split()]
    mapping = {"ins":"INS","del":"DEL","inv":"INV","mei":"MEI","snv":"SNV","indel":"INDEL"}
    return {mapping[t] for t in tokens if t in mapping}


def variant_y_from_zero(log10_density):
    """
    Anchor variant density at the same baseline (0) as coverage
    and force it downward with increasing density:
        y = -(log10_density + 2)
    Because with eps=0.01: count=0 -> log10_density=-2 -> y=0.
    """
    return -(log10_density + 2.0)

# ---------------- Main ----------------

def main():
    parser = argparse.ArgumentParser(description="Plot coverage + variant densities (SV/MEI/SNV/INDEL) on CHM13.")
    parser.add_argument("-i", "--input", default="bin_coverage_spread.tsv",
                        help="Input TSV from the spreading script (default: bin_coverage_spread.tsv)")
    parser.add_argument("--bin-size", type=int, default=100_000,
                        help="Bin size for coverage and variant densities (default: 100000)")
    parser.add_argument("--title", default="HPRCv2: Minigraph Sequence vs CHM13",
                        help="Plot title")
    parser.add_argument("--out1", default="HPRCv2_minigraph_vs_CHM13_log10.pdf",
                        help="Output PDF with Total / Non-ref / Nested + densities")
    parser.add_argument("--out2", default="HPRCv2_minigraph_vs_CHM13_simple.pdf",
                        help="Output PDF with Total / Non-ref + densities")

    # SV sources (can extend this to more sources of variation; as we need)
    # https://www.internationalgenome.org/data-portal/data-collection/hgsvc3
        # SV sources
    parser.add_argument("--vcf-insdel",
                        default="variants_T2T-CHM13_sv_insdel_sym_HGSVC2024v1.0.vcf",
                        help="VCF with large INS/DEL (HGSVC3, CHM13 coords)")
    parser.add_argument("--vcf-inv",
                        default="variants_T2T-CHM13_sv_inv_sym_HGSVC2024v1.0.vcf",
                        help="VCF with large INV (HGSVC3, CHM13 coords)")
    # MEI source
    parser.add_argument("--mei-csv",
                        default="MEI_Callset_T2T-CHM13.ALL.20241211.csv",
                        help="CSV with MEIs (CHM13 coords)")
    # SNV / small INDEL sources
    parser.add_argument("--vcf-snv",
                        default="variants_T2T-CHM13_snv_snv_alt_HGSVC2024v1.0.vcf",
                        help="VCF with SNVs (CHM13 coords)")
    parser.add_argument("--vcf-indel",
                        default="variants_T2T-CHM13_indel_insdel_alt_HGSVC2024v1.0.vcf",
                        help="VCF with small INDELs (CHM13 coords)")
                        

    # Individual toggles (kept for compatibility)
    parser.add_argument("--show-ins",   action="store_true", help="Show large insertion density")
    parser.add_argument("--show-del",   action="store_true", help="Show large deletion density")
    parser.add_argument("--show-inv",   action="store_true", help="Show inversion density")
    parser.add_argument("--show-mei",   action="store_true", help="Show MEI density")
    parser.add_argument("--show-snv",   action="store_true", help="Show SNV density")
    parser.add_argument("--show-indel", action="store_true", help="Show small INDEL (aggregated) density")
    parser.add_argument("--no-sv", action="store_true", help="Disable all variant density tracks")

    # Convenience selectors
    parser.add_argument(
        "--include-tracks",
        type=str, default=None,
        help="Comma/space-separated list of tracks to include (choices: ins, del, inv, mei, snv, indel). "
             "Overrides the individual --show-* flags if provided."
    )
    parser.add_argument(
        "--exclude-tracks",
        type=str, default=None,
        help="Comma/space-separated list of tracks to exclude (choices: ins, del, inv, mei, snv, indel)."
    )

    args = parser.parse_args()

    # Decide which tracks are enabled
    track_names = ["INS","DEL","INV","MEI","SNV","INDEL"]
    enabled = {t: False for t in track_names}

    if args.no_sv:
        pass  # all remain False
    else:
        include = _parse_track_list(args.include_tracks)
        exclude = _parse_track_list(args.exclude_tracks)

        if include:
            for t in include:
                enabled[t] = True
            if any([args.show_ins, args.show_del, args.show_inv, args.show_mei, args.show_snv, args.show_indel]):
                print("[note] --include-tracks provided; ignoring individual --show-* flags.", file=sys.stderr)
        elif any([args.show_ins, args.show_del, args.show_inv, args.show_mei, args.show_snv, args.show_indel]):
            enabled["INS"]   = args.show_ins
            enabled["DEL"]   = args.show_del
            enabled["INV"]   = args.show_inv
            enabled["MEI"]   = args.show_mei
            enabled["SNV"]   = args.show_snv
            enabled["INDEL"] = args.show_indel
        else:
            # default = include all
            for t in track_names:
                enabled[t] = True

        # apply excludes last
        for t in exclude:
            enabled[t] = False

    print("Tracks enabled:", ", ".join([t for t in track_names if enabled[t]]) or "none")

    print("Reading spread bin coverage data...")
    df = pd.read_csv(args.input, sep="\t")

    # Clean CHM13 prefix
    df["chromosome"] = df["chromosome"].astype(str).str.replace(r"^CHM13#", "", regex=True)

    # Canonical chromosomes
    all_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

    # Coverage
    data = df[df["chromosome"].isin(all_chroms)].copy()
    bs  = args.bin_size
    eps = 0.01

    data["total_coverage"]  = data["total_bp"]  / bs
    data["ref_coverage"]    = data["ref_bp"]    / bs
    data["nonref_coverage"] = data["nonref_bp"] / bs
    data["nested_coverage"] = data["nested_bp"] / bs

    data["log10_total"]   = np.log10(data["total_coverage"]  + eps)
    data["log10_nonref"]  = np.log10(data["nonref_coverage"] + eps)
    data["log10_nested"]  = np.log10(data["nested_coverage"] + eps)
    data["bin_mb"]        = data["bin_start"] / 1e6

    data["chromosome"] = pd.Categorical(data["chromosome"], categories=all_chroms, ordered=True)
    data = data.sort_values(["chromosome", "bin_start"])

    # Summary
    total_gb    = data["total_bp"].sum()   / 1e9
    nonref_gb   = data["nonref_bp"].sum()  / 1e9
    nonref_pct  = 100.0 * data["nonref_bp"].sum() / max(1.0, data["total_bp"].sum())
    max_cov     = data["total_coverage"].max()
    n_bins_high = (data["total_coverage"] > 10).sum()

    print("\nSummary for visualization:")
    print(f"  Total sequence: {total_gb:.2f} Gb")
    print(f"  Non-reference:  {nonref_gb:.2f} Gb ({nonref_pct:.1f}%)")
    print(f"  Maximum coverage: {max_cov:.1f}x")
    print(f"  Bins with >10x coverage: {n_bins_high}")

    # ---------- Variant densities ----------
    sv_tracks = {}  # key -> tidy df with ['chromosome','bin_start','svclass','count','log10_density']
    if not args.no_sv:
        print("\nReading variant files and computing densities...")

        # Large INS/DEL (from sv_insdel_sym)
        if enabled["INS"] or enabled["DEL"]:
            insdel_df = read_vcf_positions(args.vcf_insdel, all_chroms) if args.vcf_insdel else pd.DataFrame(columns=["chromosome","pos","svtype"])
            if enabled["INS"]:
                sv_tracks["INS"] = sv_density_by_bin(insdel_df, bs, sv_classes=["INS"])
                print("  INS variants:", sv_tracks["INS"]["count"].sum())
            if enabled["DEL"]:
                sv_tracks["DEL"] = sv_density_by_bin(insdel_df, bs, sv_classes=["DEL"])
                print("  DEL variants:", sv_tracks["DEL"]["count"].sum())

        # Inversions
        if enabled["INV"]:
            inv_df = read_vcf_positions(args.vcf_inv, all_chroms) if args.vcf_inv else pd.DataFrame(columns=["chromosome","pos","svtype"])
            sv_tracks["INV"] = sv_density_by_bin(inv_df, bs, sv_classes=["INV"])
            print("  INV variants:", sv_tracks["INV"]["count"].sum())

        # MEI
        if enabled["MEI"]:
            mei_df = read_mei_positions(args.mei_csv, all_chroms)
            sv_tracks["MEI"] = sv_density_by_bin(mei_df, bs, sv_classes=["MEI"])
            print("  MEI variants:", sv_tracks["MEI"]["count"].sum())

        # SNV
        if enabled["SNV"]:
            snv_df = read_vcf_positions(args.vcf_snv, all_chroms) if args.vcf_snv else pd.DataFrame(columns=["chromosome","pos","svtype"])
            sv_tracks["SNV"] = sv_density_by_bin(snv_df, bs, sv_classes=["SNV"])
            print("  SNV variants:", sv_tracks["SNV"]["count"].sum())

        # Small INDELs (aggregate INS+DEL from indel_alt)
        if enabled["INDEL"]:
            indel_df = read_vcf_positions(args.vcf_indel, all_chroms) if args.vcf_indel else pd.DataFrame(columns=["chromosome","pos","svtype"])
            sv_tracks["INDEL"] = sv_density_by_bin(indel_df, bs, sv_classes=["INS","DEL"], aggregate_label="INDEL")
            print("  INDEL variants:", sv_tracks["INDEL"]["count"].sum())

        # Order & category for plotting
        for k, dfk in list(sv_tracks.items()):
            if dfk is None or dfk.empty:
                continue
            dfk["chromosome"] = pd.Categorical(dfk["chromosome"], categories=all_chroms, ordered=True)
            sv_tracks[k] = dfk.sort_values(["chromosome", "bin_start"])

    # ---------- Plot helpers ----------
    def setup_facets(nrows, height_per=0.55, width=14, top_space=0.07, bottom_space=0.06, hspace=0.15):
        fig_height = max(6, nrows * height_per)
        fig, axes = plt.subplots(nrows=nrows, ncols=1, sharex=True, sharey=True, figsize=(width, fig_height))
        if nrows == 1:
            axes = [axes]
        plt.subplots_adjust(top=1 - top_space, bottom=bottom_space, hspace=hspace)
        return fig, axes

    # Determine lower bound so variant tracks don't clip.
    # Find the largest (most positive) log10 density across enabled tracks
    # then map it to the down-anchored y and extend a bit.
    min_variant_y = -2.5  # default if no tracks
    if not args.no_sv and any((k in sv_tracks and not sv_tracks[k].empty) for k in track_names):
        max_sv_log = 0.0
        for k, dfk in sv_tracks.items():
            if dfk is not None and not dfk.empty:
                max_sv_log = max(max_sv_log, float(dfk["log10_density"].max()))
        # y_variant = -(max_sv_log + 2); show a hair more
        min_variant_y = -(max_sv_log + 2.0) - 0.1
    lower_limit = min(-2.5, min_variant_y)  # at least show to -2.5
    y_limits = (lower_limit, 2.5)

    y_ticks      = [-2, -1, 0, 1, 2]
    y_ticklabels = ["0.01x", "0.1x", "1x", "10x", "100x"]  # 0 ≡ 1x

    present = [c for c in all_chroms if c in data["chromosome"].unique()]
    chrom_groups = [data[data["chromosome"] == c] for c in present]

    # Colors for variant tracks (valid matplotlib names/hex)
    sv_colors = {
        "INS":   "forestgreen",
        "DEL":   "crimson",
        "INV":    "purple",
        "MEI":   "teal",
        "SNV":   "dodgerblue",
        "INDEL": "sienna",
    }

    # ---------- Plot 1: with nested + densities anchored at 0 and going down ----------
    print("\nCreating final log10 overlay visualization (variant densities anchored at baseline and downward)...")
    fig1, axes1 = setup_facets(nrows=len(chrom_groups))
    for ax, chrom, sub in zip(axes1, present, chrom_groups):
        # Coverage (positive above baseline)
        ax.plot(sub["bin_mb"], sub["log10_total"],  linewidth=0.6, alpha=0.85, color="darkblue",   label="Total sequence")
        ax.plot(sub["bin_mb"], sub["log10_nonref"], linewidth=0.6, alpha=0.85, color="darkorange", label="Non-reference sequence")
        ax.plot(sub["bin_mb"], sub["log10_nested"], linewidth=0.6, alpha=0.80, color="darkviolet", label="Nested sequence")

        # Baseline (1x)
        ax.axhline(0.0, color="black", linewidth=0.3, alpha=0.5)

        # Variant densities (anchored at 0, descending with density)
        if not args.no_sv:
            for sv_name in ["INS","DEL","INV","MEI","SNV","INDEL"]:
                if sv_name not in sv_tracks or sv_tracks[sv_name].empty or not enabled[sv_name]:
                    continue
                sub_sv = sv_tracks[sv_name][sv_tracks[sv_name]["chromosome"] == chrom]
                if sub_sv.empty:
                    continue
                y_down = variant_y_from_zero(sub_sv["log10_density"].to_numpy())
                ax.plot(sub_sv["bin_start"] / 1e6, y_down,
                        linewidth=0.55, alpha=0.85, color=sv_colors.get(sv_name, "black"),
                        label=f"{sv_name} density")

        ax.set_ylim(y_limits)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_ticklabels, fontsize=6)
        ax.set_ylabel(chrom, rotation=0, ha="right", va="center", fontsize=7, labelpad=20)
        ax.grid(axis="y", linewidth=0.3, color="gray", alpha=0.3)
        ax.grid(axis="x", which="both", alpha=0.0)

    axes1[-1].set_xlabel("Position (Mb)", fontsize=8)
    axes1[-1].set_xlim(left=0)
    for ax in axes1:
        ax.set_xticks(np.arange(0, 251, 50))
        ax.tick_params(axis="x", labelsize=6)

    subtitle = (f"Coverage of minigraph nodes in 100kb bins across CHM13 reference | "
                f"Total: {total_gb:.2f} Gb | Non-reference: {nonref_gb:.2f} Gb ({nonref_pct:.1f}%) | "
                f"Max: {max_cov:.0f}x")

    # Dynamic color key for shown tracks
    shown_keys = [k for k in ["INS","DEL","INV","MEI","SNV","INDEL"]
                  if enabled[k] and k in sv_tracks and not sv_tracks[k].empty]
    color_key  = ", ".join([f"{k}={sv_colors[k]}" for k in shown_keys]) if shown_keys else "none"


    fig1.suptitle(args.title, fontsize=12, fontweight="bold", y=0.995)
    fig1.text(0.5, 0.965, subtitle, ha="center", va="top", fontsize=9)

    # Legend
    handles = [
        plt.Line2D([0], [0], color="darkblue",   lw=1.5, label="Total sequence"),
        plt.Line2D([0], [0], color="darkorange", lw=1.5, label="Non-reference sequence"),
        plt.Line2D([0], [0], color="darkviolet", lw=1.5, label="Nested sequence"),
    ]
    for k in ["INS","DEL","INV","MEI","SNV","INDEL"]:
        if enabled[k] and k in sv_tracks and sv_tracks[k] is not None and not sv_tracks[k].empty:
            handles.append(plt.Line2D([0], [0], color=sv_colors[k], lw=1.5, label=f"{k} density"))
    fig1.legend(handles=handles, loc="upper center", ncol=min(6, len(handles)), fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.985))

    fig1.savefig(args.out1, dpi=300, bbox_inches="tight")
    plt.close(fig1)

    # ---------- Plot 2: simplified (no nested) ----------
    print("\nCreating simplified version (no nested; densities anchored at baseline and downward)...")
    fig2, axes2 = setup_facets(nrows=len(chrom_groups))
    for ax, chrom, sub in zip(axes2, present, chrom_groups):
        ax.plot(sub["bin_mb"], sub["log10_total"],  linewidth=0.75, alpha=0.9, color="darkblue",   label="Total sequence")
        ax.plot(sub["bin_mb"], sub["log10_nonref"], linewidth=0.75, alpha=0.9, color="darkorange", label="Non-reference sequence")

        ax.axhline(0.0, color="black", linewidth=0.3, alpha=0.5)

        if not args.no_sv:
            for sv_name in ["INS","DEL","INV","MEI","SNV","INDEL"]:
                if sv_name not in sv_tracks or sv_tracks[sv_name].empty or not enabled[sv_name]:
                    continue
                sub_sv = sv_tracks[sv_name][sv_tracks[sv_name]["chromosome"] == chrom]
                if sub_sv.empty:
                    continue
                y_down = variant_y_from_zero(sub_sv["log10_density"].to_numpy())
                ax.plot(sub_sv["bin_start"] / 1e6, y_down,
                        linewidth=0.55, alpha=0.85, color=sv_colors.get(sv_name, "black"),
                        label=f"{sv_name} density")

        ax.set_ylim(y_limits)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_ticklabels, fontsize=6)
        ax.set_ylabel(chrom, rotation=0, ha="right", va="center", fontsize=7, labelpad=20)
        ax.grid(axis="y", linewidth=0.3, color="gray", alpha=0.3)
        ax.grid(axis="x", which="both", alpha=0.0)

    axes2[-1].set_xlabel("Position (Mb)", fontsize=8)
    axes2[-1].set_xlim(left=0)
    for ax in axes2:
        ax.set_xticks(np.arange(0, 251, 50))
        ax.tick_params(axis="x", labelsize=6)

    subtitle2 = (f"Coverage of minigraph nodes in 100kb bins across CHM13 reference | "
                 f"Total: {total_gb:.2f} Gb | Non-reference: {nonref_gb:.2f} Gb ({nonref_pct:.1f}%)")
    
    fig2.suptitle(args.title, fontsize=12, fontweight="bold", y=0.995)
    fig2.text(0.5, 0.965, subtitle2, ha="center", va="top", fontsize=9)

    handles2 = [
        plt.Line2D([0], [0], color="darkblue",   lw=1.5, label="Total sequence"),
        plt.Line2D([0], [0], color="darkorange", lw=1.5, label="Non-reference sequence"),
    ]
    for k in ["INS","DEL","INV","MEI","SNV","INDEL"]:
        if enabled[k] and k in sv_tracks and sv_tracks[k] is not None and not sv_tracks[k].empty:
            handles2.append(plt.Line2D([0], [0], color=sv_colors[k], lw=1.5, label=f"{k} density"))
    fig2.legend(handles=handles2, loc="upper center", ncol=min(6, len(handles2)), fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.985))

    fig2.savefig(args.out2, dpi=300, bbox_inches="tight")
    plt.close(fig2)

    print("\nDone!")
    print(f"  - {args.out1}")
    print(f"  - {args.out2}")

if __name__ == "__main__":
    main()