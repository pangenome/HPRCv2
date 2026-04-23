# Supplementary Figure: Inter-chromosomal matches in HPRC vs CHM13

## (a) Caption

**Figure SXXX. Inter-chromosomal matches across the CHM13 genome.**
Left: karyogram of CHM13 v2.0. Each 100 kb window is coloured by the number of other chromosomes from which at least one HPRC haplotype contig reaches that window: grey, none; blue to dark red, 1 to 11+. Overlaid tracks mark centromeres (CEN, black), pseudoautosomal regions (PAR, orange), acrocentric pseudo-homologous regions (PHR, green) and the X-transposed region (XTR, blue). Black boxes on the karyogram mark the three regions shown to the right as grey-scale heat maps: chr1p tip (0–0.4 Mbp, top right), chr13 short arm (0–5 Mbp, middle) and chr4q tip (193.2–193.6 Mbp, bottom right). In each heat map, cells show the number of haplotypes supporting the other-chromosome matches (row) in each window (column); darker cells mean more haplotypes.

## (b) Results

Inter-chromosomal matches are rare and concentrate at a small set of well-known shared regions. The most widespread signal is on the acrocentric short arms (chr13, 14, 15, 21, 22), which share sequence almost exclusively among themselves (chr13 p-arm inset). Smaller blocks appear at the pseudoautosomal regions (PAR), the X-transposed region (XTR), and around centromeres. The most extreme windows are subtelomeric: the chr1p tip matches many other chromosomes (top-right inset), while the chr4q tip is a near-exclusive match to chr10q through the shared D4Z4 macrosatellite, which is linked to facioscapulohumeral muscular dystrophy (bottom-right inset).

## (c) Methods

**Inter-chromosomal match quantification.**
We used WFMASH v0.23.0 (`-p 95 -P inf`) to produce a random subset of all-vs-all pairwise haplotype alignments (25,272 pairs spanning the 466 HPRC haplotypes from 234 samples). CHM13 was then divided into non-overlapping 100 kb windows. For each window *W* on chromosome *C*, we queried the implicit pangenome graph formed by the pairwise alignments with IMPG v0.4.1 (`-x -m 5 --min-identity 0.98 --min-transitive-len 50000 --min-output-length 50000`) to obtain every haplotype contig reachable from *W* within 5 transitive hops, with each hop ≥ 50 kb and ≥ 0.98 estimated identity. Each contig was assigned to its primary chromosome *P′*; chrUn and chrM were ignored. For each *W* we counted the number of distinct chromosomes (*P′* ≠ *C*) from which at least one contig was reachable.
