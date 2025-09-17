#!/usr/bin/env python3
"""
Extract pangenome nodes and spread long sequences across all bins they overlap.
Instead of anchoring sequences at a single position, distribute them across
all reference bins they span.
"""

import gzip
import sys
from collections import defaultdict, deque
import re

def parse_pan_sn(sn_field):
    """Parse the SN field to extract sample, haplotype, and contig."""
    parts = sn_field.split('#')
    if len(parts) == 3:
        sample_with_hap = parts[0]
        haplotype = parts[1]
        contig = parts[2]

        # Check if sample name has haplotype suffix (e.g., NA20809.1)
        if '.' in sample_with_hap and sample_with_hap.split('.')[-1].isdigit():
            sample_base = '.'.join(sample_with_hap.split('.')[:-1])
            sample_hap = sample_with_hap.split('.')[-1]
            # Reconstruct pan-SN without redundant haplotype
            pan_sn = f"{sample_base}#{sample_hap}#{contig}"
        else:
            # Keep original format
            pan_sn = sn_field

        return pan_sn, sample_with_hap, haplotype, contig
    return sn_field, None, None, None

def main():
    gfa_file = "hprc-v2.0-mc-chm13.sv.gfa.gz"

    print("Phase 1: Reading graph structure...")
    sys.stdout.flush()

    nodes = {}  # node_id -> {pan_sn, start, end, length}
    edges = defaultdict(set)  # node -> set of connected nodes
    chm13_nodes = set()  # Track CHM13 reference nodes
    ref_positions = {}  # node_id -> (ref_chrom, ref_pos)

    with gzip.open(gfa_file, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 100000 == 0:
                print(f"  Processed {line_num:,} lines...")
                sys.stdout.flush()

            if line.startswith('S'):
                fields = line.strip().split('\t')
                node_id = fields[1]

                # Parse attributes
                attrs = {}
                for field in fields[3:]:
                    if ':' in field:
                        key, type_tag, value = field.split(':', 2)
                        attrs[key] = value

                if 'SN' in attrs and 'SO' in attrs and 'LN' in attrs:
                    sn_raw = attrs['SN']
                    start = int(attrs['SO'])
                    length = int(attrs['LN'])
                    end = start + length

                    # Parse pan-SN
                    pan_sn, sample, hap, contig = parse_pan_sn(sn_raw)

                    nodes[node_id] = {
                        'pan_sn': pan_sn,
                        'start': start,
                        'end': end,
                        'length': length,
                        'sample': sample,
                        'contig': contig
                    }

                    # Check if CHM13 reference
                    if sample and 'CHM13' in sample:
                        chm13_nodes.add(node_id)
                        # For CHM13, contig is the chromosome
                        if contig:
                            ref_positions[node_id] = (contig, start)

            elif line.startswith('L'):
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    from_node = fields[1]
                    to_node = fields[3]
                    edges[from_node].add(to_node)
                    edges[to_node].add(from_node)

    print(f"\nTotal nodes: {len(nodes):,}")
    print(f"CHM13 reference nodes: {len(chm13_nodes):,}")

    # Phase 2: BFS to find reference positions for all nodes
    print("\nPhase 2: Finding reference positions via graph traversal...")
    sys.stdout.flush()

    visited = set()
    node_distance = {}  # Track distance from reference
    nested_nodes = set()  # Nodes at distance >= 2

    # BFS from CHM13 nodes
    queue = deque()
    for node_id in chm13_nodes:
        queue.append((node_id, ref_positions[node_id], 0))
        visited.add(node_id)
        node_distance[node_id] = 0

    nodes_processed = 0
    while queue:
        node_id, (ref_chrom, ref_pos), distance = queue.popleft()
        nodes_processed += 1

        if nodes_processed % 10000 == 0:
            print(f"  Processed {nodes_processed:,} nodes, queue size: {len(queue):,}")
            sys.stdout.flush()

        # Set reference position for current node if not already set
        if node_id not in ref_positions and node_id in nodes:
            ref_positions[node_id] = (ref_chrom, ref_pos)
            node_distance[node_id] = distance
            if distance >= 2:
                nested_nodes.add(node_id)

        # Add unvisited neighbors to queue
        for neighbor in edges.get(node_id, []):
            if neighbor not in visited and neighbor in nodes:
                visited.add(neighbor)
                # Neighbor inherits the reference position
                queue.append((neighbor, (ref_chrom, ref_pos), distance + 1))

    print(f"\nNodes with reference positions: {len(ref_positions):,}")
    print(f"Nested nodes (distance >= 2): {len(nested_nodes):,}")

    # Phase 3: Write BED with bin spreading
    print("\nPhase 3: Writing BED with sequences spread across bins...")

    bin_size = 100000
    bin_coverage = defaultdict(lambda: defaultdict(lambda: {
        'total_bp': 0,
        'ref_bp': 0,
        'nonref_bp': 0,
        'nested_bp': 0,
        'node_count': 0,
        'ref_nodes': 0,
        'nonref_nodes': 0,
        'nested_nodes': 0
    }))

    # Process each node and spread it across bins
    nodes_written = 0
    for node_id, node_info in nodes.items():
        if node_id in ref_positions:
            ref_chrom, ref_start = ref_positions[node_id]
            ref_chrom_clean = f"CHM13#{ref_chrom}"

            pan_sn = node_info['pan_sn']
            node_start = node_info['start']
            node_end = node_info['end']
            node_length = node_info['length']

            # Determine node properties
            is_reference = 1 if node_id in chm13_nodes else 0
            is_nested = 1 if node_id in nested_nodes else 0
            distance = node_distance.get(node_id, -1)

            # CRITICAL: Spread this node across all bins it overlaps
            # The node spans from ref_start to ref_start + node_length
            ref_end = ref_start + node_length

            # Find all bins this node overlaps
            start_bin = (ref_start // bin_size) * bin_size
            end_bin = (ref_end // bin_size) * bin_size

            # If node is entirely within one bin
            if start_bin == end_bin:
                bin_key = (ref_chrom_clean, start_bin)
                bin_coverage[ref_chrom_clean][start_bin]['total_bp'] += node_length
                bin_coverage[ref_chrom_clean][start_bin]['node_count'] += 1

                if is_reference:
                    bin_coverage[ref_chrom_clean][start_bin]['ref_bp'] += node_length
                    bin_coverage[ref_chrom_clean][start_bin]['ref_nodes'] += 1
                else:
                    bin_coverage[ref_chrom_clean][start_bin]['nonref_bp'] += node_length
                    bin_coverage[ref_chrom_clean][start_bin]['nonref_nodes'] += 1

                if is_nested:
                    bin_coverage[ref_chrom_clean][start_bin]['nested_bp'] += node_length
                    bin_coverage[ref_chrom_clean][start_bin]['nested_nodes'] += 1

            else:
                # Node spans multiple bins - split proportionally
                for bin_start_pos in range(start_bin, end_bin + bin_size, bin_size):
                    bin_end_pos = bin_start_pos + bin_size

                    # Calculate overlap with this bin
                    overlap_start = max(ref_start, bin_start_pos)
                    overlap_end = min(ref_end, bin_end_pos)

                    if overlap_end > overlap_start:
                        overlap_length = overlap_end - overlap_start

                        # Add this portion to the bin
                        bin_coverage[ref_chrom_clean][bin_start_pos]['total_bp'] += overlap_length
                        # For node count, we could either count fractionally or count 1 if any overlap
                        # Let's count 1 if there's any overlap (each bin knows this node touches it)
                        bin_coverage[ref_chrom_clean][bin_start_pos]['node_count'] += 1

                        if is_reference:
                            bin_coverage[ref_chrom_clean][bin_start_pos]['ref_bp'] += overlap_length
                            bin_coverage[ref_chrom_clean][bin_start_pos]['ref_nodes'] += 1
                        else:
                            bin_coverage[ref_chrom_clean][bin_start_pos]['nonref_bp'] += overlap_length
                            bin_coverage[ref_chrom_clean][bin_start_pos]['nonref_nodes'] += 1

                        if is_nested:
                            bin_coverage[ref_chrom_clean][bin_start_pos]['nested_bp'] += overlap_length
                            bin_coverage[ref_chrom_clean][bin_start_pos]['nested_nodes'] += 1

            nodes_written += 1
            if nodes_written % 10000 == 0:
                print(f"  Processed {nodes_written:,} nodes...")

    # Write BED file with original format for compatibility
    print("\nWriting BED file...")
    with open("hprc-v2.0-mc-chm13.pansn.spread.bed", 'w') as out:
        # Write header
        out.write("#pan_sn\tstart\tend\tnode_id\tref_chrom\tref_pos\tis_reference\tdistance_from_ref\tis_nested\n")

        for node_id, node_info in nodes.items():
            if node_id in ref_positions:
                ref_chrom, ref_pos = ref_positions[node_id]
                ref_chrom_clean = f"CHM13#{ref_chrom}"

                pan_sn = node_info['pan_sn']
                node_start = node_info['start']
                node_end = node_info['end']

                is_reference = 1 if node_id in chm13_nodes else 0
                is_nested = 1 if node_id in nested_nodes else 0
                distance = node_distance.get(node_id, -1)

                out.write(f"{pan_sn}\t{node_start}\t{node_end}\t{node_id}\t"
                         f"{ref_chrom_clean}\t{ref_pos}\t{is_reference}\t{distance}\t{is_nested}\n")

    # Write bin statistics file
    print("\nWriting bin statistics with proper spreading...")
    with open("bin_coverage_spread.tsv", 'w') as out:
        out.write("chromosome\tbin_start\ttotal_bp\tref_bp\tnonref_bp\tnested_bp\t"
                  "coverage_depth\tnode_count\tref_nodes\tnonref_nodes\tnested_nodes\n")

        for chrom in sorted(bin_coverage.keys()):
            for bin_start in sorted(bin_coverage[chrom].keys()):
                stats = bin_coverage[chrom][bin_start]
                coverage_depth = stats['total_bp'] / bin_size

                out.write(f"{chrom}\t{bin_start}\t{stats['total_bp']}\t{stats['ref_bp']}\t"
                         f"{stats['nonref_bp']}\t{stats['nested_bp']}\t{coverage_depth:.2f}\t"
                         f"{stats['node_count']}\t{stats['ref_nodes']}\t"
                         f"{stats['nonref_nodes']}\t{stats['nested_nodes']}\n")

    print(f"\nComplete! Wrote {nodes_written:,} nodes to BED")
    print("Created files:")
    print("  - hprc-v2.0-mc-chm13.pansn.spread.bed (original format)")
    print("  - bin_coverage_spread.tsv (bin statistics with proper spreading)")

    # Print some statistics
    print("\nBin coverage statistics:")
    max_coverage = 0
    total_bins = 0
    high_coverage_bins = 0

    for chrom in bin_coverage:
        for bin_start, stats in bin_coverage[chrom].items():
            coverage = stats['total_bp'] / bin_size
            total_bins += 1
            if coverage > max_coverage:
                max_coverage = coverage
                max_chrom = chrom
                max_bin = bin_start
            if coverage > 2:
                high_coverage_bins += 1

    print(f"  Total bins with data: {total_bins:,}")
    print(f"  Bins with >2x coverage: {high_coverage_bins:,}")
    print(f"  Maximum coverage: {max_coverage:.1f}x at {max_chrom}:{max_bin/1e6:.1f}Mb")

if __name__ == "__main__":
    main()