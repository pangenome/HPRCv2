# HPRCv2 Minigraph Coverage Analysis

This analysis visualizes the coverage depth of HPRCv2 minigraph sequences when mapped to the CHM13 reference genome.

## Overview

The HPRCv2 pangenome represents genetic variation from 44 diverse human genomes plus the CHM13 reference. This analysis:
1. Extracts sequence nodes from the minigraph structure
2. Maps them to CHM13 reference coordinates
3. Distributes long sequences across all overlapping 100kb bins
4. Visualizes coverage patterns to reveal structural variation hotspots

## Key Files

### Data Processing
- `extract_with_bin_spreading.py` - Extracts nodes from the RGFA file and properly distributes long sequences across bins
  - Input: `hprc-v2.0-mc-chm13.sv.gfa.gz` (minigraph RGFA)
  - Output: `bin_coverage_spread.tsv` (coverage statistics per 100kb bin)

### Visualization
- `create_final_log10_overlay.R` - Creates log10-scaled coverage visualization
  - Input: `bin_coverage_spread.tsv`
  - Output: `HPRCv2_minigraph_vs_CHM13_log10.pdf`

## Usage

1. Extract and spread sequences:
```bash
python3 extract_with_bin_spreading.py
```

2. Create visualization:
```bash
Rscript create_final_log10_overlay.R
```

## Key Findings

- **Total sequence**: 3.83 Gb across the pangenome
- **Non-reference**: 0.72 Gb (18.7%) represents variation not in CHM13
- **Coverage distribution**: Most bins show ~1x coverage (CHM13 only), with 165 bins exceeding 10x
- **Maximum coverage**: 128x at chr17:22Mb (centromeric region)
- **Nested sequences**: ~60% of non-reference sequences are >1 hop from the reference backbone

## Visualization Interpretation

The log10-scaled plot shows three overlaid metrics:
- **Blue line**: Total sequence coverage (reference + non-reference)
- **Orange line**: Non-reference sequence coverage only
- **Purple line**: Nested sequences (complex variations)

Most of the genome appears as a flat line near 1x coverage because the CHM13 reference provides continuous baseline coverage. Peaks indicate structural variation hotspots, particularly in:
- Centromeric regions
- Acrocentric short arms (chr13, 14, 15, 21, 22)
- Heterochromatic regions

## Requirements

- Python 3 with pandas
- R with ggplot2, dplyr, tidyr
- Input RGFA file: `hprc-v2.0-mc-chm13.sv.gfa.gz`