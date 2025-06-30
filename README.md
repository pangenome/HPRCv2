# HPRCv2
pangenome alignment, implicit/explicit graph, and variants for human pangenome project release 2

## Pangenome alignment

Get the alignment at ...................

## Implicit pangenome graph

Get the IMPG index at ...................

### What you can do with an implicit pangenome graph

Query the graph for a specific region:

```bash
```

Make a region-of-interest (ROI) explicit pangenome graph with PGGB:

```bash
impg C4
pggb
visualize
```

Make a ROI explicit pangenome graph with SPOA:

```bash
impg C4
spoa
visualize
```

Compute pairwise similarity in a ROI pangenome:

```bash
```

Perform principal component analysis (PCA) on a ROI pangenome:

```bash
```

Partition the pangenome:

```
```

The latter command was used to partition the whole HPRCv2 pangenome, build explicit pangenome graphs for each partition with PGGB, and lace all partition-specific graphs into a single "explicit" pangenome graph with GFALACE.

## Explicit pangenome graph

Get the explicit pangenome graph build with PGGB at ...................

## Variants

Get the PGGB-graph variants called with VG DECONSTRUCT at ...................
