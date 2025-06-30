# HPRCv2
pangenome alignment, implicit/explicit graph, and variants for human pangenome project release 2

## Pangenome alignment

Get the pangenome alignment - generated with [WFMASH](https://github.com/waveygang/wfmash) - at https://garrisonlab.s3.amazonaws.com/hprcv2/pafs/hprc25272.aln.paf.gz.

## Implicit pangenome graph

Get the IMPG index - generated with [IMPG](https://github.com/pangenome/impg) - at https://garrisonlab.s3.amazonaws.com/hprcv2/impg/hprc25272.impg.idx.

We built the IMPG index with the following command:

```bash
impg index -p hprc25272.aln.paf.gz 
```

which requires the bgzipped PAF file (representing the pangenome alignment) and its bgzip index (`.gzi` file). The latter can be found at https://garrisonlab.s3.amazonaws.com/hprcv2/pafs/hprc25272.aln.paf.gz.gzi or generated with the following command:

```bash
bgzip -r hprc25272.aln.paf.gz 
```

Put the IMPG index in the same directory as the PAF file, and then you can query the pangenome alignment with [IMPG](https://github.com/pangenome/impg).

### What you can do with an implicit pangenome graph

Query the pangenome for a specific region-of-interest (ROI):

```bash
impg query -p hprc25272.paf.gz -r GRCh38#0#chr8:5748405-13676927 --merge-distance 1000000 > hprcv2.human8p23-1.bed

awk '$3-$2>=2000000' hprcv2.human8p23-1.bed | sort | head | column -t
    CHM13#0#chr8                 7491000  11605998  .  .  -
    GRCh38#0#chr8                5748405  13676927  .  .  +
    HG00097#1#CM094064.1         7660033  11772000  .  .  -
    HG00097#2#CM094079.1         7506003  11613398  .  .  -
    HG00099#1#CM087320.1         7529012  11754000  .  .  -
    HG00099#2#CM087363.1         7598012  11827000  .  .  -
    HG00126#1#JBHIKU010000039.1  7663003  11762000  .  .  -
    HG00126#2#CM090126.1         7753017  11898998  .  .  +
    HG00128#1#CM090082.1         7696008  11797000  .  .  -
    HG00128#2#JBHIKT010000047.1  5598009  9704000   .  .  -
```

Make a ROI explicit pangenome graph with [PGGB](https://github.com/pangenome/pggb):

```bash
ls ls /lizardfs/guarracino/pangenomes/HPRCv2/*.fa.gz > hprcv2.fasta-paths.txt
impg query -p hprc25272.aln.paf.gz -r GRCh38#0#chr6:31972057-32055418 -o fasta --fasta-list hprcv2.fasta-paths.txt | bgzip -l 9 -@ 16 > hprc25272.C4.fa.gz
samtools faidx hprc25272.C4.fa.gz
pggb -i hprc25272.C4.fa.gz -o pggb.hprc25272.C4
```

This is the result graph visualized with [ODGI](https://github.com/pangenome/odgi):

![C4 pangenome graph layout](./images/hprc25272.C4.fa.gz.a65af12.11fba48.3bf8f48.smooth.final.og.lay.draw.png)


Compute pairwise similarity in a ROI pangenome:

```bash
impg similarity -p hprc25272.aln.paf.gz -r GRCh38#0#chr6:31972057-32055418 --fasta-list hprcv2.fasta-paths.txt > hprc25272.C4.similarity.tsv
```

Perform principal component analysis (PCA) on a ROI pangenome:

```bash
impg similarity -p hprc25272.aln.paf.gz -r GRCh38#0#chr6:31972057-32055418 --pca
```

######## CHECK REGION IN WINDOWEDPCA PAPER FOR A PC1 PLOT

Partition the pangenome:

```bash
impg ppartition -p hprc25272.aln.paf.gz -w 1000000 -m 5 -f 10000 -d 1000000 -l 1000 --chm13-total
```

The latter command was used to partition the whole HPRCv2 pangenome, build explicit pangenome graphs for each partition with PGGB, and lace all partition-specific graphs into a single "explicit" pangenome graph with GFALACE.

## Explicit pangenome graph

Get the explicit whole-pangenome PGGB-graph at https://garrisonlab.s3.amazonaws.com/hprcv2/gfas/w1000000-m5-f10000-d1000000-l1000-chm13-total.tmp.fix.gfa.gz.

## Variants

Get the whole-pangenome PGGB-graph variants called with VG DECONSTRUCT with respect to:
- GRCh38 at https://garrisonlab.s3.amazonaws.com/hprcv2/vcfs/hprc8424.GRCh38.merged.norm.vcf.gz
- CHM13 at https://garrisonlab.s3.amazonaws.com/hprcv2/vcfs/hprc8424.CHM13.merged.norm.vcf.gz
