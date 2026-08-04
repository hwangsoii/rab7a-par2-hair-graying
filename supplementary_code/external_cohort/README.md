# External cohort and Figure S2B

Figure S2B uses the human scalp hair-follicle scRNA-seq cohort reported by Wu
et al. (2022):

- Article: <https://doi.org/10.1038/s41421-022-00394-2>
- Public data: NODE/BioSino accession OEP002321,
  <https://www.biosino.org/node/project/detail/OEP002321>
- Source authors' processing code: <https://github.com/zhendejuzi/scRNA_HF>

The six public samples used were F18 (`ryg029`), F31 black (`ryg047`), F31
white (`ryg048`), F59 (`ryg035`), F62 black (`black1_F62B`) and F62 white
(`white2_F62W`). F18 was classified as black hair according to the source
article and its Supplementary Table S1.

For the downstream analysis, keratinocyte-lineage cells were re-clustered and
subclusters were ranked by mean normalized KRT35/KRT85 expression. The top
quartile (6 of 21 subclusters) was retained. The released source-data table
contains 2,932 cells: 2,393 from black-hair samples and 539 from white-hair
samples. RAB7A distributions were compared using a two-sided Wilcoxon rank-sum
test at the cell level.

Run:

```sh
Rscript supplementary_code/external_cohort/01_make_Figure_S2B.R
```

This writes `figures/Figure_S2B.pdf` and
`results/Figure_S2B_statistics.tsv`. `data/external/` also includes the six-row
sample summary so that the pooled cell-level display is not mistaken for six
independent, balanced biological replicates.

The repository does not redistribute the third-party FASTQ files or the full
expression object. It provides the exact normalized values used in Figure S2B,
the cluster-selection summary and the downstream code needed to reproduce the
panel and reported cell-level statistic.
