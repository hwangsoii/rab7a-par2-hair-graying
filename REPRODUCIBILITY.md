# Reproducibility notes

## Environment

The reported analyses used R 4.3.3 and Seurat 5.2.1. The saved Seurat objects
report object version 5.1.0. `sessionInfo.txt` records the package environment.
Seurat 5 or later is required because the code calls `JoinLayers()`.

## Main-cohort data flow

```text
01_analysis_pipeline.R
  filtered Matrix Exchange files -> integrated and annotated Seurat objects

02_keratinocyte_analysis.R
  integrated object -> keratinocyte subclusters and differential-expression results

03_make_Fig1.R
  processed objects -> Figure 1B-E and Figure S1C-E source panels

04_make_SuppFig1.R
  processed objects + released Figure S2B table -> Figure S1A-B and Figure S2A-B panels
```

Main parameters:

- QC: `nFeature_RNA > 500`, `nFeature_RNA < 5000`, `percent.mt < 10`,
  `nCount_RNA < 40000`
- Normalization: Seurat log normalization and 2,000 variable features per
  library
- Integration: canonical correlation analysis, dimensions 1-20
- PCA: 30 components; UMAP, t-SNE and neighbors: dimensions 1-20
- Clustering resolution: 0.5 for all cells and 0.6 for keratinocytes
- Internal KRT35/KRT85 selection: subcluster mean composite score greater than
  0.5
- Differential expression: Wilcoxon rank-sum test; Seurat `p_val_adj` is
  Bonferroni adjusted across genes

The four libraries come from two paired donors. Cell-level tests describe
differences among profiled cells and do not increase the number of independent
donors.

## External cohort and Figure S2B

Figure S2B uses the public human hair-follicle cohort from Wu et al. 2022
([article](https://doi.org/10.1038/s41421-022-00394-2),
[OEP002321](https://www.biosino.org/node/project/detail/OEP002321),
[source code](https://github.com/zhendejuzi/scRNA_HF)). Six samples were used:
F18, F31 black, F31 white, F59, F62 black and F62 white. F18 was classified as
black hair according to the source article and Supplementary Table S1.

Keratinocyte-lineage cells were re-clustered into 21 subclusters. Subclusters
were ranked by their mean normalized KRT35/KRT85 composite expression, and the
top quartile (six subclusters) was retained. The released source table contains
2,932 selected cells, including 2,393 from black-hair samples and 539 from
white-hair samples. Figure S2B reports a two-sided cell-level Wilcoxon rank-sum
test. The six-row sample summary is provided to make the cohort structure
explicit.

Run `supplementary_code/external_cohort/01_make_Figure_S2B.R` to recreate the
panel and statistics from the exact released values. The full third-party
expression object and raw reads are not redistributed.

## Stochastic steps

Seurat defaults used in the released workflow are `FindClusters(random.seed =
0)`, `RunUMAP(seed.use = 42)` and `RunTSNE(seed.use = 1)`. Cosmetic point jitter
is unseeded and does not affect any statistic. Results remain sensitive to
package versions; compare reruns with `sessionInfo.txt`.
