# Reproducibility and interpretation

## Internal cohort (unchanged analysis)

Four libraries represent two paired donors, not four independent people.
QC uses >500 and <5,000 detected genes, <10% mitochondrial UMIs and <40,000
UMIs. Per-library LogNormalize and 2,000 variable genes precede CCA integration
(dimensions 1-20). PCA uses 30 components; UMAP/neighbors use dimensions 1-20.
Clustering resolutions are 0.5 for all cells and 0.6 for keratinocytes.
The internal KRT35/KRT85 cluster-mean threshold is >0.5. The existing selected
object comprises all 618 cells of cluster 6 (382 black, 236 gray).

`01`-`03` and `config.R` are copied unchanged from the prior code release.
`04` retains the internal calculations unchanged, but no longer builds the
obsolete external panel. Its melanocyte analysis tests 10 specified genes,
not 19, with a 10-test Bonferroni family. Other internal correction families
are 49 RAB genes and nine cell types. Seurat `FindMarkers` uses its own
Bonferroni-adjusted `p_val_adj`; it is not a Benjamini-Hochberg result.
The standalone `05` and `06` scripts replot existing internal objects.

## External counts and identity

Input is the six Cell Ranger-filtered integer UMI matrices deposited by
Audrey Onfroy at https://doi.org/10.5281/zenodo.15103193, derived from
https://www.biosino.org/node/project/detail/OEP002321. All 43 downloaded files
matched the Zenodo MD5 records. The 42 payload files also matched the author's
42-entry checksum list; the 43rd file is that checksum list itself.
Filtered columns exactly matched the raw-droplet
matrices, and retained cell counts/median UMIs matched the web summaries.
The deposit description says Cell Ranger 3.0.1/hg19; all six web summaries say
3.1.0 with an hg19ens91 reference. That documentation discrepancy remains
recorded rather than silently resolved.

| Sample | Official FASTQ prefix | Phenotype | Donor |
| --- | --- | --- | --- |
| F18 | ryg035 | black | F18 |
| F59 | ryg029 | black | F59 |
| F31B | black-1 | black | F31 |
| F31W | white-2 | white | F31 |
| F62B | ryg047 | black | F62 |
| F62W | ryg048 | white | F62 |

F18 and F59 contribute unpaired black samples. F31 and F62 contribute paired
black/white samples. The former local donor assignments were swapped; this
release uses the official mapping. The old h5ad input remains unauthenticated
and is not an input to the new analysis.

## Processing and population definition

The reconstruction uses Seurat 5.2.1 on R 4.3.3. Each sample uses
`CreateSeuratObject(min.cells=3, min.features=200)`, the internal QC thresholds,
LogNormalize (scale factor 10,000), 2,000 variable genes and CCA dimensions
1-20. A fixed lineage-marker module-score rule annotates integrated clusters
at resolution 0.5. No doublet-removal step was performed. The mitochondrial
filter removes approximately 53-55% of the F31 libraries, versus 11-31% in
the other libraries. Differing chemistry and QC retention remain limitations.

The first frozen >0.5 KRT35/KRT85 rule selected all 20 keratinocyte subclusters:
15,214 cells. This is the all-keratinocyte context, not an enriched subset.
Amendment 01 instead selected the integrated clusters assigned the
`KC_cortex_cuticle` marker lineage (clusters 3, 6, 9, 14), retaining all
3,315 cells. The lineage marker set is KRT31, KRT35, KRT85, KRT32 and KRT33A.
Marker review supports the display label "KRT35/KRT85-enriched keratinocytes"
but the subset includes broader differentiating/matrix-TAC features; it is
not a pure cortex population or the same selection rule as the internal cohort.

## Timing and statistics

The original count-based plan was frozen before new RAB7A outcome extraction,
with prior exploratory results disclosed. The amendment followed inspection
of lineage/marker distributions, before the new target-gene results. Analysts
were not blind to the earlier h5ad-based exploratory results. The old top-quartile
rule was selected after outcome inspection and is not reused here.

The frozen outcome runner produced descriptive per-sample and paired summaries
once. Nominal pooled-cell Wilcoxon tests were added subsequently, after viewing
those descriptive results, at the authors' request. They use two-sided,
unpaired, asymptotic rank-sum tests (`exact=FALSE`, `correct=TRUE`) with ties
and all zero values retained. No multiple-testing or donor-clustering adjustment
is applied. These P values characterize sampled-cell distributions, not
independent biological replication. The selected population is nested within
the all-keratinocyte population.

For the selected population, the paired mean-lognormalized changes have
opposite signs (F31 lower in white, F62 slightly higher), although both paired
pseudobulk log2(CPM+1) changes are lower. Do not summarize these findings as
uniform donor-level validation. Exact sample/paired summaries remain in the
release and aggregate provenance records.

## What was verified for this release

The small-table export is checked against the frozen objects, approved
membership and recorded sample summaries. Standalone plotting and its nominal
statistics are checked independently of access to the full expression objects.
R syntax, public-package contents and original-source hashes are checked.
The full internal and external preprocessing pipelines were NOT rerun for
this packaging update. Historical freeze gates are preserved and must not be
bypassed to advertise a fresh prospective analysis. See the provenance README
for files/paths required to study or reconstruct the historical workflow.
