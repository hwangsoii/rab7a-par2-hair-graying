# Processed data package

This package contains the processed expression data needed to rerun the main
single-cell analysis. It is deposited separately from the code repository at
[Mendeley Data](https://doi.org/10.17632/nk9cmvsf6k.1).

Raw sequencing reads are not included. Direct participant identifiers, exact
ages, sex, collection dates and institutional identifiers are also excluded.

## Study design

Four 10x Genomics libraries were generated from two donors. Each donor
contributed one black-hair and one gray-hair sample from the same scalp.

| sample_id | donor_id | phenotype | Matrix cells | Cells after QC |
| --- | --- | --- | --- | --- |
| Black1 | Donor1 | Black | 10,615 | 9,750 |
| Grey1 | Donor1 | Grey | 7,629 | 7,070 |
| Black2 | Donor2 | Black | 5,302 | 4,602 |
| Grey2 | Donor2 | Grey | 3,049 | 2,571 |
| | | **Total** | **26,595** | **23,993** |

## Package contents

| File | Description |
| --- | --- |
| `<sample>_matrix.mtx.gz` | Filtered sparse UMI matrix for one library |
| `<sample>_features.tsv.gz` | Ensembl identifier, gene symbol and feature type |
| `<sample>_barcodes.tsv.gz` | Barcodes in matrix-column order |
| `cell_metadata.tsv.gz` | Coded metadata for 23,993 QC-retained cells |
| `sample_manifest.tsv` | Four coded libraries and paired-donor mapping |
| `DATA_DICTIONARY.md` | Field definitions |
| `MANIFEST.tsv` | File sizes and SHA-256 checksums |

The package does not contain the Wu et al. external cohort. Revised external-cohort source values belong with the corrected analysis
code, not in this internal-cohort deposit. The original matrix and metadata
payloads are unchanged in the 2 October 2026 code update.

## Matrix generation

- Platform: 10x Genomics Chromium single-cell 5-prime
- Alignment and counting: Cell Ranger 3.1.0
- Reference: GRCh38-3.0.0
- Shared matrices: Cell Ranger `filtered_feature_bc_matrix` outputs, with only
  the filenames changed to include the coded library identifier

## Manuscript analysis

The analysis code is available at
<https://github.com/hwangsoii/rab7a-par2-hair-graying>. The main
steps were:

1. Cells were retained with `nFeature_RNA > 500`, `nFeature_RNA < 5000`,
   `percent.mt < 10` and `nCount_RNA < 40000`.
2. Each library underwent log normalization and variable-feature selection
   (`vst`, 2,000 features).
3. Libraries were integrated by canonical correlation analysis using dimensions
   1-20.
4. PCA used 30 components; UMAP, t-SNE and neighbor finding used dimensions
   1-20.
5. All cells were clustered at resolution 0.5 and annotated into nine major
   cell populations.
6. Keratinocytes were re-clustered at resolution 0.6. KRT35/KRT85-enriched
   subclusters were selected using a mean composite-expression threshold of
   greater than 0.5.
7. Differential expression used Seurat's Wilcoxon test; `p_val_adj` is the
   Bonferroni-adjusted value returned by `FindMarkers`.

The cell metadata provide the final cell-type assignment and UMAP coordinates.
The `sample` and `matrix_barcode` fields map each retained cell to its deposited
library matrix.


## Availability status checked on 2 October 2026

Mendeley Data version 1 is deposited under embargo until 6 August 2027
(00:00 UTC). Earlier public release is planned upon publication. The corrected
code destination is https://github.com/hwangsoii/rab7a-par2-hair-graying.
The corrected code and external source data are supplied in this repository.
This code revision does not modify the Mendeley files or embargo. Direct
comparison with the authenticated Mendeley upload remains pending. No new
FASTQ or participant data are added.
