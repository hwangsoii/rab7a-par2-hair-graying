# RAB7A-PAR2 trafficking in hair graying

Analysis code and source data accompanying the manuscript **"RAB7A regulates
PAR2 surface trafficking and melanosome uptake in keratinocytes."**

The repository reproduces the single-cell RNA-sequencing analyses in Figure 1
and Supplementary Figures S1-S2. It contains no raw sequencing reads, full
single-cell objects, direct participant identifiers or clinical demographics.

## Data availability

- **This study:** processed filtered count matrices and coded cell metadata for
  four libraries from two paired donors are available at
  [Mendeley Data](https://doi.org/10.17632/nk9cmvsf6k.1). Raw FASTQ files are
  not included.
- **External cohort:** Wu et al., *Cell Discovery* 2022
  ([doi:10.1038/s41421-022-00394-2](https://doi.org/10.1038/s41421-022-00394-2));
  NODE/BioSino accession
  [OEP002321](https://www.biosino.org/node/project/detail/OEP002321); source
  authors' code at [scRNA_HF](https://github.com/zhendejuzi/scRNA_HF).
  `data/external/` provides the exact normalized values and subcluster summary
  used for the downstream Figure S2B analysis. Third-party FASTQ files and the
  full expression object are not redistributed.

## Study design

Four single-cell libraries were generated from two donors. Each donor
contributed one black-hair and one gray-hair sample from the same scalp. The
four libraries are paired samples from two biological donors, not four
independent individuals.

| sample_id | donor_id | phenotype | paired_design |
| --- | --- | --- | --- |
| Black1 | Donor1 | Black | yes |
| Grey1 | Donor1 | Grey | yes |
| Black2 | Donor2 | Black | yes |
| Grey2 | Donor2 | Grey | yes |

## Repository layout

```text
supplementary_code/                    R analysis scripts and instructions
supplementary_code/external_cohort/    Figure S2B downstream analysis
data/sample_manifest.tsv               coded four-library manifest
data/external/                          exact Figure S2B source-data tables
DATA_DICTIONARY.md                      field definitions
PROCESSED_DATA_PACKAGE_README.md        deposited-data contents and processing
REPRODUCIBILITY.md                      versions, parameters and data flow
sessionInfo.txt                         analysis environment
```

## Quick start

The main cohort can be rerun from the four deposited Matrix Exchange triplets:

```sh
export HAIR_BASE_DIR=/path/to/project_root
export HAIR_RAW_DIR=/path/to/cellranger_samples
export HAIR_DATA_DIR=/path/to/working_data

cd supplementary_code
Rscript 01_analysis_pipeline.R
Rscript 02_keratinocyte_analysis.R
Rscript 03_make_Fig1.R
Rscript 04_make_SuppFig1.R
```

Figure S2B can be reproduced independently from the small source-data tables:

```sh
Rscript supplementary_code/external_cohort/01_make_Figure_S2B.R
```

The code requires Seurat 5 or later. The reported analyses used R 4.3.3 and
Seurat 5.2.1. See `supplementary_code/README.md` and `REPRODUCIBILITY.md` for
the full panel map and analysis parameters.

## License

Code is released under the MIT License. The external-cohort source-data tables
are provided solely to reproduce the reported downstream analysis and remain
subject to the terms of the original public dataset.
