# RAB7A-PAR2 trafficking in hair graying

Corrected analysis-code release, 2 October 2026. This revision replaces the
previous external-cohort Figure S2B analysis. It does not replace or reprocess
the internal human cohort. The previous release remains in the Git history;
the changes are documented in `CHANGELOG.md`.

## Data and code locations

- **This study:** four processed count matrices and coded cell metadata from
  two paired donors are deposited in [Mendeley Data](https://doi.org/10.17632/nk9cmvsf6k.1).
  Version 1 is under embargo until 6 August 2027 (00:00 UTC); earlier release
  is planned upon article publication. It is not currently an openly
  downloadable dataset. No raw FASTQs are included.
- **Analysis code and external source data:** this repository,
  [rab7a-par2-hair-graying](https://github.com/hwangsoii/rab7a-par2-hair-graying).
  `CODE_MANIFEST.json` records the files and SHA-256 hashes for this revision.
- **External input:** Wu et al., *Cell Discovery* (2022),
  [article](https://doi.org/10.1038/s41421-022-00394-2),
  [NODE OEP002321](https://www.biosino.org/node/project/detail/OEP002321).
  The reconstructed analysis uses Audrey Onfroy's documented integer count
  matrices, [Zenodo 15103193](https://doi.org/10.5281/zenodo.15103193).
  This is a third-party deposit, not our study's data repository.

## Reproduce the corrected external panels

From the repository root:

```sh
Rscript supplementary_code/external_cohort/01_make_Figure_S2B.R
```

See `supplementary_code/external_cohort/README.md` for dependencies, optional
output paths, exact source-table fields and checks. The released small tables
reproduce the cell-level plots and nominal Wilcoxon statistics without any
private Seurat object, patient metadata or original cell barcode.

The selected population contains **3,315 cells: 2,684 black and 631 white**.
The all-keratinocyte context contains **15,214 cells: 10,382 black and 4,832 white**.
Both retain six samples from four donors, including two paired donors. These
overlapping populations are not independent replications. Nominal cell-level
P values do not account for donor clustering and do not establish donor-level
validation. Sample summaries are supplied even though the plots show cells.

## Internal cohort

Four libraries came from two biological donors, each contributing black and
gray scalp follicles. Coded library and donor assignments are in
`data/sample_manifest.tsv`. Internal matrices and coded metadata are unchanged
by this update. For the existing full pipeline:

```sh
export HAIR_BASE_DIR=/path/to/code_repository
export HAIR_RAW_DIR=/path/to/cellranger_samples
export HAIR_DATA_DIR=/path/to/working_data
cd supplementary_code
Rscript 01_analysis_pipeline.R
Rscript 02_keratinocyte_analysis.R
Rscript 03_make_Fig1.R
Rscript 04_make_SuppFig1.R
Rscript 05_replot_Fig1B_RAB7A_violin.R --nominal-cell-p
Rscript 06_replot_Fig1A_cluster6_umap.R
```

The last two scripts require the existing processed objects or regenerated
equivalents and assert the expected cell counts and cluster identity. They do
not perform a new analysis. Montages use historical internal panel letters;
see `supplementary_code/README.md` before assembling manuscript figures.

## Provenance and scope

`provenance/external_stage2/` preserves dated analysis plans, preprocessing
scripts and selected aggregate records. It is an auditable historical archive,
not a claim that the entire count-to-figure pipeline has been rerun in a clean
environment. The standalone table-to-plot workflow is separately verified.
`REPRODUCIBILITY.md` describes selection timing, sample mapping and limitations.

No raw reads, complete expression matrices, full single-cell objects, direct
participant identifiers or clinical demographics are included in this code
package. Code is MIT licensed; external derived data retain their source
attribution and applicable source-license terms.
