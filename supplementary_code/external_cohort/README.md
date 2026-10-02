# External Cohort: Frozen Count-Derived Replot

This downstream recipe replaces the legacy external S2B recipe **in this new
release only**. No historical figure, selection, object, or script is changed.
The current panel is **KRT35/KRT85-enriched keratinocytes** (v2); **All
keratinocytes** (v1) is retained as a secondary descriptive comparison. Both use
the count-derived populations frozen on 2026-10-01, not the old exploratory
top-quartile selection. This is a technical reproduction, not independent
biological validation or donor-level replication.

## Run

From the `code_repository` directory, with R, ggplot2, digest, jsonlite, and
Cairo graphics support already available:

```sh
Rscript supplementary_code/external_cohort/01_make_Figure_S2B.R \
  --output-dir external_cohort_replot
```

The script locates `data/external` relative to itself, so it also runs from an
unrelated working directory. `--data-dir DIR`, `--output-dir NEW_DIR`, and
`--font-family FONT` are supported. Output directories must not already exist.
SHA-256 checksums for all nine source tables and this script are checked against
`data/external/export_qa.json` before data are plotted. That manifest is a
local integrity record, not a cryptographic proof of external authorship.
No installation, network request, Seurat dependency, absolute source path,
upstream reconstruction, or outcome-script execution is needed. The local
`export_jid_external_release.R` helper is not a downstream dependency.

Two PDF/PNG panels (3.3 x 4.1 inches; PNG 300 dpi) and three verification TSVs
are written only to the chosen output directory. The script reproduces the
current width-scaled violins, all-cell x-only jitter (seed 20261001), white
median/IQR boxes, common y scale, colors, nominal P labels, and Arial typography.
Use `--font-family sans` if Arial is unavailable; font and ggplot2 versions can
change rendering without changing values or statistics. The reference
environment is R 4.3.3 and ggplot2 3.5.2. No cells, including zero values, are
downsampled or removed.

## Populations And Inference

| Population | Role | Total | black | white | W | Nominal P |
|---|---|---:|---:|---:|---:|---:|
| `v2_KC_cortex_cuticle` | Current S2B | 3,315 | 2,684 | 631 | 952884 | 8.0794986333324072e-07 |
| `v1_all_keratinocytes` | Secondary descriptive | 15,214 | 10,382 | 4,832 | 26497275.5 | 9.4861419954486193e-09 |

Both nominal tests use black as x and white as y:

```r
stats::wilcox.test(black, white, alternative = "two.sided",
                   paired = FALSE, exact = FALSE, correct = TRUE)
```

These are **nominal pooled-cell P values**, with the asymptotic tie correction
and continuity correction, not donor-adjusted or multiple-testing-adjusted
inference. Cells are treated as independent despite within-donor dependence.
The P annotations were requested **after the frozen descriptive outcomes had
been examined**; the amended descriptive outcome plan itself did not authorize
inferential P values. The populations are nested (v2 is a subset of v1) and
their results are **not independent replications**. Small pooled-cell P values
do not establish donor-level replication.

The old exploratory selection predates count reconstruction and was chosen
after exploratory outcomes were visible. The count-derived v2 amendment was
written after inspection of marker profiles, before the new count-derived
RAB7A outcomes, with historical results already known. Neither the amended
population nor the later P annotation is represented as an originally
prespecified, outcome-blind independent validation. This release does not
reproduce or endorse the legacy S2B top-quartile result.

## Sample Mapping

Six article samples represent four donors, with two paired donors. The official
article sample IDs and FASTQ prefixes below are deliberately retained. Donor
groups D01-D04 are new release-only labels; they are not additional subjects.
No ages, sex, undocumented local aliases, original cell barcodes, or
participant-level clinical metadata are released.

| Article sample | Public FASTQ prefix | phenotype | Donor group | Paired | v2 cells | v1 cells |
|---|---|---|---|---|---:|---:|
| F18 | ryg035 | black | D01 | no | 179 | 1,436 |
| F59 | ryg029 | black | D02 | no | 820 | 2,539 |
| F31B | black-1 | black | D03 | yes | 547 | 3,012 |
| F31W | white-2 | white | D03 | yes | 247 | 2,487 |
| F62B | ryg047 | black | D04 | yes | 1,138 | 3,395 |
| F62W | ryg048 | white | D04 | yes | 384 | 2,345 |

Both pairs and both descriptive metrics are retained. In v2, F31 has lower
white values on both metrics; F62 has lower white pseudobulk but slightly
higher white mean LogNormalize (metric disagreement). In v1, F31 has higher
white values on both metrics and F62 lower white values on both. No paired
inferential test is performed. Chemistry differs between F31 and the other
libraries; neither the cell-level test nor the RNA expression is made
donor-adjusted by the upstream integration used for clustering.

## Meaning Of The V2 Label

Membership is the entire existing `lineage == "KC_cortex_cuticle"` group in
the frozen integrated object. Original integrated resolution-0.5 clusters
**3, 6, 9, and 14** follow from that label; they are not new hard-coded selection
targets or the later keratinocyte-subcluster IDs. No cluster was excluded,
reclustered, rescored, or retuned for this release. V1 is the original frozen
non-discriminating selection, which retained all keratinocytes, and must not
be called KRT35/KRT85-enriched.

The current short v2 title denotes a **heterogeneous differentiating
cortex/cuticle-lineage compartment**, provisionally described in the frozen
marker review as hair-shaft differentiating keratinocytes with cortex/cuticle
marker enrichment. It is not a pure mature cortex/cuticle population or an
exact equivalent of the internal-cohort population:

| Integrated cluster | Cells | Frozen descriptive interpretation |
|---|---:|---|
| 3 | 1,430 | KRT35/KRT85-dominant, intermediate KRT31/KRT32 |
| 6 | 926 | Cortex-like, high KRT31 |
| 9 | 749 | Proliferating matrix/TAC-like, early hair-keratin expression; retained |
| 14 | 210 | Cuticle-like; elevated KRT25 is mixed evidence, not reassignment |

The annotation-consistency decision was `provisionally_consistent`, not
independent biological validation. Cluster 9 has a small cortex-versus-matrix
module-score margin (0.1244). KRT33A offers little support, hair-keratin/KRT14
detection is high even in comparators, and sample/chemistry background differs.
Mean expression and detection must therefore be considered together. Small
sample-by-cluster groups remain present. The aggregate marker tables include
competing basal/ORS, IRS, matrix/progenitor, and KRT1/KRT10 profiles rather than
only supportive markers. These interpretations preserve the dated review;
the short plot title does not relabel the biology.

## Released Tables

All tables are tab-delimited under `data/external`. Floating-point values are
serialized with 17 significant decimal digits. The cell table also stores an
exact hexadecimal expression representation because some R decimal readers
shift the last binary digit. The script checks decimal/hex agreement to 1e-14
and uses the hex representation to preserve the source doubles and Wilcoxon
ties exactly. Public phenotype labels are lowercase.

| File | Content |
|---|---|
| `rab7a_cells.tsv` | 18,529 population rows: opaque sequential `row_id`, `population`, article `sample`, `phenotype`, single-gene integer `rab7a_umi`, decimal `rab7a_lognorm`, exact `rab7a_lognorm_hex` |
| `sample_manifest.tsv` | Six official mappings, release-only donor groups, pairing flags |
| `sample_summaries.tsv` | 12 population/sample rows: cell and zero counts, summed RAB7A/denominator UMI, pseudobulk, mean expression, detection |
| `paired_descriptive.tsv` | Four population/pair rows: both white-minus-black differences and metric-disagreement flag; no P values |
| `nominal_cell_tests.tsv` | Reference counts, zero counts, W, nominal P, method settings, inference flags |
| `cluster_selection.tsv` | Original integrated cluster counts, labels, runner-up labels, margins, v1/v2 membership flags |
| `v2_cluster_by_sample.tsv` | All 24 candidate cluster/sample counts, including sparse groups |
| `markers_by_cluster.tsv` | Frozen marker means/detection by integrated keratinocyte cluster, candidate and comparator clusters |
| `markers_by_lineage.tsv` | Frozen marker means/detection by keratinocyte lineage group |
| `export_qa.json` | Technical export checks, relative source paths, before/after SHA-256 hashes, table hashes, and replot-script hash |

### Field Dictionary

Keys are stated explicitly; omitted fields are not present in these files.

| File | Fields and definitions |
|---|---|
| `rab7a_cells.tsv` | `row_id`: unique E000001-style release row; `population`: frozen population ID; `sample`: article sample ID; `phenotype`: black/white; `rab7a_umi`: nonnegative integer RAB7A count; `rab7a_lognorm`: nonnegative natural-log LogNormalize expression; `rab7a_lognorm_hex`: exact hexadecimal floating-point representation of the same expression, readable with R `as.numeric()` |
| `sample_manifest.tsv` | Key `sample`; `public_fastq_prefix`: official mapping; `phenotype`: black/white; `donor_group`: D01-D04 release-only grouping; `paired`: whether the donor contributes both colors |
| `sample_summaries.tsv` | Key `population`, `sample`; `phenotype`, `donor_group`: manifest mapping; `n_cells`, `n_zero`: all population cells and zero-RAB7A cells; `rab7a_umi`, `denominator_umi`: summed counts; `pseudobulk_log2cpm1`: log2(CPM + 1); `mean_lognorm`: arithmetic mean including zeros; `detection`: fraction with RAB7A count > 0 |
| `paired_descriptive.tsv` | Key `population`, `donor_group`; `black_sample`, `white_sample`: article IDs; `n_black`, `n_white`: cell counts; `diff_pseudobulk_log2cpm1`, `diff_mean_lognorm`: white minus black; `metric_disagreement`: opposite signs between the two metrics |
| `nominal_cell_tests.tsv` | Key `population`; `role`: current_S2B or secondary_descriptive; `n_black`, `n_white`, `zeros_black`, `zeros_white`: cell/zero counts; `statistic_w`: black-versus-white W; `p_nominal`: unadjusted P; `alternative`: two.sided; `paired`, `exact`: FALSE; `correct`: TRUE; `donor_adjusted`, `multiple_testing_adjusted`: FALSE |
| `cluster_selection.tsv` | Key `cluster`: original integrated resolution-0.5 ID; `n_cells`: original cluster size; `lineage`, `runner_up_lineage`: top two existing labels; `score_margin`: frozen top-minus-runner-up mean module score; `in_v1`, `in_v2`: frozen population membership flags |
| `v2_cluster_by_sample.tsv` | Key `cluster`, `sample`; `n_cells`: fixed v2 membership count, not a filtering threshold |
| `markers_by_cluster.tsv` | Key `cluster`, `gene`; `lineage`: frozen label; `candidate`: v2 membership; `n_cells`: cluster size; `gene`: marker symbol; `marker_set`: frozen review group; `mean_lognorm`, `detection`: aggregate expression and fraction with count > 0, rounded to four decimals in the source |
| `markers_by_lineage.tsv` | Key `lineage`, `gene`; `n_cells`, `gene`, `marker_set`, `mean_lognorm`, `detection`: same definitions, aggregated by keratinocyte lineage |
| `export_qa.json` | `release`, `export_script`, `export_script_sha256`, `replot_script`, `source_manifest_sha256`: version anchors; `sources`: relative paths and before/after SHA-256; `files`: relative paths, row counts, SHA-256; `populations`: reference cell tests; `validation`: technical checks; `limitations`: interpretation limits; `software`: extraction versions |

`row_id` is a new release-row identifier with no biological meaning, not a
barcode or a reversible encoding of a barcode. Rows in the two populations do
not expose a cross-population cell-identity key. Nesting was checked privately
against both frozen membership lists. It cannot be re-established from these
anonymized row IDs alone.

`rab7a_lognorm` is taken directly from the relevant frozen RNA `data` layer:
v2 uses `stage2_integrated_annotated.rds`; v1 uses
`stage2_keratinocytes_subclustered.rds`. The RNA `counts` layers provide the
single-gene counts and denominator summaries. Normalized values are not
substituted between objects, reconstructed from rounded summaries, or taken
from an integrated assay. The objects and full matrices are not distributed.

Detection is the fraction of population cells with RAB7A UMI > 0. Pseudobulk
is `log2(sum(RAB7A UMI) / sum(denominator UMI) * 1e6 + 1)`. Denominators are
RNA-count column sums, checked equal to `nCount_RNA`, over the object's 18,491
features (the union retained after per-sample `min.cells = 3` filtering), not
Cell Ranger total UMI. These are relative expression summaries, not absolute
transcript abundance. Per-cell denominator vectors are not released, so their
sample sums are verified against the frozen source at export, not independently
reconstructed by the downstream script. Marker summaries preserve the frozen
four-decimal rounding; they are explanatory, not a fresh selection rule.

## Validation And Boundaries

The local export checks C10 and every bound dependency, reads the actual
frozen objects and membership lists, checks all six sample counts, RNA-count
denominators, RAB7A sums/detection/means against frozen E1, and checks paired
differences against E2. It matches the current plotting QA's W and P values,
asserts an exact expression round trip, and rehashes its sources afterward.
The source paths in `export_qa.json` are informational repository-relative
provenance pointers, not downstream dependencies.

The downstream script checks package-source SHA-256 hashes and recomputes the *technical* pooled-cell
statistics from the small table, verifies sample/paired summaries to 1e-12,
asserts that every cell reaches the plot, and checks that input tables remain
unchanged. It does not independently validate biology or reproduce upstream
count reconstruction, annotation, clustering, or raw-data provenance. Frozen
outcome script 08 must not be rerun for this release. The main release's
provenance archive and top-level documentation are maintained separately.

### Checked On 2026-10-02

- Both R scripts parsed, and the existing frozen-source/figure smoke suite
  passed all five tests before export.
- All 34 source-file SHA-256 hashes matched before and after extraction and
  after isolated downstream testing; the nine released tables stayed unchanged.
- Both exact W/P results, all 12 sample summaries, all four pair summaries,
  and zero counts matched. V2 retains 602 black and 158 white zero values;
  v1 retains 3,634 black and 1,791 white zero values.
- A copy containing only these downstream files ran from an isolated directory,
  both with defaults and with explicit data/output options. No upstream object
  or script was available in that copy.
- Both 990 x 1230 PNGs were nonblank, visually inspected, pixel-identical to
  the frozen current PNGs, and identical on a repeated seeded run. Both PDF
  files had valid headers. PDF metadata scrubbing/distribution packaging was
  not performed here and remains part of the release packaging step.
- Deliberately altered table and script fixtures were rejected by SHA-256
  validation before output creation. A pre-existing output directory was
  rejected without modifying its contents.
- No upstream reconstruction, clustering, outcome-script-08 rerun, new
  biological validation, or donor-level inferential analysis was performed.

These are local technical checks, not an additional independent validation
cohort. The original frozen release remains unchanged.
