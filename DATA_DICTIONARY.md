# Data dictionary

## Main-cohort processed data

### `sample_manifest.tsv`

| Field | Description |
| --- | --- |
| `sample_id` | Coded library: `Black1`, `Grey1`, `Black2` or `Grey2` |
| `donor_id` | Coded donor: `Donor1` or `Donor2` |
| `phenotype` | `Black` or `Grey` |
| `paired_design` | `yes`; each donor contributed both phenotypes |
| `tissue` | `scalp_hair_follicle` |

### `cell_metadata.tsv.gz`

| Field | Description |
| --- | --- |
| `cell_barcode` | Cell identifier in the merged Seurat object |
| `matrix_barcode` | Original 10x barcode; maps to the matrix with `sample` |
| `sample` | Coded library identifier |
| `donor_id` | Coded donor identifier |
| `phenotype` | `Black` or `Grey` |
| `nCount_RNA` | Total UMI count |
| `nFeature_RNA` | Number of detected genes |
| `percent.mt` | Percentage of UMIs assigned to mitochondrial genes |
| `seurat_clusters` | Integrated-object cluster at resolution 0.5 |
| `cell_type` | Final major cell-type annotation |
| `umap_1`, `umap_2` | Integrated UMAP coordinates |

### Per-library Matrix Exchange files

| File | Description |
| --- | --- |
| `<sample>_matrix.mtx.gz` | Sparse genes-by-cells UMI matrix |
| `<sample>_features.tsv.gz` | Ensembl identifier, gene symbol and feature type |
| `<sample>_barcodes.tsv.gz` | Barcodes in matrix-column order |


## Revised external-cohort tables

The fields and validation rules for every revised table are documented in
`supplementary_code/external_cohort/README.md`. These are derived from public
third-party data. Opaque IDs replace original cell barcodes. Public sample and
donor labels are retained to make the four-donor/six-sample structure explicit.
The old 2,932-cell table and top-quartile cluster ranking are superseded, not
combined with the revised 3,315-cell analysis. Internal metadata definitions
above are unchanged.
