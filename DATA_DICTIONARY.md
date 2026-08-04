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

## External-cohort Figure S2B source data

These tables contain values derived from the public Wu et al. cohort, not from
participants enrolled in this study.

### `data/external/Figure_S2B_source_data.tsv`

| Field | Description |
| --- | --- |
| `public_cell_id` | New sequential identifier used only in this repository |
| `source_sample` | Sample identifier in the working public dataset |
| `donor` | Public donor label from the source study |
| `phenotype` | `black` or `white` hair |
| `krt_cluster` | Re-clustered keratinocyte subcluster |
| `KRT35`, `KRT85` | Normalized expression used for subcluster ranking |
| `KRT35_KRT85_score` | Per-cell mean of normalized KRT35 and KRT85 expression |
| `RAB7A` | Normalized expression plotted in Figure S2B |

### `data/external/Figure_S2B_cluster_scores.tsv`

| Field | Description |
| --- | --- |
| `krt_cluster` | Public-cohort keratinocyte subcluster |
| `mean_krt35_85` | Subcluster mean of the composite KRT35/KRT85 score |
| `total_cells` | Number of keratinocyte-lineage cells in the subcluster |
| `rank` | Descending rank by `mean_krt35_85` |
| `selected` | Whether the subcluster was in the top quartile |

### `data/external/Figure_S2B_sample_summary.tsv`

| Field | Description |
| --- | --- |
| `source_sample`, `donor`, `phenotype` | Public-cohort sample descriptors |
| `n_cells` | Selected cells in the sample |
| `mean_RAB7A`, `median_RAB7A` | Sample-level normalized-expression summaries |
| `pct_expressing_RAB7A` | Percentage of selected cells with RAB7A above zero |
