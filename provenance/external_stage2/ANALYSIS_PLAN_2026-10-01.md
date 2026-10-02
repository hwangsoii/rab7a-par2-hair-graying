# External cohort (OEP002321): dated prospective analysis plan

**Date written:** 2026-10-01, before any RAB7A value was computed from the reconstructed counts.
**Status:** the plan is frozen once its SHA-256 is recorded in `outputs/B0_plan_freeze.txt`. Later changes require a new dated version, and the reason must be stated there. The frozen version is never edited.
**Supersedes:** nothing. The current Figure S2B stays on HOLD and is not repaired by this plan. This is a **new, separately documented analysis**. It does not reuse any of the old h5ad or Seurat objects, or the legacy cluster or cell selections, as inputs.

## 0. Disclosure of prior exploratory results (analyst not blind)

The analyst (Claude Code), the user and the reviewers have all seen the following earlier results from the undocumented local h5ad input. They are disclosed here so that the reader can judge the risk of outcome-driven choices.
- **Legacy rule:** keratinocyte subclusters with mean KRT31/35/85 score > 0.5 selected all 21 clusters. Cell-level p 6.66e-14; sample-mean p 0.487.
- **2026-08-04 candidate sweep** in a Codex session (top-k under the legacy score; "Cortex"; KRT35or85+). It included a sample-mean p of 0.82 for the legacy-score top-6 candidate.
- **Post hoc top-quartile rule (the current S2B):** cell-level p 6.3e-4. It was adopted after candidate outcomes were visible.
- **Within-pair directions in the current S2B selection** (sample means): ryg047 → ryg048, white higher; black1_F62B → white2_F62W, white lower. Opposite directions.

None of the choices below were tuned against these numbers. Every threshold is copied unchanged from the internal-cohort scripts (`supplementary_code/01_analysis_pipeline.R`, `02_keratinocyte_analysis.R`). The one new element, the marker-based annotation rule (§4), does not use RAB7A or F2RL1.

## 1. Input (contingent on checkpoint A)

- **Source:** Zenodo record 15103193 (DOI 10.5281/zenodo.15103193; A. Onfroy, 2025; CC-BY-4.0), "Count matrices of the OEP00002321 scRNA-Seq dataset". It provides Cell Ranger outputs per article sample ID; code is in GitHub `audrey-onfroy/Onfroy_12tips_ACM_REP_25`.
- **Analysis input:** the **Cell Ranger-filtered** matrices (`<sample>_filtered_{matrix.mtx,features.tsv,barcodes.tsv}.gz`).
  - The raw-droplet matrices are used only to verify the filtered ones: subset, identical values, and cell calling.
  - Choosing filtered matrices mirrors the internal cohort, which starts from `filtered_feature_bc_matrix`.
  - Onfroy's own empty-droplet calling, which needs `aquarius` and is not installed, is **not** used.
- **Checkpoint A must pass before §3 runs.** It requires:
  - all MD5s match;
  - all values are nonnegative integers;
  - filtered barcodes are a subset of the raw barcodes, with identical values;
  - the per-sample web_summary cell numbers match the filtered matrices;
  - the sample labels are consistent with the official NODE mapping.
  If any check fails, the analysis stops and the failure is reported.

## 2. Samples, phenotype, donor (official mapping; all six retained)

| Article sample | Official FASTQ prefix (id2FastqID.xlsx) | Hair (Wu et al. Table S1) | Donor | Role |
|---|---|---|---|---|
| F18 | ryg035 | Black | F18 | unpaired black |
| F59 | ryg029 | Black | F59 | unpaired black |
| F31B | black-1 | Black | F31 | paired |
| F31W | white-2 | White | F31 | paired |
| F62B | ryg047 | Black | F62 | paired |
| F62W | ryg048 | White | F62 | paired |

- All six samples and **both** paired donors are retained whatever the direction of their results.
- No sample is relabelled by expression.
- Chemistry differs by donor: F31 is 3′ v3; the others are 3′ v2. This is handled by integration (§3) and stated as a limitation. It is not a reason to drop anyone.

## 3. Processing (internal-cohort parameters, unchanged)

All steps use Seurat 5.2.1 with sparse matrices. Samples are loaded one at a time.
1. Per sample, `CreateSeuratObject(min.cells = 3, min.features = 200)`, then `percent.mt` (`^MT-`).
2. **QC:** `nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA < 40000`.
   - There is no doublet removal, mirroring the internal pipeline. scDblFinder is not installed and will not be installed. This is a stated limitation.
3. Per sample, `NormalizeData` (LogNormalize, scale factor 1e4) and `FindVariableFeatures` (vst, 2000).
4. **CCA integration:** `FindIntegrationAnchors(dims = 1:20)` and `IntegrateData(dims = 1:20)`.
5. On the integrated assay: `ScaleData`, `RunPCA(npcs = 30)`, `FindNeighbors(dims = 1:20)`, `FindClusters(resolution = 0.5)`, `RunUMAP(dims = 1:20)`. Seurat default seeds are used.

## 4. Cell-type annotation rule (new; deterministic; no manual overrides)

- On the RNA assay (joined layers, LogNormalize data), compute `AddModuleScore` per cell for each marker set below, with Seurat defaults (`ctrl = 100`, `seed = 1`).
- Genes absent from the feature list are dropped, and the drop is reported.
- Each resolution-0.5 cluster gets the lineage with the highest **cluster-mean** module score.
- A cluster is **keratinocyte** if its top lineage is any `KC_*` set.
- The score margin between the top two lineages is reported for every cluster. A small margin does **not** trigger a manual override.

| Lineage | Markers |
|---|---|
| KC_basal_ORS | KRT5, KRT14, KRT15, KRT17, KRT6A, KRT16, COL17A1 |
| KC_IFE_suprabasal | KRT1, KRT10, KRTDAP, SBSN, DMKN |
| KC_IRS | KRT25, KRT27, KRT71, KRT73, TCHH |
| KC_cortex_cuticle | KRT31, KRT35, KRT85, KRT32, KRT33A |
| KC_matrix | MSX2, LEF1, HOXC13 |
| Fibroblast | COL1A1, COL1A2, DCN, LUM, PDGFRA |
| Endothelial | PECAM1, VWF, CDH5, CLDN5 |
| Lymphatic | LYVE1, PROX1, CCL21 |
| Mural | RGS5, ACTA2, TAGLN, MYH11, PDGFRB |
| Melanocyte | PMEL, MLANA, DCT, TYRP1, TYR |
| T_NK | PTPRC, CD3D, CD3E, CD2, NKG7 |
| Myeloid | CD68, CD14, LYZ, AIF1, CD207 |
| Mast | TPSAB1, CPA3, MS4A2 |
| Gland | DCD, SCGB2A2, PIP, MUCL1 |
| Schwann | MPZ, PLP1, PRX |

## 5. Keratinocyte subclustering and the single population rule (frozen before RAB7A)

1. **Subclustering.** Mirrors `02_keratinocyte_analysis.R` STEP 1. On the keratinocyte cells, using the RNA assay with joined layers:
   - `NormalizeData`, `FindVariableFeatures` (vst, 2000), `ScaleData`, `RunPCA(npcs = 30)`, `FindNeighbors(dims = 1:20)`, `FindClusters(resolution = 0.6)`.
   - UMAP for display only.
2. **Population rule** (internal threshold, unchanged):
   - Per cell, `KRT35_KRT85_score = (KRT35 + KRT85) / 2` from the RNA `data` layer.
   - Clusters whose mean score is **> 0.5** are selected, and every cell in a selected cluster is included.
   - **No fallback rule.**
     - If no cluster qualifies, there is no population; the RAB7A analysis is not run and the result is reported as inconclusive.
     - If every cluster qualifies, the population is all keratinocytes, reported as "rule non-discriminating".
   - No top-k, top-quartile or threshold sweep.
3. **Freeze.** Script 03 writes the following, and never reads RAB7A or F2RL1:
   - the cluster score table;
   - the selected clusters;
   - per-sample cell counts;
   - the SHA-256 of the sorted selected-cell IDs.
   The cell IDs themselves stay in the private data folder.
4. **Gate.** The RAB7A script (04) refuses to run unless both are true:
   - the freeze file exists;
   - the recorded SHA-256 of this plan matches the current file.

## 6. Outcome analysis (after independent review of the freeze)

**Primary (descriptive, donor level).**
- For each paired donor (F31, F62), compute RAB7A **pseudobulk log2 CPM** in the frozen population: the sum of RAB7A counts over the population cells of a sample, divided by the summed total counts of those cells, × 1e6, then log2(x + 1).
- Report white − black per donor. Both donors are reported regardless of direction.
- No P value is computed: two pairs cannot support an inferential test.

**Co-primary descriptive metric:** the per-sample mean of the RAB7A LogNormalize values over population cells. White − black is reported per donor.

**Secondary (exploratory, sample level).**
- Two-sided exact Wilcoxon rank-sum on the pseudobulk log2 CPM, 4 black vs 2 white samples.
- The minimum attainable two-sided p is 2/15 ≈ 0.133, stated in advance. The test cannot reach 0.05.

**Descriptive only.** Cell-level Wilcoxon (Black vs White), detection fraction and violin plots. These are labelled pseudoreplicated and are not used for inference.

**Interpretation rule** (internal cohort direction: lower in grey):
- "Direction consistent in both paired donors" only if white < black in **both** F31 and F62 on the primary metric.
- "Inconsistent across donors" if the two donors differ.
- "Opposite" if both are higher.
- No claim of validation rests on the cell-level P value.

**Not done in this plan:** DE across genes, F2RL1 analyses, threshold variants and alternative populations. Any of these would need a new dated plan.

## 7. Checkpoints

| Checkpoint | Content | Who clears it |
|---|---|---|
| A | Input verification report (checksums, integer counts, raw vs filtered, identity) | Author check; available for independent review |
| B | This plan, SHA-256 recorded before script 03 runs | Recorded automatically |
| C | Population freeze: annotation table, subcluster score table, selected clusters, per-sample counts | **Independent review (Reviewers A/B per `REVIEW_PROTOCOL.md`) before script 04 is run** |
| D | RAB7A outcome (script 04), reported in full | Independent review before any figure or manuscript change |

## 8. Outputs, privacy and controls

- Per-cell objects and IDs live only under `OneDrive_2025-09-11/Data/cell discovery/onfroy_zenodo_15103193/derived/`, which git ignores. Reports contain aggregates only and no barcodes.
- No existing figure, table, manuscript, EndNote file or released script is overwritten or edited. New outputs go only into this `stage2_reconstructed_counts_2026-10-01/` folder.
- No package installation, commit, push or upload.
- `sessionInfo()` and input hashes are recorded by each script.
