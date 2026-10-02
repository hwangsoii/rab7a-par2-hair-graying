# External cohort (OEP002321), Amendment 01 to the dated analysis plan

**Date written:** 2026-10-01, about 08:30 UTC.

**Amends:** `ANALYSIS_PLAN_2026-10-01.md`
- SHA-256 `9be9163a6661798c9626ee4e52fe340b0c34e93f5e7be0c4f30d0a89ce0a4c33`, frozen 2026-10-01T05:45:31Z.
- The original plan is **not edited**, and nor are its scripts (00–04), freeze records (B0, C9), outputs (A0–A3, C0–C9), reports or derived objects.

**Status:** frozen once its SHA-256 is recorded in `outputs/B1_amendment01_freeze.txt`. It is frozen **before** any amended marker-review output (`outputs/V2_*`) is created.

**Prompted by:**
- `CODEX_CHECKPOINT_C_REVIEW_2026-10-01.txt` (Reviewers A/Kepler and B/Aquinas);
- the user's plan review of 2026-10-01, which replaced a proposed pass/fail marker rule with the annotation-consistency review in §3.

**Change log:** `AMENDMENT_01_CHANGE_LOG.md`.

This amendment does **not** approve any population, and it does **not** authorise computing RAB7A or F2RL1 outcomes.

## 1. Timing and prior knowledge (disclosure)

**What was known when this was written**
- This is an **amendment written after inspection**. It is not an originally prespecified rule. It was written after the analysts had seen:
  - the KRT35/KRT85 subcluster score distribution (C5; all 20 subclusters > 0.5);
  - the original lineage annotation and its module-score margins (C2, C3);
  - Codex's descriptive group-mean check of 14 lineage markers on the saved keratinocyte object.
- It was written **before any new count-based RAB7A outcome** had been inspected.
- The analysts are not blind to historical results. The plan's §0 disclosures are retained unchanged:
  - **Legacy rule on the undocumented local h5ad.** The mean KRT31/35/85 score > 0.5 selected all 21 clusters: cell-level p 6.66e-14, sample-mean p 0.487.
  - **2026-08-04 candidate sweep** in a Codex session, including a sample-mean p of 0.82 for the legacy-score top-6 candidate.
  - **Post hoc top-quartile rule (the current S2B):** cell-level p 6.3e-4, adopted after the candidate outcomes were visible.
  - **Opposite within-pair directions in the current S2B selection.** By the Stage 2 barcode evidence, "white higher" belongs to donor F62 and "white lower" to donor F31.

**Attestation versus verification**
- **The analyst (Claude Code) attests** that no outcome-specific RAB7A or F2RL1 comparison has been examined on the reconstructed counts. RAB7A is **not** claimed to be excluded from general processing: normalisation, variable-feature selection, scaling and PCA used all genes.
- **Reviewer B verified something narrower:** the inspected pre-outcome scripts (02, 03 and the helper) contain no target-specific outcome extraction. Reviewers did not independently establish the broader attestation.

## 2. Populations

| Population | Definition | Role | Expression source |
|---|---|---|---|
| **v1** | Frozen in C9 (status `rule_non_discriminating_all_keratinocytes`): all 15,214 keratinocytes. Cell-list SHA-256 `688dbc3d…97cbc` | Contextual analysis. Unchanged; it must **not** be called "KRT35/KRT85-enriched" | `stage2_keratinocytes_subclustered.rds` (SHA-256 `8b1593f2…74b4fa`): its own RNA `counts` and `data` layers |
| **v2 (candidate)** | Every cell in `stage2_integrated_annotated.rds` whose existing `lineage` equals `"KC_cortex_cuticle"` | Targeted amended analysis, **if** §3 concludes `provisionally_consistent` | `stage2_integrated_annotated.rds`: RNA `counts` and `data` layers |

**How v2 membership is defined**
- Membership is derived at run time from the existing label, which script 02 assigned before this amendment.
- It is **not** derived from:
  - a hard-coded cluster or cell list;
  - a ranking of KRT35/KRT85 scores;
  - top-k or top-quartile rules;
  - any new numeric threshold.
- Cluster identifiers are the **original integrated resolution-0.5 clusters** (`seurat_clusters` in the integrated object). They are not the later RNA-assay keratinocyte-subcluster IDs of script 03.

**Audit values from C2/C3 (checks, not targets)**
- Clusters 3, 6, 9 and 14; 3,315 cells.

| Sample | Cells |
|---|---:|
| F18 | 179 |
| F59 | 820 |
| F31B | 547 |
| F31W | 247 |
| F62B | 1,138 |
| F62W | 384 |

- If the source label disagrees with these values, the source label is used and the discrepancy is reported.

**Naming and relation between the populations**
- **Provisional name**, used only if supported: "hair-shaft differentiating keratinocytes with cortex/cuticle-marker enrichment". No exact equivalence to the internal-cohort KRT35/KRT85 population is claimed.
- v2 is a subset of v1. The two populations overlap, so they must **not** be presented as independent replications of any biological association.

## 3. Annotation-consistency review of the v2 candidate (marker-only)

This review checks whether the fixed, already-assigned label is consistent with its marker profiles.
- It is **not** independent biological validation.
- Module-score labels and UMAP proximity are not used as proof.
- Mean LogNormalize expression and detection fraction (counts > 0) are reported side by side as complementary descriptive evidence.
- No significance test is used.
- The context is Wu et al. 2022 (Cell Discovery 8:49), Results "Single-cell RNA-seq reveals lineage trajectory of hsHF" and Fig. 1:
  - **Cortex (CO):** KRT31.
  - **Fibre cuticle (FC):** "expresses KRT35 but lacks KRT31".
  - **IRS:** GATA3/TCHH double positive.
  - **Matrix:** MSX2-enriched; matrix TACs LEF1-enriched.
  - **KRT10:** expressed in upper/middle ORS of the isolated follicle.
  - **Companion layer:** none discrete.

**Markers**
- Outcome genes RAB7A and F2RL1 are excluded by an assertion in code.
  - **Cortex/cuticle:** KRT31, KRT35, KRT85, KRT32, KRT33A.
  - **Basal/ORS:** KRT5, KRT14, KRT15, KRT17, KRT6A, KRT16, COL17A1.
  - **IRS:** KRT25, KRT27, KRT71, KRT73, TCHH, plus GATA3. GATA3 is a Wu et al. context marker that was not in plan §4.
  - **Matrix/progenitor:** MSX2, LEF1, HOXC13, plus MKI67 (context, not in plan §4).
  - **KRT1/KRT10.**

**Review items**

**S1, membership.**
- Derive membership from the existing label and reconcile it with C2's per-cluster `top_lineage`.
- Report the clusters, the total and the per-sample counts against the audit values, without forcing them.
- All candidate clusters are retained, **including cluster 9**, whose cortex-versus-matrix margin is small (0.1244).

**S2, hair-keratin profiles.**
- Review KRT31, KRT35, KRT85, KRT32 and KRT33A pooled, by cluster and by sample, against the other keratinocyte label groups.
- Complementary profiles are acceptable:
  - **cortex-like:** KRT31-associated;
  - **cuticle-like:** KRT35-expressing with little KRT31.
- There is **no** requirement that KRT31, KRT35 and KRT85 each exceed every comparator on both metrics in every cluster.
- KRT33A is expected to be nearly undetected. It is reported, not counted as support.

**S3, competing profiles, assessed as patterns rather than as universal marker purity.**
- **Basal/ORS profile:** the basal/ORS keratins and COL17A1.
- **IRS profile:**
  - TCHH alone is **not** treated as IRS-specific;
  - it is interpreted together with GATA3 and the IRS keratins KRT25, KRT27, KRT71 and KRT73.
- **Matrix/progenitor profile:**
  - MSX2-, LEF1-, HOXC13- or MKI67-associated signal is **not** automatically incompatible with cortex/cuticle differentiation, because hair-shaft lineages arise from matrix TACs.
- **KRT1/KRT10:** interpreted in the upper-ORS context. They are not by themselves evidence of interfollicular-epidermis contamination.
- **Mixed-lineage co-detection** is reported per cluster as fractions:
  - KRT31+;
  - KRT35+ with KRT31−;
  - (KRT35 or KRT85) together with KRT14;
  - together with GATA3 and TCHH;
  - together with MSX2 or LEF1.
- Mixed or contradictory evidence is recorded. **No cell or cluster is deleted.**

**S4, sample representation.**
- All six samples must contribute cells.
- Report sample-specific profiles, sparse groups (for example, sample × cluster cells with few cells) and discordance between samples.
- There is no per-sample pass rule, majority vote, score threshold or minimum-effect cutoff.

**Conclusion for the WHOLE fixed candidate**
- Exactly one of:
  - `provisionally_consistent` with the broad label;
  - `unresolved`;
  - `inconsistent`.
- It is recorded in `CHECKPOINT_C2_AMENDMENT01_MARKER_REVIEW.md` and, in structured form, in `outputs/V2_D0_annotation_decision.json`. That file binds the conclusion to the hashes of:
  - the marker-review evidence;
  - the source object;
  - the report;
  - the review script.
- The conclusion is the analyst's. It is not independent approval.
- If the conclusion is `unresolved` (material unresolved contradictions) or `inconsistent`, the work **stops for review**:
  - no v2 freeze;
  - no adjustment of membership;
  - no search of other subsets.

## 4. Freezing v2 (only if `provisionally_consistent`)

Script 07 writes the private sorted cell list `derived/stage2_v2_population_cells.txt` and `outputs/C10_v2_population_freeze.json`.

C10 records:
- the operational definition, the derived clusters, the total and per-sample counts, and the v2 ⊂ v1 check;
- the SHA-256 of each dependency:
  - the original plan, B0 and C9;
  - this amendment and B1;
  - the v1 and v2 cell lists;
  - the keratinocyte and integrated objects;
  - A4, V2_M0, V2_D0 and the C2 report;
  - scripts 05–08, `stage2_common.R`, `stage2_validate.R` and `stage2_amend_common.R`.

C10 does **not** contain its own checksum. Once C10 is final, the independent review record stores C10's SHA-256 (see §6). This review-record template is created with `approved: false`, and no approval is fabricated.

Both populations, v1 and v2, are frozen before any RAB7A outcome in either population is computed.

## 5. Amended descriptive outcome plan (supersedes plan §6 for the amended execution)

This is prospective. It is recorded before any outcome is viewed.

**For each population, v1 and v2, report:**
- **all six samples:** cell count, summed RAB7A UMI, summed denominator UMI, pseudobulk log2(CPM + 1) and mean LogNormalize;
- **both paired donors, F31 and F62:** White − Black on **both** metrics;
- **per-cell distribution and detection plots,** descriptive and without P values.

**Pseudobulk metric**
- log2(sum of RAB7A counts ÷ sum of denominator counts × 10⁶ + 1) over the population cells of a sample.
- The **denominator feature universe** is all features in that object's RNA `counts` layer: the union of features kept by the per-sample `CreateSeuratObject(min.cells = 3)` filter. The per-cell denominator is `colSums(counts)`, asserted equal to `nCount_RNA`. It is therefore not the Cell Ranger total UMI.
- These are relative expression summaries. They are not absolute transcript abundance per cell, and they depend on RNA composition.

**Removed from the amended execution** (changes from plan §6)
- The secondary exact Wilcoxon rank-sum test, 4 Black versus 2 White samples. The two paired donors (F31, F62) make those samples non-independent.
- The cell-level Wilcoxon, which is pseudoreplicated.
- Any P value on the violins or distribution plots and on the publication figure.
- No other inferential model is added.
- Two pairs give very limited inferential resolution, so the paired comparison stays descriptive, with no significance claim.

**Reporting categories (replace the plan §6 interpretation rule)**
- For each donor and metric: "white lower", "white higher" or "zero difference" (exact numerical equality).
- If the two metrics have different signs for a donor, the result is reported as **metric disagreement**.
- No overall concordant/discordant verdict is forced.
- "Lower in white in both paired donors" may be stated only if all four donor × metric entries are "white lower".

**Failure handling**
- The frozen definitions require all six samples. Any of the following makes the run **fail with diagnostics**:
  - a missing sample;
  - a zero-cell sample;
  - a non-finite or zero denominator;
  - a donor pair without exactly one Black and one White sample.
- No sample or donor is dropped, no undefined value is replaced with 0, and no incomplete pair is reported.

**Expression sources**
- **v1** uses the frozen keratinocyte object's `counts` and `data` layers, so the original plan's v1 definitions are preserved.
- **v2** uses the integrated object's `counts` and `data` layers.
- Equality of counts alone would not justify substituting one object's normalised `data` layer for the other's.

**Batch correction**
- CCA integration was used only to build the integrated clustering and annotation.
- Script 03's RNA subclustering and all RNA-based outcome summaries are **not** batch-corrected by that fact.
- Several RNA subclusters are dominated by F31 (3′ v3). Their composition alone does not show whether such differences are technical or biological.

**Reporting commitment:** both populations, both paired donors and all six samples are reported, regardless of the direction of the results.

## 6. Execution gate (new runner `08_outcome_amendment01.R`; the original 04 is not used)

The checks run in this order, all before any object is loaded:
1. **The review record.**
   - `--review-record` must name an existing file;
   - it must have `approved == true`, at least two named reviewers and an ISO date.
2. **C10 against the record:** the current C10's SHA-256 must equal the record's `c10_sha256`.
3. **Every C10 dependency:** each must match both the current file and the record's `approved_versions` (shared validator, `stage2_validate.R`).
4. **A4** must read `PASS`.

**Bound artifacts:**
- the original plan, B0 and C9;
- this amendment and B1;
- both population manifests;
- both source objects;
- A4;
- the executable and helper scripts.

**Gate-only mode:** `--gate-only` stops after these checks, before any `readRDS`, gene extraction or output-directory creation.

**Output folder:** `outputs/E_outcome_amendment01`, never overwritten.

**Checkpoint A enforcement**
- Before any marker review, `05_enforce_checkpoint_A.R` enforces the recorded checkpoint-A results and fails on any missing or false value. It covers:
  - every integer/nonnegative check;
  - raw-filtered equality;
  - feature identity;
  - checksum and size;
  - the web-summary checks;
  - a 42/42 hash comparison against the author manifest.
- It enforces **recorded** results and records the checksums of its inputs and of the verifier. It does not independently repeat the full matrix checks.
- Reviewer A's separate receipt (43/43 checksum and size; 42/42 author-manifest hashes; 36 barcode pairings; 12 matrix headers; the XLSX mapping) is preserved as a distinct record.

**Testing**
- Gate tests use isolated, test-only fixtures. They never use fabricated production approval records.
- The successful production-approval path remains untested until a genuine approval exists.

## 7. Stop point

This task stops after:
1. the marker-review report; and
2. **only if** the conclusion is `provisionally_consistent`, the v2 freeze (C10) and the unapproved review template.

No RAB7A or F2RL1 outcome is calculated until the independent review records approval of the exact versions. No commit, push, upload, or manuscript, figure, EndNote or release change is made.
