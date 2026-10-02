# Amendment 01 change log (2026-10-01)

This log lists every change relative to the frozen record. The frozen record is not edited:
- `ANALYSIS_PLAN_2026-10-01.md` (SHA-256 `9be9163a…0a4c33`);
- scripts 00–04 (SHA-256s in B0; script 02 post-freeze sha `bc054d7f…`);
- B0 and C9;
- A0–A3 and C0–C9;
- `CHECKPOINT_A_INPUT_VERIFICATION.md`;
- `CHECKPOINT_C_POPULATION_FREEZE.md`;
- the derived objects.

Times are UTC. Entries after the amendment freeze are appended below the line "Post-freeze entries".

## A. Analysis-plan changes (prospective; no RAB7A/F2RL1 outcome examined)

| # | Time | Change | Original | Reason |
|---|---|---|---|---|
| A1 | 08:25 | Pre-amendment snapshot of SHA-256 for every existing file in this folder and `derived/` (28 files); `git status` recorded | — | To show afterwards that nothing frozen changed |
| A2 | 08:30 | Added the **v2 candidate population**: existing `lineage == "KC_cortex_cuticle"` label in `stage2_integrated_annotated.rds` | Plan §5: single rule, no alternative populations | v1 rule non-discriminating (C9); review recommendation (Codex review C1, C5, B5) |
| A3 | 08:30 | v1 kept as the **contextual** population and named "all keratinocytes" | Plan §5.2 | Codex review C1 |
| A4 | 08:30 | **Annotation-consistency review** (S1–S4; three-level conclusion) replaces the proposed universal S1–S4 pass rule; cortex-like and cuticle-like profiles allowed; TCHH read with GATA3; matrix signal not automatically incompatible; no per-sample pass rule | Initial proposed rule (plan file only, never frozen or run) | User plan review 2026-10-01; Wu et al. 2022 cortex (KRT31) vs fibre-cuticle (KRT35 without KRT31) distinction |
| A5 | 08:30 | Context markers **GATA3** and **MKI67** added to the marker review only (not to the annotation rule) | Plan §4 marker sets | Wu et al. IRS (GATA3/TCHH) and matrix-TAC descriptions |
| A6 | 08:30 | **Removed** the 4 Black vs 2 White exact Wilcoxon and the cell-level Wilcoxon; no P value on plots; no inferential model added | Plan §6 "Secondary" and "Descriptive only" | Paired donors violate independence; pseudoreplication (Codex review C2; Reviewer B, B1) |
| A7 | 08:30 | The plan §6 interpretation rule is replaced by per-donor × metric categories (white lower / white higher / zero difference) with an explicit metric-disagreement flag; no forced verdict | Plan §6 "consistent / inconsistent / opposite" from pseudobulk sign only | Reviewer B, B4 |
| A8 | 08:30 | A missing or zero-cell sample, a non-finite or zero denominator, or an incomplete pair **fails with diagnostics** (no "missing" category) | Not specified | User plan review |
| A9 | 08:30 | v1 expression uses the **keratinocyte object's own counts and data**; v2 uses the **integrated object** | 04 used the keratinocyte object (v1 only) | User plan review: preserve v1 definitions |
| A10 | 08:30 | The denominator feature universe is stated: the union of per-sample `min.cells = 3` features; `colSums(counts) == nCount_RNA` asserted | Implicit in 04 | Reviewer B, B4 |
| A11 | 08:30 | New gate (runner 08): dated review record with `approved: true`; non-circular chain review record → C10 SHA → C10 dependencies; `--gate-only`; shared side-effect-free validator | 04: `--reviewed` flag; script hash recorded but not enforced | Codex review C4; Reviewer B, B2; user plan review §2–3 |
| A12 | 08:30 | New wrapper 05 enforces the checkpoint-A recorded results (fail-closed) plus a 42/42 author-manifest hash join; requires an A4 PASS bound to current hashes before 06 | 01 enforced only the MD5 match | Reviewer A, A3 |

## B. Documentation-scope corrections (the original reports are not rewritten)

Each row quotes or locates the original wording and gives the corrected scope. These corrections apply wherever the originals are cited.

| # | Original (location) | Corrected scope | Source |
|---|---|---|---|
| B1 | "The local object therefore came from a **different count run** of the same libraries (different reference or software)" (`CHECKPOINT_A_INPUT_VERIFICATION.md` §4c) | The legacy `metadata.tsv` totals **differ** from the Zenodo column sums for matched barcodes (0 equal; median ratio 1.16–1.25). Script 01 compares stored totals and performs no inverse normalisation. The cause has not been established. Neither the legacy field's derivation nor a particular count run, reference or software is shown | Reviewer A, A1 |
| B2 | "The web summaries are the machine record and are taken as authoritative" (Checkpoint A §3) | The machine records report Cell Ranger **3.1.0** and the label **hg19ens91-**; the deposit prose says **3.0.1 / hg19**. The discrepancy is preserved and the exact reference contents are unknown. Its effect on quantification has **not** been tested. Checksums establish that the downloaded files are unchanged; they do not establish equivalence between quantification pipelines | Reviewer A, A2; user item 7 |
| B3 | "**Name link now verified.** The Stage 1 identity table … is confirmed by data" (Checkpoint A §4b) | Barcode overlap (93–99%, against at most 24 off-target) is **strong local-to-Zenodo library-correspondence evidence**. FASTQ size and chemistry concordance is **corroboration**, not a cryptographic or reproduced FASTQ-to-count reconstruction. The original FASTQ-to-output donor authentication remains unverified | Reviewer A; user item 7 |
| B4 | "No RAB7A value has been computed" (`CHECKPOINT_C_POPULATION_FREEZE.md` line 3; B0) | After genome-wide preprocessing, the precise statement is: **no outcome-specific RAB7A comparison has been examined**. Normalisation and PCA used all genes. The B0 statement at plan-freeze time is not shown to be false | Reviewer B, B3 |
| B5 | "No P value is computed: two pairs cannot support an inferential test" (plan §6) | Two pairs permit some exact calculations but give very limited inferential resolution; the paired comparison stays descriptive | Codex review C2 |
| B6 | Plan §2: "Chemistry differs by donor … This is handled by integration (§3)" | CCA integration supports only the integrated clustering and annotation. RNA subclustering (script 03) and RNA outcome summaries are not batch-corrected by it | Codex review C3; Reviewer B, B4 |
| B7 | Author manifest "all match" (Checkpoint A §2) | The author's `md5sum.txt` lists 42 files; the Zenodo API lists 43. A row count is not a checksum comparison. The 42/42 hash equality is enforced by 05 (A4) and was separately confirmed in Reviewer A's receipt | Reviewer A, A3 |
| B8 | Annotation label "KC_IFE_suprabasal" (C2) | This is a module-score label. KRT10 marks upper/middle ORS in Wu et al.'s isolated-follicle data, so the label does not by itself show interfollicular-epidermis contamination | Codex review C5 |
| B9 | Checkpoint C §3 "possibly ambient hair-keratin RNA" | A possible explanation, not a tested finding | Codex review C1 |

## Post-freeze entries

All post-freeze entries are UTC on 2026-10-01. The A2–A12 times above are approximate (written 08:25–08:27). The actual amendment freeze time is the one in B1.

| Time | Entry |
|---|---|
| 08:27:34 | Amendment frozen. SHA-256 `0d8f368a…af84` recorded in `outputs/B1_amendment01_freeze.txt`; file made read-only. No `V2_*` output existed at freeze. |
| ~08:29 | `05_enforce_checkpoint_A.R` (sha `b6a8bd2c…`) → `A4_checkpointA_enforced.json`: **PASS**, 120/120 recorded criteria, 42/42 author-manifest hash join, 43 downloads re-hashed unchanged. A first attempt stopped on a helper-name clash (`get` → `field`) before writing any output. |
| ~08:30 | `06_marker_review_v2.R` (sha `c7ff28cb…`) → `V2_M0`–`V2_M7`. Membership was derived from the label and reconciled 21/21 with C2. Clusters 3/6/9/14, 3,315 cells, matching the audit values. Marker-only; the integrated object's hash was unchanged (`63ca4077…ffb59`). |
| ~08:32 | `CHECKPOINT_C2_AMENDMENT01_MARKER_REVIEW.md` and `V2_D0_annotation_decision.json` (hash-bound): analyst conclusion **provisionally_consistent**, with caveats (cluster 9 matrix/TAC-like; saturated detection; KRT33A; KRT25; F31 background; composition). |
| 08:33 | **Post-freeze helper fix.** `stage2_validate.R` prefixed absolute paths with `./`. The first run of 07 therefore stopped with `missing_file` **before writing anything**. Fix: absolute paths are used as given. SHA changed from `68f8ebd8…` (B1) to `ed902acf…`. This is technical and does not affect analysis logic. The new hash is bound in C10. `stage2_amend_common.R` and `stage2_common.R` are unchanged from B1. |
| 08:33:48 | `07_freeze_v2_population.R` (sha `c8407803…`) → `derived/stage2_v2_population_cells.txt` (private) and `outputs/C10_v2_population_freeze.json` (SHA-256 `9a63d03d…41c5`). 20 dependencies; **no self-checksum**. Also wrote `REVIEW_RECORD_AMENDMENT01_TEMPLATE.json` with `approved: false` and no reviewers; it stores C10's SHA-256 and 20 `approved_versions`. |
| 08:34 | `scripts/tests/test_amendment01_gate.R` → `outputs/T1_gate_tests.txt`: **12/12 passed**. Coverage: validator unchanged/tampered/deleted/missing-field on isolated scratchpad fixtures; production `--gate-only` with no record (`review_record_missing`) and with the template (`review_not_approved`), each exit 2 at stage `approval_record`, never `GATE_OK`, no outcome folder. A first test run failed because Rscript mangles spaces in absolute `--file` paths; fixed by calling with a relative path. **Untested until genuine approval exists:** the successful approval path, plus the `c10_mismatch` and `dependency_not_approved` branches. |
| 08:35 | Repository smoke spec `npx --offline playwright test seed.spec.ts`: **could not run**, "Cannot find module '@playwright/test'" (not installed; not installed here). Its stray `test-results/.last-run.json` was removed so that `git status` matches the start snapshot. `seed.spec.ts` is an empty placeholder in any case. |
| 08:35 | Verification: all 28 pre-existing files in this folder and `derived/` keep their SHA-256; `git status` is unchanged; no `D_rab7a`/`E_*` folder exists; no barcode-like strings in the new outputs or reports. **No RAB7A or F2RL1 outcome was computed.** |
| 09:40 | After independent AI-agent review (`CODEX_AMENDMENT01_CHECKPOINT_C2_REVIEW_2026-10-01.txt`), the unchanged 08 was run with `--gate-only` and `REVIEW_RECORD_AMENDMENT01_OUTCOME_APPROVED_2026-10-01.json` (sha `54862e48…`). Result: `GATE_OK`, exit 0. All 56 existing stage/`derived/` hashes were snapshotted first. The smoke spec was retried via the configured Playwright test runner and is **still blocked** (`Cannot find module '@playwright/test'`); nothing was installed and the suite is not passed. Outputs therefore stay local and are not for export. |
| 09:43:28 | The unchanged 08 (sha `3511fcf2…`) was run **once** with the same record, without `--gate-only`. Exit 0; this is the first exercise of the successful production-approval path. Wrote `outputs/E_outcome_amendment01/` E0–E3. All 56 existing files are unchanged afterwards and `git status` is unchanged. Results are in `CHECKPOINT_E_AMENDMENT01_OUTCOME_REPORT.md`, pending independent review. |
