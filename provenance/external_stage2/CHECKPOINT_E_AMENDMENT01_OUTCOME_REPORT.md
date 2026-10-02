# Checkpoint E: Amendment 01 descriptive RAB7A outcome (2026-10-01)

Status: **computed, NOT yet independently reviewed.** These are descriptive values for the next independent review.
- This is not external validation.
- It is not an approval of the historical S2B, which remains on hold.
- It is not a manuscript or figure result.

## 1. Execution record
| Item | Value |
|---|---|
| Runner | `scripts/08_outcome_amendment01.R`, unchanged, SHA-256 `3511fcf2ef2f17577fa15ebf8c199bfc8a90ea9711a8adfee2d1ed3c1353792f` |
| Review record used | `REVIEW_RECORD_AMENDMENT01_OUTCOME_APPROVED_2026-10-01.json`, SHA-256 `54862e483fb9d1bf5f25acff4cbcb371afd5850d450b81d7aa1d2f691264f031` |
| Review type | Independent AI agents, not human scientific reviewers |
| C10 | `9a63d03d0e6197e65352e29cac72d4ec5411af4985377760732880756a2e41c5` |
| Pre-run gate (`--gate-only`, same record) | stages approval_record, c10, dependencies, gate_ok; `GATE_OK`; exit 0; no outcome folder created |
| Outcome run (same record, no `--gate-only`) | executed once, 2026-10-01 09:43:28 UTC; gate passed again; exit 0 |
| Outputs (`outputs/E_outcome_amendment01/`) | `E0_outcome.json` `85878ba5…3976`<br>`E1_per_sample.tsv` `ce1c9d3a…f22c`<br>`E2_paired_donors.tsv` `be193670…5a55`<br>`E3_distribution_detection_DESCRIPTIVE.pdf` `698342bb…d2de` |
| Feature universe (denominator) | 18,491 RNA features in both objects (union of per-sample `min.cells = 3`). The denominator is colSums(counts), asserted equal to nCount_RNA. |
| Runtime assertions | All passed: nonnegative integer counts, colSums = nCount_RNA, six intended samples per population, valid denominators, and one Black plus one White sample per paired donor. |

**Integrity**
- All 56 files that existed before the run in the stage folder and in `derived/` kept the same SHA-256. This includes C10, all 20 dependencies, the template and both review records.
- `git status` is unchanged.
- The E outputs contain no cell identifiers and no absolute paths.

**Smoke test (operational, separate from scientific clearance).** `npx playwright test seed.spec.ts` is **still blocked**.
- The configured Playwright test runner fails with `Cannot find module '@playwright/test'`.
- The only cached npx packages are `playwright`/`playwright-core` and `@playwright/cli`. The `@playwright/test` package that the spec imports is not among them.
- Nothing was installed and nothing was bypassed. The smoke test is **not** passed.
- AGENTS.md requires it before outputs are shared or exported. The E outputs are therefore local, git-ignored review material only and must not be exported until the suite runs green.

## 2. Populations (frozen, unchanged)
**v1** (contextual): all 15,214 keratinocytes.
- Expression comes from the keratinocyte object's own counts and data layers.

**v2** (targeted): the complete existing `KC_cortex_cuticle` label, 3,315 cells.
- It comprises integrated clusters 3, 6, 9 and 14, including the matrix/TAC-like, proliferating cluster 9.
- Expression comes from the integrated object.
- v2 is a broad, marker-supported, differentiation-associated hair-follicle keratinocyte population with a matrix/TAC-like component. It is not pure mature cortex or cuticle, and its equivalence to the internal population or to Wu et al.'s clusters is not established.

**Relationship between v1 and v2**
- v2 is nested in v1: every v2 cell is also a v1 cell. The two results are therefore **not independent replications**.
- The lower candidate share in the White samples is a composition observation (F31 18.2% → 9.9%; F62 33.5% → 16.4%). It is a different quantity from RAB7A expression within cells.

## 3. Per-sample values (E1)
"Pseudobulk" = log2(RAB7A UMI / denominator UMI × 10^6 + 1). "Mean LN" = mean RAB7A LogNormalize per cell. "Detection" = fraction of cells with RAB7A UMI > 0.

| Population | Sample | Phenotype | Donor | Cells | RAB7A UMI | Denominator UMI | Pseudobulk | Mean LN | Detection |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| v1 | F18 | Black | F18 | 1,436 | 2,568 | 16,368,346 | 7.303 | 0.755 | 0.699 |
| v1 | F59 | Black | F59 | 2,539 | 4,101 | 26,763,216 | 7.269 | 0.730 | 0.681 |
| v1 | F31B | Black | F31 | 3,012 | 4,061 | 38,667,392 | 6.728 | 0.528 | 0.605 |
| v1 | F31W | White | F31 | 2,487 | 3,457 | 31,347,753 | 6.798 | 0.549 | 0.593 |
| v1 | F62B | Black | F62 | 3,395 | 4,751 | 34,335,908 | 7.123 | 0.681 | 0.646 |
| v1 | F62W | White | F62 | 2,345 | 3,400 | 25,808,204 | 7.052 | 0.668 | 0.667 |
| v2 | F18 | Black | F18 | 179 | 557 | 3,459,084 | 7.340 | 0.827 | 0.832 |
| v2 | F59 | Black | F59 | 820 | 2,123 | 13,185,772 | 7.340 | 0.851 | 0.861 |
| v2 | F31B | Black | F31 | 547 | 1,013 | 9,921,780 | 6.688 | 0.545 | 0.691 |
| v2 | F31W | White | F31 | 247 | 387 | 4,417,088 | 6.469 | 0.488 | 0.644 |
| v2 | F62B | Black | F62 | 1,138 | 2,179 | 16,241,481 | 7.079 | 0.705 | 0.746 |
| v2 | F62W | White | F62 | 384 | 851 | 6,897,587 | 6.959 | 0.712 | 0.818 |

**Notes on the per-sample values**
- F18 and F59 are Black-only donors, so they provide context, not a within-donor contrast.
- The paired libraries (F31, F62) sit lower on both metrics than F18 and F59. This is a between-sample observation and is not interpreted further.

## 4. Paired donors: White minus Black (E2; categories are the script's)
| Population | Donor | Cells (B / W) | Δ pseudobulk | Category | Δ mean LN | Category | Metric disagreement | Δ detection* |
|---|---|---|---:|---|---:|---|---|---:|
| v1 | F31 | 3,012 / 2,487 | +0.070 | white higher | +0.022 | white higher | no | −0.011 |
| v1 | F62 | 3,395 / 2,345 | −0.070 | white lower | −0.012 | white lower | no | +0.022 |
| v2 | F31 | 547 / 247 | −0.218 | white lower | −0.057 | white lower | no | −0.047 |
| v2 | F62 | 1,138 / 384 | −0.120 | white lower | +0.007 | white higher | **yes** | +0.072 |

\*The Δ detection column is simple arithmetic on the E1 detection fractions, done for this report. The script does not compute it or categorise it.

**Observations stated without inference**
1. **v1.** The two donors go in opposite directions, consistently on both metrics: F31 is white higher and F62 is white lower. The magnitudes are small on both metrics.
2. **v2, F31.** White is lower on both metrics.
3. **v2, F62.** White is lower on pseudobulk but higher on mean LogNormalize. This is a metric disagreement.
   - Detection is also higher in F62W.
   - From E1, mean denominator UMI per cell is higher in F62W (about 17,960) than in F62B (about 14,270). This is a descriptive observation; no cause is inferred.
4. **E0 flag.** `all_four_white_lower` is FALSE for both v1 and v2. No population shows White lower in both donors on both metrics.
5. **Zero differences.** None.

**Not reported, by design:** no cell-level P value, no 4-versus-2 unpaired P value, no inferential model and no significance claim. Two paired donors give limited resolution. CCA integration did not batch-correct these RNA-layer values. The values are relative expression, not absolute abundance.

## 5. Figure status
**E3 contents**
- Page 1: per-cell RAB7A LogNormalize violins (width-scaled) by sample, faceted by population.
- Page 2: detection fraction by sample.

**Status:** E3 is **diagnostic only, not submission-ready.**
- The subtitle reads "Descriptive only; no P values".
- It carries no phenotype labelling or project colour scheme.
- Its file name says DESCRIPTIVE rather than DIAGNOSTIC; the frozen script was not edited to change this.
- Page 1 was rendered and inspected. Page 2 was not visually rendered here because no local PDF rasteriser handled it.

## 6. Not done (out of scope)
- No changes to the manuscript, EndNote, publication figures or legends, released scripts, or the public dataset.
- The old S2B remains on hold.
- No claim of external validation.
- No package installation, commit, push, upload or project-memory update.
- Script 04 was not run, and script 08 was run once only.
- The next step is independent review of these computed values.
