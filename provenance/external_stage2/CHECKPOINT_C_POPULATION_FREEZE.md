# Checkpoint C: population freeze, awaiting independent review (2026-10-01)

**Status:** the population is frozen. **No RAB7A value has been computed.** Script 04 (outcome) is written and its SHA-256 is recorded in the freeze. It refuses to run without `--reviewed`, and it has not been run.

**Plan:** `ANALYSIS_PLAN_2026-10-01.md`, SHA-256 `9be9163a…0a4c33`, frozen at 2026-10-01T05:45:31Z (`outputs/B0_plan_freeze.txt`). The plan was not changed afterwards. Two post-freeze script edits are logged in B0, both made before any clustering output existed:
- a technical memory limit for `future`;
- descriptive per-criterion QC counts.

**Freeze record:** `outputs/C9_population_freeze.json` (05:54:43 UTC). Population cell-ID SHA-256 `688dbc3d…97cbc`; keratinocyte object SHA-256 `8b1593f2…74b4fa`; outcome script SHA-256 `01d73890…b3e760`.

## 1. QC (internal thresholds unchanged; `outputs/C1`)

| Sample | Cell Ranger cells | ≥200 genes | mt ≥ 10% | genes ≤ 500 | After QC |
|---|---:|---:|---:|---:|---:|
| F18 | 2,945 | 2,417 | 745 | 443 | 1,462 |
| F59 | 3,140 | 3,080 | 344 | 233 | 2,583 |
| F31B | 8,729 | 8,724 | 4,639 | 1,817 | 3,127 |
| F31W | 7,799 | 7,781 | 4,291 | 1,404 | 2,606 |
| F62B | 4,973 | 4,865 | 890 | 764 | 3,467 |
| F62W | 4,165 | 3,679 | 1,087 | 501 | 2,376 |

- The criteria overlap, so the columns do not sum to the removals.
- **For reviewers:** `percent.mt < 10` removes about 53–55% of cells in the two v3 (F31) libraries, against 11–31% elsewhere. This follows the frozen internal threshold. It is not changed here, but it reduces the F31 pair to about 3,100 and 2,600 cells.

## 2. Annotation (plan §4; `outputs/C2`, `C3`, `C4`)

- Total after QC: 15,621 cells in 21 integrated clusters.
- The deterministic rule assigns 18 clusters (15,214 cells) to keratinocyte lineages and 3 to non-keratinocyte lineages: myeloid (cluster 16), melanocyte (17), T/NK (20).
- No fibroblast, endothelial or mural cluster arises. This fits Wu et al.'s isolated-hair-follicle preparation.
- Marker genes absent after `min.cells = 3`: CLDN5; LYVE1 and CCL21; PDGFRB; TPSAB1 and CPA3; DCD and PIP. Each set kept at least one gene.
- Smallest top-vs-second margins: cluster 9 (0.12, cortex vs matrix), 2 (0.23), 19 (0.39). All three are keratinocyte-vs-keratinocyte or keratinocyte-vs-gland, so keratinocyte membership is unaffected.
- No manual override was made.

## 3. Population rule result (plan §5.2; `outputs/C5`–`C8`)

**Status: `rule_non_discriminating_all_keratinocytes`.** All 20 keratinocyte subclusters (resolution 0.6) have a cluster-mean (KRT35 + KRT85)/2 above 0.5:

| Cluster means | Clusters (dominant lineage from §2) |
|---|---|
| 3.1 – 5.0 | 9, 2, 17, 14, 12 (cortex/cuticle) |
| 1.36 – 1.79 | 3, 15, 10 (mixed basal/IFE; IRS/cortex; IRS) |
| 0.66 – 1.02 | the remaining 12 (basal/ORS, IFE) |

- As the plan requires, the population is **all keratinocytes, 15,214 cells**. No fallback, top-k or alternative threshold was applied.

| Sample | Phenotype | Donor | Population cells |
|---|---|---|---:|
| F18 | Black | F18 | 1,436 |
| F59 | Black | F59 | 2,539 |
| F31B | Black | F31 | 3,012 |
| F31W | White | F31 | 2,487 |
| F62B | Black | F62 | 3,395 |
| F62W | White | F62 | 2,345 |

**Interpretation for reviewers (no outcome involved).**
- The legacy finding that the > 0.5 rule "selects every cluster" (Stage 1: all 21 clusters on the undocumented input) is **reproduced on genuine counts**. It was therefore not only an input-scale artefact.
- Even basal/ORS and IFE clusters show a mean of 0.66–1.0. This suggests KRT35/KRT85 background across all keratinocytes in this dataset, possibly ambient hair-keratin RNA, which was not tested here.
- The internal 0.5 threshold does not transfer as an "enriched" selector. Consequently:
  - Any result run under this frozen plan is a result for **all keratinocytes**.
  - It must not be described as a "KRT35/KRT85-enriched keratinocyte" result.

## 4. Decision needed at this checkpoint (not taken by the author)

1. **Run script 04 under the frozen plan.** The outcome would be RAB7A in all keratinocytes; both paired donors are reported; primary per-donor pseudobulk; no P value for the pairs.
2. **Write a new dated plan, v2,** before any RAB7A value, with a different population rule.
   - Its rule would be chosen after seeing the KRT35/85 score distribution above, but still blind to RAB7A.
   - It must disclose this, cite this checkpoint and be frozen and reviewed in the same way.
   - The v1 result (option 1) should then also be reported, so that the choice between rules cannot be outcome-driven.
3. **Stop and exclude** the external cohort from the Letter with the Stage 1 explanation. This Stage 2 evidence can be added: identity confirmed, input documented, internal rule non-discriminating in this dataset.

Under any option, both paired donors are retained whatever the result direction. No figure, manuscript, EndNote or released file has been changed.
