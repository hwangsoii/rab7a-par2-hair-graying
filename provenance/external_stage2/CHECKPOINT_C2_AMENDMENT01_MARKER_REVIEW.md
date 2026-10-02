# Checkpoint C2 (Amendment 01): marker-only annotation-consistency review of the v2 candidate (2026-10-01)

**Status.** This is a marker review only.
- **No RAB7A or F2RL1 value was extracted or compared.** Script 06 extracts only the 25 predefined review markers, and an assertion in code excludes the outcome genes.
- This is an annotation-consistency review under amendment §3. It is not independent biological validation, and it is not approval.
- No RAB7A or F2RL1 outcome may be computed until the independent review record approves the exact versions (amendment §6).

**Analyst conclusion: `provisionally_consistent`, with recorded caveats (§5).**
- The whole, fixed candidate is consistent with the broad label "hair-shaft differentiating keratinocytes with cortex/cuticle-marker enrichment".
- It is **not** a pure cortex/cuticle population: original cluster 9 is a proliferating, matrix/TAC-like component with early hair-keratin expression.
- Reviewers may weigh cluster 9 differently. If they judge it a material contradiction, the conclusion becomes `unresolved` and the work stops, with no membership change.

**Inputs and integrity**

| Item | Value |
|---|---|
| Amendment (frozen 08:27:34Z, before any V2 output) | SHA-256 `0d8f368a…af84` (B1) |
| Checkpoint-A enforcement | A4 PASS: 120/120 recorded criteria, 42/42 author-manifest hashes, 43 downloads re-hashed unchanged |
| Integrated object | SHA-256 `63ca4077…ffb59`; identical before and after the run |

**Evidence files** (aggregates only, `outputs/`):
- `V2_M0_marker_review.json`: hashes and log;
- `V2_M1`: markers by integrated cluster;
- `V2_M2`: markers by lineage group;
- `V2_M3`: candidate vs other keratinocytes, per sample;
- `V2_M4`: candidate cluster × sample;
- `V2_M5`: co-detection profiles;
- `V2_M6`: dot plots;
- `V2_M7`: S1 membership.

In the tables below, values are mean LogNormalize / detection fraction.

## S1. Membership (from the existing label)

**Reconciliation**
- The object's `lineage` agrees with C2 `top_lineage` and with the cluster sizes for **21/21** integrated clusters.
- `lineage == "KC_cortex_cuticle"` gives original integrated clusters **3, 6, 9 and 14**: **3,315 cells**.
- The candidate is a subset of the frozen v1 cell list (C9 hash verified) and contains keratinocytes only.
- Every per-sample count equals the audit value. Nothing was forced; the values were derived and then compared.

| Sample | Candidate cells | v1 keratinocytes | Candidate share |
|---|---:|---:|---:|
| F18 | 179 | 1,436 | 12.5% |
| F59 | 820 | 2,539 | 32.3% |
| F31B | 547 | 3,012 | 18.2% |
| F31W | 247 | 2,487 | 9.9% |
| F62B | 1,138 | 3,395 | 33.5% |
| F62W | 384 | 2,345 | 16.4% |

All four clusters are retained, including cluster 9.

## S2. Hair-keratin profiles

**Pooled by lineage label** (`V2_M2`):

| Marker | Cortex/cuticle candidate (3,315) | Basal/ORS (8,135) | IFE-suprabasal label (2,939) | IRS (825) |
|---|---|---|---|---|
| KRT31 | **2.42 / 0.86** | 0.66 / 0.56 | 0.72 / 0.60 | 0.85 / 0.68 |
| KRT35 | **4.07 / 0.98** | 0.78 / 0.64 | 0.82 / 0.70 | 1.09 / 0.77 |
| KRT85 | **4.60 / 1.00** | 1.11 / 0.78 | 1.15 / 0.81 | 1.60 / 0.90 |
| KRT32 | **1.19 / 0.70** | 0.12 / 0.13 | 0.13 / 0.14 | 0.27 / 0.27 |
| KRT33A | 0.00 / 0.01 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 |

**By cluster** (`V2_M1`, `V2_M4`). These are complementary profiles, matching Wu et al.'s cortex (KRT31) versus fibre-cuticle (KRT35 without KRT31) distinction:

| Cluster | Cells | Profile | Key values |
|---|---:|---|---|
| **6** | 926 | **Cortex-like** | KRT31 5.13/0.99, KRT85 5.68, KRT35 4.28, KRT32 0.59 |
| **14** | 210 | **Cuticle-like** | KRT32 3.84/1.00, KRT35 3.85, KRT85 3.82, KRT31 0.92 (near background); KRT25 1.92/0.93 |
| **3** | 1,430 | KRT35/KRT85-dominant, intermediate KRT31 and KRT32 | KRT35 4.62, KRT85 4.83, KRT31 1.71, KRT32 1.45 |
| **9** | 749 | **Intermediate hair-keratin level** | KRT35 2.85/0.95, KRT85 3.06/0.99, KRT31 0.85 (background) |

- In cluster 9, KRT35 and KRT85 are about 2–4× the non-candidate cluster means (0.63–1.25 and 0.82–1.97) but below clusters 3, 6 and 14 (see S3).
- Each cluster's profile repeats in every sample, including F18's 10 cells in cluster 9 and its 12 in cluster 14 (`V2_M4`).
- **KRT33A** is essentially undetected everywhere (≤ 0.02 detection) and provides no support. Not every gene of the five-marker set supports the annotation.

**Important limitation.** KRT35 and KRT85 are *detected* in 54–78% (KRT35) and 59–92% (KRT85) of cells in every non-candidate keratinocyte cluster, at low mean levels (0.6–2.0). Detection fractions are therefore nearly saturated, and the hair-keratin contrast lies mainly in **mean expression**.
- This background matches the v1 finding that all keratinocyte subclusters exceed 0.5.
- Ambient hair-keratin RNA is a possible explanation, not a tested one.

## S3. Competing profiles (patterns, not purity)

**Basal/ORS profile.** The candidate is low in the basal/ORS keratins KRT5, KRT14, KRT16, KRT17 and KRT6A:
- candidate means 0.4–1.7, against 2.6–5.4 in the basal/ORS group;
- detection stays high: KRT14 0.94–0.98 in every candidate cluster. This is consistent with the same ambient-like background.

**IRS profile.** TCHH is read together with GATA3 and the IRS keratins.
- The candidate is low in KRT27, KRT71, KRT73, TCHH and GATA3: TCHH 0.26/0.32 and GATA3 0.12/0.16, against 3.09/0.83 and 1.22/0.65 in IRS.
- Cells co-detecting GATA3 and TCHH: 0.04–0.14 in candidate clusters, against 0.65 in IRS cluster 10.
- **Elevated KRT25 is the exception:** 1.92/0.93 in cluster 14; moderate (0.62–0.85) in clusters 3 and 9.
  - KRT25 is listed under IRS in plan §4, but cluster 14 lacks the other IRS markers.
  - This is recorded as mixed evidence, not reassignment.

**Matrix/progenitor profile.** MSX2, LEF1 and HOXC13 are elevated in **all four** candidate clusters:
- MSX2 0.64–0.97 (detection 0.67–0.86), against ≤ 0.11 in the basal/ORS and IFE-label clusters;
- also present in the IRS clusters 10 and 15.
- Per Wu et al., cortex, cuticle and IRS lineages derive from MSX2+ matrix TACs, so this pattern is compatible with matrix-derived hair-shaft differentiation, not a contradiction.
- **Cluster 9 is proliferating:** MKI67 0.68/0.58, against ≤ 0.02 in clusters 3, 6 and 14. The only comparable cluster is basal/ORS cluster 12, at 0.62/0.53.
  - Cluster 9 also has slightly higher KRT15 (1.02/0.65), KRT5 (1.42) and GATA3 (0.30/0.30) than the other candidate clusters.
  - It is best described as a **matrix/TAC-like, early hair-keratin-expressing component**. It may include precursors of more than one matrix lineage, which was not resolved here.
  - Its cortex-versus-matrix module margin is small (0.1244). The marker profiles are consistent with that ambiguity.

**KRT1/KRT10.**
- KRT10 is broadly detected, including 1.34–1.50 / 0.90–0.97 in the candidate. This is consistent with Wu et al.'s upper/middle-ORS KRT10, and it is not evidence of IFE contamination.
- KRT1 is low in the candidate (0.07/0.09).

**Co-detection (`V2_M5`).**
- (KRT35 or KRT85) together with KRT14 is 0.80–0.98 in **every** keratinocyte cluster, candidate or not.
- Because both genes are near-ubiquitously detected, binary co-detection cannot distinguish genuine mixed-lineage cells from background, so it is recorded but not interpreted.
- (KRT35 or KRT85) together with (MSX2 or LEF1) is 0.74–0.96 in candidate clusters, against 0.02–0.11 in basal/ORS and IFE-label clusters (0.70–0.78 in IRS).

No cell or cluster was deleted.

## S4. Sample representation and discordance

**Representation**
- All six samples contribute: 179–1,138 cells.
- **Sparse groups:**
  - F18 has 10 cells in cluster 9 and 12 in cluster 14;
  - F62W has 17 cells in cluster 14;
  - F31W and F31B have 26 and 28 cells in cluster 14.
- These small groups show the same qualitative profiles but have unstable means.

**Within-sample contrast** (`V2_M3`). In every sample, candidate KRT35 and KRT85 means exceed that sample's other keratinocytes:

| Sample | KRT35, candidate vs other | KRT85, candidate vs other |
|---|---|---|
| F18 | 4.52 vs 0.52 | 5.35 vs 0.73 |
| F59 | 4.34 vs 0.96 | 4.78 vs 1.09 |
| F31B | 3.48 vs 1.09 | 4.29 vs 1.77 |
| F31W | 3.78 vs 0.56 | 4.28 vs 1.04 |
| F62B | 4.12 vs 0.98 | 4.57 vs 1.24 |
| F62W | 4.18 vs 0.61 | 4.64 vs 0.74 |

- No per-sample pass rule was applied.

**Discordance**
- **F31B (3′ v3)** has the highest background: other keratinocytes KRT35 0.83 detection, KRT85 0.98, KRT31 1.10. It also has the smallest candidate-versus-other contrast.
- **F31B and F31W** have higher KRT5 and KRT14 within the candidate (1.6–1.8 and 2.0–2.1) than the v2 libraries (0.7–0.9 and 1.4–1.7).
- These are recorded as possible chemistry or ambient differences, not interpreted.

**Composition**
- The candidate share of keratinocytes is lower in the White than in the Black sample of both pairs: F31 18.2% → 9.9%; F62 33.5% → 16.4%.
- This is a **population-size observation**, not an expression outcome. It is relevant context for later interpretation, and it was not used to define membership.

## 4. Naming and relation to other populations

- **Provisional name, if approved:** "hair-shaft differentiating keratinocytes with cortex/cuticle-marker enrichment". The population comprises:
  - a cortex-like component (cluster 6);
  - a cuticle-like component (cluster 14);
  - a KRT35/KRT85-dominant component (cluster 3);
  - a proliferating matrix/TAC-like component with early hair-keratin expression (cluster 9).
- **Not claimed:**
  - exact equivalence to the internal-cohort KRT35/KRT85 population;
  - equivalence to Wu et al.'s CO or FC clusters;
  - that every marker supports the label.
- The labels come from module scores. UMAP proximity was not used.
- v2 ⊂ v1, so they are not independent replications.

## 5. Conclusion and caveats (analyst; not independent approval)

**Conclusion: `provisionally_consistent`**
- S1 is reconciled exactly.
- S2 shows complementary cortex-like and cuticle-like hair-keratin enrichment by mean expression, in every cluster and every sample.
- S3 shows the basal/ORS and IRS profiles depleted, and a matrix-associated pattern compatible with hair-shaft differentiation.
- S4 shows all samples represented, with sparse subgroups and F31-specific background recorded.

**Caveats for the reviewers**
1. **Cluster 9** (749 cells, 22.6% of the candidate) is proliferating and matrix/TAC-like. It has intermediate hair-keratin levels and background-level KRT31, and its identity is mixed. It is retained, because membership is fixed. Any result on v2 is therefore for the broad hair-shaft differentiation compartment, not for mature cortex/cuticle alone.
2. **Detection fractions** of hair keratins and KRT14 are near-saturated across keratinocytes. Mean-expression contrasts carry the evidence, and binary co-detection is uninformative.
3. **KRT33A** gives no support. **KRT25** in cluster 14 is mixed evidence.
4. **F31 libraries (3′ v3)** show higher hair-keratin background and higher basal keratins within the candidate.
5. **The candidate share** of keratinocytes is lower in the White samples of both pairs. That is composition, not expression.

On this conclusion, `07_freeze_v2_population.R` may freeze v2 (C10) together with an **unapproved** review template. No outcome is computed.
