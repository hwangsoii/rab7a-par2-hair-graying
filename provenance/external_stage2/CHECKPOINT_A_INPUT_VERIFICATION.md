# Checkpoint A: verification of the third-party reconstructed counts (2026-10-01)

**Scope.** This report checks whether the Zenodo 15103193 count matrices are a suitable, documented input. It contains **no RAB7A result**. The FASTQs were not downloaded and Cell Ranger was not run. The old h5ad and Seurat objects are not inputs; the local `metadata.tsv` barcodes were read only for the identity cross-check in §4.

**Result: PASS.** All criteria in plan §1 are met. Two documentation discrepancies (§3) are recorded; neither affects the count values.

Evidence: `scripts/00_node_manifest_by_prefix.py`, `scripts/01_verify_zenodo_counts.R` → `outputs/A0_node_manifest_by_prefix.tsv`, `A1_checksums.tsv`, `A2_local_vs_zenodo_barcode_overlap.tsv`, `A3_verification.json`. Data: `OneDrive_2025-09-11/Data/cell discovery/onfroy_zenodo_15103193/download/` (git-ignored; 776 MB; public CC-BY-4.0 deposit).

## 1. Source and documentation

| Item | Record |
|---|---|
| Deposit | Zenodo 15103193, DOI 10.5281/zenodo.15103193. "Count matrices of the OEP00002321 scRNA-Seq dataset", A. Onfroy, v1, 2025-03-28, CC-BY-4.0, linked `isSupplementTo` Wu et al. 2022 |
| Stated method | 48 NODE FASTQs (8 per sample) downloaded and run through `cellranger count`, "version 3.0.1", "hg19" |
| Code | GitHub `audrey-onfroy/Onfroy_12tips_ACM_REP_25`, HEAD `eac531f2d1` (2025-11-12); notebooks start from `raw_feature_bc_matrix`. Paper: Onfroy et al., ACM REP '25, doi 10.1145/3736731.3746138. Code licence CC-BY-NC-4.0 (code not reused here) |
| Per-sample Cell Ranger record | `<sample>_web_summary.html` (Sample ID, chemistry, transcriptome, pipeline version, reads, cells) |

## 2. Integrity and count scale (FULL checks, all six samples)

| Check | Result |
|---|---|
| MD5, 43 files vs Zenodo API | 43/43 match; sizes match |
| MD5 vs the author's `md5sum.txt` (42 listed files) | all match |
| Matrix Market header | `coordinate integer general` (raw and filtered) |
| Every stored value a nonnegative integer | TRUE for every raw and filtered matrix (min nonzero 1) |
| Features | 27,955, all "Gene Expression", identical across raw/filtered and all samples; RAB7A, KRT35, KRT85, KRT31, F2RL1 each present once |
| Raw matrices | Full barcode whitelist: 737,280 columns (v2) or 6,794,880 (v3). These are **raw-droplet** matrices |
| Filtered matrices | Cell Ranger cell calls. Every filtered barcode is in the raw matrix, with **identical** column values. Column count and median UMI per cell equal the web_summary "Estimated Number of Cells" and "Median UMI Counts per Cell" exactly |

| Sample | Chemistry | Reads | Filtered cells | Median UMI | Median genes |
|---|---|---:|---:|---:|---:|
| F18 | 3′ v2 | 144,897,014 | 2,945 | 3,518 | 1,209 |
| F59 | 3′ v2 | 272,290,543 | 3,140 | 7,648 | 2,071 |
| F31B | 3′ v3 | 343,661,650 | 8,729 | 5,188 | 1,409 |
| F31W | 3′ v3 | 302,006,188 | 7,799 | 5,479 | 1,484 |
| F62B | 3′ v2 | 246,368,459 | 4,973 | 5,953 | 1,734 |
| F62W | 3′ v2 | 219,720,015 | 4,165 | 5,278 | 1,619 |

## 3. Documentation discrepancies (recorded, not resolved)

- **Software and reference.** The deposit text says Cell Ranger 3.0.1 and hg19. All six web summaries say pipeline **3.1.0** and transcriptome **"hg19ens91-"**, which looks like a custom hg19 reference with Ensembl 91 annotation. The web summaries are the machine record and are taken as authoritative. Exact reference build files are not deposited.
- **Pipeline used here.** Onfroy's notebooks call cells from the raw matrices with the `aquarius` package, which is not installed. Following the plan, the analysis uses Cell Ranger's own filtered calls, as the internal cohort does.

## 4. Sample identity: three links, each checked independently

**(a) Zenodo label ↔ official FASTQ prefix.** Onfroy labelled the samples by article ID, and the web summaries carry "Sample ID" = F18, F59 and so on. Two independent checks against the NODE manifest (`A0`) agree with `id2FastqID.xlsx`:

| Zenodo sample | Official prefix | FASTQ total (GB) | Reads / GB | Chemistry, consistent with run |
|---|---|---:|---:|---|
| F18 | ryg035 | 22.88 | 6.33 M | v2; 2017-era `ryg` lane-style names |
| F59 | ryg029 | 44.36 | 6.14 M | v2; `ryg` |
| F62B | ryg047 | 39.81 | 6.19 M | v2; BKDL**1707** |
| F62W | ryg048 | 35.56 | 6.18 M | v2; BKDL**1707** |
| F31B | black-1 | 45.37 | 7.57 M | v3; BKDL**2100**, Illumina `_S_L_R_001` names |
| F31W | white-2 | 39.76 | 7.60 M | v3; BKDL**2100** |

- Reads per GB is uniform within each chemistry.
- If F18 and F59 were swapped, the ratios would be 3.3 M/GB and 11.9 M/GB, which are clear outliers.
- The two v3 libraries can only be the 2021-tagged `black-1`/`white-2` runs, since v3 chemistry postdates the 2017 runs.
- This is consistency evidence, not a cryptographic link, but it agrees with the official mapping on every sample.

**(b) Local h5ad sample name ↔ Zenodo sample** (barcode overlap; `A2`). Each local sample overlaps exactly one Zenodo sample, at 93–99% of its cells. Every other pairing overlaps by at most 24 barcodes, which is chance level.

| Local `sample` | Matches Zenodo | Overlap | Local donor label (metadata.tsv / Supplementary Methods) |
|---|---|---|---|
| ryg035 | **F18** | 2,023 / 2,051 | F59 |
| ryg029 | **F59** | 2,931 / 3,006 | F18 |
| ryg047 | **F62B** | 4,543 / 4,677 | F31 black |
| ryg048 | **F62W** | 3,282 / 3,343 | F31 white |
| black1_F62B | **F31B** | 5,483 / 5,902 | F62 black |
| white2_F62W | **F31W** | 5,576 / 5,763 | F62 white |

**Consequences for the Stage 1 findings:**
- **Name link now verified.** The Stage 1 identity table (conditional on the name link) is confirmed by data.
- **Donor labels swapped.** The local donor labels are swapped (F18↔F59, F31↔F62), so the local names `black1_F62B`/`white2_F62W` carry the wrong donor tag.
- **Phenotypes unaffected.** Phenotype labels in the v2 object were correct per sample.
- **Within-pair observation re-attributed.** In the current S2B, "white higher" belongs to donor **F62** (ryg047→ryg048) and "white lower" to donor **F31**.

**(c) The local h5ad was not built from these count matrices.**
- On matched barcodes, the local `total_counts` never equals the Zenodo UMI total (0 of ~23,800 cells), and neither does the gene count.
- The local value is a median 1.16–1.25× higher in every sample.
- The local object therefore came from a **different count run** of the same libraries (different reference or software). This is consistent with the unknown builder recorded in Stage 1 and does not authenticate the old `raw/X`.

## 5. Decision

- The Zenodo filtered matrices are **suitable** as the input for the dated plan. They are genuine nonnegative integer UMI counts, Cell Ranger-called, checksum-verified, documented per sample, and labelled consistently with the official mapping.
- The plan was frozen after this checkpoint: `outputs/B0_plan_freeze.txt`, plan SHA-256 `9be9163a…0a4c33`, 2026-10-01T05:45:31Z.
