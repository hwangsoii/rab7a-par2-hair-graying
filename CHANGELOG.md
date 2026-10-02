# Corrected release, 2 October 2026

- Replaced the undocumented-input, 2,932-cell external S2B with a separately
  reconstructed count-based analysis: 3,315 selected cells plus all-KC context.
- Corrected all six external sample/donor links using the official mapping.
  Black/white phenotype assignments did not change.
- Replaced external source values, sample summaries, selection descriptions
  and plotting code; retained zero values and all samples.
- Added nominal cell-level Wilcoxon statistics and clearly separated them
  from donor-level replication. Preserved contrasting paired summaries.
- Recorded selection chronology, the marker-only amendment and the subsequent
  addition of nominal P values after descriptive outcome inspection.
- Preserved dated count-reconstruction code/plans as provenance. No old h5ad
  object or legacy top-quartile selection is used in the corrected panel.
- Removed the obsolete external branch from script 04 while leaving internal
  panel calculations unchanged. Added current internal Letter replot scripts.
- Made script 05's QA file paths portable for externally configured data and
  output folders. Its cell values, statistical calculation and plot are unchanged.
- Corrected public repository and embargo wording. The internal matrices,
  coded participant metadata, original manuscripts and historical release
  are unchanged. This GitHub revision changes only the code repository and
  derived external source-data tables, not the Mendeley deposit or embargo.
- Retained the previous release in Git history. Superseded the three legacy
  `Figure_S2B_*` TSV files at the current revision with the documented tables
  under `data/external/`; they should not be mixed across release versions.
