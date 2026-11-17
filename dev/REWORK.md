# RTCGAToolbox rework

Working notes for the ground-up rework. Not part of the built package
(`^dev$` is in `.Rbuildignore`).

## Ground rules

1. **Tests first.** `tests/testthat/` pins the behaviour of the code as of
   2.43.1 before anything is changed. The suite runs fully offline
   (~20 s): network calls are mocked and `getFirehoseData()` runs end to
   end against a fake Firehose mirror built from real-format archives
   (`helper-fixtures.R`).
2. **Expected values come from the fixture text or the raw slots, never
   from the code under test.**
3. **Known bugs are pinned, not hidden.** Tests tagged `[BUG-nn]` assert the
   current (wrong) behaviour so nothing changes silently. Fixing a bug means
   flipping its test(s) to the correct expectation in the same commit.
4. **One function (or one bug) per commit.** Pure refactors must leave every
   test unchanged and green; behaviour changes must say so in NEWS.
5. **Before removing or renaming any function** (exported or internal),
   query the Bioconductor usage index (see the user's global instructions)
   and record the result here.
6. Prefer base R over heavy dependencies (e.g. drop `stringr`).

## Bug register

Severity: **D** = silent data corruption/loss, **R** = wrong result or
crash on valid input, **M** = misleading/cosmetic.

| ID | Sev | Where | Problem |
|----|-----|-------|---------|
| BUG-03 | D | `.mergePlatforms` | Same-named platforms (e.g. Agilent 07_1/07_2/07_3 -> `TX_g4502a`) are `cbind`-ed **by position**; only `nrow` is checked, not row names. Values land on the wrong gene if row order differs. |
| BUG-04 | D | `.mergePlatforms` / `.unNestList` | Same-named platforms with different row counts are silently dropped: `biocExtract()` returns an empty list. |
| BUG-27 | D | `.makeRangedSummarizedExperimentFromDataFrame` / `.hasConsistentRanges` | Consistency is checked with `setequal()` (order- and multiplicity-blind) but values are `cbind`-ed by row position: segments listed in a different order per sample, or rows with NA start/end (dropped from ranges only), put values on the wrong range. Found while fixing BUG-24. |
| BUG-28 | R | `makeSummarizedExperimentFromGISTIC` (by gene) | A single-sample table fails: `gist[, rel.cols]` drops to a vector, losing the sample name. |
| BUG-29 | R | `.makeRangedSummarizedExperimentFromDataFrame` | With `names.field` (as cBioPortalData passes via `.convertRangeBioc`), the names column is dropped before the ranges are built: "'names.field' is not found". Same in 2.43.1. Fixing it changes what cBioPortalData imports (the name column would otherwise also become an assay), so coordinate with cBioPortalData. |
| BUG-02 | R | `makeSummarizedExperimentFromGISTIC`, `biocExtract(GISTIC*)` | All GISTIC assays (AllByGene, ThresholdedByGene, Peaks) are **character** matrices. Arithmetic fails; sorting/comparison is lexicographic. |
| BUG-05 | D | `.makeExprMat` | Duplicated gene symbols with different Entrez ids (e.g. `SLC35E2|728661`, `SLC35E2|9906`) keep only the first; the other gene's data are dropped without a message. |
| BUG-26 | D | `.makeExprMat` (`RNAseq2Norm != "normalized_counts"`) | `mRNAseq_Preprocess` files have one header row but two rows are removed, so the first data row is lost. Masked in practice because the first row is usually a `?|` row. |
| BUG-06 | R | `.makeExprMat` | Fails with a single gene row or a single sample (`apply()` / `[` drop dimensions). |
| BUG-01 | R | `.stdIDs` | Uses `sample()` to decide whether to convert dotted barcodes; on mixed input the result depends on the RNG state. |
| BUG-24 | R | `.makeRangedSummarizedExperimentFromDataFrame` | Segments with identical ranges and a single value column error (`invalid 'length' argument`). |
| BUG-17 | M | `biocExtract(CNASNP/CNVSNP)` | Data are hg19 (from the file name) but no genome is recorded on the ranges. |
| BUG-10 | M | `getGISTICPeaks`, Peaks SE | The empty trailing column (`V100`/`V12`, values `""`) survives the all-NA filter; Peaks SE row names are leftover data.frame indices instead of peak names. |
| BUG-11 | M | `biocExtract` on segment/MAF data | Output class depends on the data: one sample (or identical ranges) -> RSE, otherwise RaggedExperiment. For single-sample MAFs, all 47 MAF columns become (character) assays. Fixed (decision): `biocExtract()` always gives a RaggedExperiment for CN/Mutation types; getters stay tabular; `.convertRangeBioc(ragged = FALSE)` default kept for cBioPortalData. |
| BUG-07 | M | `getFirehoseData` | `gistic2Date` check uses `&&`: any single string (e.g. `"20991231"`) passes. |
| BUG-09 | M | `getFirehoseData(clinical)` | Header row leaks in as a constant `Composite Element REF` = `"value"` column (fixed); patient ids are now standard barcodes (`TCGA-OR-A5J1`); `updateObject()` converts old objects and `accmini` was regenerated (decision: follow other software). |
| BUG-20 | D | `getFirehoseData(Mutation)` | `merge()` of per-sample MAFs with different columns is a join on all shared columns: rows are re-sorted, and rows identical in the shared columns are collapsed with values copied between them (3 rows -> 2 in the test). Fixed: row union via `rbindlist(fill = TRUE)`. |
| BUG-23 | R | `biocExtract` | Errors on a bare matrix (`.removeShell`) and on a non-ranged plain `data.frame` (`metadata()` on data.frame). |
| BUG-25 | M | `.getBuild` | Mixed builds fail inside `TCGAutils::uniformBuilds()` with a cryptic vapply error; `.validateNCBI()` is unreachable. |
| BUG-12 | M | `getData(..., platform = 0)` | Cryptic subscript error instead of the "valid platform" message. |
| BUG-13 | M | `getMutationRate` | Returns `NULL` silently with no mutation data; samples read from column **16** by position; denominator is mutated samples only (undocumented). |
| BUG-14 | R | `.makeSummarizedExperimentFromDataFrame` | **Not dead**: cBioPortalData calls it via `:::` with a `DataFrame` and `colnames = <sample ids>`, which works. Broken otherwise: the default `colnames = c("Hugo", "Entrez")` treats the annotation columns as samples (assay and rowData swapped), and a plain `data.frame` fails (`metadat` undefined). |
| BUG-15 | — | `.barcodeUUID` / `getUUIDs = TRUE` | Calls `tcga-data.nci.nih.gov`, retired in 2016. Cannot work. Deprecated with `BiocBaseUtils::lifeCycle()` via `.lifeCycleArg()`; `TRUE` warns and is ignored. `.barcodeUUID` is now unreachable (no external users in bioc324-index-20261104); removing it drops RJSONIO and `RCurl::getURL` (`RCurl::url.exists` is still used by `getBroadSubtypes`). |
| BUG-30 | D | `getFirehoseData` (RNASeq2GeneNorm, Methylation, mRNAArray, miRNAArray, RPPAArray) | Cache files are named by position (`-Methylation-<n>.txt`) counting only archives under `fileSizeLimit`. A later call with a smaller limit, or a changed listing, reads another platform's cached file under the wrong `Filename`. Same fix as CNACGH (BUG-08): name by platform. Open. |
| BUG-16 | M | `show,FirehoseData` | `dim(Mutation)[2] & dim(Mutation)[2]` typo (should check rows). |
| BUG-08 | ? | `getFirehoseData(CNACGH)` | Link pattern matches `hg_cgh_244a`/`hg_cgh_415k` (HMS) archives but not MSKCC `cgh_1x1m ... segmentation_data_computation__seg`. Fixed (decision: include): MSKCC `cgh_1x1m_g4447a` (OV, LUSC) and MSKCC `hg_cgh_244a` (GBM, GBMLGG) added; cache files named by platform since positional `CNACGH-<n>.txt` names would map a cached HMS file onto the MSKCC archive once the link list changed. |

Dead code (no callers in this package): `.getGISTIC`, `.fileSelect`,
`.runOnDupElements` (unused since BUG-03/04), `.validateNCBI` (unused
since BUG-25). None is used elsewhere in Bioconductor
(bioc324-index-20261104).

## Internals used by other Bioconductor packages

From bioc324-index-20261104 (all via `:::`, all in cBioPortalData). Their
signatures and semantics are effectively public API; change them only
in coordination with cBioPortalData.

`.convertRangeBioc`, `.findCol`, `.getBuild`, `.hasExperimentData`,
`.hasInfo`, `.hasRangeNames`, `.makeGRangesFromDataFrame`,
`.makeSummarizedExperimentFromDataFrame`, `.standardizeBC`.

Changes so far that reach cBioPortalData: `.standardizeBC` (BUG-01,
deterministic; non-TCGA ids unaffected), `.getBuild` (BUG-25, clear error
on mixed builds; all-NA column now returns NA instead of erroring),
`.convertRangeBioc` (new `defaultBuild` argument before `...`; BUG-27
alignment).

## Proposed order of work

### Phase 1 — characterisation tests (done)

### Phase 2 — refactors with no behaviour change
- [x] `.findBuild`: replace `stringr::str_extract` with base R; drop `stringr`
      from Imports.
- [ ] Factor the Firehose archive patterns into a single table shared by
      `getFirehoseData()` and `getLinks()` (they are duplicated today and
      can drift).

### Phase 3 — correctness fixes (highest risk first; one per commit)
- [x] BUG-03 / BUG-04: merge platforms by row name, never by position;
      never drop data silently.
- [x] BUG-02: numeric GISTIC assays.
- [x] BUG-05 / BUG-26 / BUG-06: rewrite `.makeExprMat` (explicit header
      handling per file type; keep or disambiguate duplicated symbols).
- [x] BUG-01: deterministic `.stdIDs`.
- [x] BUG-24 + BUG-27: align consistent-range segments by range.
- [x] BUG-23: `biocExtract` on matrix / data.frame / FirehoseGISTIC.
- [x] BUG-17: record GRCh37 for CNASNP/CNVSNP from FirehoseData.
- [x] BUG-10: GISTIC Peaks column filter, row names, no dedup by range.
- [x] BUG-07: gistic2Date validation.
- [x] BUG-09 (junk column, patient barcodes), BUG-12, BUG-13, BUG-16, BUG-20, BUG-25, BUG-28.

### Phase 4 — removals (Bioc index check first) and API decisions
- [ ] Remove dead internals.
- [x] Deprecate `getUUIDs` (BUG-15).
- [ ] Remove `.barcodeUUID` and the RJSONIO import (after the deprecation cycle).
- [ ] BUG-30: platform-named cache files for the other list slots.
- [x] BUG-11: RaggedExperiment from biocExtract for CN/Mutation.
- [x] BUG-08: MSKCC CGH archives included.

### Phase 5 — decompose `getFirehoseData()` into per-type readers
