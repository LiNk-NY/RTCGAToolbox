## Changes in version 2.44.0

### Significant user-visible changes

* Clinical data from `getFirehoseData()` now use standard patient
barcodes as row names (`TCGA-OR-A5J1`) instead of the lower-case,
dot-separated names made by `read.delim()` (`tcga.or.a5j1`). They now match
the sample barcodes used elsewhere, e.g. in `biocExtract()` output and
`TCGAutils::TCGAbarcode()`. Function arguments and the `clinical` slot are
unchanged. `updateObject()` converts objects saved with earlier versions,
and the `accmini` example data were updated this way; only the clinical
row names and the `Composite Element REF` column changed. Code that looks
patients up by the old names can use `toupper(gsub(".", "-", x, fixed =
TRUE))` to convert them.

* `biocExtract()` now always returns a `RaggedExperiment` for copy number
(`CNASNP`, `CNVSNP`, `CNASeq`, `CNACGH`) and mutation (`Mutation`) data.
Previously the class depended on the data: a single sample, or samples
sharing the same ranges, gave a `RangedSummarizedExperiment` instead. Use
`RaggedExperiment::compactAssay()` for a matrix when the ranges are shared.
`getFirehoseData()` and `getData()` still return these data as tables.

* `getFirehoseData(CNACGH = TRUE)` now also downloads the MSKCC CGH
segmentations (`cgh_1x1m_g4447a` for OV and LUSC; `hg_cgh_244a` from MSKCC
for GBM and GBMLGG), which were previously skipped. They come back as
additional `FirehoseCGHArray` elements; `biocExtract()` names them
`CGH_1x1m_g4447a_mskcc` and `CGH_hg_244a_mskcc`. The HMS platforms keep
their names. CNACGH files are now cached as
`<runDate>-<dataset>-CNACGH-<platform>.txt` instead of numbering them by
position; files cached by earlier versions are downloaded again.

* The `getUUIDs` argument of `getFirehoseData()` is deprecated. The
barcode to UUID service it used (`tcga-data.nci.nih.gov`) was retired in
2016, so `getUUIDs = TRUE` could only fail, after all downloads had
finished. It now gives a deprecation warning up front and is otherwise
ignored. Use `TCGAutils::barcodeToUUID()` instead; see
`help("RTCGAToolbox-deprecated")`.

### Bug fixes and minor improvements

* `getFirehoseData()` now names the cached files of the Methylation,
mRNAArray, miRNAArray, RPPAArray and RNASeq2GeneNorm data types after
their archive (e.g., `<runDate>-<dataset>-Methylation-humanmethylation450__jhu_usc_edu.txt`)
instead of numbering them by position. The numbers counted only archives
under `fileSizeLimit`, so a later call with a larger limit could read one
platform's cached file under another platform's name (e.g., U133A values
returned as the Agilent array). Files cached by earlier versions are
downloaded again.

* Internal converters used by cBioPortalData: ranged data with a
`names.field` (e.g., mutations) and the same ranges in every sample no
longer fail with "'names.field' is not found"; the names label the ranges
and are not an assay. `.convertRangeBioc()` gains `ragged = TRUE` to always
return a `RaggedExperiment`. `.makeSummarizedExperimentFromDataFrame()`
works with a single sample column, a plain `data.frame`, and without a
`colnames` argument (sample columns starting with "TCGA").

* `biocExtract()` now aligns rows by name when it combines array files
that share a platform name (e.g., Agilent G4502A 07_1, 07_2 and 07_3).
Previously the files were combined by position, so values could be
assigned to the wrong gene when the row order differed between files.
* `biocExtract()` no longer drops array files that share a platform name
but have different rows (or overlapping samples); it returns them as
separate elements (e.g., `TX_g4502a`, `TX_g4502a_1`) with a message.
Previously such data were silently lost.
* `makeSummarizedExperimentFromGISTIC()` and `biocExtract()` on GISTIC
data now return numeric assays: double for `AllByGene`, integer for
`ThresholdedByGene` and `Peaks` (double if any value is not whole).
Previously the assays were character matrices. Values that are not
numbers now raise an error instead of being coerced.
* `getFirehoseData()` no longer drops genes that share a gene symbol
(e.g., `SLC35E2|728661` and `SLC35E2|9906` in RNA-seq data). Such rows
keep their full `SYMBOL|ENTREZID` identifier as the row name; all other
row names are unchanged. Previously only the first gene was kept.
* `getFirehoseData(RNAseq2Norm = ...)` with a non-default value (files from
`mRNAseq_Preprocess`) no longer drops the first gene of the matrix.
* Expression data with a single gene or a single sample are now read
correctly instead of failing.
* Missing values written as `null` in array data become `NA` without a
coercion warning; any other non-numeric value is now an error.
* Standardisation of dotted TCGA barcodes (e.g., `TCGA.OR.A5J1.01A`) no
longer depends on the random number generator; previously, vectors mixing
dotted and dashed barcodes were converted or not at random, and the
check consumed random numbers.
* `biocExtract()` on segmented copy-number data in which every sample has
the same ranges now aligns each sample's values by range instead of by row
position. Previously, values were assigned to the wrong segment when
samples listed their segments in different orders, or when a row with a
missing start or end was present. Data with a single value column (e.g.,
only `Segment_Mean`) no longer fail.
* `biocExtract()` accepts a bare matrix, a plain `data.frame`, or a
`FirehoseGISTIC` object, not only a `FirehoseData` object. Non-ranged
tables with non-numeric columns now raise an error instead of being
converted to factor codes.
* `biocExtract(x, "CNASNP")` and `biocExtract(x, "CNVSNP")` on a
`FirehoseData` object now record the genome build (`GRCh37`) on the ranges.
`getFirehoseData()` only downloads the hg19 SNP6 segmentation files.
* GISTIC `Peaks`: the empty, unnamed trailing column (e.g., `V100`) is
removed by `getGISTICPeaks()`; the `SummarizedExperiment` row names are now
the peak names (e.g., `Amplification Peak  1`) instead of row indices; and
peaks that share the same limits are no longer dropped.
* Clinical data with a single patient now come back as one row; the
values were previously spread over one row per attribute.
* `getFirehoseData(GISTIC = TRUE)` now rejects a `gistic2Date` that is
not one of `getFirehoseAnalyzeDates()` before downloading anything.
* `makeSummarizedExperimentFromGISTIC()` works with single-sample data.
* `getData()` gives the "valid 'platform' index" error (and lists the
available platforms) for a zero, negative or fractional `platform`.
* Data with more than one genome build in the `NCBI_Build` column (e.g.,
`36` and `37`) now fail with "Inconsistent genome builds found" instead of
an obscure internal error. Equivalent labels (`37`, `GRCh37`, `hg19`) are
treated as one build.
* `getFirehoseData(clinical = TRUE)` no longer returns the
`Composite Element REF` column, which only repeated the word "value".
* `getFirehoseData(Mutation = TRUE)` stacks per-sample MAF files with
different columns by column name (missing columns are `NA`) and keeps file
order. Previously they were combined with `merge()`, a join on all shared
columns, which re-sorted rows and could collapse rows or copy values
between rows when the shared columns were identical.
* `getMutationRate()` finds genes and samples by the `Hugo_Symbol` and
`Tumor_Sample_Barcode` column names instead of column positions 1 and 16,
returns an empty `data.frame` with a warning (instead of `NULL`) when there
are no mutation data, and documents that the denominator is the number of
samples with at least one mutation.

## Changes in version 2.34.0

### Significant user-visible changes

* `getReport` has been removed from the package.

### Bug fixes and minor improvements

* Remove unnecessary imports e.g., from `RCircos`.

## Changes in version 2.32.0

### New features

* The functions `getCNGECorrelation`, `getDiffExpressedGenes`, and `getSurvival`
have been removed from the package.

## Changes in version 2.30.0

### New features

* Include additional options for `miRNASeqGeneType` and `RNAseq2Norm` inputs,
see `?getFirehoseData`.
* The functions `getCNGECorrelation`, `getDiffExpressedGenes`, and `getSurvival`
have been defunct and will be removed from the package (see
`?'RTCGAToolbox-defunct'`).
* The `RNAseq2Norm` argument in `getFirehoseData` allows additional
options: "RSEM_normalized_log2", "raw_counts", "scaled_estimate" from the
'preprocessed' tarballs in Firehose. The "normalized_counts" default remains.

### Bug fixes and minor improvements

* `getFirehoseData` when used with the `miRNASeqGene` argument was
downloading read counts data rather than RPM. This has been fixed with the
`miRNASeqGeneType` argument. "read_count" and "cross-mapped" data are still
available but must be entered explicitly in `miRNASeqGeneType`.

## Changes in version 2.28.0

### New features

* Resolve disparate columns in mutation files from `GBMLGG` (@andreagrioni, #45)
* Update codebase to download https resources from newer layout in
<gdac.broadinstitute.org> (@biostars-nyc, #44).
* Update `makeSummarizedExperimentFromGISTIC` interface with `rownameCol` input,
type checking, and set row names when the are not duplicated.

### Bug fixes and minor improvements

* Use cache directory obtained from `tools::R_user_dir`
* Clean up class membership checks, e.g., with `is(x, "classname")`
* Set missing rows to `""` for downstream compatibility with
`SummarizedExperiment`.

## Changes in version 2.24.0

### Bug fixes and minor improvements

* The deprecated functionality vignette moved to the `vignettes/analysis/`
folder

## Changes in version 2.22.0

### New features

* `makeSummarizedExperimentFromGISTIC` transferred from `TCGAutils`

### Bug fixes and minor improvements

* `biocExtract` now merges datasets from the same platforms
* Re-worked and simplified `getData` method for `FirehoseData` and
`FirehoseGISTIC`
* Remove missing ranges when creating `GRanges` and `GRangesList` from
`DataFrame`

## Changes in version 2.20.0

### New features

* Added the `RNASeq2Gene` slot to the `FirehoseData` class. This data type
mainly obtains RNASeq v2 `raw_counts` from the pipeline (`scaled_estimates`
also available; @mherberg #39)
* Added an `accmini` example dataset as obtained from `getFirehoseData`
* `getLinks` function shows the user some file provenance based on data
requested
* Newly deprecated functions: `getDiffExpressedGenes`, `getCNGECorrelation`,
`getSurvival` and `getReport`. It is no longer possible for the maintainer to
update these functions in a way that would benefit users. A transfer of
responsibility would be required, i.e. to another package.
* Vignettes are updated to reflect changes in the codebase.

### Bug fixes and minor improvements

* Improvements to internal functions for converting tabular data to
Bioconductor classes
* Missing (NA) `seqnames` are removed when converting to `RaggedExperiment`
* Remove static file dependencies from GitHub and use text inside function
(@DavisWeaver, #34)
* Added default values to helper for making `SummarizedExperiment` datasets
* Coerce sample names to character when in (the rare) case they're numeric
* Added an ellipsis argument to `biocExtract` for specifying the `names.field`
in tabular data that will correspond to the row names of a
`SummarizedExperiment`

## Changes in version 2.18.0

### New features

* Warning for Windows users added when file paths are too long
* `getGISTICPeaks` now requires a `FirehoseGISTIC` data object obtained from
`getFirehoseData`

### Bug fixes and minor improvements

* Consolidate GISTIC data download methods in `getFirehoseData` and
`getGISTICPeaks`
* Increase robustness of internal helper functions that work with strands
* 'TCGA' sample column identification is less strict

## Changes in version 2.16.0

### New features

* `RNASeq2GeneNorm` slot in the `FirehoseData` class is a `list` now (from
`matrix`)
* Use `tempdir()` as the default directory for downloading data in
`getFirehoseData`

### Bug fixes and minor improvements

* Save all `RNASeq2GeneNorm` datasets within the output object as a list.
Previously, only the last dataset would get returned (#30)
* Read files from the appropriate download location in `getFirehoseData`
* Move static text file references from 'canevolve.org' to GitHub hosted
locations
* Check file sizes using `httr` instead of 'canevolve.org' query
(@mksamur, #32)

## Changes in version 2.14.0

### Bug fixes and minor improvements

* GISTIC data for `SKCM` and `LAML` is correctly returned by the
`getFirehoseData` function
* Helpers now correctly assign row names in presence of ranged data when
using `biocExtract`

## Changes in version 2.12.0

### New features

* `biocExtract` can now return GISTIC peak data as well with the 'GISTICP'
option
* Methylation data is now converted to `DelayedArray` and contained in
`SummarizedExperiment` objects via `biocExtract`
* Added an `isEmpty` method to the `FirehoseGISTIC` data class

### Bug fixes and minor improvements

* Set default 'peak' argument to 'wide'
* Functionality from the `TCGAutils` package is now employed such as
`findGRangesCols` and `uniformBuilds`

## Changes in version 2.9.41

### Bug fixes and minor improvements

* Fixed bug where file rename was not working due to incorrect path location
in `getFirehoseData`

## Changes in version 2.9.40

### New features

* `getBroadSubtypes` and `getGISTICPeaks` now included in the package.
Thanks to @lgeistlinger!
* Updated import directives
* Added peaks data slot to the `FirehoseGISTIC` class
* Numerous improvements to internal helpers
* Added a `GISTIC` argument to `getFirehoseData`

## Changes in version 2.7.44

### New features

* `selectType` function replaces all previous data type extractor functions
* `getData` improved, removed `CN` argument
* Cleaner documentation for `biocExtract`

### Bug fixes and minor improvements

* Included helper functions for enabling full functionality of `biocExtract`
* Added appropriate import directives
* Re-organized documentation for readability

## Changes in version 2.7.21

### New features

* `biocExtract` function now available - convert data objects to Bioconductor
friendly data classes

## Changes in version 2.7.21

### New features

* New maintainer for package
* Updated API - see `FirehoseData` slots for changes
* New extractor functions for slots of `FirehoseData`
* Standardize argument names and slots
* Export functions in package
* Updated NAMESPACE

### Bug fixes and minor improvements

* Cleaner documentation
* Import directives appropriate
* Example dataset updated to new API
