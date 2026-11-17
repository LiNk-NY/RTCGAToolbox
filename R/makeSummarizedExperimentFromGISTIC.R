#' Create a SummarizedExperiment from FireHose GISTIC
#'
#' @description Use the output of `getFirehoseData` to create a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class].
#' This can be done for three types of data, G-scores threshold by gene, copy
#' number by gene, and copy number by peak regions.
#'
#' @param gistic A [RTCGAToolbox::FirehoseGISTIC-class] object
#' @param dataType character(1) One of "ThresholdedByGene" (default),
#'   "AllByGene", or "Peaks"
#' @param rownameCol character(1) The name of the column in the data to use as
#'   rownames in the data matrix (default: 'Gene.Symbol'). The row names are
#'   only set when the column name is found in the data and all values are
#'   unique.
#' @param ... Additional arguments passed to 'getGISTICPeaks'.
#'
#' @author L. Geistlinger, M. Ramos
#' @examples
#'
#' co <- getFirehoseData("COAD", clinical = FALSE, GISTIC = TRUE,
#'     destdir = tempdir())
#' makeSummarizedExperimentFromGISTIC(co, "AllByGene")
#'
#' @return A `SummarizedExperiment` object
#' @export
makeSummarizedExperimentFromGISTIC <-
    function(
        gistic,
        dataType = c("AllByGene", "ThresholdedByGene", "Peaks"),
        rownameCol = "Gene.Symbol", ...
    )
{
    if (!is(gistic, "FirehoseData") && !is(gistic, "FirehoseGISTIC"))
        stop("'gistic' must be a 'FirehoseData' or 'FirehoseGISTIC' object")
    dataType <- match.arg(dataType)
    if (identical(dataType, "Peaks")) {
        gist <- getGISTICPeaks(gistic, ...)
        if (!length(gist)) return(list())
        rel.cols <- grepl("^TCGA", colnames(gist))

        gistData <- as.matrix(gist[, rel.cols])
        # get the peak type (amplification / deletion)
        peak.type <- vapply(strsplit(gist[["Unique.Name"]], " "),
            function(x) x[[1L]], character(1L))
        rowdata <- cbind.data.frame(gist[, !rel.cols], type = peak.type,
            stringsAsFactors = FALSE)
        feats <- grepl("gene|ranges", names(rowdata), ignore.case = TRUE)
        rows <- rowdata[, feats]
        if (length(rows)) {
            if (as.logical(anyDuplicated(rows))) {
                uniq <- !duplicated(rows)
                rows <- rows[uniq]
                rowdata <- rowdata[uniq, ]
                gistData <- gistData[uniq, ]
                gist <- gist[uniq, , drop = FALSE]
            }
            rownames(rowdata) <- rows
        }
        rowranges <- gist[["rowRanges"]]
        rowranges <- as(rowranges, "GRanges")
        # get the peak type (amplification / deletion)
        peak.type <- vapply(strsplit(gist[["Unique.Name"]], " "),
            function(x) x[[1L]], character(1L))
        # create the SE
        gistData <- .gisticToNumeric(gistData, integer = TRUE)
        gisticSE <- SummarizedExperiment(gistData,
            rowRanges = rowranges)
        rowData(gisticSE) <- rowdata
    } else {
        gist <- getData(gistic, "GISTIC", dataType)
        if (!length(gist)) return(list())
        rel.cols <- grepl("^TCGA", colnames(gist))
        gistData <- as.matrix(gist[, rel.cols])
        annoteRowDF <- gist[, !rel.cols, drop = FALSE]
        anyduprows <- anyDuplicated(annoteRowDF[[rownameCol]])
        if (rownameCol %in% names(annoteRowDF) && !anyduprows)
            rownames(gistData) <- annoteRowDF[[rownameCol]]
        colnames(gistData) <- .stdIDs(colnames(gistData))
        gistData <- .gisticToNumeric(gistData,
            integer = identical(dataType, "ThresholdedByGene"))
        gisticSE <- SummarizedExperiment(gistData, rowData = annoteRowDF)
    }
    gisticSE
}

## GISTIC tables are read as character; convert the sample columns of the
## assay without losing information. Empty strings and "NA" become NA; any
## other value that is not a number is an error rather than a silent NA.
## With 'integer = TRUE' (thresholded calls, peak calls) the result is
## integer when every value is a whole number, otherwise double.
.gisticToNumeric <- function(m, integer = FALSE) {
    if (is.numeric(m))
        return(m)
    vals <- trimws(m)
    vals[vals %in% c("", "NA")] <- NA_character_
    num <- suppressWarnings(as.numeric(vals))
    bad <- is.na(num) & !is.na(vals)
    if (any(bad))
        stop("Non-numeric values in GISTIC data: ",
            paste(utils::head(unique(vals[bad]), 5L), collapse = ", "),
            call. = FALSE)
    fin <- num[is.finite(num)]
    if (integer && all(fin == round(fin)) &&
            all(abs(fin) <= .Machine[["integer.max"]]))
        num <- as.integer(num)
    dim(num) <- dim(m)
    dimnames(num) <- dimnames(m)
    num
}
