#' @include utils.R
NULL

#' Extract and convert data from a `FirehoseData` object to a
#' `Bioconductor` object
#'
#' This function processes data from a
#' [FirehoseData][FirehoseData-class] object. Raw data is
#' converted to a conventional Bioconductor object. Copy number segments
#' (`CNASNP`, `CNVSNP`, `CNASeq`, `CNACGH`) and mutations (`Mutation`) are
#' always returned as a
#' [RaggedExperiment][RaggedExperiment::RaggedExperiment-class], whatever the
#' number of samples; use [RaggedExperiment::compactAssay()] to get a matrix
#' when all samples share the same ranges. The other data types are returned
#' as a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class].
#' The raw tables are available with [getData()]. In cases
#' where there are multiple platforms in a data type, an attempt to consolidate
#' datasets will be made based on matching dimension names. For ranged data,
#' this functionality is provided with more control as part of the
#' `RaggedExperiment` features. See the
#' [RaggedExperiment-class][RaggedExperiment::RaggedExperiment-class] for
#' more details.
#'
#' @details A typical additional argument for this function passed down to
#' lower level functions is the `names.field` which indicates the row names
#' in the data. By default, it is the "Hugo_Symbol" column in the internal
#' code that converts `data.frame`s to `SummarizedExperiment` representations
#' (via the `.makeSummarizedExperimentFromDataFrame` internal function).
#'
#' @section type:
#' Choices include the following:
#' * clinical: Get the clinical data slot
#' * RNASeqGene: RNASeqGene, RNASeq v1
#' * RNASeqGene: RNASeq2Gene, RNASeq v2
#' * RNASeq2GeneNorm: RNASeq v2 Normalized
#' * miRNASeqGene: micro RNA SeqGene
#' * CNASNP: Copy Number Alteration
#' * CNVSNP: Copy Number Alteration (germline CNVs removed)
#' * CNASeq: Copy Number Alteration
#' * CNACGH: Copy Number Alteration
#' * Methylation: Methylation
#' * mRNAArray: Messenger RNA
#' * miRNAArray: micro RNA
#' * RPPAArray: Reverse Phase Protein Array
#' * Mutation: Mutations
#' * GISTICA: GISTIC v2 ('AllByGene' only)
#' * GISTICT: GISTIC v2 ('ThresholdedByGene' only)
#' * GISTICP: GISTIC v2 ('Peaks' only)
#' * GISTIC: GISTIC v2 scores, probabilities, and peaks
#'
#' @param object A `FirehoseData` object from which to extract data.
#' @param type The type of data to extract from the "FirehoseData" object,
#' see type section.
#' @param ... Additional arguments passed to lower level functions that
#' convert tabular data into Bioconductor object such as
#' `.makeRangedSummarizedExperimentFromDataFrame` or
#' `.makeRaggedExperimentFromDataFrame`
#'
#' @return A [RaggedExperiment][RaggedExperiment::RaggedExperiment-class]
#' for copy number and mutation data (a `GRanges` if the table has no sample
#' column), a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class]
#' for other data types, a list of these for data types with several
#' platforms, or the `data.frame` for `clinical`.
#'
#' @author Marcel Ramos \email{marcel.ramos@@sph.cuny.edu}
#'
#' @examples
#'
#' data(accmini)
#' biocExtract(accmini, "RNASeq2Gene")
#' biocExtract(accmini, "miRNASeqGene")
#' biocExtract(accmini, "RNASeq2GeneNorm")
#' biocExtract(accmini, "CNASNP")
#' biocExtract(accmini, "CNVSNP")
#' biocExtract(accmini, "Methylation")
#' biocExtract(accmini, "Mutation")
#' biocExtract(accmini, "RPPAArray")
#' biocExtract(accmini, "GISTIC")
#'
#' @export biocExtract
biocExtract <- function(object, type = c("clinical", "RNASeqGene",
    "RNASeq2Gene", "miRNASeqGene", "RNASeq2GeneNorm", "CNASNP", "CNVSNP",
    "CNASeq", "CNACGH", "Methylation", "Mutation", "mRNAArray", "miRNAArray",
    "RPPAArray", "GISTIC", "GISTICA", "GISTICT", "GISTICP"), ...) {
    if (!identical(length(type), 1L))
        stop("Please specify a single data type")
    message("working on: ", type)
    fromShell <- is(object, "FirehoseData")
    if (fromShell)
        object <- .removeShell(object, type)
    if (!length(object)) { return(object) }
    if (is.list(object) && !is.data.frame(object)) {
        object <- .unNestList(object)
    }
    if (type == "clinical") { return(object) }
    if (is.matrix(object)) {
        return(SummarizedExperiment(assays = SimpleList(object)))
    }
    if (is(object, "list") && !is(object, "DataFrame") &&
        type != "Methylation") {
        return(.extractList(object, type = type))
    }
    if (is(object, "SummarizedExperiment")) { return(object) }

    gisticType <- grepl("^GISTIC", type, ignore.case = TRUE)
    if (gisticType) {
        slotreq <- switch(type,
            GISTICA = "AllByGene", GISTICT = "ThresholdedByGene",
            GISTICP = "Peaks",
            GISTIC = c("AllByGene", "ThresholdedByGene", "Peaks")
        )
        result <-
        if (type == "GISTIC") {
            names(slotreq) <- slotreq
            Filter(length,
                lapply(slotreq, function (x) {
                    makeSummarizedExperimentFromGISTIC(
                        gistic = object, dataType = x
                    )
                })
            )
        } else {
            makeSummarizedExperimentFromGISTIC(object, slotreq)
        }
        if (length(result) == 1L)
            result <- result[[1L]]
        return(result)
    }

    if (type == "Methylation") {
        if (is.list(object)) {
            result <- lapply(object, .getMethyl)
        } else {
            result <- .getMethyl(object)
        }
        return(result)
    }

    hasRanged <- .hasRangeNames(object)
    if (hasRanged) {
        ## getFirehoseData() only downloads the '*_hg19__seg' SNP6 archives
        hg19seg <- fromShell && type %in% c("CNASNP", "CNVSNP")
        ## segments and mutations are always RaggedExperiment (BUG-11)
        ragged <- type %in% c("CNASNP", "CNVSNP", "CNASeq", "CNACGH",
            "Mutation")
        object <- .convertRangeBioc(object,
            defaultBuild = if (hg19seg) "hg19" else NA_character_,
            ragged = ragged)
    } else {
        object <- .standardizeBC(object)
        metadat <- if (is(object, "Annotated")) metadata(object) else list()
        ## data.matrix() would turn character columns into factor codes
        isNum <- vapply(seq_len(ncol(object)), function(j)
            is.numeric(object[[j]]) || is.logical(object[[j]]), logical(1L))
        if (!all(isNum))
            stop("Cannot convert non-numeric columns to an assay: ",
                paste(utils::head(names(object)[!isNum], 5L), collapse = ", "))
        object <- data.matrix(object)
        object <- SummarizedExperiment(assays = SimpleList(object))
        metadata(object) <- metadat
    }
    return(object)
}

## 'defaultBuild' is used when the data have no NCBI_Build column.
## With 'ragged = TRUE', data with a sample column always become a
## RaggedExperiment; otherwise identical ranges across samples give a
## RangedSummarizedExperiment (cBioPortalData relies on this default).
.convertRangeBioc <- function(object, defaultBuild = NA_character_,
    ragged = FALSE, ...) {
    build <- unique(.getBuild(object))
    if (all(is.na(build)))
        build <- defaultBuild
    build <- TCGAutils::correctBuild(build, "NCBI")
    if (!ragged && .hasConsistentRanges(object)) {
        object <- .makeRangedSummarizedExperimentFromDataFrame(object,
            build = build, ...)
    } else {
        split.field <- .findSampleCol(object)
        if (is.na(split.field) || !length(split.field))
            object <- .makeGRangesFromDataFrame(object,
                build = build)
        else
            object <- .makeRaggedExperimentFromDataFrame(
                object, split.field = split.field, build = build, ...)
    }
    object
}
