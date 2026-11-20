#' Get resource links from inputs
#'
#' This function provides a reference to the resources downloaded from the
#' GDAC Firehose pipeline. Based on the input, the function returns a
#' URL location to the resource if there exists one.
#'
#' @inheritParams getFirehoseData
#'
#' @param data_date Either a runDate or analysisDate typically entered in
#' `getFirehoseData`
#'
#' @return A character URL to a dataset location
#' @examples
#'
#' getLinks("BRCA", CNASeq = TRUE)
#'
#' @export
getLinks <- function(dataset, data_date="20160128",
    RNASeqGene=FALSE, RNASeq2Gene=FALSE, clinical=FALSE, miRNASeqGene=FALSE,
    RNASeq2GeneNorm=FALSE,
    RNAseq2Norm = c(
      "normalized_counts", "RSEM_normalized_log2", "raw_counts", "scaled_estimate"
    ),
    CNASNP=FALSE, CNVSNP=FALSE, CNASeq=FALSE,
    CNACGH=FALSE, Methylation=FALSE, Mutation=FALSE, mRNAArray=FALSE,
    miRNAArray=FALSE, RPPAArray=FALSE, GISTIC=FALSE)
{
    fh_url <- "https://gdac.broadinstitute.org/runs/stddata__"

    if (GISTIC)
        fh_url <- "https://gdac.broadinstitute.org/runs/analyses__"

    fh_url <- paste0(fh_url, substr(data_date, 1, 4), "_",
        substr(data_date, 5, 6), "_", substr(data_date, 7, 8), "/data/")
    fh_url <- paste0(fh_url, dataset, "/", data_date, "/")

    doc <- .files_from_html_table(fh_url)

    ## the first selected data type, in this order, is used
    ## (each flag is tested with if (), as before)
    types <- vapply(list(RNASeqGene = RNASeqGene, RNASeq2Gene = RNASeq2Gene,
        clinical = clinical, miRNASeqGene = miRNASeqGene,
        RNASeq2GeneNorm = RNASeq2GeneNorm, CNASNP = CNASNP, CNVSNP = CNVSNP,
        CNASeq = CNASeq, CNACGH = CNACGH, Methylation = Methylation,
        Mutation = Mutation, mRNAArray = mRNAArray, miRNAArray = miRNAArray,
        RPPAArray = RPPAArray, GISTIC = GISTIC),
        function(x) if (x) TRUE else FALSE, logical(1L))
    type <- names(types)[types][1L]
    if (is.na(type))
        return(fh_url)
    if (identical(type, "RNASeq2GeneNorm") &&
        !identical(match.arg(RNAseq2Norm), "normalized_counts"))
        type <- "mRNAseq_Preprocess"
    if (identical(type, "GISTIC")) {
        tag <- switch(dataset, SKCM = "TM", LAML = "TB", "TP")
        dataset <- paste0(dataset, "-", tag)
    }
    paste0(fh_url, .archiveLinks(type, dataset, doc))
}
