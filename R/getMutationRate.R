#' Make a table for mutation rate of each gene in the cohort
#'
#' For each gene, the proportion of samples with at least one mutation in
#' that gene. The denominator is the number of samples that appear in the
#' mutation data (i.e., samples with at least one mutation), not all samples
#' in the cohort.
#'
#' @param dataObject This must be `FirehoseData` object.
#' @return A `data.frame` with columns `Genes` (from `Hugo_Symbol`) and
#'   `MutationRatio`, one row per gene, with genes as row names. Samples are
#'   identified by the `Tumor_Sample_Barcode` column. With no mutation data,
#'   a zero-row `data.frame` is returned with a warning.
#' @examples
#' data(accmini)
#' mutRate <- getMutationRate(dataObject=accmini)
#' mutRate <- mutRate[order(mutRate[,2],decreasing = TRUE),]
#' head(mutRate)
#' @export getMutationRate
getMutationRate <- function(dataObject)
{
    if (!is(dataObject, "FirehoseData"))
        stop("'dataObject' must be 'FirehoseData' class!")
    mutAll <- dataObject@Mutation
    if (!nrow(mutAll) || !ncol(mutAll)) {
        warning("No mutation data in 'dataObject'")
        return(data.frame(Genes = character(0L), MutationRatio = numeric(0L)))
    }
    geneCol <- .findCol(mutAll, "Hugo_Symbol")
    sampleCol <- .findCol(mutAll, "Tumor_Sample_Barcode")
    if (!length(geneCol) || !length(sampleCol))
        stop("Mutation data must have 'Hugo_Symbol' and",
            " 'Tumor_Sample_Barcode' columns")
    genes <- as.character(mutAll[[geneCol]])
    samples <- as.character(mutAll[[sampleCol]])
    uniqueGenes <- unique(genes)
    uniqueSamples <- unique(samples)
    countsMut <- matrix(0, length(uniqueGenes), length(uniqueSamples),
        dimnames = list(uniqueGenes, uniqueSamples))
    countsMut[cbind(match(genes, uniqueGenes),
        match(samples, uniqueSamples))] <- 1
    mutRatio <- rowMeans(countsMut)
    data.frame(Genes = uniqueGenes, MutationRatio = mutRatio,
        row.names = uniqueGenes)
}
