#' Deprecated and defunct features in RTCGAToolbox
#'
#' These features are provided for compatibility with older versions of
#' `RTCGAToolbox` only, and will be removed in the next release.
#'
#' @section Deprecated arguments:
#' * `getFirehoseData(getUUIDs = TRUE)`: the barcode-to-UUID lookup used
#'   the TCGA Data Coordinating Center web service
#'   (`tcga-data.nci.nih.gov`), which was retired in 2016, so it can no
#'   longer work. The argument is ignored with a warning and the
#'   `BarcodeUUID` slot is left empty. `TCGAutils::barcodeToUUID()`
#'   translates barcodes with the NCI Genomic Data Commons API.
#'
#' @name RTCGAToolbox-deprecated
#' @keywords internal
NULL

## BiocBaseUtils::lifeCycle() names the deprecated feature after the
## function that calls it. Calling it from a function named, e.g.,
## 'getFirehoseData(getUUIDs = TRUE)' makes the standard message refer to
## the argument rather than to the whole function.
.lifeCycleArg <- function(fun, arg, cycle = c("deprecated", "defunct")) {
    cycle <- match.arg(cycle)
    what <- paste0(fun, "(", arg, ")")
    env <- new.env(parent = environment())
    assign(what, function() BiocBaseUtils::lifeCycle(
        package = "RTCGAToolbox", cycle = cycle, title = "RTCGAToolbox"),
        envir = env)
    eval(call(what), env)
}
