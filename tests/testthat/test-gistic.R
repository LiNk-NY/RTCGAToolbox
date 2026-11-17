## GISTIC2 helpers: getGISTICPeaks, makeSummarizedExperimentFromGISTIC,
## biocExtract(type = "GISTIC*")

suppressPackageStartupMessages(library(SummarizedExperiment))
data(accmini, package = "RTCGAToolbox", envir = environment())

## a GISTIC object built from the fixture files (same parsing as
## getFirehoseData: fread(colClasses = "character") then data.frame())
fixture_gistic <- function() {
    rd <- function(lines) {
        f <- tempfile()
        writeLines(lines, f)
        data.frame(data.table::fread(f, header = TRUE,
            colClasses = "character", data.table = FALSE))
    }
    new("FirehoseGISTIC", Dataset = "ACC",
        AllByGene = rd(fx_gistic_all),
        ThresholdedByGene = rd(fx_gistic_thr),
        Peaks = rd(fx_gistic_peaks))
}

test_that("getGISTICPeaks keeps only 'Peak N' rows and parses limits", {
    g <- fixture_gistic()
    p <- getGISTICPeaks(g, "wide", rm.chrX = FALSE)
    ## the '- CN values' row is dropped
    expect_identical(p[["Unique.Name"]], c("Amplification Peak  1",
        "Deletion Peak  1", "Amplification Peak  2"))
    ## sorted by genomic position (chr1, chr3, chrX)
    expect_identical(p[["rowRanges"]], c("chr1:120523956-155582840",
        "chr3:59000000-61000000", "chrX:1-2496016"))
    expect_identical(names(p)[1L], "rowRanges")
    ## barcodes standardised
    expect_true(all(GBC %in% names(p)))
    ## BUG-10: the empty trailing column (values "") is dropped
    expect_false(any(grepl("^V[0-9]+$", names(p))))
    expect_identical(names(p), c("rowRanges", "Unique.Name", "Descriptor",
        "Wide.Peak.Limits", "Peak.Limits", "Region.Limits", "q.values",
        "Residual.q.values.after.removing.segments.shared.with.higher.peaks",
        "Broad.or.Focal", "Amplitude.Threshold", GBC))
    ## sample values untouched (still character)
    expect_identical(p[[GBC[1L]]], c("2", "0", "1"))
    expect_identical(p[[GBC[2L]]], c("0", "2", "2"))
})

test_that("getGISTICPeaks peak types select the right limits column", {
    g <- fixture_gistic()
    expect_identical(getGISTICPeaks(g, "narrow")[["rowRanges"]],
        c("chr1:155024310-155083455", "chr3:60000000-60500000"))
    expect_identical(getGISTICPeaks(g, "full")[["rowRanges"]],
        c("chr1:155024310-155095653", "chr3:58000000-62000000"))
    ## rm.chrX = TRUE (default) drops the chrX peak
    expect_false(any(grepl("chrX", getGISTICPeaks(g)[["rowRanges"]])))
    expect_error(getGISTICPeaks(g, "bogus"))
    expect_error(getGISTICPeaks(data.frame()))
})

test_that("getGISTICPeaks works on FirehoseData (accmini)", {
    p <- getGISTICPeaks(accmini, "wide")
    ## rowRanges + 9 annotation columns + 90 samples; V100 (empty) dropped
    expect_identical(dim(p), c(6L, 100L))
    expect_false("V100" %in% names(p))
    expect_identical(p[["rowRanges"]][1:2],
        c("chr1:120523956-155582840", "chr4:1-2496016"))
})

test_that("makeSummarizedExperimentFromGISTIC by-gene values and names", {
    g <- fixture_gistic()
    for (dt in c("AllByGene", "ThresholdedByGene")) {
        se <- makeSummarizedExperimentFromGISTIC(g, dt)
        raw <- slot(g, dt)
        expect_s4_class(se, "SummarizedExperiment")
        expect_identical(rownames(se), c("ACAP3", "AGRN", "TP53"))
        expect_identical(colnames(se), GBC)
        expect_identical(names(rowData(se)),
            c("Gene.Symbol", "Locus.ID", "Cytoband"))
        expect_equal(unname(assay(se)),
            unname(apply(as.matrix(raw[, 4:5]), 2L, as.numeric)), info = dt)
    }
    se <- makeSummarizedExperimentFromGISTIC(g, "AllByGene")
    expect_identical(unname(assay(se)["AGRN", ]), c(0.125, -1.25))
})

test_that("GISTIC assays are numeric and match the fixture text", {
    ## BUG-02: these used to be character matrices
    g <- fixture_gistic()
    a <- assay(makeSummarizedExperimentFromGISTIC(g, "AllByGene"))
    expect_identical(storage.mode(a), "double")
    expect_identical(unname(a), matrix(c(0.03, 0.125, -0.5, -0.07, -1.25, 0.9),
        3L))
    expect_identical(dimnames(a), list(c("ACAP3", "AGRN", "TP53"), GBC))
    t <- assay(makeSummarizedExperimentFromGISTIC(g, "ThresholdedByGene"))
    expect_identical(unname(t), matrix(c(0L, 1L, -2L, 0L, -1L, 2L), 3L))
    p <- assay(makeSummarizedExperimentFromGISTIC(g, "Peaks"))
    expect_identical(unname(p), matrix(c(2L, 0L, 0L, 2L), 2L))
    expect_identical(colnames(p), GBC)
    for (ty in c("GISTICA", "GISTICT", "GISTICP"))
        expect_true(is.numeric(assay(quiet(biocExtract(accmini, ty)))),
            info = ty)
    ## values agree with the character slot
    ta <- assay(makeSummarizedExperimentFromGISTIC(accmini, "ThresholdedByGene"))
    raw <- accmini@GISTIC@ThresholdedByGene
    raw <- as.matrix(raw[grepl("^TCGA", names(raw))])
    expect_identical(unname(ta), unname(array(as.integer(raw), dim(raw))))
})

test_that(".toNumericMatrix is strict and lossless", {
    m <- matrix(c("1", " 2", "", "NA"), 2, dimnames = list(c("a", "b"), NULL))
    expect_identical(ns$.toNumericMatrix(m, integer = TRUE),
        matrix(c(1L, 2L, NA, NA), 2, dimnames = list(c("a", "b"), NULL)))
    expect_identical(ns$.toNumericMatrix(m), matrix(c(1, 2, NA, NA), 2,
        dimnames = list(c("a", "b"), NULL)))
    ## Firehose arrays write missing values as 'null'
    expect_identical(ns$.toNumericMatrix(matrix(c("null", "-0.5"))),
        matrix(c(NA, -0.5)))
    ## non-whole values are never truncated
    expect_identical(ns$.toNumericMatrix(matrix(c("1", "1.5")), integer = TRUE),
        matrix(c(1, 1.5)))
    expect_error(
        ns$.toNumericMatrix(matrix(c("1", "amp")), context = "GISTIC"),
        "Non-numeric values in GISTIC data: amp")
    expect_identical(ns$.toNumericMatrix(matrix(1:2)), matrix(1:2))
    expect_identical(ns$.toNumericMatrix(matrix(character(0L), 0L, 2L)),
        matrix(numeric(0L), 0L, 2L))
})

test_that("makeSummarizedExperimentFromGISTIC duplicated symbols keep rows", {
    g <- fixture_gistic()
    g@AllByGene[["Gene.Symbol"]][2L] <- "ACAP3"
    se <- makeSummarizedExperimentFromGISTIC(g, "AllByGene")
    expect_identical(nrow(se), 3L)
    expect_null(rownames(se))
})

test_that("makeSummarizedExperimentFromGISTIC Peaks: ranges and rowData", {
    g <- fixture_gistic()
    se <- makeSummarizedExperimentFromGISTIC(g, "Peaks")
    expect_s4_class(se, "RangedSummarizedExperiment")
    expect_identical(unname(as.character(rowRanges(se))),
        c("chr1:120523956-155582840", "chr3:59000000-61000000"))
    expect_identical(rowData(se)[["type"]], c("Amplification", "Deletion"))
    expect_identical(colnames(se), GBC)
    expect_identical(unname(assay(se)[, GBC[1L]]), c(2L, 0L))
    expect_identical(unname(assay(se)[, GBC[2L]]), c(0L, 2L))
    ## BUG-10: row names are the peak names, not data.frame row indices
    expect_identical(rownames(se),
        c("Amplification Peak  1", "Deletion Peak  1"))
    expect_identical(rowData(se)[["Unique.Name"]], rownames(se))
})

test_that("Peaks with identical limits are all kept", {
    g <- fixture_gistic()
    p <- g@Peaks
    del <- p[["Unique.Name"]] == "Deletion Peak  1"
    amp <- p[["Unique.Name"]] == "Amplification Peak  1"
    p[del, "Wide.Peak.Limits"] <- p[amp, "Wide.Peak.Limits"]
    g@Peaks <- p
    se <- makeSummarizedExperimentFromGISTIC(g, "Peaks")
    expect_identical(nrow(se), 2L)
    expect_identical(rowData(se)[["type"]], c("Amplification", "Deletion"))
    expect_identical(unname(assay(se)["Deletion Peak  1", ]), c(0L, 2L))
})

test_that("accmini Peaks rowData has no empty trailing column", {
    ## BUG-10
    se <- makeSummarizedExperimentFromGISTIC(accmini, "Peaks")
    expect_false("V100" %in% names(rowData(se)))
    ## a named column that happens to be blank is kept
    expect_true("Broad.or.Focal" %in% names(rowData(se)))
    expect_false(anyDuplicated(rownames(se)) > 0L)
    expect_true(all(grepl("Peak +[0-9]+$", rownames(se))))
})

test_that("makeSummarizedExperimentFromGISTIC input validation", {
    expect_error(makeSummarizedExperimentFromGISTIC(data.frame()),
        "must be a 'FirehoseData' or 'FirehoseGISTIC'")
    expect_error(makeSummarizedExperimentFromGISTIC(accmini, "bogus"))
    expect_identical(
        makeSummarizedExperimentFromGISTIC(new("FirehoseGISTIC"), "AllByGene"),
        list())
})

test_that("biocExtract GISTIC types", {
    all3 <- quiet(biocExtract(accmini, "GISTIC"))
    expect_type(all3, "list")
    expect_identical(names(all3), c("AllByGene", "ThresholdedByGene", "Peaks"))
    expect_identical(quiet(biocExtract(accmini, "GISTICA")),
        makeSummarizedExperimentFromGISTIC(accmini, "AllByGene"))
    expect_identical(quiet(biocExtract(accmini, "GISTICT")),
        makeSummarizedExperimentFromGISTIC(accmini, "ThresholdedByGene"))
    expect_identical(quiet(biocExtract(accmini, "GISTICP")),
        makeSummarizedExperimentFromGISTIC(accmini, "Peaks"))
    expect_identical(all3[["AllByGene"]],
        makeSummarizedExperimentFromGISTIC(accmini, "AllByGene"))
})
