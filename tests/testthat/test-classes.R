## S4 classes, accessors, show and updateObject methods

data(accmini, package = "RTCGAToolbox", envir = environment())

test_that("FirehoseData has the expected slots and accmini is valid", {
    expect_true(validObject(accmini))
    expect_identical(slotNames("FirehoseData"), c("Dataset", "runDate",
        "gistic2Date", "clinical", "RNASeqGene", "RNASeq2Gene",
        "RNASeq2GeneNorm", "miRNASeqGene", "CNASNP", "CNVSNP", "CNASeq",
        "CNACGH", "Methylation", "mRNAArray", "miRNAArray", "RPPAArray",
        "Mutation", "GISTIC", "BarcodeUUID"))
    expect_identical(slotNames("FirehoseGISTIC"),
        c("Dataset", "AllByGene", "ThresholdedByGene", "Peaks"))
})

test_that("getData returns raw slots", {
    expect_identical(getData(accmini, "clinical"), accmini@clinical)
    expect_identical(getData(accmini, "RNASeq2Gene"), accmini@RNASeq2Gene)
    expect_identical(getData(accmini, "CNASNP"), accmini@CNASNP)
    ## RNASeq2GeneNorm is not in the platform list: the list is returned
    expect_identical(getData(accmini, "RNASeq2GeneNorm"),
        accmini@RNASeq2GeneNorm)
    expect_error(getData(accmini, "foo"), "no slot")
    expect_error(getData(accmini, c("a", "b")))
    expect_error(getData(accmini, NA_character_))
})

test_that("getData with platform index for list slots", {
    expect_identical(getData(accmini, "Methylation", 1),
        accmini@Methylation[[1L]]@DataMatrix)
    expect_identical(getData(accmini, "RPPAArray", 1L),
        accmini@RPPAArray[[1L]]@DataMatrix)
    expect_message(
        expect_error(getData(accmini, "Methylation", 2),
            "valid 'platform' index"),
        "Accessible platforms")
    expect_error(quiet(getData(accmini, "Methylation", "a")),
        "valid 'platform' index")
    expect_error(getData(accmini, "mRNAArray", 1), "No data available")
    ## BUG-12: invalid indices get the "valid platform" message
    for (bad in list(0, -1, 1.5))
        expect_error(quiet(getData(accmini, "Methylation", bad)),
            "valid 'platform' index", info = bad)
})

test_that("getData for GISTIC platforms", {
    for (p in c("AllByGene", "ThresholdedByGene", "Peaks"))
        expect_identical(getData(accmini, "GISTIC", p), slot(accmini@GISTIC, p))
    expect_error(getData(accmini, "GISTIC", "foo"), "GISTIC platforms")
    expect_identical(getData(accmini@GISTIC, "GISTIC", "Peaks"),
        accmini@GISTIC@Peaks)
})

test_that("selectType returns the slot", {
    expect_identical(selectType(accmini, "CNVSNP"), accmini@CNVSNP)
    expect_identical(selectType(accmini, "GISTIC"), accmini@GISTIC)
})

test_that("isEmpty for FirehoseGISTIC", {
    expect_false(isEmpty(accmini@GISTIC))
    expect_true(isEmpty(new("FirehoseGISTIC")))
})

test_that("show methods print a summary", {
    out <- capture.output(show(accmini))
    expect_match(out[1L], "^ACC FirehoseData object")
    expect_true(any(grepl("RNASeq2Gene: A matrix .* dim:  6 x 79", out)))
    expect_true(any(grepl("CNASNP: A data.frame, dim:  6 x 6", out)))
    expect_true(any(grepl("Mutation: A data.frame, dim:  6 x 52", out)))
    expect_true(any(grepl("GISTIC: A FirehoseGISTIC", out)))
    expect_false(any(grepl("^  RNASeqGene:", out)))

    ## BUG-16: a Mutation table with columns but no rows is not listed
    x <- accmini
    x@Mutation <- x@Mutation[0L, ]
    expect_false(any(grepl("Mutation:", capture.output(show(x)))))

    empty <- capture.output(show(new("FirehoseData", Dataset = "X")))
    expect_identical(empty[length(empty)],
        "To export data, use the 'getData' function.")

    msgs <- function(x) paste(capture_messages(show(x)), collapse = "")
    expect_match(msgs(accmini@GISTIC), "Dataset:ACC")
    expect_match(msgs(accmini@GISTIC), "FirehoseGISTIC object, dim: 6\t93")
    expect_match(msgs(accmini@Methylation[[1L]]),
        "FirehoseMethylationArray object, dim: 6\t83")
    expect_match(msgs(accmini@RPPAArray[[1L]]), "dim: 6\t46")
    expect_match(msgs(new("FirehoseCGHArray", Filename = "f",
        DataMatrix = data.frame(a = 1))), "Platform:f")
})

test_that("DGEResult / CorResult showResults return the tables", {
    tt <- data.frame(gene = letters[1:8], p = seq(0.01, 0.08, by = 0.01))
    d <- new("DGEResult", Dataset = "ACC", Toptable = tt)
    expect_message(out <- capture.output(res <- showResults(d)), "Dataset: ACC")
    expect_match(out[1L], "gene")
    expect_identical(res, tt)
    cr <- new("CorResult", Dataset = "ACC", Correlations = tt)
    expect_message(out <- capture.output(res <- showResults(cr)), "Dataset: ACC")
    expect_identical(res, tt)
    expect_identical(capture_messages(show(d)),
        c("Dataset:ACC\n", "DGEResult object, dim: 8\t2\n"))
    expect_identical(capture_messages(show(cr)),
        c("Dataset:ACC\n", "CorResult object, dim: 8\t2\n"))
})

test_that("updateObject is a no-op on current objects", {
    expect_identical(updateObject(accmini), accmini)
    expect_identical(updateObject(accmini@GISTIC), accmini@GISTIC)
    expect_match(paste(capture_messages(updateObject(accmini, verbose = TRUE)),
        collapse = ""), "updateObject\\(object = 'FirehoseData'\\)")
    expect_false(ns$.hasOldAPI(accmini))
    expect_false(ns$.hasOldGISTIC(accmini@GISTIC))
})
