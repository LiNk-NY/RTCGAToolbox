## Conversion of tabular data to GRanges / RaggedExperiment / RSE objects.

suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(RaggedExperiment)
})

seg_df <- function() data.frame(
    Sample = c("TCGA-A1-0001-01A", "TCGA-A1-0001-01A", "TCGA-A1-0002-01A"),
    Chromosome = c(1, 2, 1),
    Start = c(100, 500, 150),
    End = c(200, 900, 300),
    Num_Probes = c(10, 20, 30),
    Segment_Mean = c(0.5, -1.2, 0.1)
)

maf_df <- function() data.frame(
    Hugo_Symbol = c("TP53", "KRAS", "TP53"),
    NCBI_Build = "37",
    Chromosome = c("17", "12", "17"),
    Start_position = c(7577120, 25398284, 7577121),
    End_position = c(7577120, 25398284, 7577121),
    Strand = c("+", "+", "-"),
    Variant_Classification = c("Missense", "Missense", "Nonsense"),
    Tumor_Sample_Barcode = c("TCGA-A1-0001-01A", "TCGA-A1-0001-01A",
        "TCGA-A1-0002-01A")
)

test_that(".makeGRangesListFromDataFrame splits rows by sample", {
    grl <- ns$.makeGRangesListFromDataFrame(seg_df(), split.field = "Sample",
        seqnames.field = "Chromosome", start.field = "Start",
        end.field = "End", ignore.strand = TRUE)
    expect_s4_class(grl, "GRangesList")
    expect_identical(names(grl), c("TCGA-A1-0001-01A", "TCGA-A1-0002-01A"))
    expect_identical(lengths(grl), c("TCGA-A1-0001-01A" = 2L,
        "TCGA-A1-0002-01A" = 1L))
    g1 <- grl[["TCGA-A1-0001-01A"]]
    expect_identical(as.character(seqnames(g1)), c("1", "2"))
    expect_identical(start(g1), c(100L, 500L))
    expect_identical(end(g1), c(200L, 900L))
    expect_identical(mcols(g1)[["Segment_Mean"]], c(0.5, -1.2))

    expect_error(ns$.makeGRangesListFromDataFrame(seg_df(), split.field = "x"),
        "not in 'names\\(df\\)'")
    expect_error(ns$.makeGRangesListFromDataFrame(seg_df(),
        split.field = c("a", "b")), "single string")
})

test_that(".makeGRangesListFromDataFrame drops rows with missing start/end", {
    df <- seg_df()
    df[["Start"]][2L] <- NA
    grl <- ns$.makeGRangesListFromDataFrame(df, split.field = "Sample",
        seqnames.field = "Chromosome", start.field = "Start",
        end.field = "End", ignore.strand = TRUE)
    expect_identical(sum(lengths(grl)), 2L)
    expect_identical(mcols(grl[[1L]])[["Segment_Mean"]], 0.5)
})

test_that(".makeRaggedExperimentFromDataFrame keeps every value in place", {
    re <- quiet(ns$.makeRaggedExperimentFromDataFrame(seg_df(),
        split.field = "Sample", build = "hg19"))
    expect_s4_class(re, "RaggedExperiment")
    expect_identical(colnames(re), c("TCGA-A1-0001-01A", "TCGA-A1-0002-01A"))
    expect_identical(assayNames(re), c("Num_Probes", "Segment_Mean"))
    rr <- rowRanges(re)
    expect_identical(as.character(seqnames(rr)), c("1", "2", "1"))
    expect_identical(start(rr), c(100L, 500L, 150L))
    expect_identical(end(rr), c(200L, 900L, 300L))
    sm <- assay(re, "Segment_Mean")
    expect_identical(unname(sm[, 1L]), c(0.5, -1.2, NA))
    expect_identical(unname(sm[, 2L]), c(NA, NA, 0.1))
    expect_identical(unname(genome(re)), c("hg19", "hg19"))
})

test_that(".makeRaggedExperimentFromDataFrame standardises strand", {
    df <- maf_df()
    df[["Strand"]] <- c("1", "-1", "null")
    re <- quiet(ns$.makeRaggedExperimentFromDataFrame(df,
        split.field = "Tumor_Sample_Barcode"))
    expect_identical(as.character(strand(rowRanges(re))), c("+", "-", "*"))
})

consistent_df <- function() data.frame(
    Sample = rep(c("TCGA-A1-0001-01A", "TCGA-A1-0002-01A"), each = 2),
    Chromosome = 1, Start = c(1, 100, 1, 100), End = c(50, 200, 50, 200),
    Num_Probes = c(5, 6, 7, 8),
    Segment_Mean = c(0.1, 0.2, 0.3, 0.4))

test_that(".makeRangedSummarizedExperimentFromDataFrame for identical ranges", {
    df <- consistent_df()
    expect_true(ns$.hasConsistentRanges(df))
    rse <- ns$.makeRangedSummarizedExperimentFromDataFrame(df, build = "hg19")
    expect_s4_class(rse, "RangedSummarizedExperiment")
    expect_identical(dim(rse), c(2L, 2L))
    expect_identical(assayNames(rse), c("Num_Probes", "Segment_Mean"))
    expect_identical(colnames(rse), c("TCGA-A1-0001-01A", "TCGA-A1-0002-01A"))
    expect_identical(unname(assay(rse, "Num_Probes")),
        matrix(c(5, 6, 7, 8), 2))
    m <- assay(rse, "Segment_Mean")
    expect_identical(unname(m[, 1L]), c(0.1, 0.2))
    expect_identical(unname(m[, 2L]), c(0.3, 0.4))
    expect_identical(start(rse), c(1L, 100L))
    expect_identical(unname(genome(rse)), "hg19")
})

test_that("[BUG-24] identical ranges with a single value column error", {
    df <- consistent_df()[, -5L]
    expect_true(ns$.hasConsistentRanges(df))
    expect_error(ns$.makeRangedSummarizedExperimentFromDataFrame(df),
        "invalid 'length' argument")
    expect_error(quiet(biocExtract(df, "CNASNP")),
        "invalid 'length' argument")
})

test_that(".hasConsistentRanges is FALSE when samples have different ranges", {
    expect_false(ns$.hasConsistentRanges(seg_df()))
    expect_false(ns$.hasConsistentRanges(data.frame(a = 1)))
})

test_that(".makeGRangesFromDataFrame without sample column", {
    df <- seg_df()[, -1L]
    df[["Start"]][3L] <- NA
    gr <- ns$.makeGRangesFromDataFrame(df, build = "hg19")
    expect_s4_class(gr, "GRanges")
    expect_identical(length(gr), 2L)
    expect_identical(mcols(gr)[["Segment_Mean"]], c(0.5, -1.2))
    expect_identical(unname(genome(gr)), c("hg19", "hg19"))
})

test_that(".unNestList converts platform objects to annotated DataFrames", {
    data(accmini, package = "RTCGAToolbox", envir = environment())
    dm <- accmini@RNASeq2GeneNorm[[1L]]@DataMatrix
    res <- ns$.unNestList(accmini@RNASeq2GeneNorm)
    expect_s4_class(res, "DataFrame")
    expect_identical(metadata(res), list(
        filename = accmini@RNASeq2GeneNorm[[1L]]@Filename,
        build = NA_character_, platform = "illuminahiseq"))
    expect_identical(rownames(res), rownames(dm))
    ## DataFrame() mangles the barcodes; they are restored downstream
    expect_identical(names(res), make.names(colnames(dm)))
    expect_identical(unname(as.matrix(res)), unname(dm))
    ## lists of anything else pass through
    expect_identical(ns$.unNestList(list(1, 2)), list(1, 2))
})

mk_arr <- function(fname, m) new("FirehosemRNAArray", Filename = fname,
    DataMatrix = m)
agil1 <- "gdac_ACC.Merge_transcriptome__agilentg4502a_07_1__unc"
agil2 <- "gdac_ACC.Merge_transcriptome__agilentg4502a_07_2__unc"
u133 <- "gdac_ACC.Merge_transcriptome__ht_hg_u133a__broad"
m12 <- matrix(c(1, 2, 3, 4), 2,
    dimnames = list(c("g1", "g2"), c("TCGA-A1-0001-01A", "TCGA-A1-0002-01A")))

test_that(".unNestList merges same-named platforms when rows align", {
    m34 <- matrix(c(5, 6, 7, 8), 2,
        dimnames = list(c("g1", "g2"), c("TCGA-A1-0003-01A", "TCGA-A1-0004-01A")))
    res <- ns$.unNestList(list(mk_arr(agil1, m12), mk_arr(agil2, m34)))
    expect_identical(names(res), "TX_g4502a")
    expect_identical(dim(res[[1L]]), c(2L, 4L))
    expect_identical(unname(as.matrix(res[[1L]])),
        unname(cbind(m12, m34)))
    expect_identical(metadata(res[[1L]])[["filename"]], c(agil1, agil2))
})

test_that("[BUG-03] same-named platforms are cbind-ed by position, not name", {
    ## the second platform lists g2 before g1; values end up on the wrong gene
    m_rev <- matrix(c(60, 50, 80, 70), 2,
        dimnames = list(c("g2", "g1"), c("TCGA-A1-0003-01A", "TCGA-A1-0004-01A")))
    obj <- new("FirehoseData", Dataset = "ACC",
        mRNAArray = list(mk_arr(agil1, m12), mk_arr(agil2, m_rev)))
    res <- quiet(biocExtract(obj, "mRNAArray"))
    a <- assay(res[[1L]])
    ## true value of g1 in sample 3 is 50
    expect_identical(a["g1", "TCGA-A1-0003-01A"], 60)
})

test_that("[BUG-04] same-named platforms with different row counts vanish", {
    m3 <- matrix(9:14 + 0, 3, dimnames = list(c("g1", "g2", "g3"),
        c("TCGA-A1-0005-01A", "TCGA-A1-0006-01A")))
    obj <- new("FirehoseData", Dataset = "ACC",
        mRNAArray = list(mk_arr(agil1, m12), mk_arr(agil2, m3)))
    res <- quiet(biocExtract(obj, "mRNAArray"))
    expect_length(res, 0L)
})

test_that("distinct platforms are returned as a named list of SEs", {
    m34 <- matrix(c(5, 6, 7, 8), 2,
        dimnames = list(c("g1", "g2"), c("TCGA-A1-0003-01A", "TCGA-A1-0004-01A")))
    obj <- new("FirehoseData", Dataset = "ACC",
        mRNAArray = list(mk_arr(agil1, m12), mk_arr(u133, m34)))
    res <- quiet(biocExtract(obj, "mRNAArray"))
    expect_type(res, "list")
    expect_identical(names(res), c("TX_g4502a", "TX_ht_hg_u133a"))
    expect_identical(unname(assay(res[["TX_g4502a"]])), unname(m12))
    expect_identical(unname(assay(res[["TX_ht_hg_u133a"]])), unname(m34))
    expect_identical(colnames(res[["TX_ht_hg_u133a"]]), colnames(m34))
    expect_identical(metadata(res[["TX_ht_hg_u133a"]])[["filename"]], u133)
})

test_that("[BUG-14] .makeSummarizedExperimentFromDataFrame is broken (unused)", {
    df <- data.frame(Hugo_Symbol = c("A", "B"),
        "TCGA-A1-0001-01A" = c(1, 2), "TCGA-A1-0002-01A" = c(3, 4),
        check.names = FALSE)
    expect_error(ns$.makeSummarizedExperimentFromDataFrame(df))
    expect_error(ns$.makeSummarizedExperimentFromDataFrame(
        S4Vectors::DataFrame(df, check.names = FALSE)))
})
