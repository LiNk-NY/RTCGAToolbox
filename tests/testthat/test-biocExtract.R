## biocExtract() on the bundled real ACC subset. Every test compares the
## Bioconductor object against the raw slot it was built from, value by value.

suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(RaggedExperiment)
})
data(accmini, package = "RTCGAToolbox", envir = environment())

test_that("biocExtract requires a single type", {
    expect_error(biocExtract(accmini, c("RNASeq2Gene", "CNASNP")),
        "single data type")
    expect_error(biocExtract(accmini))
})

test_that("clinical is returned unchanged", {
    expect_identical(quiet(biocExtract(accmini, "clinical")), accmini@clinical)
})

test_that("empty slots return the empty slot", {
    expect_identical(quiet(biocExtract(accmini, "RNASeqGene")), accmini@RNASeqGene)
    expect_identical(quiet(biocExtract(accmini, "CNASeq")), accmini@CNASeq)
    expect_identical(quiet(biocExtract(accmini, "mRNAArray")), list())
    expect_identical(quiet(biocExtract(accmini, "CNACGH")), list())
})

test_that("matrix slots become SummarizedExperiments with identical values", {
    for (type in c("RNASeq2Gene", "miRNASeqGene")) {
        raw <- slot(accmini, type)
        se <- quiet(biocExtract(accmini, type))
        expect_s4_class(se, "SummarizedExperiment")
        expect_identical(assay(se), raw, info = type)
        expect_identical(dimnames(se), dimnames(raw), info = type)
    }
})

test_that("list-of-platform slots carry values and metadata", {
    for (type in c("RNASeq2GeneNorm", "RPPAArray")) {
        plat <- slot(accmini, type)[[1L]]
        raw <- plat@DataMatrix
        se <- quiet(biocExtract(accmini, type))
        expect_s4_class(se, "SummarizedExperiment")
        a <- assay(se)
        expect_identical(storage.mode(a), "double", info = type)
        expect_identical(unname(a), unname(raw), info = type)
        expect_identical(rownames(se), rownames(raw), info = type)
        expect_identical(colnames(se), colnames(raw), info = type)
        expect_identical(metadata(se)[["filename"]], plat@Filename, info = type)
    }
})

test_that("Methylation: beta values numeric, probe ids rows, annotation rowData", {
    raw <- accmini@Methylation[[1L]]@DataMatrix
    se <- quiet(biocExtract(accmini, "Methylation"))
    expect_s4_class(se, "SummarizedExperiment")
    expect_s4_class(assay(se), "DelayedMatrix")
    a <- as.matrix(assay(se))
    expect_identical(rownames(se), rownames(raw))
    samples <- grep("^TCGA", names(raw), value = TRUE)
    expect_identical(colnames(se), samples)
    expect_identical(unname(a),
        unname(apply(as.matrix(raw[, samples]), 2L, as.numeric)))
    expect_identical(names(rowData(se)),
        c("Gene_Symbol", "Chromosome", "Genomic_Coordinate"))
    expect_identical(rowData(se)[["Gene_Symbol"]], raw[["Gene_Symbol"]])
})

test_that("CNASNP/CNVSNP (single sample) become RangedSummarizedExperiment", {
    for (type in c("CNASNP", "CNVSNP")) {
        raw <- slot(accmini, type)
        rse <- quiet(biocExtract(accmini, type))
        expect_s4_class(rse, "RangedSummarizedExperiment")
        expect_identical(colnames(rse), unique(raw[["Sample"]]), info = type)
        expect_identical(assayNames(rse), c("Num_Probes", "Segment_Mean"))
        expect_identical(unname(assay(rse, "Segment_Mean")[, 1L]),
            raw[["Segment_Mean"]], info = type)
        expect_identical(unname(assay(rse, "Num_Probes")[, 1L]),
            raw[["Num_Probes"]], info = type)
        expect_identical(start(rse), as.integer(raw[["Start"]]), info = type)
        expect_identical(end(rse), as.integer(raw[["End"]]), info = type)
        expect_identical(as.character(seqnames(rse)),
            as.character(raw[["Chromosome"]]), info = type)
        ## [BUG-17] the data are hg19 but no genome is recorded
        expect_true(all(is.na(genome(rse))), info = type)
    }
})

test_that("multi-sample segments become a RaggedExperiment", {
    seg <- data.frame(
        Sample = c("TCGA-A1-0001-01A", "TCGA-A1-0001-01A", "TCGA-A1-0002-01A"),
        Chromosome = c(1, 2, 1), Start = c(100, 500, 150),
        End = c(200, 900, 300), Num_Probes = c(10, 20, 30),
        Segment_Mean = c(0.5, -1.2, 0.1))
    re <- quiet(biocExtract(seg, "CNASNP"))
    expect_s4_class(re, "RaggedExperiment")
    sm <- assay(re, "Segment_Mean")
    expect_identical(unname(sm[, "TCGA-A1-0001-01A"]), c(0.5, -1.2, NA))
    expect_identical(unname(sm[, "TCGA-A1-0002-01A"]), c(NA, NA, 0.1))
})

test_that("[BUG-11] single-sample Mutation becomes an RSE with MAF columns as assays", {
    raw <- accmini@Mutation
    rse <- quiet(biocExtract(accmini, "Mutation"))
    expect_s4_class(rse, "RangedSummarizedExperiment")
    expect_identical(colnames(rse), unique(raw[["Tumor_Sample_Barcode"]]))
    expect_identical(unname(genome(rse)), "GRCh37")
    expect_identical(unname(assay(rse, "Hugo_Symbol")[, 1L]),
        raw[["Hugo_Symbol"]])
    expect_identical(start(rse), as.integer(raw[["Start_position"]]))
    expect_identical(as.character(strand(rse)), raw[["Strand"]])
    expect_length(assayNames(rse), 47L)
})

test_that("multi-sample MAF becomes a RaggedExperiment with GRCh37", {
    maf <- data.frame(Hugo_Symbol = c("TP53", "KRAS", "TP53"),
        NCBI_Build = "37", Chromosome = c("17", "12", "17"),
        Start_position = c(7577120, 25398284, 7577121),
        End_position = c(7577120, 25398284, 7577121),
        Strand = "+", Variant_Classification = c("Missense", "Missense",
            "Nonsense"),
        Tumor_Sample_Barcode = c("TCGA-A1-0001-01A", "TCGA-A1-0001-01A",
            "TCGA-A1-0002-01A"))
    re <- quiet(biocExtract(maf, "Mutation"))
    expect_s4_class(re, "RaggedExperiment")
    expect_identical(unique(unname(genome(re))), "GRCh37")
    vc <- assay(re, "Variant_Classification")
    expect_identical(unname(vc[, "TCGA-A1-0001-01A"]),
        c("Missense", "Missense", NA))
    expect_identical(unname(vc[, "TCGA-A1-0002-01A"]),
        c(NA, NA, "Nonsense"))
})

test_that("[BUG-23] biocExtract errors on bare matrix / non-ranged data.frame", {
    m <- matrix(1:4, 2, dimnames = list(c("a", "b"), c("x", "y")))
    expect_error(quiet(biocExtract(m, "RNASeq2Gene")), "subscript out of bounds")
    df <- data.frame("TCGA.A1.0001.01A" = c(1, 2), "TCGA.A1.0002.01A" = c(3, 4),
        row.names = c("g1", "g2"), check.names = FALSE)
    expect_error(quiet(biocExtract(df, "RNASeq2Gene")), "metadata")
})
