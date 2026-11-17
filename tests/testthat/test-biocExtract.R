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

test_that("CNASNP/CNVSNP (single sample) become a RaggedExperiment", {
    ## BUG-11: one sample used to give a RangedSummarizedExperiment
    for (type in c("CNASNP", "CNVSNP")) {
        raw <- slot(accmini, type)
        re <- quiet(biocExtract(accmini, type))
        expect_s4_class(re, "RaggedExperiment")
        expect_identical(colnames(re), unique(raw[["Sample"]]), info = type)
        expect_identical(assayNames(re), c("Num_Probes", "Segment_Mean"))
        expect_identical(unname(assay(re, "Segment_Mean")[, 1L]),
            raw[["Segment_Mean"]], info = type)
        expect_identical(unname(assay(re, "Num_Probes")[, 1L]),
            raw[["Num_Probes"]], info = type)
        gr <- rowRanges(re)
        expect_identical(start(gr), as.integer(raw[["Start"]]), info = type)
        expect_identical(end(gr), as.integer(raw[["End"]]), info = type)
        expect_identical(as.character(seqnames(gr)),
            as.character(raw[["Chromosome"]]), info = type)
        ## BUG-17: the SNP6 segment archives are hg19 (seqnames '1', '2', ...)
        expect_identical(unique(unname(genome(re))), "GRCh37", info = type)
    }
})

test_that("the hg19 default applies only to CNASNP/CNVSNP from FirehoseData", {
    ## the same table outside a FirehoseData object has no known build
    re <- quiet(biocExtract(accmini@CNASNP, "CNASNP"))
    expect_true(all(is.na(genome(re))))
    ## an NCBI_Build column takes precedence over the default
    x <- accmini
    x@CNASNP[["NCBI_Build"]] <- "36"
    re <- quiet(biocExtract(x, "CNASNP"))
    expect_identical(unique(unname(genome(re))), "NCBI36")
})

test_that("segments with identical ranges across samples stay ragged", {
    seg <- data.frame(
        Sample = rep(c("TCGA-A1-0001-01A", "TCGA-A1-0002-01A"), each = 2L),
        Chromosome = 1, Start = c(1, 100, 1, 100), End = c(50, 200, 50, 200),
        Num_Probes = 1:4, Segment_Mean = c(0.1, 0.2, 0.3, 0.4))
    for (type in c("CNASNP", "CNVSNP", "CNASeq", "CNACGH")) {
        re <- quiet(biocExtract(seg, type))
        expect_s4_class(re, "RaggedExperiment")
        expect_identical(dim(re), c(4L, 2L), info = type)
        ## compactAssay() recovers the matrix form, values in place
        ca <- compactAssay(re, "Segment_Mean")
        expect_identical(unname(ca), matrix(c(0.1, 0.2, 0.3, 0.4), 2L),
            info = type)
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

test_that("single-sample Mutation becomes a RaggedExperiment", {
    ## BUG-11: used to be an RSE with every MAF column as an assay
    raw <- accmini@Mutation
    re <- quiet(biocExtract(accmini, "Mutation"))
    expect_s4_class(re, "RaggedExperiment")
    expect_identical(colnames(re), unique(raw[["Tumor_Sample_Barcode"]]))
    expect_identical(unique(unname(genome(re))), "GRCh37")
    expect_identical(unname(assay(re, "Hugo_Symbol")[, 1L]),
        raw[["Hugo_Symbol"]])
    gr <- rowRanges(re)
    expect_identical(start(gr), as.integer(raw[["Start_position"]]))
    expect_identical(as.character(strand(gr)), raw[["Strand"]])
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

test_that("biocExtract accepts a bare matrix or a plain data.frame", {
    ## BUG-23: both used to error
    m <- matrix(1:4, 2, dimnames = list(c("a", "b"), c("x", "y")))
    se <- quiet(biocExtract(m, "RNASeq2Gene"))
    expect_s4_class(se, "SummarizedExperiment")
    expect_identical(assay(se), m)
    df <- data.frame("TCGA.A1.0001.01A" = c(1, 2), "TCGA.A1.0002.01A" = c(3, 4),
        row.names = c("g1", "g2"), check.names = FALSE)
    se <- quiet(biocExtract(df, "RNASeq2Gene"))
    expect_identical(assay(se), matrix(c(1, 2, 3, 4), 2, dimnames = list(
        c("g1", "g2"), c("TCGA-A1-0001-01A", "TCGA-A1-0002-01A"))))
    expect_identical(metadata(se), list())
})

test_that("biocExtract never turns character columns into factor codes", {
    df <- data.frame(TCGA.A1.0001.01A = c("high", "low"),
        TCGA.A1.0002.01A = c(3, 4), row.names = c("g1", "g2"))
    expect_error(quiet(biocExtract(df, "RNASeq2Gene")),
        "non-numeric columns to an assay: TCGA-A1-0001-01A")
})

test_that("biocExtract on a FirehoseGISTIC object", {
    data(accmini, package = "RTCGAToolbox", envir = environment())
    expect_identical(quiet(biocExtract(accmini@GISTIC, "GISTICT")),
        quiet(biocExtract(accmini, "GISTICT")))
})
