## getFirehoseData() end to end against a fake, offline Firehose mirror.
## Expected values below are taken directly from the fixture text in
## helper-fixtures.R, not from the code under test.

fh <- build_fake_firehose()

all_types <- c("clinical", "RNASeq2Gene", "RNASeq2GeneNorm", "miRNASeqGene",
    "CNASNP", "CNVSNP", "CNASeq", "CNACGH", "Methylation", "mRNAArray",
    "miRNAArray", "RPPAArray", "Mutation", "GISTIC")

fetch <- function(types, destdir, ...) {
    args <- list(dataset = FH_DATASET, clinical = FALSE, destdir = destdir, ...)
    for (t in types) args[[t]] <- TRUE
    suppressWarnings(quiet(do.call(getFirehoseData, args)))
}

## ---- argument validation ----------------------------------------------

test_that("getFirehoseData validates dataset and runDate", {
    expect_error(getFirehoseData(), "dataset")
    expect_error(getFirehoseData(1), "dataset")
    expect_error(getFirehoseData(c("ACC", "BRCA")), "dataset")
    expect_error(getFirehoseData("A"), "dataset")
    expect_error(getFirehoseData("NOTACOHORT", destdir = tempdir()),
        "valid dataset name")
    expect_error(getFirehoseData("ACC", runDate = "2016-01-28",
        destdir = tempdir()), "runDate")
    expect_error(getFirehoseData("ACC", runDate = "20991231",
        destdir = tempdir()), "valid run date")
})

test_that("an invalid gistic2Date is rejected before any download", {
    ## BUG-07: '&&' let any single string through
    cnt <- local_fake_firehose(fh)
    dest <- withr::local_tempdir()
    for (bad in list("20991231", "2016-01-28", c("20160128", "20150821"),
            NA_character_, 20160128))
        expect_error(fetch("GISTIC", dest, gistic2Date = bad),
            'valid "gistic2Date"', info = paste(bad, collapse = ","))
    expect_length(cnt$urls_listed, 0L)
    expect_length(cnt$downloads, 0L)
    ## not checked when GISTIC = FALSE
    expect_no_error(fetch("clinical", dest, gistic2Date = "20991231"))
})

## ---- URLs and file selection -----------------------------------------

test_that("listing URLs follow the Firehose layout", {
    cnt <- local_fake_firehose(fh)
    dest <- withr::local_tempdir()
    fetch(c("clinical", "GISTIC"), dest)
    expect_identical(cnt$urls_listed, c(
        "https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ACC/20160128/",
        "https://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/ACC-TP/20160128/"
    ))
})

test_that("each data type downloads exactly its archive (no decoys)", {
    cnt <- local_fake_firehose(fh)
    dest <- withr::local_tempdir()
    fetch(all_types, dest)
    got <- basename(cnt$downloads)
    expect_length(got, 14L)
    expect_false(any(grepl("\\.aux\\.|mage-tab|\\.md5$|hg18|isoform|exon_quant|Preprocess|Oncotated|Coverage|LowPass", got)))
    ## everything on the mirror except the mRNAseq_Preprocess archive,
    ## which is only used for non-default RNAseq2Norm values
    expect_setequal(got, c(
        grep("mRNAseq_Preprocess", list.files(fh$std, pattern = "tar.gz$"),
            value = TRUE, invert = TRUE),
        list.files(fh$ana, pattern = "tar.gz$")))
    ## extracted directories and tarballs are cleaned up
    expect_setequal(list.files(dest), c(
        "20160128-ACC-Clinical.txt",
        "20160128-ACC-RNAseq2Gene.txt",
        "20160128-ACC-RNAseq2GeneNorm-normalized_counts-1.txt",
        "20160128-ACC-read_count-miRNAseqGene.txt",
        "20160128-ACC-CNASNPHg19.txt",
        "20160128-ACC-CNVSNPHg19.txt",
        "20160128-ACC-CNAseq.txt",
        "20160128-ACC-CNACGH-1.txt",
        "20160128-ACC-Methylation-1.txt",
        "20160128-ACC-mRNAArray-1.txt",
        "20160128-ACC-miRNAArray-1.txt",
        "20160128-ACC-RPPAArray-1.txt",
        "20160128-ACC-Mutations-AllSamples.txt",
        "20160128-ACC-all_data_by_genes.txt",
        "20160128-ACC-all_thresholded.by_genes.txt",
        "20160128-ACC-all_lesions.conf_99.txt"))
})

## ---- parsed contents ---------------------------------------------------

local({
    cnt <- local_fake_firehose(fh)
    dest <- withr::local_tempdir()
    x <- fetch(all_types, dest)

    test_that("object metadata", {
        expect_s4_class(x, "FirehoseData")
        expect_identical(x@Dataset, "ACC")
        expect_identical(x@runDate, "20160128")
        expect_identical(x@gistic2Date, "20160128")
        expect_true(validObject(x))
    })

    test_that("clinical: one row per patient, attributes as columns", {
        cl <- x@clinical
        expect_identical(rownames(cl), c("tcga.or.a5j1", "tcga.or.a5j2"))
        expect_identical(cl[["years_to_birth"]], c("58", "44"))
        expect_identical(cl[["vital_status"]], c("1", "0"))
        expect_identical(cl[["gender"]], c("male", "female"))
        ## BUG-09: the 'Composite Element REF' header row is not a column
        expect_identical(names(cl), c("years_to_birth", "vital_status",
            "gender"))
        ## patient ids stay as read.delim() names them (lower case, dots)
    })

    test_that("RNASeq2Gene: raw_count columns", {
        m <- x@RNASeq2Gene
        rn <- c("A1BG", "A2M", "SLC35E2|728661", "SLC35E2|9906")
        expect_identical(dimnames(m), list(rn, unname(BC[c("S1", "S2")])))
        expect_identical(m[, BC[["S1"]]],
            setNames(c(197, 5000.5, 11, 99), rn))
        expect_identical(m[, BC[["S2"]]],
            setNames(c(23.84, 10, 12, 98), rn))
    })

    test_that("RNASeq2GeneNorm: list of FirehosemRNAArray", {
        expect_length(x@RNASeq2GeneNorm, 1L)
        a <- x@RNASeq2GeneNorm[[1L]]
        expect_s4_class(a, "FirehosemRNAArray")
        expect_match(a@Filename, "RSEM_genes_normalized__data.Level_3")
        expect_identical(a@DataMatrix["A1BG", ],
            setNames(c(16.3305, 9.5987, 20.7377), BC))
        expect_identical(a@DataMatrix["A2M", ],
            setNames(c(1000.1, 2000.2, 3000.3), BC))
    })

    test_that("miRNASeqGene: read_count by default", {
        m <- x@miRNASeqGene
        expect_identical(m[, BC[["S1"]]],
            c("hsa-let-7a-1" = 76213, "hsa-let-7a-2" = 151321))
        expect_identical(m[, BC[["S2"]]],
            c("hsa-let-7a-1" = 45441, "hsa-let-7a-2" = 90039))
    })

    test_that("segmented copy-number slots", {
        expected <- data.frame(
            Sample = rep(c("TCGA-OR-A5J1-01A-11D-A29I-01",
                "TCGA-OR-A5J2-01A-11D-A29I-01"), each = 2),
            Chromosome = c(1, 2, 1, 23),
            Start = c(61735, 100, 61735, 1000),
            End = c(12834253, 50000, 20000000, 2000),
            Num_Probes = c(6520, 48, 9000, 8),
            Segment_Mean = c(0.0339, -0.3543, 0.5, -1.25))
        expect_identical(x@CNASNP, expected)
        expect_identical(x@CNVSNP, expected)
        expect_identical(x@CNASeq, expected)
        expect_length(x@CNACGH, 1L)
        expect_s4_class(x@CNACGH[[1L]], "FirehoseCGHArray")
        expect_identical(x@CNACGH[[1L]]@DataMatrix, expected)
        expect_match(x@CNACGH[[1L]]@Filename, "hg_cgh_415k_g4124a")
    })

    test_that("Methylation: beta values with probe annotation", {
        expect_length(x@Methylation, 1L)
        dm <- x@Methylation[[1L]]@DataMatrix
        expect_identical(rownames(dm),
            c("cg00000029", "cg00000108", "cg00000109"))
        expect_identical(names(dm), c("Gene_Symbol", "Chromosome",
            "Genomic_Coordinate", "TCGA-OR-A5J1-01A-11D-A29J-05",
            "TCGA-OR-A5J2-01A-11D-A29J-05"))
        expect_identical(dm[["Gene_Symbol"]], c("RBL2", "C3orf35", "FNDC3B"))
        expect_identical(dm[["TCGA-OR-A5J1-01A-11D-A29J-05"]],
            c("0.119877013723081", NA, "0.9"))
        expect_identical(dm[["TCGA-OR-A5J2-01A-11D-A29J-05"]],
            c("0.107120474727399", "0.5", NA))
    })

    test_that("array platforms", {
        m <- x@mRNAArray[[1L]]@DataMatrix
        expect_identical(m["A1BG", ], setNames(c(-0.5, 0.25), BC[1:2]))
        expect_identical(m["A2M", ], setNames(c(NA, 1.5), BC[1:2]))
        m <- x@miRNAArray[[1L]]@DataMatrix
        expect_identical(m["hsa-let-7b", ], setNames(c(8.125, -2), BC[1:2]))
        m <- x@RPPAArray[[1L]]@DataMatrix
        expect_identical(dimnames(m), list(c("14-3-3_beta", "4E-BP1_pS65"),
            c("TCGA-OR-A5J2-01A-21-A39K-20", "TCGA-OR-A5J3-01A-21-A39K-20")))
        expect_identical(unname(m["4E-BP1_pS65", ]), c(-1.5, 2.75))
    })

    test_that("Mutation: per-sample MAFs stacked by column name", {
        mu <- x@Mutation
        expect_identical(nrow(mu), 3L)
        expect_true("Extra_Col" %in% names(mu))
        key <- paste(mu[["Hugo_Symbol"]], mu[["Start_position"]],
            mu[["Tumor_Sample_Barcode"]])
        expect_setequal(key, c(
            "TP53 7577120 TCGA-OR-A5J1-01A-11D-A29I-10",
            "KRAS 25398284 TCGA-OR-A5J1-01A-11D-A29I-10",
            "TP53 7578406 TCGA-OR-A5J2-01A-11D-A29I-10"))
        expect_identical(
            mu[["Variant_Classification"]][mu[["Start_position"]] == "7578406"],
            "Nonsense_Mutation")
        expect_identical(mu[["Extra_Col"]][mu[["Start_position"]] == "7578406"],
            "x")
        expect_true(all(is.na(
            mu[["Extra_Col"]][mu[["Start_position"]] != "7578406"])))
        ## BUG-20: rows stay in file order (merge() used to re-sort them)
        expect_identical(mu[["Hugo_Symbol"]], c("TP53", "KRAS", "TP53"))
        expect_identical(names(mu)[length(mu)], "Extra_Col")
    })

    test_that("GISTIC: three tables, read as character", {
        g <- x@GISTIC
        expect_s4_class(g, "FirehoseGISTIC")
        expect_identical(g@Dataset, "ACC")
        expect_identical(names(g@AllByGene), c("Gene.Symbol", "Locus.ID",
            "Cytoband", "TCGA.OR.A5J1.01A.11D.A29H.01",
            "TCGA.OR.A5J2.01A.11D.A29H.01"))
        expect_identical(g@AllByGene[["TCGA.OR.A5J2.01A.11D.A29H.01"]],
            c("-0.070", "-1.250", "0.900"))
        expect_identical(g@ThresholdedByGene[["TCGA.OR.A5J1.01A.11D.A29H.01"]],
            c("0", "1", "-2"))
        expect_identical(nrow(g@Peaks), 4L)
        expect_identical(g@Peaks[["Unique.Name"]][4L],
            "Amplification Peak  1 - CN values")
    })

    test_that("a second call reuses the cache and gives identical data", {
        n <- length(cnt$downloads)
        y <- fetch(all_types, dest)
        expect_length(cnt$downloads, n)
        for (s in slotNames(x))
            expect_identical(slot(y, s), slot(x, s), info = s)
    })

    test_that("forceDownload = TRUE downloads again, same data", {
        n <- length(cnt$downloads)
        y <- fetch(c("Mutation", "GISTIC", "RNASeq2Gene"), dest,
            forceDownload = TRUE)
        expect_length(cnt$downloads, n + 3L)
        expect_identical(y@Mutation, x@Mutation)
        expect_identical(y@GISTIC, x@GISTIC)
        expect_identical(y@RNASeq2Gene, x@RNASeq2Gene)
    })
})

## ---- other options ----------------------------------------------------

test_that("miRNASeqGeneType selects the column and the cache file name", {
    cnt <- local_fake_firehose(fh)
    dest <- withr::local_tempdir()
    x <- fetch("miRNASeqGene", dest,
        miRNASeqGeneType = "reads_per_million_miRNA_mapped")
    expect_identical(unname(x@miRNASeqGene[, BC[["S1"]]]), c(15000.5, 30000.75))
    expect_true(file.exists(file.path(dest,
        "20160128-ACC-reads_per_million_miRNA_mapped-miRNAseqGene.txt")))
    expect_error(fetch("miRNASeqGene", dest, miRNASeqGeneType = "bogus"))
})

test_that("RNAseq2Norm = RSEM_normalized_log2 reads mRNAseq_Preprocess", {
    cnt <- local_fake_firehose(fh)
    dest <- withr::local_tempdir()
    x <- fetch("RNASeq2GeneNorm", dest, RNAseq2Norm = "RSEM_normalized_log2")
    expect_identical(basename(cnt$downloads),
        "gdac.broadinstitute.org_ACC.mRNAseq_Preprocess.Level_3.2016012800.0.0.tar.gz")
    expect_true(file.exists(file.path(dest,
        "20160128-ACC-RNAseq2GeneNorm-RSEM_normalized_log2-1.txt")))
    m <- x@RNASeq2GeneNorm[[1L]]@DataMatrix
    expect_identical(colnames(m), PBC)
    ## PARADIGM / Z_Score files (value -999) are not selected
    expect_false(any(m == -999, na.rm = TRUE))
    expect_identical(m["TP53", ], setNames(c(10.25, 11.5), PBC))
    expect_identical(m["A2M", ], setNames(c(2.5, NA), PBC))
    ## BUG-26: these files have one header row; the first data row (A1BG)
    ## used to be dropped
    expect_identical(m["A1BG", ],
        setNames(c(1.7075689633722, 1.42304025336151), PBC))
    expect_identical(rownames(m), c("A1BG", "A2M", "TP53"))

    y <- fetch("RNASeq2GeneNorm", dest, RNAseq2Norm = "raw_counts")
    m <- y@RNASeq2GeneNorm[[1L]]@DataMatrix
    expect_identical(m["TP53", ], setNames(c(300, 400), PBC))
    expect_identical(m["A1BG", ], setNames(c(6.8, 6.66), PBC))
    expect_identical(rownames(m), c("A1BG", "A2M", "TP53"))
    expect_error(fetch("RNASeq2GeneNorm", dest, RNAseq2Norm = "bogus"))
})

test_that("files above fileSizeLimit are skipped", {
    cnt <- local_fake_firehose(fh)
    local_mocked_bindings(.checkFileSize = function(dataURL, fileSizeLimit)
        FALSE, .package = "RTCGAToolbox")
    dest <- withr::local_tempdir()
    x <- fetch(all_types, dest)
    expect_length(cnt$downloads, 0L)
    expect_identical(dim(x@RNASeq2Gene), c(0L, 0L))
    expect_identical(x@CNACGH, list())
    expect_identical(x@RNASeq2GeneNorm, list())
    expect_true(isEmpty(x@GISTIC))
})

test_that("missing destdir is an error", {
    cnt <- local_fake_firehose(fh)
    expect_error(fetch("RNASeq2Gene", file.path(tempdir(), "does-not-exist")),
        "Directory does not exist")
})

test_that("runDate = NULL skips standard data; GISTIC uses cohort tags", {
    cnt <- local_fake_firehose(fh)
    dest <- withr::local_tempdir()
    x <- fetch(c("RNASeq2Gene", "GISTIC"), dest, runDate = NULL)
    expect_identical(x@runDate, NA_character_)
    expect_identical(dim(x@RNASeq2Gene), c(0L, 0L))
    expect_false(isEmpty(x@GISTIC))
    expect_false(any(grepl("stddata", cnt$urls_listed)))
})

test_that("GISTIC listing tag is TM for SKCM and TB for LAML", {
    seen <- character(0L)
    local_mocked_bindings(
        .files_from_html_table = function(url) {
            seen <<- c(seen, url)
            character(0L)
        },
        .package = "RTCGAToolbox")
    for (ds in c("SKCM", "LAML", "BRCA"))
        quiet(getFirehoseData(ds, runDate = NULL, clinical = FALSE,
            GISTIC = TRUE, destdir = tempdir()))
    expect_identical(sub(".*/data/([^/]+)/.*", "\\1", seen),
        c("SKCM-TM", "LAML-TB", "BRCA-TP"))
})

test_that("gistic2Date slot is NA when GISTIC not requested", {
    cnt <- local_fake_firehose(fh)
    x <- fetch(character(0L), withr::local_tempdir())
    expect_identical(x@gistic2Date, NA_character_)
    expect_length(cnt$downloads, 0L)
})
