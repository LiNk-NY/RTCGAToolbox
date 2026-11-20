## Column-finding and range-annotation helpers in utils.R

test_that(".findCol matches ignoring case and punctuation", {
    df <- data.frame(Hugo_Symbol = 1, NCBIbuild = 2, Tumor.Sample = 3)
    expect_identical(ns$.findCol(df, "NCBI_Build"), "NCBIbuild")
    expect_identical(ns$.findCol(df, "hugo symbol"), "Hugo_Symbol")
    expect_identical(ns$.findCol(df, "tumor_sample"), "Tumor.Sample")
    expect_identical(ns$.findCol(df, "absent"), character(0L))
    expect_error(ns$.findCol(df, 1L), "colname is not character")
    ## multiple matches: the first one is returned
    df2 <- data.frame(NCBI_Build = 1, ncbi.build = 2)
    expect_identical(ns$.findCol(df2, "NCBI_Build"), "NCBI_Build")
})

test_that(".hasInfo reports presence of a column", {
    expect_true(ns$.hasInfo(data.frame(NCBI_Build = 37)))
    expect_false(ns$.hasInfo(data.frame(x = 37)))
    expect_true(ns$.hasInfo(data.frame(Hugo_Symbol = "A"), "Hugo_Symbol"))
})

test_that(".findUniqueAnnoCol / .setAnnoRows use the first unique annotation", {
    df <- data.frame(Hugo_Symbol = c("A", "B"), Entrez_Gene_Id = c(1, 2))
    expect_identical(ns$.findUniqueAnnoCol(df), "Hugo_Symbol")
    expect_identical(rownames(ns$.setAnnoRows(df)), c("A", "B"))

    ## duplicated symbols fall back to Entrez ids
    df <- data.frame(Hugo_Symbol = c("A", "A"), Entrez_Gene_Id = c(1, 2))
    expect_identical(ns$.findUniqueAnnoCol(df), "Entrez_Gene_Id")
    expect_identical(rownames(ns$.setAnnoRows(df)), c("1", "2"))

    ## nothing unique: row names untouched
    df <- data.frame(Hugo_Symbol = c("A", "A"), x = 1:2)
    expect_identical(ns$.findUniqueAnnoCol(df), character(0L))
    expect_identical(rownames(ns$.setAnnoRows(df)), c("1", "2"))
})

test_that(".findSampleCol finds the sample column by priority", {
    expect_identical(
        ns$.findSampleCol(data.frame(ID = 1, Sample = 2, Tumor_Sample_Barcode = 3)),
        "Tumor_Sample_Barcode")
    expect_identical(ns$.findSampleCol(data.frame(ID = 1, Sample = 2)), "Sample")
    expect_identical(ns$.findSampleCol(data.frame(id = 1)), "id")
    expect_identical(ns$.findSampleCol(data.frame(a = 1)), NA_character_)
})

test_that(".standardstrand maps Firehose strand codes to +, -, *", {
    expect_identical(
        ns$.standardstrand(c("1", "-1", "null", "NULL", NA, "+", "-", "*")),
        c("+", "-", "*", "*", "*", "+", "-", "*"))
    df <- data.frame(Strand = c("1", "-1"))
    expect_identical(ns$.standardizeStrand(df, "Strand")[["Strand"]], c("+", "-"))
})

test_that(".getBuild reads NCBI_Build columns", {
    expect_identical(ns$.getBuild(data.frame(NCBI_Build = c("37", "37"))), "37")
    expect_identical(
        ns$.getBuild(data.frame(NCBI_Build = c("GRCh37", "GRCh37"))), "GRCh37")
    expect_identical(ns$.getBuild(data.frame(x = 1)), NA_character_)
    ## BUG-25: mixed builds give a clear error (was a cryptic vapply error)
    expect_error(ns$.getBuild(data.frame(NCBI_Build = c("36", "37"))),
        "Inconsistent genome builds found: 36, 37")
    expect_error(ns$.getBuild(data.frame(NCBI_Build = c(rep("37", 9), "36"))),
        "Inconsistent genome builds found: 37, 36")
    expect_error(ns$.getBuild(data.frame(NCBI_Build = c("hg19", "hg18"))),
        "Inconsistent genome builds")
    ## equivalent labels are one build
    expect_identical(ns$.getBuild(data.frame(NCBI_Build = c("37", "GRCh37"))),
        "GRCh37")
    expect_identical(ns$.getBuild(data.frame(NCBI_Build = c("hg19", "37"))),
        "GRCh37")
    ## missing values: tolerated below 20%, otherwise an error
    expect_identical(
        ns$.getBuild(data.frame(NCBI_Build = c(rep("37", 9), NA))), "37")
    expect_error(ns$.getBuild(data.frame(NCBI_Build = c("37", NA))),
        "missing in 1 of 2 rows")
    expect_identical(ns$.getBuild(data.frame(NCBI_Build = c(NA, ""))),
        NA_character_)
    expect_identical(ns$.getBuild(data.frame(NCBI_Build = 37L)), "37")
})

test_that(".ansRangeNames maps TCGA column names to GRanges fields", {
    seg <- data.frame(Sample = "s", Chromosome = 1, Start = 1, End = 2,
        Num_Probes = 1, Segment_Mean = 0)
    expect_identical(ns$.ansRangeNames(seg), list(
        seqnames.field = "Chromosome", start.field = "Start",
        end.field = "End", ignore.strand = TRUE))

    maf <- data.frame(Chromosome = 1, Start_position = 1, End_position = 2,
        Strand = "+")
    expect_identical(ns$.ansRangeNames(maf), list(
        seqnames.field = "Chromosome", start.field = "Start_position",
        end.field = "End_position", strand.field = "Strand",
        ignore.strand = FALSE))

    ## strand column that is all NA is ignored
    maf[["Strand"]] <- NA
    expect_true(ns$.ansRangeNames(maf)[["ignore.strand"]])

    expect_identical(ns$.ansRangeNames(list()), list())
})

test_that(".hasRangeNames detects genomic coordinate columns", {
    expect_true(ns$.hasRangeNames(data.frame(Chromosome = 1, Start = 1, End = 2)))
    expect_false(ns$.hasRangeNames(data.frame(Chromosome = 1, Start = 1)))
    expect_false(ns$.hasRangeNames(data.frame(TCGA.1 = 1, TCGA.2 = 2)))
    expect_false(ns$.hasRangeNames(list()))
    expect_error(ns$.hasRangeNames(c(a = 1, b = 2)), "must be rectangular")
    ## unnamed input short-circuits: all(grepl(., NULL)) is TRUE
    expect_false(ns$.hasRangeNames(1:3))
})

test_that(".samplesAsCols flags TCGA-named columns", {
    df <- data.frame(Hugo_Symbol = "A", "TCGA-1" = 1, "tcga-2" = 2,
        check.names = FALSE)
    expect_identical(ns$.samplesAsCols(df), c(FALSE, TRUE, TRUE))
    ## with explicit prefixes, columns *starting with* those prefixes are
    ## flagged instead
    expect_identical(unname(ns$.samplesAsCols(df, "Hugo")), c(TRUE, FALSE, FALSE))
})

test_that(".omitAdditionalIdx drops extra coordinate-like columns", {
    df <- data.frame(Chromosome = 1, Start = 1, End = 2, chr = 1, width = 3,
        x = 1)
    rn <- list(seqnames.field = "Chromosome", start.field = "Start",
        end.field = "End", ignore.strand = TRUE)
    expect_identical(ns$.omitAdditionalIdx(df, rn), c(4L, 5L))
})

test_that(".missingRanges / .rmNAse / .removeNASeq drop incomplete rows", {
    df <- data.frame(Chromosome = c(1, NA, 1), Start = c(1, 2, NA),
        End = c(5, 6, 7))
    rn <- list(start.field = "Start", end.field = "End")
    expect_identical(ns$.missingRanges(df, rn), c(FALSE, FALSE, TRUE))
    expect_identical(nrow(ns$.rmNAse(df, rn)), 2L)
    expect_message(res <- ns$.removeNASeq(df, "Chromosome"),
        "Removing 1 rows")
    expect_identical(res[["Start"]], c(1, NA))
})

test_that(".removeShell maps GISTIC* types to the GISTIC slot", {
    data(accmini, package = "RTCGAToolbox", envir = environment())
    expect_identical(ns$.removeShell(accmini, "GISTICA"), accmini@GISTIC)
    expect_identical(ns$.removeShell(accmini, "CNASNP"), accmini@CNASNP)
})
