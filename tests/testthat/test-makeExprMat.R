## .makeExprMat parses Firehose expression matrices (one or two header rows).

write_fx <- function(lines) {
    f <- tempfile(fileext = ".txt")
    writeLines(lines, f)
    f
}

test_that(".makeExprMat selects columns by normalisation prefix", {
    f <- write_fx(fx_rnaseq2)
    raw <- quiet(ns$.makeExprMat(f, "raw_count", "RNAseq"))
    expect_true(is.matrix(raw))
    expect_identical(storage.mode(raw), "double")
    expect_identical(colnames(raw), unname(BC[c("S1", "S2")]))
    expect_identical(rownames(raw),
        c("A1BG", "A2M", "SLC35E2|728661", "SLC35E2|9906"))
    expect_identical(unname(raw[, 1]), c(197, 5000.5, 11, 99))
    expect_identical(unname(raw[, 2]), c(23.84, 10, 12, 98))

    se <- quiet(ns$.makeExprMat(f, "scaled_estimate", "RNAseq"))
    expect_identical(unname(se[1:3, 1]), c(6.6e-06, 1.5e-04, 4e-07))
    expect_identical(unname(se[1:3, 2]), c(1.2e-06, 3.0e-07, 5e-07))
})

test_that(".makeExprMat errors when no column matches normMethod", {
    f <- write_fx(fx_rnaseq2)
    expect_error(quiet(ns$.makeExprMat(f, "RPKM", "RNAseq")),
        "No data columns found for 'RPKM'")
    expect_error(quiet(ns$.makeExprMat(f, "raw", "RNAseq", headerRows = 1L)),
        "requires two header rows")
})

test_that(".makeExprMat drops '?|' rows and strips Entrez ids from names", {
    f <- write_fx(fx_rnaseq2)
    raw <- quiet(ns$.makeExprMat(f, "raw_count", "RNAseq"))
    expect_false(any(grepl("?", rownames(raw), fixed = TRUE)))
    expect_identical(rownames(raw)[1:2], c("A1BG", "A2M"))
})

test_that(".makeExprMat keeps every gene that shares a symbol", {
    ## BUG-05: SLC35E2|728661 and SLC35E2|9906 are distinct genes; the
    ## second (99, 98) used to be dropped. Both keep their full identifier.
    f <- write_fx(fx_rnaseq2)
    raw <- quiet(ns$.makeExprMat(f, "raw_count", "RNAseq"))
    expect_identical(rownames(raw),
        c("A1BG", "A2M", "SLC35E2|728661", "SLC35E2|9906"))
    expect_identical(unname(raw["SLC35E2|728661", ]), c(11, 12))
    expect_identical(unname(raw["SLC35E2|9906", ]), c(99, 98))
})

test_that(".geneRowNames disambiguates only duplicated symbols", {
    expect_identical(ns$.geneRowNames(c("A|1", "B|2", "A|3", "C")),
        c("A|1", "B", "A|3", "C"))
    expect_identical(ns$.geneRowNames(character(0L)), character(0L))
    expect_warning(res <- ns$.geneRowNames(c("A|1", "A|1")),
        "Duplicated row identifiers: A\\|1")
    expect_identical(res, c("A|1", "A|1"))
})

test_that(".makeExprMat with a prefix matching no column keeps only ids", {
    ## 'raw' is a prefix of 'raw_count', so it matches like raw_count
    f <- write_fx(fx_rnaseq2)
    expect_identical(
        quiet(ns$.makeExprMat(f, "raw", "RNAseq")),
        quiet(ns$.makeExprMat(f, "raw_count", "RNAseq")))
})

test_that(".makeExprMat handles miRNA-seq read_count columns", {
    f <- write_fx(fx_mirna)
    m <- quiet(ns$.makeExprMat(f, "read_count", "miRNAseq"))
    expect_identical(rownames(m), c("hsa-let-7a-1", "hsa-let-7a-2"))
    expect_identical(unname(m[, 1]), c(76213, 151321))
    expect_identical(unname(m[, 2]), c(45441, 90039))
    rpm <- quiet(ns$.makeExprMat(f, "reads_per_million_miRNA_mapped", "miRNAseq"))
    expect_identical(unname(rpm[, 1]), c(15000.5, 30000.75))
})

test_that(".makeExprMat arrayData = TRUE takes every data column", {
    f <- write_fx(fx_mrnaarray)
    expect_no_warning(
        m <- quiet(ns$.makeExprMat(f, "", "mRNAArray", arrayData = TRUE)))
    expect_identical(colnames(m), unname(BC[c("S1", "S2")]))
    expect_identical(rownames(m), c("A1BG", "A2M"))
    ## 'null' strings become NA
    expect_identical(unname(m[, 1]), c(-0.5, NA))
    expect_identical(unname(m[, 2]), c(0.25, 1.5))

    f <- write_fx(fx_rnaseq2norm)
    m <- quiet(ns$.makeExprMat(f, "normalized_counts", "RNAseq2",
        arrayData = TRUE))
    expect_identical(colnames(m), unname(BC))
    expect_identical(rownames(m), c("A1BG", "A2M"))
    expect_identical(unname(m["A1BG", ]), c(16.3305, 9.5987, 20.7377))
})

test_that(".makeExprMat works with one gene row or one sample", {
    ## BUG-06: both used to error
    one_row <- write_fx(c(
        .tsv("Hybridization REF", BC[["S1"]], BC[["S2"]]),
        .tsv("gene_id", "raw_count", "raw_count"),
        .tsv("A1BG|1", "1", "2")))
    expect_identical(quiet(ns$.makeExprMat(one_row, "raw_count", "RNAseq")),
        matrix(c(1, 2), 1L, dimnames = list("A1BG", unname(BC[c("S1", "S2")]))))
    one_col <- write_fx(c(
        .tsv("Hybridization REF", BC[["S1"]]),
        .tsv("gene_id", "raw_count"),
        .tsv("A1BG|1", "1"),
        .tsv("A2M|2", "3")))
    expect_identical(quiet(ns$.makeExprMat(one_col, "raw_count", "RNAseq")),
        matrix(c(1, 3), 2L, dimnames = list(c("A1BG", "A2M"), BC[["S1"]])))
})

test_that(".makeExprMat reads one-header-row Preprocess files in full", {
    ## BUG-26: the first data row used to be dropped
    f <- write_fx(fx_pre_raw)
    m <- quiet(ns$.makeExprMat(f, "raw_counts", "RNAseq2", arrayData = TRUE,
        headerRows = 1L))
    expect_identical(m, matrix(c(6.8, 100, 300, 6.66, 200, 400), 3L,
        dimnames = list(c("A1BG", "A2M", "TP53"), PBC)))
    f <- write_fx(fx_pre_log2)
    m <- quiet(ns$.makeExprMat(f, "RSEM_normalized_log2", "RNAseq2",
        arrayData = TRUE, headerRows = 1L))
    expect_identical(rownames(m), c("A1BG", "A2M", "TP53"))
    expect_identical(unname(m["A1BG", ]), c(1.7075689633722, 1.42304025336151))
    expect_identical(unname(m["A2M", ]), c(2.5, NA))
})

test_that(".makeExprMat errors on non-numeric values", {
    bad <- write_fx(c(
        .tsv("Hybridization REF", BC[["S1"]]),
        .tsv("Composite Element REF", "x"),
        .tsv("A1BG", "1.5abc")))
    expect_error(quiet(ns$.makeExprMat(bad, "", "mRNAArray", arrayData = TRUE)),
        "Non-numeric values in mRNAArray data: 1.5abc")
})
