## getLinks(), listing parser, file size check, cache dir, getBroadSubtypes.
## All network access is mocked.

fh <- build_fake_firehose()
std_url <- "https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ACC/20160128/"

test_that(".files_from_html_table parses a real Firehose index page", {
    res <- ns$.files_from_html_table(test_path("fixtures", "firehose_listing.html"))
    expect_identical(unname(res), c(
        "gdac.broadinstitute.org_ACC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz",
        "gdac.broadinstitute.org_ACC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz.md5",
        "gdac.broadinstitute.org_ACC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz",
        "gdac.broadinstitute.org_ACC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz.md5",
        "gdac.broadinstitute.org_ACC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg18__seg.Level_3.2016012800.0.0.tar.gz",
        "gdac.broadinstitute.org_ACC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg18__seg.Level_3.2016012800.0.0.tar.gz.md5",
        "gdac.broadinstitute.org_ACC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3.2016012800.0.0.tar.gz",
        "gdac.broadinstitute.org_ACC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3.2016012800.0.0.tar.gz.md5"
    ))
})

test_that(".getLinks filters on both keywords and the dataset prefix", {
    doc <- c("x.ACC.Merge_a__b.tar.gz", "x.ACC.Merge_a__b.tar.gz.md5",
        "x.BRCA.Merge_a__b.tar.gz", "x.ACC.Other.tar.gz")
    expect_identical(ns$.getLinks("Merge_a", "tar[.]gz$", NULL, doc),
        c("x.ACC.Merge_a__b.tar.gz", "x.BRCA.Merge_a__b.tar.gz"))
    expect_identical(ns$.getLinks("Merge_a", "[.]Merge_a.*tar[.]gz$", "ACC", doc),
        "x.ACC.Merge_a__b.tar.gz")
    expect_identical(ns$.getLinks("nothing", "x", NULL, doc), character(0L))
})

test_that("getLinks returns the archive URL for each data type", {
    local_mocked_bindings(.files_from_html_table = function(url) {
        if (grepl("analyses__", url)) fh$ana_listing else fh$std_listing
    }, .package = "RTCGAToolbox")
    lk <- function(...) getLinks("ACC", ...)
    pat <- c(
        clinical = "Clinical_Pick_Tier1.Level_4",
        RNASeq2Gene = "RSEM_genes__data.Level_3",
        miRNASeqGene = "miR_gene_expression__data.Level_3",
        RNASeq2GeneNorm = "RSEM_genes_normalized__data.Level_3",
        CNASNP = "segmented_scna_hg19__seg.Level_3",
        CNVSNP = "segmented_scna_minus_germline_cnv_hg19__seg.Level_3",
        CNASeq = "dnaseqc__hms_harvard_edu__Level_3__segmentation__seg.Level_3",
        Methylation = "humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3",
        Mutation = "Mutation_Packager_Calls.Level_3",
        mRNAArray = "agilentg4502a_07_3__unc_edu__Level_3",
        miRNAArray = "h_mirna_8x15k__unc_edu__Level_3",
        RPPAArray = "mda_rppa_core__mdanderson_org__Level_3"
    )
    for (type in names(pat)) {
        args <- list(TRUE)
        names(args) <- type
        res <- do.call(lk, args)
        expect_length(res, 1L)
        expect_true(startsWith(res, std_url), info = type)
        expect_true(grepl(pat[[type]], res, fixed = TRUE), info = type)
        expect_true(endsWith(res, ".tar.gz"), info = type)
    }
    ## BUG-08: CNACGH includes the MSKCC segmentations
    cgh <- lk(CNACGH = TRUE)
    expect_identical(cgh, paste0(std_url, FH_PREFIX, "Merge_cna__", c(
        "cgh_1x1m_g4447a__mskcc_org__Level_3__segmentation_data_computation__seg",
        "hg_cgh_415k_g4124a__hms_harvard_edu__Level_3__segmentation__seg"),
        FH_SUFFIX, ".tar.gz"))
    g <- lk(GISTIC = TRUE)
    expect_identical(g, paste0(
        "https://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/ACC/20160128/",
        "gdac.broadinstitute.org_ACC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz"))
    pre <- lk(RNASeq2GeneNorm = TRUE, RNAseq2Norm = "raw_counts")
    expect_identical(pre, paste0(std_url,
        "gdac.broadinstitute.org_ACC.mRNAseq_Preprocess.Level_3.2016012800.0.0.tar.gz"))
    ## nothing requested: only the base URL
    expect_identical(lk(), std_url)
})

test_that(".checkFileSize compares content-length (bytes) to a MB limit", {
    local_mocked_bindings(HEAD = function(url, ...)
        list(headers = list(`content-length` = "2500000")), .package = "httr")
    expect_true(ns$.checkFileSize("http://x", 3))
    msgs <- capture_messages(res <- ns$.checkFileSize("http://x", 2))
    expect_false(res)
    expect_match(paste(msgs, collapse = ""), "File Size: ~2.5MB")
    expect_match(paste(msgs, collapse = ""), "won't be downloaded")
    ## exactly at the limit is allowed
    expect_true(ns$.checkFileSize("http://x", 2.5))
})

test_that(".setCache creates and returns the cache directory", {
    d <- file.path(withr::local_tempdir(), "cache", "nested")
    expect_message(res <- ns$.setCache(d, verbose = TRUE, ask = FALSE),
        "cache directory set to")
    expect_identical(res, d)
    expect_true(dir.exists(d))
    expect_silent(ns$.setCache(d, verbose = FALSE))
    expect_error(ns$.setCache(c("a", "b")))
    expect_error(ns$.setCache(NA_character_))
    ## non-interactive 'ask' answers "n"
    expect_identical(ns$.getAnswer("q", c("y", "n")), "n")
    expect_error(ns$.setCache(file.path(d, "new"), ask = TRUE),
        "not created")
})

test_that("getBroadSubtypes builds the URL and falls back to mRNAseq", {
    seen <- character(0L)
    local_mocked_bindings(url.exists = function(url, ...) FALSE,
        .package = "RCurl")
    local_mocked_bindings(read.delim = function(file, ...) {
        seen <<- c(seen, file)
        data.frame(SampleName = c("TCGA-A1", "TCGA-A2"), cluster = c(1L, 2L),
            silhouetteValue = c(0.5, 0.25))
    }, .package = "RTCGAToolbox")
    expect_warning(res <- getBroadSubtypes("SKCM", "CNMF"),
        "mRNAseq clustering is taken instead")
    expect_identical(seen, paste0(
        "https://gdac.broadinstitute.org/runs/analyses__latest/reports/cancer/",
        "SKCM-TM/mRNAseq_Clustering_CNMF/SKCM-TM.bestclus.txt"))
    expect_identical(rownames(res), c("TCGA-A1", "TCGA-A2"))
    expect_identical(res[["cluster"]], c(1L, 2L))
    expect_error(getBroadSubtypes("SKCM"), "valid clustering")
    expect_error(getBroadSubtypes(c("A", "B"), "CNMF"), "valid cancer code")
})
