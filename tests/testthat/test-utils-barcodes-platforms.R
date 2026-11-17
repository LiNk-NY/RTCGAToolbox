## Barcode standardisation, platform-name and genome-build parsing.
##
## Tests tagged [BUG-nn] pin *current* behaviour that is known to be wrong;
## see dev/REWORK.md. They are flipped to the correct expectation in the
## commit that fixes the bug.

test_that(".stdIDs converts dotted TCGA barcodes to dashes and upper case", {
    dotted <- c("TCGA.OR.A5J1.01A.11D.A29H.01", "TCGA.OR.A5J2.01A.11D.A29H.01")
    expect_identical(
        ns$.stdIDs(dotted),
        c("TCGA-OR-A5J1-01A-11D-A29H-01", "TCGA-OR-A5J2-01A-11D-A29H-01")
    )
    ## already standard: unchanged
    std <- c("TCGA-OR-A5J1-01A", "TCGA-OR-A5J2-01A")
    expect_identical(ns$.stdIDs(std), std)
    ## upper-casing applies when every id starts with 'TCGA'
    expect_identical(ns$.stdIDs("TCGA-or-a5j1"), "TCGA-OR-A5J1")
})

test_that(".stdIDs leaves non-TCGA identifiers untouched", {
    expect_identical(ns$.stdIDs(c("foo.bar", "baz")), c("foo.bar", "baz"))
    ## lower-case 'tcga' prefix is not recognised (case-sensitive startsWith)
    expect_identical(ns$.stdIDs("tcga.or.a5j1"), "tcga.or.a5j1")
    ## a single non-TCGA id disables conversion for the whole vector
    expect_identical(
        ns$.stdIDs(c("TCGA.OR.A5J1", "Hugo_Symbol")),
        c("TCGA.OR.A5J1", "Hugo_Symbol")
    )
})

test_that("[BUG-01] .stdIDs result on mixed input depends on the RNG", {
    ## One dotted and one dashed barcode: whether the dots are replaced
    ## depends on which elements sample() draws.
    mixed <- c("TCGA-OR-A5J1-01A", "TCGA.OR.A5J2.01A")
    outs <- vapply(1:2000, function(s)
        withr::with_seed(s, paste(ns$.stdIDs(mixed), collapse = " ")),
        character(1L))
    expect_setequal(unique(outs), c(
        "TCGA-OR-A5J1-01A TCGA.OR.A5J2.01A",
        "TCGA-OR-A5J1-01A TCGA-OR-A5J2-01A"
    ))
})

test_that(".standardizeBC standardises column names only", {
    m <- matrix(1:4, 2, dimnames = list(c("a.b", "c.d"),
        c("TCGA.A1.0001.01A", "TCGA.A1.0002.01A")))
    res <- ns$.standardizeBC(m)
    expect_identical(colnames(res), c("TCGA-A1-0001-01A", "TCGA-A1-0002-01A"))
    expect_identical(rownames(res), c("a.b", "c.d"))
    expect_identical(unname(res), unname(m))
})

test_that(".findBuild extracts UCSC and NCBI builds from file names", {
    expect_identical(
        ns$.findBuild("BRCA.Merge_cna__hms_harvard_edu__hg19__Level_3"), "hg19")
    expect_identical(ns$.findBuild("foo_HG18_bar"), "HG18")
    expect_identical(ns$.findBuild("x_GRCh38_y", "NCBI"), "GRCh38")
    expect_identical(ns$.findBuild("no build here"), NA_character_)
    expect_identical(ns$.findBuild(character(0L)), NA_character_)
    ## first match only
    expect_identical(ns$.findBuild("hg18_and_hg19"), "hg18")
    ## vectorised
    expect_identical(ns$.findBuild(c("a_hg19", "b")), c("hg19", NA))
    expect_identical(ns$.findBuild(c(NA, "x_hg38")), c(NA, "hg38"))
    expect_identical(ns$.findBuild(NA_character_), NA_character_)
    expect_error(ns$.findBuild("hg19", "Ensembl"))
})

test_that(".nameClean shortens platform words", {
    expect_identical(ns$.nameClean("humanmethylation450"), "methyl450")
    expect_identical(ns$.nameClean("agilentg4502a"), "g4502a")
    expect_identical(ns$.nameClean("Transcriptome"), "tx")
    expect_identical(ns$.nameClean(character(0L)), character(0L))
})

test_that(".searchPlatform derives short platform names from Firehose files", {
    fn <- function(core)
        paste0("gdac.broadinstitute.org_ACC.", core,
            ".Level_3.2016012800.0.0.tar.gz")
    expect_identical(ns$.searchPlatform(fn(
        "Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data")),
        "methyl450")
    expect_identical(ns$.searchPlatform(fn(
        "Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data")),
        "methyl27")
    expect_identical(ns$.searchPlatform(fn(
        "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data")),
        "illuminahiseq")
    expect_identical(ns$.searchPlatform(fn(
        "Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data")),
        "TX_g4502a")
    expect_identical(ns$.searchPlatform(fn(
        "Merge_transcriptome__ht_hg_u133a__broad_mit_edu__Level_3__gene_rma__data")),
        "TX_ht_hg_u133a")
    expect_identical(ns$.searchPlatform(fn(
        "Merge_exon__huex_1_0_st_v2__lbl_gov__Level_3__quantile_normalization_gene__data")),
        "huex")
    expect_identical(ns$.searchPlatform(fn(
        "Merge_mirna__h_mirna_8x15k__unc_edu__Level_3__unc_DWD_Batch_adjusted__data")),
        "MIRNA_8x15k")
    expect_identical(ns$.searchPlatform(fn(
        "Merge_cna__hms_harvard_edu__hg19__Level_3__segmentation__seg")),
        "_hg19")
    ## RPPA has no recognised platform token
    expect_identical(ns$.searchPlatform(fn(
        "Merge_protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data")),
        character(0L))
    ## all three Agilent 07_x platforms collapse to the same name; this is
    ## what triggers .mergePlatforms() downstream
    expect_identical(
        ns$.searchPlatform(fn("Merge_transcriptome__agilentg4502a_07_1__unc_edu__Level_3__x__data")),
        ns$.searchPlatform(fn("Merge_transcriptome__agilentg4502a_07_2__unc_edu__Level_3__x__data"))
    )
})

test_that(".mergeNames combines platform and version tokens", {
    expect_identical(ns$.mergeNames(character(0L), character(0L)), character(0L))
    expect_identical(ns$.mergeNames("cgh", character(0L)), "cgh")
    expect_identical(ns$.mergeNames(character(0L), "hg19"), "_hg19")
    expect_identical(ns$.mergeNames("cgh", "415k"), "CGH_415k")
    expect_identical(ns$.mergeNames("cgh", c("415k", "hg19")), "CGH_415k_hg19")
    ## shortest platform name wins
    expect_identical(ns$.mergeNames(c("methyl450", "methyl"), character(0L)),
        "methyl")
})
