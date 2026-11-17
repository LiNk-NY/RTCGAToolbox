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

test_that(".stdIDs is deterministic on mixed dotted/dashed input", {
    ## BUG-01: the result used to depend on which elements sample() drew
    mixed <- c("TCGA-OR-A5J1-01A", "TCGA.OR.A5J2.01A")
    outs <- vapply(1:200, function(s)
        withr::with_seed(s, paste(ns$.stdIDs(mixed), collapse = " ")),
        character(1L))
    expect_identical(unique(outs), "TCGA-OR-A5J1-01A TCGA-OR-A5J2-01A")
    ## does not touch the RNG
    withr::with_seed(1, {
        before <- .Random.seed
        ns$.stdIDs(rep("TCGA.A", 20))
        expect_identical(.Random.seed, before)
    })
})

test_that(".stdIDs handles NA and empty input", {
    expect_identical(ns$.stdIDs(character(0L)), character(0L))
    expect_identical(ns$.stdIDs(c("TCGA.A1", NA)), c("TCGA.A1", NA))
})

test_that(".clinicalIDs gives standard patient barcodes", {
    expect_identical(ns$.clinicalIDs(c("tcga-or-a5j1", "tcga.or.a5j2",
        "TCGA-OR-A5J3")), c("TCGA-OR-A5J1", "TCGA-OR-A5J2", "TCGA-OR-A5J3"))
    ## only TCGA ids change
    expect_identical(ns$.clinicalIDs(c("p.1", "tcga.or.a5j1")),
        c("p.1", "TCGA-OR-A5J1"))
    expect_identical(ns$.clinicalIDs(character(0L)), character(0L))
})

test_that(".clinicalFrame transposes clinical tables", {
    raw <- data.frame(`Hybridization REF` = c("years_to_birth", "gender"),
        `tcga-or-a5j1` = c("58", "male"), `tcga-or-a5j2` = c("44", NA),
        check.names = FALSE)
    expect_identical(ns$.clinicalFrame(raw),
        data.frame(years_to_birth = c("58", "44"), gender = c("male", NA),
            row.names = c("TCGA-OR-A5J1", "TCGA-OR-A5J2")))
    ## a single patient or a single attribute keeps its shape
    one <- ns$.clinicalFrame(raw[, 1:2])
    expect_identical(dim(one), c(1L, 2L))
    expect_identical(rownames(one), "TCGA-OR-A5J1")
    expect_identical(ns$.clinicalFrame(raw[1L, ]),
        data.frame(years_to_birth = c("58", "44"),
            row.names = c("TCGA-OR-A5J1", "TCGA-OR-A5J2")))
    ## attribute names are kept as they are
    odd <- raw
    odd[[1L]] <- c("a b", "1x")
    expect_identical(names(ns$.clinicalFrame(odd)), c("a b", "1x"))
    ## ids that collide after standardisation are kept, with a warning
    dup <- raw
    names(dup)[3L] <- "tcga.or.a5j1"
    expect_warning(res <- ns$.clinicalFrame(dup), "TCGA-OR-A5J1")
    expect_identical(rownames(res), c("TCGA-OR-A5J1", "TCGA-OR-A5J1.1"))
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
    ## BUG-08: MSKCC CGH platforms get their own names; HMS names unchanged
    cgh <- function(x) ns$.searchPlatform(fn(paste0("Merge_cna__", x)))
    expect_identical(cgh(
        "cgh_1x1m_g4447a__mskcc_org__Level_3__segmentation_data_computation__seg"),
        "CGH_1x1m_g4447a_mskcc")
    expect_identical(cgh(
        "hg_cgh_244a__mskcc_org__Level_3__segmentation_data_computation__seg"),
        "CGH_hg_244a_mskcc")
    expect_identical(cgh("hg_cgh_244a__hms_harvard_edu__Level_3__segmentation__seg"),
        "CGH_hg_244a")
    expect_identical(cgh(
        "hg_cgh_415k_g4124a__hms_harvard_edu__Level_3__segmentation__seg"),
        "CGH_hg_415k_g4124a")
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
