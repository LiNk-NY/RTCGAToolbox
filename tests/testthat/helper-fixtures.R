## Test fixtures that mimic the GDAC Firehose file layout.
##
## Every archive is a .tar.gz with a single top-level directory named after
## the archive (as Firehose does) holding the data file(s) and a MANIFEST.txt.
## File contents follow the real Firehose formats (two header rows for
## expression/methylation matrices, 'gene|entrez' row ids for RNA-seq, etc.)
## so that the parsing code in getFirehoseData() runs exactly as it would on
## downloaded data.

ns <- asNamespace("RTCGAToolbox")

FH_RUNDATE <- "20160128"
FH_DATASET <- "ACC"
FH_PREFIX <- "gdac.broadinstitute.org_ACC."
FH_SUFFIX <- ".Level_3.2016012800.0.0"

## sample barcodes used across fixtures
BC <- c(
    S1 = "TCGA-OR-A5J1-01A-11R-A29S-07",
    S2 = "TCGA-OR-A5J2-01A-11R-A29S-07",
    S3 = "TCGA-OR-A5J3-01A-11R-A29S-07"
)

.tsv <- function(...) paste(c(...), collapse = "\t")

## --- file contents -------------------------------------------------------

fx_clinical <- c(
    .tsv("Hybridization REF", "tcga-or-a5j1", "tcga-or-a5j2"),
    .tsv("Composite Element REF", "value", "value"),
    .tsv("years_to_birth", "58", "44"),
    .tsv("vital_status", "1", "0"),
    .tsv("gender", "male", "female")
)

## RNASeq2Gene: three columns per sample
fx_rnaseq2 <- c(
    .tsv("Hybridization REF",
        rep(BC[["S1"]], 3), rep(BC[["S2"]], 3)),
    .tsv("gene_id",
        "raw_count", "scaled_estimate", "transcript_id",
        "raw_count", "scaled_estimate", "transcript_id"),
    .tsv("?|100130426", "0.00", "0", "uc011lsn.1", "3.00", "1e-07", "uc011lsn.1"),
    .tsv("A1BG|1", "197.00", "6.6e-06", "uc002qsd.2", "23.84", "1.2e-06", "uc002qsd.2"),
    .tsv("A2M|2", "5000.50", "1.5e-04", "uc001qvk.1", "10.00", "3.0e-07", "uc001qvk.1"),
    ## duplicated symbol with a different Entrez id (real case in TCGA)
    .tsv("SLC35E2|728661", "11.00", "4e-07", "uc001agw.2", "12.00", "5e-07", "uc001agw.2"),
    .tsv("SLC35E2|9906", "99.00", "9e-06", "uc001agx.2", "98.00", "8e-06", "uc001agx.2")
)

## RNASeq2GeneNorm (normalized_counts): one column per sample
fx_rnaseq2norm <- c(
    .tsv("Hybridization REF", BC[["S1"]], BC[["S2"]], BC[["S3"]]),
    .tsv("gene_id", "normalized_count", "normalized_count", "normalized_count"),
    .tsv("?|100130426", "0.0000", "0.5000", "1.0000"),
    .tsv("A1BG|1", "16.3305", "9.5987", "20.7377"),
    .tsv("A2M|2", "1000.1", "2000.2", "3000.3")
)

## mRNAseq_Preprocess files (used for RNAseq2Norm other than
## "normalized_counts"): ONE header row, then data. The first data row is a
## named gene on purpose (in real files it is usually a '?|' row).
PBC <- c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01")
fx_pre_log2 <- c(
    .tsv("gene", PBC),
    .tsv("A1BG|1", "1.7075689633722", "1.42304025336151"),
    .tsv("?|100134869", "1.9", "3.1"),
    .tsv("A2M|2", "2.5", "NA"),
    .tsv("TP53|7157", "10.25", "11.5")
)
fx_pre_raw <- c(
    .tsv("HYBRIDIZATION R", PBC),
    .tsv("A1BG|1", "6.8", "6.66"),
    .tsv("A2M|2", "100", "200"),
    .tsv("TP53|7157", "300", "400")
)
fx_pre_decoy <- c(.tsv("gene", PBC), .tsv("A1BG|1", "-999", "-999"))

## miRNASeqGene: three columns per sample
fx_mirna <- c(
    .tsv("Hybridization REF",
        rep(BC[["S1"]], 3), rep(BC[["S2"]], 3)),
    .tsv("miRNA_ID",
        "read_count", "reads_per_million_miRNA_mapped", "cross-mapped",
        "read_count", "reads_per_million_miRNA_mapped", "cross-mapped"),
    .tsv("hsa-let-7a-1", "76213", "15000.5", "N", "45441", "12000.25", "N"),
    .tsv("hsa-let-7a-2", "151321", "30000.75", "Y", "90039", "24000.5", "N")
)

## segmented copy number (CNASNP / CNVSNP / CNASeq / CNACGH)
fx_seg <- c(
    .tsv("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean"),
    .tsv("TCGA-OR-A5J1-01A-11D-A29I-01", "1", "61735", "12834253", "6520", "0.0339"),
    .tsv("TCGA-OR-A5J1-01A-11D-A29I-01", "2", "100", "50000", "48", "-0.3543"),
    .tsv("TCGA-OR-A5J2-01A-11D-A29I-01", "1", "61735", "20000000", "9000", "0.5"),
    .tsv("TCGA-OR-A5J2-01A-11D-A29I-01", "23", "1000", "2000", "8", "-1.25")
)

## methylation: four columns per sample
fx_methyl <- c(
    .tsv("Hybridization REF",
        rep("TCGA-OR-A5J1-01A-11D-A29J-05", 4),
        rep("TCGA-OR-A5J2-01A-11D-A29J-05", 4)),
    .tsv("Composite Element REF",
        "Beta_value", "Gene_Symbol", "Chromosome", "Genomic_Coordinate",
        "Beta_value", "Gene_Symbol", "Chromosome", "Genomic_Coordinate"),
    .tsv("cg00000029", "0.119877013723081", "RBL2", "16", "53468112",
        "0.107120474727399", "RBL2", "16", "53468112"),
    .tsv("cg00000108", "NA", "C3orf35", "3", "37459206",
        "0.5", "C3orf35", "3", "37459206"),
    .tsv("cg00000109", "0.9", "FNDC3B", "3", "171916037",
        "NA", "FNDC3B", "3", "171916037")
)

## array data (mRNA / miRNA / RPPA): one column per sample
fx_mrnaarray <- c(
    .tsv("Hybridization REF", BC[["S1"]], BC[["S2"]]),
    .tsv("Composite Element REF",
        "log2 lowess normalized (cy5/cy3) collapsed by gene symbol",
        "log2 lowess normalized (cy5/cy3) collapsed by gene symbol"),
    .tsv("A1BG", "-0.5", "0.25"),
    .tsv("A2M", "null", "1.5")
)
fx_mirnaarray <- c(
    .tsv("Hybridization REF", BC[["S1"]], BC[["S2"]]),
    .tsv("Composite Element REF", "unc_DWD_Batch_adjusted", "unc_DWD_Batch_adjusted"),
    .tsv("hsa-let-7a", "10.5", "11.25"),
    .tsv("hsa-let-7b", "8.125", "-2")
)
fx_rppa <- c(
    .tsv("Sample REF", "TCGA-OR-A5J2-01A-21-A39K-20", "TCGA-OR-A5J3-01A-21-A39K-20"),
    .tsv("Composite Element REF", "Protein Expression", "Protein Expression"),
    .tsv("14-3-3_beta", "0.2233428", "-0.14206456"),
    .tsv("4E-BP1_pS65", "-1.5", "2.75")
)

## per-sample MAF files; the second one has an extra column so the
## merge(..., all = TRUE) branch is exercised
fx_maf1 <- c(
    .tsv("Hugo_Symbol", "Entrez_Gene_Id", "NCBI_Build", "Chromosome",
        "Start_position", "End_position", "Strand", "Variant_Classification",
        "Tumor_Sample_Barcode"),
    .tsv("TP53", "7157", "37", "17", "7577120", "7577120", "+",
        "Missense_Mutation", "TCGA-OR-A5J1-01A-11D-A29I-10"),
    .tsv("KRAS", "3845", "37", "12", "25398284", "25398284", "+",
        "Missense_Mutation", "TCGA-OR-A5J1-01A-11D-A29I-10")
)
fx_maf2 <- c(
    .tsv("Hugo_Symbol", "Entrez_Gene_Id", "NCBI_Build", "Chromosome",
        "Start_position", "End_position", "Strand", "Variant_Classification",
        "Tumor_Sample_Barcode", "Extra_Col"),
    .tsv("TP53", "7157", "37", "17", "7578406", "7578406", "+",
        "Nonsense_Mutation", "TCGA-OR-A5J2-01A-11D-A29I-10", "x")
)

## GISTIC2 level 4 files
GBC <- c("TCGA-OR-A5J1-01A-11D-A29H-01", "TCGA-OR-A5J2-01A-11D-A29H-01")
fx_gistic_all <- c(
    .tsv("Gene Symbol", "Locus ID", "Cytoband", GBC),
    .tsv("ACAP3", "116983", "1p36.33", "0.030", "-0.070"),
    .tsv("AGRN", "375790", "1p36.33", "0.125", "-1.250"),
    .tsv("TP53", "7157", "17p13.1", "-0.500", "0.900")
)
fx_gistic_thr <- c(
    .tsv("Gene Symbol", "Locus ID", "Cytoband", GBC),
    .tsv("ACAP3", "116983", "1p36.33", "0", "0"),
    .tsv("AGRN", "375790", "1p36.33", "1", "-1"),
    .tsv("TP53", "7157", "17p13.1", "-2", "2")
)
## note the trailing tab on every line: real Firehose lesions files have it
fx_gistic_peaks <- paste0(c(
    .tsv("Unique Name", "Descriptor", "Wide Peak Limits", "Peak Limits",
        "Region Limits", "q values",
        "Residual q values after removing segments shared with higher peaks",
        "Broad or Focal", "Amplitude Threshold", GBC),
    .tsv("Amplification Peak  1", "1q22",
        "chr1:120523956-155582840(probes 67527:70182)",
        "chr1:155024310-155083455(probes 69961:69975)",
        "chr1:155024310-155095653(probes 69961:69978)",
        "0.11282", "0.11282", "Broad", "0: t<0.1; 1: 0.1<t< 0.9; 2: t>0.9",
        "2", "0"),
    .tsv("Amplification Peak  2", "Xq28",
        "chrX:1-2496016(probes 1:20)",
        "chrX:1778802-1816689(probes 5:9)",
        "chrX:1-2709194(probes 1:25)",
        "0.038099", "0.038099", "Focal", "0: t<0.1; 1: 0.1<t< 0.9; 2: t>0.9",
        "1", "2"),
    .tsv("Deletion Peak  1", "3p14.2",
        "chr3:59000000-61000000(probes 100:200)",
        "chr3:60000000-60500000(probes 120:150)",
        "chr3:58000000-62000000(probes 90:210)",
        "1e-08", "1e-08", "Focal", "0: t>-0.1; 1: -0.1>t> -1.3; 2: t<-1.3",
        "0", "2"),
    .tsv("Amplification Peak  1 - CN values", "1q22",
        "chr1:120523956-155582840(probes 67527:70182)",
        "chr1:155024310-155083455(probes 69961:69975)",
        "chr1:155024310-155095653(probes 69961:69978)",
        "0.11282", "0.11282", "Broad", "Actual Copy Change Given",
        "0.95", "0.01")
), "\t")

## --- archive helpers -----------------------------------------------------

## Write `files` (named list: filename -> lines) into a Firehose-style
## tarball at `outdir/<archive>.tar.gz`; returns the archive filename.
make_fh_tarball <- function(outdir, archive, files, manifest = TRUE) {
    stage <- tempfile("stage")
    top <- file.path(stage, archive)
    dir.create(top, recursive = TRUE)
    for (f in names(files))
        writeLines(files[[f]], file.path(top, f))
    if (manifest)
        writeLines(paste("d41d8cd98f00b204e9800998ecf8427e", names(files)),
            file.path(top, "MANIFEST.txt"))
    tarname <- paste0(archive, ".tar.gz")
    ## internal tar warns about >100 byte paths; Firehose names are long
    withr::with_dir(stage, suppressWarnings(
        utils::tar(file.path(outdir, tarname), files = archive,
            compression = "gzip", tar = "internal")
    ))
    unlink(stage, recursive = TRUE)
    tarname
}

## Build the complete fake 'stddata' and 'analyses' mirrors.
build_fake_firehose <- function(root = tempfile("firehose")) {
    std <- file.path(root, "stddata")
    ana <- file.path(root, "analyses")
    dir.create(std, recursive = TRUE)
    dir.create(ana, recursive = TRUE)
    P <- FH_PREFIX
    S <- FH_SUFFIX
    std_files <- c(
        make_fh_tarball(std, "gdac.broadinstitute.org_ACC.Clinical_Pick_Tier1.Level_4.2016012800.0.0",
            list("All_CDEs.txt" = c("bcr_patient_barcode\ttcga-or-a5j1"),
                "ACC.clin.merged.picked.txt" = fx_clinical)),
        make_fh_tarball(std, paste0(P, "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data", S),
            list("ACC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt" = fx_rnaseq2)),
        make_fh_tarball(std, paste0(P, "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data", S),
            list("ACC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt" = fx_rnaseq2norm)),
        make_fh_tarball(std, paste0(P, "Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data", S),
            list("ACC.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt" = fx_mirna)),
        make_fh_tarball(std, paste0(P, "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg", S),
            list("ACC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt" = fx_seg)),
        make_fh_tarball(std, paste0(P, "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg", S),
            list("ACC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt" = fx_seg)),
        make_fh_tarball(std, paste0(P, "Merge_cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg", S),
            list("ACC.cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg.seg.txt" = fx_seg)),
        make_fh_tarball(std, paste0(P, "Merge_cna__hg_cgh_415k_g4124a__hms_harvard_edu__Level_3__segmentation__seg", S),
            list("ACC.cna__hg_cgh_415k_g4124a__hms_harvard_edu__Level_3__segmentation__seg.seg.txt" = fx_seg)),
        make_fh_tarball(std, paste0(P, "Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data", S),
            list("ACC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt" = fx_methyl)),
        make_fh_tarball(std, paste0(P, "Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data", S),
            list("ACC.transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.txt" = fx_mrnaarray)),
        make_fh_tarball(std, paste0(P, "Merge_mirna__h_mirna_8x15k__unc_edu__Level_3__unc_DWD_Batch_adjusted__data", S),
            list("ACC.mirna__h_mirna_8x15k__unc_edu__Level_3__unc_DWD_Batch_adjusted__data.data.txt" = fx_mirnaarray)),
        make_fh_tarball(std, paste0(P, "Merge_protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data", S),
            list("ACC.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt" = fx_rppa)),
        make_fh_tarball(std, "gdac.broadinstitute.org_ACC.mRNAseq_Preprocess.Level_3.2016012800.0.0",
            list(
                "ACC.uncv2.mRNAseq_RSEM_normalized_log2.txt" = fx_pre_log2,
                "ACC.uncv2.mRNAseq_raw_counts.txt" = fx_pre_raw,
                "ACC.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt" = fx_pre_decoy,
                "ACC.uncv2.mRNAseq_RSEM_Z_Score.txt" = fx_pre_decoy,
                "ACC.uncv2.mRNAseq_geneid_transcriptid_mapping.txt" =
                    c(.tsv("HYBRIDIZATION R", PBC[1L]), .tsv("A1BG|1", "uc002qsd.2"))
            )),
        make_fh_tarball(std, "gdac.broadinstitute.org_ACC.Mutation_Packager_Calls.Level_3.2016012800.0.0",
            list("TCGA-OR-A5J1-01.maf.txt" = fx_maf1, "TCGA-OR-A5J2-01.maf.txt" = fx_maf2))
    )
    ## Archives that appear in real listings but must never be selected.
    ## They have no tarball on the fake mirror, so selecting one makes the
    ## mocked download fail loudly.
    decoys <- c(
        sub("[.]Level_[34][.]", ".aux.", std_files),
        sub("[.]Level_[34][.]", ".mage-tab.", std_files),
        paste0(P, "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg18__seg", S, ".tar.gz"),
        paste0(P, "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg", S, ".tar.gz"),
        paste0(P, "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms__data", S, ".tar.gz"),
        paste0(P, "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data", S, ".tar.gz"),
        paste0(P, "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__exon_quantification__data", S, ".tar.gz"),
        paste0(P, "Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data", S, ".tar.gz"),
        "gdac.broadinstitute.org_ACC.Mutation_Packager_Coverage.Level_3.2016012800.0.0.tar.gz",
        "gdac.broadinstitute.org_ACC.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0.tar.gz",
        "gdac.broadinstitute.org_ACC.Methylation_Preprocess.Level_3.2016012800.0.0.tar.gz",
        "gdac.broadinstitute.org_ACC.miRseq_Preprocess.Level_3.2016012800.0.0.tar.gz",
        "gdac.broadinstitute.org_ACC.RPPA_AnnotateWithGene.Level_3.2016012800.0.0.tar.gz"
    )
    std_listing <- c(std_files, paste0(std_files, ".md5"), decoys)
    ana_files <- make_fh_tarball(ana,
        "gdac.broadinstitute.org_ACC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0",
        list(
            "all_data_by_genes.txt" = fx_gistic_all,
            "all_thresholded.by_genes.txt" = fx_gistic_thr,
            "all_lesions.conf_99.txt" = fx_gistic_peaks
        )
    )
    ana_listing <- c(ana_files, paste0(ana_files, ".md5"),
        sub("[.]Level_4[.]", ".aux.", ana_files),
        "gdac.broadinstitute.org_ACC-TP.CopyNumberLowPass_Gistic2.Level_4.2016012800.0.0.tar.gz")
    list(root = root, std = std, ana = ana,
        std_listing = std_listing, ana_listing = ana_listing)
}

## Redirect all network access in getFirehoseData()/getLinks() to the fake
## mirror. Returns an environment that counts downloads.
local_fake_firehose <- function(fh, env = parent.frame()) {
    counter <- new.env()
    counter$downloads <- character(0L)
    counter$urls_listed <- character(0L)
    testthat::local_mocked_bindings(
        .files_from_html_table = function(url) {
            counter$urls_listed <- c(counter$urls_listed, url)
            if (grepl("analyses__", url, fixed = TRUE))
                fh$ana_listing
            else
                fh$std_listing
        },
        .checkFileSize = function(dataURL, fileSizeLimit) TRUE,
        download.file = function(url, destfile, ...) {
            counter$downloads <- c(counter$downloads, url)
            src <- if (grepl("analyses__", url, fixed = TRUE)) fh$ana else fh$std
            stopifnot(file.copy(file.path(src, basename(url)), destfile,
                overwrite = TRUE))
            invisible(0L)
        },
        .package = "RTCGAToolbox",
        .env = env
    )
    counter
}

## run quietly, keep result
quiet <- function(expr) suppressMessages(expr)
