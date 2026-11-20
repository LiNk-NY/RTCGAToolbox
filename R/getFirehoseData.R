#' @importFrom data.table fread
## Read a Firehose expression matrix.
##
## 'headerRows' is the number of header lines: 2 for the 'Merge_*' files
## (sample barcodes, then the measurement type of each column), 1 for the
## mRNAseq_Preprocess files. Row 1 always holds the sample identifiers.
## With 'arrayData = FALSE' only the columns whose measurement type (row 2)
## starts with 'normMethod' are kept; otherwise every data column is.
## Rows whose identifier contains '?|' (no gene symbol) are dropped. Gene
## symbols are taken from 'SYMBOL|ENTREZ' identifiers; when a symbol occurs
## more than once, those rows keep the full identifier so no gene is lost.
.makeExprMat <- function(file, normMethod, dataType, mergeSize = 1000,
    arrayData = FALSE, headerRows = 2L)
{
    if (!arrayData && headerRows != 2L)
        stop("Column selection by 'normMethod' requires two header rows")
    hdr <- strsplit(readLines(file, n = headerRows), "\t", fixed = TRUE)
    if (length(hdr) < headerRows)
        stop("File has fewer than ", headerRows, " header lines: ", file)
    if (arrayData) {
        dataCols <- seq_along(hdr[[1L]])[-1L]
    } else {
        dataCols <- which(startsWith(hdr[[2L]], normMethod))
        dataCols <- dataCols[dataCols > 1L]
    }
    if (!length(dataCols))
        stop("No data columns found",
            if (!arrayData) paste0(" for '", normMethod, "'"),
            " in ", basename(file))
    samples <- hdr[[1L]][dataCols]

    message(dataType, " data will be imported! This may take a while!")
    message("Start: ", Sys.time())
    dat <- fread(file, sep = "\t", header = FALSE, skip = headerRows,
        select = c(1L, dataCols), colClasses = "character",
        data.table = FALSE)
    message("Done: ", Sys.time())

    ids <- dat[[1L]]
    keep <- !grepl("?|", ids, fixed = TRUE)
    mat <- .toNumericMatrix(
        as.matrix(dat[keep, -1L, drop = FALSE]), context = dataType
    )
    dimnames(mat) <- list(.geneRowNames(ids[keep]), samples)
    mat
}

## 'SYMBOL|ENTREZ' -> 'SYMBOL', keeping the full identifier for every
## member of a group of duplicated symbols.
.geneRowNames <- function(ids) {
    syms <- sub("\\|.*$", "", ids)
    dups <- syms %in% syms[duplicated(syms)]
    syms[dups] <- ids[dups]
    if (anyDuplicated(syms))
        warning("Duplicated row identifiers: ",
            paste(utils::head(unique(syms[duplicated(syms)]), 5L),
                collapse = ", "), call. = FALSE)
    syms
}

## Transpose a Firehose clin.merged.picked table (attributes in rows,
## patients in columns) to one row per patient, named by standard barcode.
.clinicalFrame <- function(raw.clin) {
    vals <- t(as.matrix(raw.clin[, -1L, drop = FALSE]))
    ids <- .clinicalIDs(rownames(vals))
    if (anyDuplicated(ids)) {
        warning("Duplicated patient identifiers in clinical data: ",
            paste(unique(ids[duplicated(ids)]), collapse = ", "),
            call. = FALSE)
        ids <- make.unique(ids)
    }
    df.clin <- data.frame(vals, stringsAsFactors = FALSE,
        check.names = FALSE)
    dimnames(df.clin) <- list(ids, raw.clin[[1L]])
    df.clin
}

## Stack per-sample MAF tables in file order. Files whose columns differ
## are bound by column name and missing columns are filled with NA (a row
## union, never a join on shared columns).
.rbindMAFs <- function(datalist) {
    cols <- lapply(datalist, names)
    if (length(unique(cols)) == 1L)
        return(do.call(rbind, datalist))
    as.data.frame(
        data.table::rbindlist(datalist, use.names = TRUE, fill = TRUE)
    )
}

## Name of the cached file for a Firehose archive. Files are named by
## archive, never by position in the link list: positions shift with
## 'fileSizeLimit' (skipped archives) and would map a cached file onto
## another platform's archive (BUG-08, BUG-30).
## 'Merge_' archives give platform and center
## ('humanmethylation450__jhu_usc_edu'); others the archive name between the
## cohort prefix and the Level suffix ('mRNAseq_Preprocess').
.archiveKey <- function(archive) {
    b <- basename(archive)
    key <- sub("^.*?[.]Merge_[[:alnum:]_]*?__(.+?)__Level_[0-9]__.*$", "\\1",
        b, perl = TRUE)
    if (identical(key, b))
        key <- sub("^gdac[.]broadinstitute[.]org_[^.]+[.](.+?)[.]Level_[0-9][.].*$",
            "\\1", b, perl = TRUE)
    if (identical(key, b))
        stop("Unrecognised Firehose archive name: ", archive)
    key
}

## Single-table RNA-seq slots (RNASeqGene, RNASeq2Gene) hold one platform.
## Some cohorts (e.g. COAD, UCEC) have Illumina GA and HiSeq archives, with
## different samples and, for aliquots run on both, different values
## (BUG-31). The first platform in 'platform' that the cohort has is used.
.selectRNAseqPlatform <- function(links, platform, dataset, type) {
    if (!length(links))
        return(links)
    plat <- sub("_.*", "", vapply(links, .archiveKey, character(1L),
        USE.NAMES = FALSE))
    use <- platform[platform %in% plat][1L]
    if (is.na(use))
        stop(dataset, " ", type, " data are available for ",
            paste0("'", unique(plat), "'", collapse = ", "),
            " only; see 'RNAseqPlatform'", call. = FALSE)
    sel <- links[plat == use]
    if (length(sel) != 1L)
        stop("<internal> ", length(sel), " '", use, "' archives for ",
            dataset, " ", type)
    if (length(unique(plat)) > 1L)
        message(dataset, " ", type, ": using '", use, "' data; also ",
            "available: ", paste0("'", setdiff(unique(plat), use), "'",
                collapse = ", "), " (see 'RNAseqPlatform')")
    sel
}

## 'archives': all archives of the data type; their keys must differ
.cacheName <- function(type, archive, archives) {
    keys <- vapply(archives, .archiveKey, character(1L), USE.NAMES = FALSE)
    if (anyDuplicated(keys))
        stop("Firehose archives share a cache name: ",
            paste(unique(keys[duplicated(keys)]), collapse = ", "))
    paste0("-", type, "-", .archiveKey(archive), ".txt")
}

## Firehose archive patterns for each data type, shared by getFirehoseData()
## and getLinks(). Each spec is passed to .getLinks(): 'filter' selects
## candidate file names, 'archive' the archive itself; with 'cohort = TRUE'
## the archive pattern must follow the cohort name (e.g. 'BRCA.Merge_...',
## which excludes 'BRCA-FFPE.Merge_...'). Specs are searched in order.
.fhSpec <- function(filter, archive, cohort = TRUE)
    list(filter = filter, archive = archive, cohort = cohort)

.FH_ARCHIVES <- list(
    clinical = list(.fhSpec(".Clinical_Pick_Tier1.Level_4", "*.tar[.]gz$",
        cohort = FALSE)),
    RNASeqGene = list(.fhSpec("Level_3__gene_expression__data.Level_3",
        "*.Merge_rnaseq__.*._rnaseq__.*.tar[.]gz$", cohort = FALSE)),
    RNASeq2Gene = list(.fhSpec("Level_3__RSEM_genes__data.Level_3",
        "*.Merge_rnaseqv2__.*._rnaseqv2__.*.tar[.]gz$", cohort = FALSE)),
    ## RNASeq2GeneNorm with RNAseq2Norm = "normalized_counts"
    RNASeq2GeneNorm = list(.fhSpec(
        "Level_3__RSEM_genes_normalized__data.Level_3",
        "*.Merge_rnaseqv2__.*._rnaseqv2__.*.tar[.]gz$", cohort = FALSE)),
    ## RNASeq2GeneNorm with any other RNAseq2Norm
    mRNAseq_Preprocess = list(.fhSpec("mRNAseq_Preprocess",
        ".*\\.Level_3\\..*\\.tar\\.gz$", cohort = FALSE)),
    miRNASeqGene = list(.fhSpec("Level_3__miR_gene_expression__data.Level_3",
        "[.]Merge_mirnaseq__.*.hiseq_mirnaseq__.*.tar[.]gz$")),
    CNASNP = list(.fhSpec("Level_3__segmented_scna_hg19__seg.Level_3",
        "[.]Merge_snp__.*.__Level_3__segmented_scna_hg19__seg.Level_3.*.tar[.]gz$")),
    CNVSNP = list(.fhSpec(
        "Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3",
        paste0("[.]Merge_snp__.*.__Level_3__segmented_scna_minus_germline",
            "_cnv_hg19__seg.Level_3.*.tar[.]gz$"))),
    CNASeq = list(.fhSpec("__Level_3__segmentation__seg.Level_3",
        "[.]Merge_cna__.*.dnaseq.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$")),
    ## HMS archives end in 'segmentation__seg', MSKCC archives (e.g.
    ## cgh_1x1m_g4447a) in 'segmentation_data_computation__seg'
    CNACGH = list(.fhSpec(
        "__Level_3__segmentation(_data_computation)?__seg[.]Level_3",
        paste0("[.]Merge_cna__.*cgh.*__Level_3__",
            "segmentation(_data_computation)?__seg[.]Level_3.*[.]tar[.]gz$"))),
    Methylation = list(.fhSpec(
        "__Level_3__within_bioassay_data_set_function__data.Level_3",
        paste0("[.]Merge_methylation__.*.methylation.*.__Level_3__within",
            "_bioassay_data_set_function__data.Level_3.*.tar[.]gz$"))),
    mRNAArray = list(
        .fhSpec("Merge_transcriptome__agilentg4502a_07",
            paste0("[.]Merge_transcriptome__agilentg4502a_.*.__Level_3__unc",
                "_lowess_normalization_gene_level__data.Level_3.*.tar[.]gz$")),
        .fhSpec("Merge_transcriptome__ht_hg_u133a",
            paste0("[.]Merge_transcriptome__ht_hg_u133a__.*.__Level_3__gene",
                "_rma__data.Level_3.*.tar[.]gz$")),
        .fhSpec("Merge_exon__huex_1_0_st_v2",
            paste0("[.]Merge_exon__huex_1_0_st_v2__.*.__Level_3__quantile",
                "_normalization_gene__data.Level_3.*.tar[.]gz$"))),
    miRNAArray = list(.fhSpec("h_mirna_8x15k",
        "[.]Merge_mirna__h_mirna_8x15k.*.data.Level_3.*.tar[.]gz$")),
    RPPAArray = list(.fhSpec("rppa_core",
        "[.]Merge_protein_exp.*.protein_normalization__data.Level_3.*.tar[.]gz$")),
    Mutation = list(.fhSpec("Mutation_Packager_Calls",
        "[.]Mutation_Packager_Calls[.]Level_3[.].*.tar[.]gz$")),
    ## 'dataset' is the analysis cohort, e.g. 'ACC-TP'
    GISTIC = list(.fhSpec("CopyNumber_Gistic2.Level_4",
        "[.]CopyNumber_Gistic2[.]Level_4.*.tar[.]gz$"))
)

## Archive names for a data type in a Firehose listing ('doc')
.archiveLinks <- function(type, dataset, doc) {
    specs <- .FH_ARCHIVES[[type]]
    if (is.null(specs))
        stop("<internal> no Firehose archive pattern for '", type, "'")
    links <- lapply(specs, function(sp)
        .getLinks(sp[["filter"]], sp[["archive"]],
            if (sp[["cohort"]]) dataset, doc))
    unlist(links)
}

#' @importFrom utils setTxtProgressBar
.getLinks <- function(keyWord1,keyWord2,datasetLink=NULL,doc) {
  plinks <- grep(keyWord1, doc, value = TRUE)
  search <-
    if (is.null(datasetLink)) keyWord2 else paste0("*.", datasetLink, keyWord2)
  grep(search,plinks,value=TRUE)
}

#' @importFrom utils txtProgressBar
#' @importFrom RJSONIO fromJSON
#' @importFrom RCurl getURL
.barcodeUUID <- function(object) {
  message("Converting barcodes to UUID")
  barcodes <- NULL
  if (dim(object@RNASeqGene)[1] > 0 & dim(object@RNASeqGene)[2] > 0) {
    barcodes <- c(barcodes,colnames(object@RNASeqGene))
  }
  if (dim(object@RNASeq2Gene)[1] > 0 & dim(object@RNASeq2Gene)[2] > 0) {
    barcodes <- c(barcodes,colnames(object@RNASeq2Gene))
  }
  if (dim(object@RNASeq2GeneNorm)[1] > 0 & dim(object@RNASeq2GeneNorm)[2] > 0) {
    barcodes <- c(barcodes,colnames(object@RNASeq2GeneNorm))
  }
  if (dim(object@miRNASeqGene)[1] > 0 & dim(object@miRNASeqGene)[2] > 0) {
    barcodes <- c(barcodes,colnames(object@miRNASeqGene))
  }
  if (dim(object@CNASNP)[1] > 0 & dim(object@CNASNP)[2] > 0) {
    barcodes <- c(barcodes,as.character(object@CNASNP[,1]))
  }
  if (dim(object@CNVSNP)[1] > 0 & dim(object@CNVSNP)[2] > 0) {
    barcodes <- c(barcodes,as.character(object@CNVSNP[,1]))
  }
  if (dim(object@CNASeq)[1] > 0 & dim(object@CNASeq)[2] > 0) {
    barcodes <- c(barcodes,as.character(object@CNASeq[,1]))
  }
  if (length(object@CNACGH) > 0 ) {
    barcodes <- c(barcodes,as.character(object@CNACGH[,1]))
  }
  if (length(object@Methylation) > 0 ) {
    for(i in 1:length(object@Methylation)) {
      barcodes <- c(barcodes,colnames(object@Methylation[[i]]@DataMatrix))
    }
  }
  if (length(object@mRNAArray) > 0 ) {
    for(i in 1:length(object@mRNAArray)) {
      barcodes <- c(barcodes,colnames(object@mRNAArray[[i]]@DataMatrix))
    }
  }
  if (length(object@miRNAArray) > 0 ) {
    for(i in 1:length(object@miRNAArray)) {
      barcodes <- c(barcodes,colnames(object@miRNAArray[[i]]@DataMatrix))
    }
  }
  if (length(object@RPPAArray) > 0 ) {
    for(i in 1:length(object@RPPAArray)) {
      barcodes <- c(barcodes,colnames(object@RPPAArray[[i]]@DataMatrix))
    }
  }
  if (length(object@GISTIC@Dataset) > 0) {
    barcodes <- c(barcodes,colnames(object@GISTIC@AllByGene)[-c(1:3)])
  }
  if (dim(object@Mutation)[1] > 0 & dim(object@Mutation)[2] > 0) {
    barcodes <- c(barcodes,unique(as.character(object@Mutation[,16])))
    barcodes <- c(barcodes,unique(as.character(object@Mutation[,17])))
  }
  barcodes <- unique(barcodes)
  barcodes <- toupper(barcodes)
  barcodes <- gsub(pattern = "\\.",replacement = "-",x = barcodes)

  pb <- txtProgressBar(min = 0, max = length(barcodes), style = 3)
  breakPoints <- seq(1,length(barcodes),20)
  for(i in breakPoints) {
    setTxtProgressBar(pb, i)
    urlToGo <- "https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/barcode/batch"
    endPoint <- (i+19)
    if (endPoint > length(barcodes)) {endPoint = length(barcodes)}
    searchBarcode <- paste(barcodes[i:endPoint], collapse=",")
    uuids = fromJSON(getURL(urlToGo, customrequest="POST",
                            httpheader=c("Content-Type: text/plain"),
                            postfields=searchBarcode))$uuidMapping
    if (is.character(uuids)) {
      if (exists("convertTable")) {
        convertTable <- rbind(convertTable,c(uuids[1],uuids[2]))
      } else {
        convertTable <- data.frame(Barcode=uuids[1],UUID=uuids[2])
      }
    } else {
      if (exists("convertTable")) {
        convertTable <- rbind(convertTable,
                              data.frame(matrix(unlist(uuids), nrow=length(uuids), byrow=TRUE),stringsAsFactors=FALSE))
      } else {
        convertTable <- data.frame(matrix(unlist(uuids), nrow=length(uuids), byrow=TRUE),stringsAsFactors=FALSE)
      }
    }

    if ((endPoint %% 400) == 0) {
      message("RTCGAToolbox will wait 185 seconds due to TCGA web service limitations")
      Sys.sleep(185)
    }

  }

  shortNames <- apply(convertTable,1,function(x) {
    as.character(paste(strsplit(x[1],split = "-")[[1]][1:3],collapse = "-"))
  })
  convertTable <- cbind(convertTable,shortNames)
  colnames(convertTable) <- c("Barcode","UUID","ShortName")
  setTxtProgressBar(pb, length(barcodes))
  return(convertTable)
}

.checkFileSize <- function(dataURL, fileSizeLimit) {
  fileSize <- httr::HEAD(dataURL)
  # this should convert file size from byte to MB
  fileSize <- as.numeric(fileSize$headers$`content-length`)/1e6

  if (fileSize > fileSizeLimit) {
    message(dataURL)
    message(paste("File Size: ~"), fileSize, "MB")
    message("File above won't be downloaded due to data size, RTCGAToolbox will skip this data!")
    FALSE
  } else {
    TRUE
  }
}

#' @importFrom utils download.file untar
.exportFiles <- function(fileLink, dataset, fileExt, searchName,subSearch=FALSE,
    exportName, manifest=FALSE, destdir=destdir, forceDownload=FALSE, runDate)
{

  if (!dir.exists(destdir))
      stop("Directory does not exist")

  tcgafile <- file.path(destdir, paste0(dataset, fileExt, sep=""))
  destfile <- file.path(destdir, paste0(runDate,"-",dataset,exportName))

  if (forceDownload || !file.exists(destfile)) {
    download.file(url=fileLink, destfile=tcgafile, method="auto",
        quiet = FALSE, mode = "wb")
    fileList <- untar(tcgafile, list = TRUE)
    fullList <- paste(destdir, fileList, sep = "/")
    if (any(nchar(fullList) > 259) && identical(.Platform$OS.type, "windows"))
        warning("File path too long, make 'destdir' shorter")
    if (!subSearch) {
      fileList = fileList[grepl(searchName,fileList)]
    } else {
      if (!manifest) {
        grepSearch = paste0("*.",dataset,searchName)
        fileList = fileList[grepl(grepSearch,fileList)]
      } else {
        fileList = fileList[!grepl("MANIFEST.txt",fileList)]
      }
    }

    untar(tcgafile, files = fileList, exdir = destdir)
    file.rename(from=file.path(destdir, fileList), to=destfile)
    file.remove(tcgafile)

    message(dirname(fileList))
    unlink(file.path(destdir, dirname(fileList)), recursive = TRUE)
  } else {
    message(sprintf('Using locally cached version of %s',destfile))
  }
  return(destfile)
}


.files_from_html_table <- function(url) {
    page_table <- rvest::html_table(rvest::read_html(url))[[1L]]
    rm.cols <- vapply(page_table, function(x) !all(is.na(x)), logical(1L))
    page_table <- page_table[, rm.cols]
    doc <- page_table[
      -which(page_table[["Name"]] %in% c("", "Parent Directory")), "Name"
    ]
    unlist(doc)
}

#' Get data from Firehose portal.
#'
#' `getFirehoseData` returns `FirehoseData` object that stores TCGA data.
#'
#' This is a main client function to download data from Firehose TCGA portal.
#'
#' @details To avoid unnecessary downloads, we use
#'   `tools::R_user_dir("RTCGAToolbox", "cache")` to set the default `destdir`
#'   parameter to the cached directory. To get the actual default directory,
#'   one can run `RTCGAToolbox:::.setCache()`.
#'
#' @param dataset A cohort disease code. TCGA cancer codes can be obtained via
#'   [getFirehoseDatasets]
#' @param runDate Standard data run dates. Date list can be accessible via
#'   [getFirehoseRunningDates]
#' @param gistic2Date Analysis run date for GISTIC obtained via
#'   [getFirehoseAnalyzeDates]
#' @param RNASeqGene Logical (default FALSE) RNAseq TPM data.
#' @param clinical Logical (default TRUE) clinical data.
#' @param RNASeq2Gene Logical (default FALSE) RNAseq v2 (RSEM processed) data;
#'   see `RNAseqNorm` argument.
#' @param RNASeq2GeneNorm Logical (default FALSE) RNAseq v2 (RSEM processed)
#'   data.
#' @param miRNASeqGene Logical (default FALSE) smallRNAseq data.
#' @param miRNASeqGeneType Character (default "read_count") Indicate which type
#'   of data should be pulled from the miRNASeqGene data. Must be one of
#'   "reads_per_million_miRNA_mapped", "read_count", or "cross-mapped".
#' @param CNASNP Logical (default FALSE) somatic copy number alterations data
#'   from SNP array, not filtered for germline CNVs
#'   (Firehose `segmented_scna_hg19`).
#' @param CNVSNP Logical (default FALSE) somatic copy number alterations data
#'   from SNP array with germline CNVs removed
#'   (Firehose `segmented_scna_minus_germline_cnv_hg19`).
#' @param CNASeq Logical (default FALSE) somatic copy number alterations data
#'   from sequencing.
#' @param CNACGH Logical (default FALSE) somatic copy number alterations data
#'   from CGH. Segmentations from both Harvard (HMS) and Memorial Sloan
#'   Kettering (MSKCC) arrays are downloaded, one `FirehoseCGHArray` per
#'   platform.
#' @param Methylation Logical (default FALSE) methylation data.
#' @param Mutation Logical (default FALSE) mutation data from sequencing.
#' @param mRNAArray Logical (default FALSE) mRNA expression data from
#'   microarray.
#' @param miRNAArray Logical (default FALSE) miRNA expression data from
#'   microarray.
#' @param RPPAArray Logical (default FALSE) RPPA data
#' @param RNAseqNorm RNAseq data normalization method. (Default raw_count)
#' @param RNAseq2Norm RNAseq v2 data normalization method. (Default
#'   normalized_count or one of RSEM_normalized_log2, raw_count,
#'   scaled_estimate)
#' @param RNAseqPlatform Character. Sequencing platform for `RNASeqGene` and
#'   `RNASeq2Gene`, in order of preference: the first one available for the
#'   cohort is used. Some cohorts (COAD, COADREAD, READ, STAD, STES, UCEC)
#'   have data from both "illuminaga" and "illuminahiseq", with different
#'   samples; the default, `c("illuminaga", "illuminahiseq")`, gives the
#'   Illumina GA data for those cohorts and the only platform available for
#'   all others. Use `"illuminahiseq"` for the HiSeq data. `RNASeq2GeneNorm`
#'   always returns all platforms.
#' @param GISTIC logical (default FALSE) processed copy number data
#' @param forceDownload Logical (default FALSE). Download the files again
#'   even when they are already in `destdir`. By default, files downloaded
#'   before are read from `destdir`. Either way, the same data are returned.
#' @param destdir Directory in which to store the resulting downloaded file.
#'   Defaults to a cache directory given by `RTCGAToolbox:::.setCache()`.
#' @param fileSizeLimit Files that are larger than set value (megabyte) won't be
#'   downloaded (Default: 500)
#' @param getUUIDs Deprecated. Logical (default FALSE). The barcode to UUID
#'   service it used was retired; `TRUE` gives a warning and is ignored. See
#'   [RTCGAToolbox-deprecated].
#' @param ... Additional arguments to pass down.
#'
#' @return A `FirehoseData` data object that stores data for selected data types.
#'
#' @seealso \link{getLinks}, \url{https://gdac.broadinstitute.org/}
#'
#' @importFrom utils write.table
#'
#' @examples
#' # Sample Dataset
#' data(accmini)
#' accmini
#' \dontrun{
#' BRCAdata <- getFirehoseData(dataset="BRCA",
#' runDate="20140416",gistic2Date="20140115",
#' RNASeqGene=TRUE,clinical=TRUE,mRNAArray=TRUE,Mutation=TRUE)
#' }
#' @export getFirehoseData
getFirehoseData <- function(dataset, runDate="20160128", gistic2Date="20160128",
    RNASeqGene=FALSE, RNASeq2Gene=FALSE, clinical=TRUE, miRNASeqGene=FALSE,
    miRNASeqGeneType =
      c("read_count", "reads_per_million_miRNA_mapped", "cross-mapped"),
    RNASeq2GeneNorm=FALSE, CNASNP=FALSE, CNVSNP=FALSE, CNASeq=FALSE,
    CNACGH=FALSE, Methylation=FALSE, Mutation=FALSE, mRNAArray=FALSE,
    miRNAArray=FALSE, RPPAArray=FALSE, GISTIC=FALSE, RNAseqNorm="raw_count",
    RNAseq2Norm = c(
      "normalized_counts", "RSEM_normalized_log2", "raw_counts", "scaled_estimate"
    ),
    RNAseqPlatform = c("illuminaga", "illuminahiseq"),
    forceDownload=FALSE, destdir=.setCache(),
    fileSizeLimit=500, getUUIDs=FALSE, ...) {
  ## the UUID web service was retired in 2016 (BUG-15)
  if (isTRUE(getUUIDs))
      .lifeCycleArg("getFirehoseData", "getUUIDs = TRUE")
  #check input parameters
  if (!is.character(RNAseqPlatform) || !length(RNAseqPlatform) ||
      !all(RNAseqPlatform %in% c("illuminaga", "illuminahiseq")))
      stop("'RNAseqPlatform' must be \"illuminaga\" and/or \"illuminahiseq\"")
  if (!is.character(dataset) || is.null(dataset) || !length(dataset) == 1 || nchar(dataset) < 2) {
      stop('Please set "dataset" parameter! You should specify one dataset name. Ex: dataset="BRCA"...')
      }
  runDatasets <- getFirehoseDatasets()
  if (!any(runDatasets==dataset)) {
      stop('Please use valid dataset name! "getFirehoseDatasets" function gives you the vector of valid dataset names!')
      }
  if (!is.null(runDate)) {
    if (!is.character(runDate) || !length(runDate) == 1 || !nchar(runDate) == 8) {
        stop('Please set "runDate" parameter! You should specify one Firehose run date. Ex: runDate="20140416"...')
        }
    runDateList <- getFirehoseRunningDates()
    if (!any(runDateList==runDate)) {
        stop('Please use valid run date! "getFirehoseRunningDates" function gives you the vector of valid dates!')
        }
  }

    if (GISTIC) {
      if (!S4Vectors::isSingleString(gistic2Date) ||
          !gistic2Date %in% getFirehoseAnalyzeDates())
          stop('Please set a valid "gistic2Date" parameter.')
    }
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)

    resultClass <- new("FirehoseData", Dataset = dataset)
    if (!is.null(runDate)) {
      resultClass@runDate <- runDate
    } else {
      resultClass@runDate <- NA_character_
    }
    if (GISTIC) {
      resultClass@gistic2Date <- gistic2Date
    } else {
      resultClass@gistic2Date <- NA_character_
    }

  if (!is.null(runDate))
  {
    ##build URL for getting file links
    fh_url <- "https://gdac.broadinstitute.org/runs/stddata__"
    fh_url <- paste0(fh_url,substr(runDate,1,4),"_",substr(runDate,5,6),"_",substr(runDate,7,8),"/data/")
    fh_url <- paste0(fh_url,dataset,"/",runDate,"/")
    doc <- .files_from_html_table(fh_url)

    #Download clinical data
    if (clinical)
    {
      #Search for links
      plinks <- .archiveLinks("clinical", dataset, doc)

      for(i in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,i),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,i),dataset,"-Clinical.tar.gz","*.clin.merged.picked.txt$",FALSE,
                       "-Clinical.txt",FALSE,destdir,forceDownload,runDate)
          raw.clin <- read.delim(export.file, colClasses = "character",
              check.names = FALSE)
          ## the 'Composite Element REF' row only holds the word 'value'
          raw.clin <- raw.clin[raw.clin[[1L]] != "Composite Element REF", ,
              drop = FALSE]
          resultClass@clinical <- .clinicalFrame(raw.clin)
          gc()
        }
      }
    }

    #Download RNAseq gene level data
    if (RNASeqGene)
    {
      #Search for links
      plinks <- trim(.archiveLinks("RNASeqGene", dataset, doc))
      ## one platform only (BUG-31); cached by archive, so neither
      ## forceDownload nor an earlier call changes which data are returned
      for(i in .selectRNAseqPlatform(plinks, RNAseqPlatform, dataset,
          "RNASeqGene"))
      {
        if (.checkFileSize(paste0(fh_url,i),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,i),dataset,"-RNAseqGene.tar.gz",
                       "[.]rnaseq__.*.__Level_3__gene_expression__data.data.txt$",
                       TRUE,.cacheName("RNAseqGene", i, plinks),FALSE,destdir,
                       forceDownload,runDate)
          #Get selected type only
          resultClass@RNASeqGene <- .makeExprMat(export.file,RNAseqNorm,"RNAseq",mergeSize=1000,arrayData=FALSE)
          gc()
        }
      }
    }

    #Download RNAseq2 gene level data (scaled estimates = TPM / 1e6 or raw_count, not rounded)
    if (RNASeq2Gene)
    {
      #Search for links
      plinks <- trim(.archiveLinks("RNASeq2Gene", dataset, doc))
      for(i in .selectRNAseqPlatform(plinks, RNAseqPlatform, dataset,
          "RNASeq2Gene"))
      {
        if (.checkFileSize(paste0(fh_url,i),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,i),dataset,"-RNAseq2Gene.tar.gz",
                                      "[.]rnaseqv2__.*.__Level_3__RSEM_genes__data.data.txt$",
                                      TRUE,.cacheName("RNAseq2Gene", i, plinks),
                                      FALSE,destdir,forceDownload,runDate)
          #Get selected type only
          resultClass@RNASeq2Gene <- .makeExprMat(export.file, RNAseqNorm, "RNAseq2",mergeSize=1000,arrayData=FALSE)
          gc()
        }
      }
    }

    #Download RNAseq2 gene level data
    if (RNASeq2GeneNorm)
    {
      #Search for links

      RNAseq2Norm <- match.arg(RNAseq2Norm)
      if (identical(RNAseq2Norm, "normalized_counts")) {
        searchName <- "[.]rnaseqv2__.*.__Level_3__RSEM_genes_normalized__data.data.txt$"
        plinks <- .archiveLinks("RNASeq2GeneNorm", dataset, doc)
      } else {
        searchName <- paste0(".*\\.mRNAseq_", RNAseq2Norm, "\\.txt$")
        plinks <- .archiveLinks("mRNAseq_Preprocess", dataset, doc)
      }
      plinks <- stats::setNames(
        vapply(plinks, function(x)
            .checkFileSize(paste0(fh_url, x), fileSizeLimit), logical(1L)
        ),
        plinks)
      plinks <- plinks[plinks]
      res <- lapply(seq_along(plinks), function(k) {
        i <- names(plinks)[k]
        export.file <- .exportFiles(
            fileLink = paste0(fh_url,i),
            dataset = dataset,
            fileExt = "-RNAseq2GeneNorm.tar.gz",
            searchName = searchName,
            subSearch = TRUE,
            exportName = .cacheName(paste0("RNAseq2GeneNorm-", RNAseq2Norm),
              i, names(plinks)),
            manifest = FALSE,
            destdir = destdir,
            forceDownload = forceDownload,
            runDate = runDate
        )
            new("FirehosemRNAArray",Filename=i,
                DataMatrix=.makeExprMat(export.file,
                    RNAseq2Norm,"RNAseq2",mergeSize=1000,arrayData=TRUE,
                    headerRows =
                        if (identical(RNAseq2Norm, "normalized_counts")) 2L
                        else 1L))
      })
      # res <- res[[which.max(vapply(res, ncol, integer(1L)))]]
      resultClass@RNASeq2GeneNorm <- res
    }

    #Download miRNAseq gene level data
    if (miRNASeqGene)
    {
      #Search for links
      plinks <- .archiveLinks("miRNASeqGene", dataset, doc)
      miRNAtype <- match.arg(miRNASeqGeneType)

      for(i in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,i),fileSizeLimit))
        {
          export.file <- .exportFiles(
              fileLink = paste0(fh_url,i),
              dataset = dataset,
              fileExt = "-miRNAseqGene.tar.gz",
              searchName = "[.]mirnaseq__.*.__Level_3__miR_gene_expression__data.data.txt$",
              subSearch = TRUE,
              exportName = paste0("-", miRNAtype, "-miRNAseqGene.txt"),
              manifest = FALSE,
              destdir = destdir,
              forceDownload = forceDownload,
              runDate = runDate
          )
          resultClass@miRNASeqGene <- .makeExprMat(export.file,miRNAtype,"miRNAseq",mergeSize=100,arrayData=FALSE)
          gc()
        }
      }
    }

    #Download CNA SNP data
    if (CNASNP)
    {
      #Search for links
      plinks <- .archiveLinks("CNASNP", dataset, doc)

      for(i in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,i),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,i),dataset,
                      "-CNASNPHg19.tar.gz",
                      "[.]snp__.*.__Level_3__segmented_scna_hg19__seg.seg.txt$",
                      TRUE,
                      "-CNASNPHg19.txt",FALSE,destdir,forceDownload,runDate)
          #Get selected type only
          tmpMat = fread(export.file,header=TRUE,colClasses=c("character","numeric","numeric",
                                                                                                "numeric","numeric","numeric"),data.table = FALSE)
          resultClass@CNASNP <- tmpMat
        }
      }
    }

    #Download CNV SNP data
    if (CNVSNP)
    {
      #Search for links
      plinks <- .archiveLinks("CNVSNP", dataset, doc)

      for(i in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,i),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,i),dataset,
                      "-CNVSNPHg19.tar.gz",
                      "[.]snp__.*.__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt$",
                      TRUE,
                      "-CNVSNPHg19.txt",FALSE,destdir,forceDownload,runDate)
          #Get selected type only
          tmpMat = fread(export.file,header=TRUE,colClasses=c("character","numeric","numeric",
                                                                                                "numeric","numeric","numeric"),data.table = FALSE)
          resultClass@CNVSNP <- tmpMat
        }
      }
    }

    #Download CNA DNAseq data
    if (CNASeq)
    {
      #Search for links
      plinks <- .archiveLinks("CNASeq", dataset, doc)

      for(i in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,i),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,i),dataset,
                      "-CNAseq.tar.gz",
                      "[.]cna__.*.__Level_3__segmentation__seg.seg.txt$",
                      TRUE,
                      "-CNAseq.txt",FALSE,destdir,forceDownload,runDate)
          #Get selected type only
          tmpMat = fread(export.file,
                         header=TRUE,colClasses=c("character","numeric","numeric","numeric","numeric","numeric"),
                         data.table = FALSE)
          #tmpMat = read.delim(paste0(runDate,"-",dataset,"-CNAseq.txt"),header=TRUE,colClasses=c("character","numeric","numeric",
          #                                                                                 "numeric","numeric"))
          resultClass@CNASeq <- tmpMat
        }
      }
    }

    #Download CNA CGH data
    if (CNACGH)
    {
      #Search for links
      ## HMS archives end in 'segmentation__seg', MSKCC archives (e.g.
      ## cgh_1x1m_g4447a) in 'segmentation_data_computation__seg'
      plinks <- .archiveLinks("CNACGH", dataset, doc)

      dataLists <- list()
      listCount = 1
      for(i in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,i),fileSizeLimit))
        {
          ## the cache file is named by platform, not by position, so a
          ## cached file is never matched to another platform's archive
          export.file <- .exportFiles(paste0(fh_url,i),dataset,
                      "-CNACGH.tar.gz",
                      "[.]cna__.*__Level_3__segmentation(_data_computation)?__seg[.]seg[.]txt$",
                      TRUE,
                      .cacheName("CNACGH", i, trim(plinks)),
                      FALSE,destdir,forceDownload,runDate)
          #Get selected type only
          tmpMat = fread(export.file,
                         header=TRUE,colClasses=c("character","numeric","numeric","numeric","numeric","numeric"),
                         data.table = FALSE)
          #tmpMat = read.delim(paste0(runDate,"-",dataset,"-CNACGH-",listCount,".txt",sep=""),header=TRUE,colClasses=c("character","numeric","numeric",
          #                                                                                               "numeric","numeric"))
          tmpReturn <- new("FirehoseCGHArray",Filename=i,DataMatrix=tmpMat)
          dataLists[[listCount]] <- tmpReturn
          listCount = listCount + 1
        }
      }
      resultClass@CNACGH <- dataLists
    }

    #Download methylation
    if (Methylation)
    {
      #Search for links
      plinks <- .archiveLinks("Methylation", dataset, doc)

      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,ii),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,ii),dataset,
                      "-Methylation.tar.gz",
                      "[.]methylation__.*.__Level_3__within_bioassay_data_set_function__data.data.txt$",
                      TRUE,
                      .cacheName("Methylation", ii, trim(plinks)),FALSE,destdir,forceDownload,runDate)

          #Get selected type only
          tmpCols <- read.delim(export.file,nrows=1,colClasses="character")
          BetaCols <- which(tmpCols[1L,] == "Beta_value")
          fixedCols <- match(c("Composite Element REF", "Gene_Symbol", "Chromosome", "Genomic_Coordinate"), tmpCols[1L, ])
          colOrder <- c(fixedCols, BetaCols)

          tmpMat <- fread(export.file,header=FALSE,colClasses = "character", select=colOrder, data.table = FALSE)
          fixedCols <- match(c("Composite Element REF", "Gene_Symbol", "Chromosome", "Genomic_Coordinate"), tmpMat[2L, ])
          BetaCols <- which(tmpMat[2L, ] == "Beta_value")
          colOrder <- c(fixedCols, BetaCols)
          tmpMat <- tmpMat[, colOrder]
          tmpMat[2L, ] <- gsub(" ", "", tmpMat[2L, ])
          firstFixed <- seq_along(fixedCols)
          BetaCols <- setdiff(seq_along(tmpMat), firstFixed)
          colnames(tmpMat) <- unlist(c(tmpMat[2L, firstFixed], tmpMat[1L, BetaCols]))

          tmpMat <- tmpMat[-c(1:2),]
          removeQM <- grepl("\\?\\|",tmpMat[,1])
          tmpMat <- tmpMat[!removeQM,]
          names1 <- tmpMat[,1]
          tmpMat <- tmpMat[,-1]
          rownames(tmpMat) <- names1
          tmpReturn <- new("FirehoseMethylationArray",Filename=ii,DataMatrix=tmpMat)
          dataLists[[listCount]] <- tmpReturn
          listCount = listCount + 1
        }
      }
      resultClass@Methylation <- dataLists
    }

    #Download mRNA array
    if (mRNAArray)
    {
      #Search for links
      plinks <- .archiveLinks("mRNAArray", dataset, doc)
      plinks = unique(plinks[plinks != ""])
      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,ii),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,ii),dataset,
                      "-mRNAArray.tar.gz",
                      "",
                      TRUE,
                      .cacheName("mRNAArray", ii, trim(plinks)),
                      TRUE,destdir,forceDownload,runDate)
          tmpReturn <- new("FirehosemRNAArray",Filename=ii,
                           DataMatrix=.makeExprMat(export.file,"","mRNAArray",1000,TRUE))
          dataLists[[listCount]] <- tmpReturn
          listCount = listCount + 1
        }
      }
      resultClass@mRNAArray <- dataLists
    }

    #Download miRNA array
    if (miRNAArray)
    {
      #Search for links
      plinks <- .archiveLinks("miRNAArray", dataset, doc)
      plinks = unique(plinks[plinks != ""])

      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,ii),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,ii),dataset,
                      "-miRNAArray.tar.gz",
                      "",
                      TRUE,
                      .cacheName("miRNAArray", ii, trim(plinks)),
                      TRUE,destdir,forceDownload,runDate)
          tmpReturn <- new("FirehosemRNAArray",Filename=ii,
                           DataMatrix=.makeExprMat(export.file,"","miRNAArray",100,TRUE))
          dataLists[[listCount]] <- tmpReturn
          listCount = listCount + 1
        }
      }
      resultClass@miRNAArray <- dataLists
    }

    #Download RPPA array
    if (RPPAArray)
    {
      #Search for links
      plinks <- .archiveLinks("RPPAArray", dataset, doc)
      plinks = unique(plinks[plinks != ""])

      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,ii),fileSizeLimit))
        {
          export.file <- .exportFiles(paste0(fh_url,ii),dataset,
                      "-RPPAArray.tar.gz",
                      "",
                      TRUE,
                      .cacheName("RPPAArray", ii, trim(plinks)),
                      TRUE,destdir,forceDownload,runDate)
          tmpReturn <- new("FirehosemRNAArray",Filename=ii,
                           DataMatrix=.makeExprMat(export.file,"","RPPAArray",100,TRUE))
          dataLists[[listCount]] <- tmpReturn
          listCount = listCount + 1
        }
      }
      resultClass@RPPAArray <- dataLists
    }

    #Download RPPA array
    if (Mutation)
    {
      #Search for links
      plinks <- .archiveLinks("Mutation", dataset, doc)
      plinks = unique(plinks[plinks != ""])

      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        if (.checkFileSize(paste0(fh_url,ii),fileSizeLimit))
        {
          allsampsfile <- paste0(destdir,"/",runDate,"-",dataset,"-Mutations-AllSamples.txt")
          if (forceDownload || !file.exists(allsampsfile)) {
            download_link = paste(fh_url,ii,sep="")
            tarfile <- file.path(destdir, paste0(dataset, "-Mutation.tar.gz"))
            download.file(url=download_link,destfile=tarfile,method="auto",quiet = FALSE, mode = "wb")
            fileList <- untar(tarfile, list = TRUE)
            grepSearch = "MANIFEST.txt"
            fileList = fileList[!grepl(grepSearch,fileList)]
            ###
            untar(tarfile, files = fileList, exdir = destdir)
            datalist <- lapply(file.path(destdir, fileList),FUN=function(files) {
              read.delim(files,header=TRUE,colClasses="character")
            })
            retMutations <- .rbindMAFs(datalist)
            delFolder <- file.path(destdir, dirname(fileList[[1L]]))
            unlink(delFolder, recursive = TRUE)
            file.remove(tarfile)
            write.table(retMutations,file=allsampsfile,sep="\t",row.names=FALSE,quote=FALSE)
          } else {
            retMutations <- fread(allsampsfile,header=TRUE,colClasses="character", data.table = FALSE)
          }
          resultClass@Mutation <- retMutations
        }
      }
    }
  }
  if (GISTIC) {
    tag <- switch(dataset, SKCM = "TM", LAML = "TB", "TP")
    dset <- paste0(dataset, "-", tag)
    ##build URL for getting file links
    fh_url <- "https://gdac.broadinstitute.org/runs/analyses__"
    fh_url <- paste(fh_url, substr(gistic2Date, 1, 4), "_",
        substr(gistic2Date, 5, 6), "_", substr(gistic2Date, 7, 8),
        "/data/", sep="")
    fh_url <- paste(fh_url, dset, "/", gistic2Date, "/", sep="")
    doc <- .files_from_html_table(fh_url)
    #Search for links
    plinks <- .archiveLinks("GISTIC", dset, doc)
    for (ii in trim(plinks))
    {
      if (.checkFileSize(paste0(fh_url,ii),fileSizeLimit)) {
        if (forceDownload || !file.exists(paste0(destdir,"/",gistic2Date,"-",dataset,"-all_thresholded.by_genes.txt"))) {
          download_link <- paste(fh_url, ii, sep = "")
          fileLoc <- file.path(destdir, paste0(dataset,"-Gistic2.tar.gz"))
          download.file(url = download_link, destfile = fileLoc,
            method = "auto", quiet = FALSE, mode = "wb")

          fileList <- untar(fileLoc, list = TRUE)
          grepSearch = "all_data_by_genes.txt"
          fileList = fileList[grepl(grepSearch,fileList)]
          untar(fileLoc, files = fileList, exdir = destdir)
          filepaths <- file.path(destdir, fileList)
          tmpCNAll = fread(filepaths, header = TRUE, colClasses = "character",
              data.table = FALSE)
          file.rename(from = filepaths,
            to = paste0(destdir, "/", gistic2Date, "-",
                dataset, "-all_data_by_genes.txt"))

          fileList <- untar(fileLoc, list = TRUE)
          grepSearch = "all_thresholded.by_genes.txt"
          fileList = fileList[grepl(grepSearch,fileList)]
          untar(fileLoc, files = fileList, exdir = destdir)
          filepaths <- file.path(destdir, fileList)
          tmpCNThreshhold = fread(filepaths, header = TRUE,
              colClasses = "character", data.table = FALSE)
          file.rename(from = filepaths,
            to = paste0(destdir, "/", gistic2Date, "-",
                dataset, "-all_thresholded.by_genes.txt"))

          fileList <- untar(fileLoc, list = TRUE)
          grepSearch = "all_lesions.conf_99.txt"
          fileList = fileList[grepl(grepSearch,fileList)]
          untar(fileLoc, files = fileList, exdir = destdir)
          filepaths <- file.path(destdir, fileList)
          tmpPeaks <- fread(filepaths, header = TRUE,
              colClasses = "character", data.table = FALSE)
          file.rename(from = filepaths,
            to = paste0(destdir, "/", gistic2Date, "-",
                dataset, "-all_lesions.conf_99.txt"))

          delFodler <- file.path(destdir, dirname(fileList))
          unlink(delFodler, recursive = TRUE)
          file.remove(fileLoc)
        } else {
          tmpCNThreshhold = fread(paste0(destdir,"/",gistic2Date,"-",dataset,"-all_thresholded.by_genes.txt"),header=TRUE,colClasses = "character", data.table = FALSE)
          tmpCNAll = fread(paste0(destdir,"/",gistic2Date,"-",dataset,"-all_data_by_genes.txt"),header=TRUE,colClasses="character", data.table = FALSE)
          tmpPeaks = fread(paste0(destdir, "/", gistic2Date, "-", dataset, "-all_lesions.conf_99.txt"), header = TRUE, colClasses = "character", data.table = FALSE)
        }
        tmpReturn <- new("FirehoseGISTIC",Dataset=dataset,AllByGene=data.frame(tmpCNAll),
            ThresholdedByGene=data.frame(tmpCNThreshhold), Peaks = data.frame(tmpPeaks))
        resultClass@GISTIC <- tmpReturn
      }
    }
  }

  return(resultClass)
}

.setCache <- function(
    directory = tools::R_user_dir("RTCGAToolbox", "cache"),
    verbose = TRUE,
    ask = interactive()
) {
    stopifnot(
        is.character(directory),
        S4Vectors::isSingleString(directory),
        !is.na(directory)
    )

    if (!dir.exists(directory)) {
        if (ask) {
            qtxt <- sprintf(
                "Create RTCGAToolbox cache at \n    %s? [y/n]: ",
                directory
            )
            answer <- .getAnswer(qtxt, allowed = c("y", "Y", "n", "N"))
            if ("n" == answer)
                stop("RTCGAToolbox directory not created")
        }
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }

    if (verbose)
        message("RTCGAToolbox cache directory set to:\n    ", directory)
    invisible(directory)
}

.getAnswer <- function(msg, allowed) {
    if (interactive()) {
        repeat {
            message(msg)
            answer <- readLines(n = 1)
            if (answer %in% allowed)
                break
        }
        tolower(answer)
    } else {
        "n"
    }
}
