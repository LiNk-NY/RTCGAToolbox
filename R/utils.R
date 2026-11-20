#' @importFrom GenomicRanges makeGRangesListFromDataFrame makeGRangesFromDataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#'     makeSummarizedExperimentFromDataFrame rowData<-
#' @importFrom DelayedArray DelayedArray
#' @importFrom RaggedExperiment RaggedExperiment
#' @importFrom S4Vectors SimpleList metadata metadata<- DataFrame mcols mcols<- splitAsList
#' @importFrom utils type.convert
#' @importFrom methods .hasSlot
## Helper functions for data extraction
.getDataMatrix <- function(object) {
    getElement(object, "DataMatrix")
}

.getFilenames <- function(object) {
    getElement(object, "Filename")
}

## Standardize barcode format
## When every identifier is a TCGA barcode, replace '.' separators (from
## make.names()) with '-' and upper-case; otherwise return the input as is.
.stdIDs <- function(sampleBarcode) {
    if (isTRUE(all(startsWith(sampleBarcode, "TCGA"))))
        sampleBarcode <- toupper(gsub(".", "-", sampleBarcode, fixed = TRUE))
    sampleBarcode
}

## Firehose clinical files name patients in lower case ('tcga-or-a5j1');
## objects made before 2.44.0 hold the read.delim() names ('tcga.or.a5j1').
## Both become standard barcodes ('TCGA-OR-A5J1'); other ids are unchanged.
.clinicalIDs <- function(ids) {
    isTCGA <- grepl("^tcga[-.]", ids, ignore.case = TRUE)
    ids[isTCGA] <- toupper(gsub(".", "-", ids[isTCGA], fixed = TRUE))
    ids
}

.standardizeBC <- function(x) {
    colnames(x) <- .stdIDs(colnames(x))
    return(x)
}

.getMethyl <- function(x) {
    headers <- names(x)
    annote <- x[, !grepl("TCGA", headers)]
    isNumRow <- all(grepl("^[0-9]*$",
        sample(rownames(x), size = 100L, replace = TRUE)))
    if (isNumRow) {
        geneSymbols <- annote[, grep("symbol", names(annote),
            ignore.case = TRUE, value = TRUE)]
        rNames <- geneSymbols
    } else { rNames <- rownames(x) }
    dm <- data.matrix(x[, grepl("TCGA", names(x))])
    mode(dm) <- "numeric"
    rownames(dm) <- rNames
    dm <- .standardizeBC(dm)
    dm <- DelayedArray::DelayedArray(dm)
    SummarizedExperiment::SummarizedExperiment(dm, rowData = annote)
}

.removeShell <- function(x, type) {
    if (startsWith(type, "GISTIC"))
        type <- "GISTIC"
    getElement(x, type)
}

.findBuild <- function(fname, type = c("UCSC", "NCBI")) {
    type <- match.arg(type)
    pattrn <- switch(type, UCSC = "[Hh][Gg][0-9]{2}",
        NCBI = "[Gg][Rr][Cc][Hh][0-9]{2}")
    if (!length(fname))
        return(NA_character_)
    bno <- rep(NA_character_, length(fname))
    pos <- regexpr(pattrn, fname)
    hit <- !is.na(pos) & pos != -1L
    bno[hit] <- regmatches(fname, pos)
    bno
}

.nameClean <- function(x) {
    x <- gsub("human|hum|agilent", "", x)
    x <- gsub("transcriptome", "tx", x, ignore.case = TRUE)
    x <- gsub("methylation", "methyl", x, ignore.case = TRUE)
    x
}

.mergeNames <- function(platform, version) {
    plat <- Filter(function(x) { !is.na(x) && length(x) }, tolower(platform))
    plat <- plat[which.min(nchar(plat))]
    if (!length(version))
        return(plat)
    ver <- tolower(version)
    logRM <- ver %in% plat
    version <- version[!logRM]
    relNames <- c(plat, version)
    if (length(plat) > 1L) {
        warning("Multiple platform names found, taking first one")
        plat <- plat[[1L]]
    }
    if (length(plat) && any(grepl(plat, tolower(version)))) {
        keep <- grepl("[0-9]{2}$", relNames, ignore.case = TRUE)
        result <- relNames[keep]
    } else if (length(version) > 1L) {
        result <- paste(toupper(plat), paste0(version, collapse = "_"),
                        sep = "_")
    } else if (length(version)) {
        result <- paste(toupper(plat), version, sep = "_")
    } else {
        result <- ""
    }
    return(result)
}

.searchPlatform <- function(x) {
    brokenUP <- unlist(strsplit(x, "_"))
    brokenUP <- Filter(function(y) nchar(y) != 0L, brokenUP)
    platNumExp <- "[0-9]k$|[0-9]a$|450$|27$|ht|hg|^1x1m$"
    namePlat <- unique(grep("cgh|mirna|meth|huex|^trans|illu", brokenUP,
        ignore.case = TRUE, value = TRUE))
    namePlat <- .nameClean(namePlat)
    vers <- grep(platNumExp, brokenUP, ignore.case = TRUE, value = TRUE)
    vers <- .nameClean(vers)
    result <- .mergeNames(namePlat, vers)
    ## MSKCC and HMS both ran the CGH 244A array; tell them apart
    if (nzchar(result) && "mskcc" %in% tolower(brokenUP))
        result <- paste(result, "mskcc", sep = "_")
    return(result)
}

#' @importFrom methods is
.unNestList <- function(x) {
    suppclasses <- all(vapply(x, function(y) {
        any(is(y, "FirehosemRNAArray"), is(y, "FirehoseCGHArray"),
            is(y, "FirehoseMethylationArray")) },
        logical(1L)))
    if (suppclasses) {
        ## DataFrame() makes duplicated row names unique ('g1', 'g1.1');
        ## record which inputs had unique row names before that happens
        uniqRows <- vapply(x, function(y) {
            rn <- rownames(.getDataMatrix(y))
            !is.null(rn) && !anyNA(rn) && !anyDuplicated(rn)
        }, logical(1L))
        x <- lapply(x, function(y) {
            fname <- .getFilenames(y)
            platform <- .searchPlatform(fname)
            build <- .findBuild(fname)
            y <- .getDataMatrix(y)
            ## Use DataFrame for metadata
            y <- DataFrame(y)
            metadata(y) <- list(filename = fname, build = build,
                platform = platform)
            return(y)
        })
        if (length(x) > 1L) {
            platNames <- vapply(x, function(y) {
                metadata(y)[["platform"]] }, character(1L))
            platNames <- gsub("human|hum|agilent", "", platNames)
            names(x) <- platNames
            if (anyDuplicated(platNames))
                x <- .mergePlatforms(x, uniqRows)
            names(x) <- make.unique(names(x), sep = "_")
        } else if (length(x) == 1L) { x <- x[[1L]] }
    }
    return(x)
}

.findCol <- function(x, colname) {
    if (!is.character(colname))
        stop("<internal> colname is not character")
    dataNames <- tolower(gsub("[^A-Za-z0-9]", "", names(x)))
    colname <- tolower(gsub("[^A-Za-z0-9]", "", colname))
    foundInData <- dataNames %in% colname
    if (sum(foundInData) > 1L)
        foundInData <- which.max(foundInData)
    if (!sum(foundInData))
        return(character(0L))
    names(x)[foundInData]
}

.hasInfo <- function(x, info = "NCBI_Build") {
    ## check "Hugo_Symbol" also possible
    buildInfo <- .findCol(x, info)
    as.logical(length(buildInfo))
}

.findUniqueAnnoCol <-
    function(df, annoCol = c("Hugo_Symbol", "Entrez_Gene_Id"))
{
    resname <- unlist(lapply(annoCol, function(anno) {
        annoName <- .findCol(df, anno)
        if (!length(annoName)) {
            character(0L)
        } else {
            annos <- df[[annoName]]
            if (identical(length(annos), length(unique(annos))))
                annoName
            else
                character(0L)
        }
    }))
    resname <- Filter(nchar, resname)
    if (length(resname) > 1L)
        resname[[1L]]
    else
        resname
}

.setAnnoRows <- function(df, rowAnnotation = c("Hugo_Symbol", "Entrez_Gene_Id"))
{
    annoName <- .findUniqueAnnoCol(df, rowAnnotation)
    if (length(annoName)) {
        annos <- df[[annoName]]
        rownames(df) <- annos
    }
    df
}

.standardstrand <- function(strandv) {
    strandv <- gsub("null", "*", strandv, ignore.case = TRUE)
    isnullna <- is.null(strandv) | is.na(strandv)
    strandv[isnullna] <- "*"
    strandv[strandv == 1] <- "+"
    strandv[strandv == -1] <- "-"
    strandv
}

.standardizeStrand <- function(x, strandcol) {
    x[[strandcol]] <- .standardstrand(x[[strandcol]])
    x
}

## Single genome build from a build column (e.g. NCBI_Build in MAFs).
## Equivalent labels ('37', 'GRCh37', 'hg19') count as one build; distinct
## builds are an error (coordinates on different assemblies must not be
## combined or relabelled). As in TCGAutils::uniformBuilds(), missing
## values are tolerated in fewer than 20% of rows.
.getBuild <- function(x, type = "NCBI_Build") {
    if (!.hasInfo(x, type))
        return(NA_character_)
    bvals <- as.character(x[[.findCol(x, type)]])
    known <- !(is.na(bvals) | bvals %in% c("", "NA"))
    ubuilds <- unique(bvals[known])
    if (!length(ubuilds))
        return(NA_character_)
    if (mean(!known) >= 0.2)
        stop("Genome build is missing in ", sum(!known), " of ",
            length(bvals), " rows", call. = FALSE)
    if (length(ubuilds) == 1L)
        return(ubuilds)
    ncbi <- unique(vapply(ubuilds, TCGAutils::correctBuild, character(1L),
        style = "NCBI", USE.NAMES = FALSE))
    if (anyNA(ncbi) || length(ncbi) > 1L)
        stop("Inconsistent genome builds found: ",
            paste(ubuilds, collapse = ", "), call. = FALSE)
    ncbi
}

.ansRangeNames <- function(x) {
    if (is(x, "list")) { return(list()) }
    granges_cols <- TCGAutils::findGRangesCols(names(x))
    fielders <- list(seqnames.field = "seqnames", start.field = "start",
        end.field = "end", strand.field = "strand")
    Fargs <- lapply(fielders, function(name) { names(x)[granges_cols[[name]]] })
    strd <- Fargs[["strand.field"]]
    allStrandNA <- if (!is.na(strd)) all(is.na(x[[strd]])) else TRUE
    Fargs[["ignore.strand"]] <- allStrandNA
    Filter(function(g) {!is.na(g)}, Fargs)
}

#' @importFrom stats na.omit
.findSampleCol <-
    function(x, sampcols = c("tumor_sample_barcode", "sample", "id"))
{
    sampcols <- tolower(sampcols)
    tsb <- na.omit(match(sampcols, tolower(names(x))))
    if (length(tsb)) {
        names(x)[tsb[[1L]]]
    } else {
        NA_character_
    }
}

.hasConsistentRanges <- function(object) {
    primary <- .findSampleCol(object)
    if (is.na(primary)) {
        return(FALSE)
    }
    ansRanges <- .ansRangeNames(object)
    # check if all ranges are of the same length
    grl <- do.call(.makeGRangesListFromDataFrame,
        c(list(df = object, split.field = primary), ansRanges))
    uniranges <- S4Vectors::isSingleInteger(unique(lengths(grl)))
    # then check if all ranges have same values
    if (!uniranges)
        return(FALSE)
    # repeated ranges within a sample cannot be aligned across samples
    if (length(grl) > 1L &&
            any(vapply(grl, function(gr) anyDuplicated(gr) > 0L, logical(1L))))
        return(FALSE)
    else
        all(vapply(grl[-1L], function(gr)
            S4Vectors::setequal(gr, grl[[1L]]), logical(1L))
        )
}

.hasRangeNames <- function(x) {
    if (is(x, "list")) { return(FALSE) }
    if (all(grepl("^TCGA", names(x)))) { return(FALSE) }
    if (!any(is.data.frame(x), is(x, "DataFrame"), is.matrix(x)))
        stop("(internal) 'x' must be rectangular")
    res <- is.na(TCGAutils::findGRangesCols(names(x)))
    if (any(res[c("seqnames", "start", "end")]))
        FALSE
    else
        !all(res)
}

.samplesAsCols <- function(x, sampleNames = character(0L)) {
    tcganames <- grepl("^TCGA", names(x), ignore.case = TRUE)
    sampleNames <- as.character(sampleNames)
    if (length(sampleNames))
        vapply(names(x), function(y) any(startsWith(y, sampleNames)),
            logical(1L))
    else
        tcganames
}

.hasExperimentData <- function(x, colnames = c("Hugo", "Entrez")) {
    anySamplesAsCols <- any(.samplesAsCols(x, colnames))
    sampcols <- na.omit(.findSampleCol(x))
    .hasRangeNames(x) || length(sampcols) || anySamplesAsCols
}

## Safe to assume equal number of ranges == equal ranges (?)
## 'colnames': prefixes of the sample columns (cBioPortalData passes the
## sample identifiers); by default, columns starting with 'TCGA'
.makeSummarizedExperimentFromDataFrame <-
    function(df, ..., colnames = character(0L))
{
    samplesAsCols <- .samplesAsCols(df, colnames)
    if (!any(samplesAsCols))
        stop("No sample columns found")
    metadat <- if (is(df, "DataFrame")) metadata(df) else list()
    rowData <- df[, !samplesAsCols, drop = FALSE]
    df <- data.matrix(df[, samplesAsCols, drop = FALSE])

    df <- .standardizeBC(df)
    args <- list(...)
    names.field <- args[["names.field"]]
    if (is.null(names.field) || !length(names.field)) {
        ## the annotation columns are in rowData, not in the matrix
        annoName <- .findUniqueAnnoCol(rowData)
        if (length(annoName))
            rownames(df) <- rowData[[annoName]]
    } else {
        rownames(df) <- rowData[[names.field]]
    }
    ## Use "" instead of missing due to changes in SummarizedExperiment
    ## constructor
    if (any(is.na(rownames(df))))
        rownames(df)[is.na(rownames(df))] <- ""
    if (length(rowData))
    object <- SummarizedExperiment(assays = SimpleList(df),
        rowData = rowData)
    else
    object <- makeSummarizedExperimentFromDataFrame(df)
    if (length(metadat))
        metadata(object) <- metadat
    return(object)
}

## Build an RSE from segment data in which every sample has the same set of
## ranges (see .hasConsistentRanges). Each sample's values are aligned to
## the ranges of the first sample by matching ranges, not by row position.
## Rows with a missing start or end are dropped from the values as well as
## from the ranges.
.makeRangedSummarizedExperimentFromDataFrame <-
    function(df, ..., seqinfo = NULL, starts.in.df.are.0based = FALSE) {
    args <- list(...)
    build <- args[["build"]]
    names.field <- args[["names.field"]]
    if (is.null(names.field) || !length(names.field)) {
        names.field <- NULL
        df <- .setAnnoRows(df)
    }
    metadat <- if (is(df, "DataFrame")) metadata(df) else list()
    split.field <- .findSampleCol(df)
    ansRanges <- .ansRangeNames(df)
    df <- df[!.missingRanges(df, ansRanges), , drop = FALSE]
    strictRanges <- Filter(function(x) !is.logical(x), ansRanges)
    RangeInfo <- c(strictRanges, list(split.field = split.field))
    ## the names column labels the ranges (BUG-29); it is not an assay
    numInfo <- df[, !(names(df) %in% c(RangeInfo, names.field)),
        drop = FALSE]
    rowRanges <- do.call(.makeGRangesListFromDataFrame,
        c(list(df = df[, c(unlist(RangeInfo), names.field), drop = FALSE],
            split.field = split.field, names.field = names.field),
            ansRanges)
    )
    rowIdx <- split(seq_len(nrow(df)), df[[split.field]])[names(rowRanges)]
    for (j in seq_along(rowIdx)) {
        ord <- match(rowRanges[[1L]], rowRanges[[j]])
        if (length(ord) != length(rowIdx[[j]]) || anyNA(ord) ||
                anyDuplicated(ord))
            stop("Ranges differ between samples; cannot build a",
                " RangedSummarizedExperiment")
        rowIdx[[j]] <- rowIdx[[j]][ord]
    }
    countList <- lapply(names(numInfo), function(nm) {
        vals <- numInfo[[nm]]
        do.call(cbind, lapply(rowIdx, function(i) vals[i]))
    })
    names(countList) <- names(numInfo)
    if (!is.null(build))
        Seqinfo::genome(rowRanges) <- build
    newSE <- SummarizedExperiment(assays = SimpleList(countList),
        rowRanges = rowRanges[[1L]])
    metadata(newSE) <- metadat
    return(newSE)
}

.removeNASeq <- function(x, colname) {
    nas <- is.na(x[[colname]])
    if (any(nas))
        message("Removing ", sum(nas), " rows where 'is.na(seqnames.field)'")
    x[!is.na(x[[colname]]), ]
}

.makeRaggedExperimentFromDataFrame <- function(df, ...) {
    args <- list(...)
    build <- args[["build"]]
    names.field <- args[["names.field"]]
    if (is.null(names.field) || !length(names.field))
        df <- .setAnnoRows(df)
    metadat <- if (is(df, "DataFrame")) { metadata(df) } else { list() }
    split.field <- args[["split.field"]]
    if (is.null(split.field))
        split.field <- .findSampleCol(df)

    ansRanges <- .ansRangeNames(df)
    rangeInfo <- c(ansRanges, list(split.field = split.field,
        names.field = names.field))

    df <- .removeNASeq(df, ansRanges[["seqnames.field"]])
    if (!is.null(ansRanges[["strand.field"]]) || length(ansRanges[["strand.field"]]))
        df <- .standardizeStrand(df, ansRanges[["strand.field"]])
    dropIdx <- .omitAdditionalIdx(df, ansRanges)
    if (length(dropIdx))
        df <- df[, -dropIdx]

    newGRL <- do.call(.makeGRangesListFromDataFrame,
        args = c(list(df = df, keep.extra.columns = TRUE), rangeInfo))
    if (!is.null(build))
        Seqinfo::genome(newGRL) <- build
    newRE <- RaggedExperiment::RaggedExperiment(newGRL)
    metadata(newRE) <- metadat
    return(newRE)
}

.rmNAse <- function(x, ansranges) {
    naRanges <- .missingRanges(x, ansranges)
    x[!naRanges, ]
}

.missingRanges <- function(x, ansranges) {
    startf <- ansranges[["start.field"]]
    endf <- ansranges[["end.field"]]
    is.na(x[[startf]]) | is.na(x[[endf]])
}

.makeGRangesFromDataFrame <- function(df, ...) {
    args <- list(...)
    build <- args[["build"]]
    metadat <- if (is(df, "DataFrame")) { metadata(df) } else { list() }
    ansRanges <- .ansRangeNames(df)
    df <- .rmNAse(df, ansRanges)
    dropIdx <- .omitAdditionalIdx(df, ansRanges)
    if (length(dropIdx))
        df <- df[, -dropIdx]
    df <- .setAnnoRows(df)
    newgr <- do.call(GenomicRanges::makeGRangesFromDataFrame,
        args = c(list(df = df, keep.extra.columns = TRUE), ansRanges))
    if (!is.null(build))
        Seqinfo::genome(newgr) <- build
    metadata(newgr) <- metadat
    return(newgr)
}

## replacing .makeGRangesFromDataFrame in GenomicRanges
.makeGRangesListFromDataFrame <-
    function(df, split.field = NULL, names.field = NULL, ...)
{
    splitIdx <- namesIdx <- integer()
    if (!is.null(split.field)) {
        if (!isSingleString(split.field))
            stop("'split.field' must be a single string")
        splitIdx <- which(names(df) %in% split.field)
        if (!length(splitIdx))
            stop("'split.field' is not in 'names(df)'")
        if (length(splitIdx) > 1L)
            stop("'split.field' matched more than one 'names(df)'")
        splitField <- df[[split.field]]
    }
    else splitField <- seq_len(nrow(df))
    if (!is.null(names.field)) {
        if (!isSingleString(names.field))
            stop("'names.field' must be a single string")
        namesIdx <- which(names(df) %in% names.field)
        if (!length(namesIdx))
            stop("'names.field' is not found in 'names(df)'")
        if (length(namesIdx) > 1L)
            stop("'names.field' matched more than one 'names(df)'")
        namesField <- df[[names.field]]
    }
    else namesField <- NULL
    if (length(c(splitIdx, namesIdx)))
        df <- df[, -c(splitIdx, namesIdx)]

    ansRanges <- .ansRangeNames(df)
    NAranges <- .missingRanges(df, ansRanges)
    df <- df[!NAranges, , drop = FALSE]
    splitField <- splitField[!NAranges]
    if (!is.null(namesField))
        namesField <- namesField[!NAranges]

    gr <- .makeGRangesFromDataFrame(df, ...)
    names(gr) <- namesField
    S4Vectors::split(gr, splitField)
}


.omitAdditionalIdx <- function(object, rangeNames) {
    rangeNames <- Filter(function(x) !is.logical(x), rangeNames)
    rangeIdx <- match(rangeNames, names(object))
    omitAdditional <- c("seqnames", "seqname", "chromosome", "chrom",
        "chromosome_name", "ranges", "seqlevels", "seqlengths", "seq_id",
        "iscircular", "start", "end", "strand", "width", "element", "chr")
    rmIdx <- which(tolower(names(object)) %in% omitAdditional)
    setdiff(rmIdx, rangeIdx)
}

## Combine DataFrames that share a platform name (e.g. Agilent 07_1, 07_2,
## 07_3 all map to 'TX_g4502a'). Columns are bound only when every
## element has the same set of unique, non-missing row names (rows are
## aligned by name, never by position) and no sample column is repeated.
## Otherwise the elements are returned unmerged, so no data are dropped.
## 'uniqRows' says whether each element's original row names were unique.
.canMergePlatforms <- function(dup, uniqRows = rep(TRUE, length(dup))) {
    rnames <- lapply(dup, rownames)
    ref <- rnames[[1L]]
    validRows <- vapply(rnames, function(rn) {
        !is.null(rn) && !anyNA(rn) && !anyDuplicated(rn) &&
            length(rn) == length(ref) && all(rn %in% ref)
    }, logical(1L))
    cnames <- unlist(lapply(dup, colnames), use.names = FALSE)
    all(uniqRows) && all(validRows) && !anyDuplicated(cnames)
}

.mergePlatforms <- function(x, uniqRows = rep(TRUE, length(x))) {
    pnames <- names(x)
    groups <- split(seq_along(x), factor(pnames, levels = unique(pnames)))
    res <- lapply(groups, function(idx) {
        dup <- unname(x[idx])
        if (length(dup) == 1L)
            return(dup)
        if (!.canMergePlatforms(dup, uniqRows[idx])) {
            fnames <- vapply(dup, function(y) {
                paste(metadata(y)[["filename"]], collapse = ", ")
            }, character(1L))
            message(
                "Platform '", pnames[idx[1L]], "' has ", length(dup),
                " files whose rows or samples do not match;",
                " returning them separately:\n  ",
                paste(fnames, collapse = "\n  ")
            )
            return(dup)
        }
        ref <- rownames(dup[[1L]])
        mets <- lapply(dup, metadata)
        meta <- split(
            unlist(mets, use.names = FALSE), names(unlist(unname(mets)))
        )
        dup <- lapply(dup, function(y) y[ref, , drop = FALSE])
        merged <- do.call(cbind, dup)
        metadata(merged) <- meta
        list(merged)
    })
    lens <- lengths(res)
    res <- unlist(unname(res), recursive = FALSE)
    names(res) <- rep(names(groups), lens)
    res
}

## Genome build from FILENAME
## RSE helper function from genome symbols to build (RNASeq, ExpSets)

.extractList <- function(object, type) {
    for (i in seq_along(object))
        object[[i]] <- biocExtract(object[[i]], type)
    return(object)
}

## Y/N flags (the miRNA-seq 'cross-mapped' columns) to a logical matrix:
## "Y" is TRUE, "N" is FALSE, "" and "NA" are NA; anything else is an error
.toLogicalMatrix <- function(m, context = "the") {
    vals <- as.vector(m)
    vals[vals %in% c("", "NA")] <- NA_character_
    bad <- !is.na(vals) & !vals %in% c("Y", "N")
    if (any(bad))
        stop("Values other than Y/N in ", context, " data: ",
            paste(utils::head(unique(vals[bad]), 5L), collapse = ", "),
            call. = FALSE)
    lgl <- vals == "Y"
    dim(lgl) <- dim(m)
    dimnames(lgl) <- dimnames(m)
    lgl
}

## Convert a character matrix read with colClasses = "character" to numeric
## without losing information. "", "NA" and "null" (used by Firehose array
## files) become NA; any other value that is not a number is an error
## rather than a silent NA. With 'integer = TRUE' the result is integer
## when every value is a whole number, otherwise double.
.toNumericMatrix <- function(m, integer = FALSE, context = "the") {
    if (is.numeric(m))
        return(m)
    vals <- as.vector(m)
    vals[vals %in% c("", "NA", "null")] <- NA_character_
    num <- suppressWarnings(as.numeric(vals))
    bad <- is.na(num) & !is.na(vals)
    if (any(bad))
        stop("Non-numeric values in ", context, " data: ",
            paste(utils::head(unique(vals[bad]), 5L), collapse = ", "),
            call. = FALSE)
    fin <- num[is.finite(num)]
    if (integer && all(fin == round(fin)) &&
            all(abs(fin) <= .Machine[["integer.max"]]))
        num <- as.integer(num)
    dim(num) <- dim(m)
    dimnames(num) <- dimnames(m)
    num
}
