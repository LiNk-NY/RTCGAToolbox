## getMutationRate(): fraction of mutated samples carrying each gene

data(accmini, package = "RTCGAToolbox", envir = environment())

maf_obj <- function() {
    maf <- read.delim(text = paste(fx_maf1, collapse = "\n"),
        colClasses = "character")
    extra <- maf[1L, ]
    extra[["Hugo_Symbol"]] <- "TP53"
    extra[["Tumor_Sample_Barcode"]] <- "TCGA-OR-A5J2-01A-11D-A29I-10"
    maf <- rbind(maf, extra)
    new("FirehoseData", Dataset = "ACC", Mutation = maf)
}

test_that("getMutationRate gives per-gene fraction over mutated samples", {
    res <- getMutationRate(maf_obj())
    ## sample 1: TP53, KRAS ; sample 2: TP53
    expect_identical(res[["Genes"]], c("TP53", "KRAS"))
    expect_identical(res[["MutationRatio"]], c(1, 0.5))
    expect_identical(rownames(res), c("TP53", "KRAS"))
})

test_that("getMutationRate on accmini (one sample)", {
    res <- getMutationRate(accmini)
    expect_identical(res[["Genes"]], unique(accmini@Mutation[["Hugo_Symbol"]]))
    expect_true(all(res[["MutationRatio"]] == 1))
})

test_that("getMutationRate input handling", {
    expect_error(getMutationRate(data.frame()), "must be 'FirehoseData'")
    ## BUG-13: no mutation data used to return NULL silently
    expect_warning(res <- getMutationRate(new("FirehoseData", Dataset = "X")),
        "No mutation data")
    expect_identical(res,
        data.frame(Genes = character(0L), MutationRatio = numeric(0L)))
    x <- maf_obj()
    x@Mutation[["Tumor_Sample_Barcode"]] <- NULL
    expect_error(getMutationRate(x), "Tumor_Sample_Barcode")
})

test_that("getMutationRate finds columns by name, not position", {
    ## BUG-13: genes and samples were read from columns 1 and 16
    x <- maf_obj()
    x@Mutation <- x@Mutation[, rev(names(x@Mutation))]
    expect_identical(getMutationRate(x), getMutationRate(maf_obj()))
    ## a gene mutated twice in one sample counts once
    y <- maf_obj()
    y@Mutation <- rbind(y@Mutation, y@Mutation[1L, ])
    expect_identical(getMutationRate(y), getMutationRate(maf_obj()))
})

test_that(".rbindMAFs is a row union, never a join", {
    a <- data.frame(Hugo_Symbol = c("TP53", "TP53"), Sample = c("s1", "s1"))
    b <- data.frame(Hugo_Symbol = "TP53", Sample = "s1", Extra = "x")
    ## identical shared columns: a join would collapse or multiply rows
    res <- ns$.rbindMAFs(list(a, b))
    expect_identical(nrow(res), 3L)
    expect_identical(names(res), c("Hugo_Symbol", "Sample", "Extra"))
    expect_identical(res[["Extra"]], c(NA, NA, "x"))
    ## same columns in a different order are matched by name
    c2 <- data.frame(Sample = "s2", Hugo_Symbol = "KRAS")
    res <- ns$.rbindMAFs(list(a, c2))
    expect_identical(res[["Hugo_Symbol"]], c("TP53", "TP53", "KRAS"))
    expect_identical(res[["Sample"]], c("s1", "s1", "s2"))
    ## identical columns: plain rbind
    expect_identical(ns$.rbindMAFs(list(a, a)), rbind(a, a))
})
