## getMutationRate(): fraction of mutated samples carrying each gene

data(accmini, package = "RTCGAToolbox", envir = environment())

maf_obj <- function() {
    maf <- read.delim(text = paste(fx_maf1, collapse = "\n"),
        colClasses = "character")
    extra <- maf[1L, ]
    extra[["Hugo_Symbol"]] <- "TP53"
    extra[["Tumor_Sample_Barcode"]] <- "TCGA-OR-A5J2-01A-11D-A29I-10"
    maf <- rbind(maf, extra)
    ## getMutationRate() reads samples from column 16 (MAF layout); pad
    ## the fixture to the full MAF column order
    full <- as.data.frame(matrix("", nrow(maf), 17L))
    full[, 1L] <- maf[["Hugo_Symbol"]]
    full[, 16L] <- maf[["Tumor_Sample_Barcode"]]
    names(full)[c(1L, 16L)] <- c("Hugo_Symbol", "Tumor_Sample_Barcode")
    new("FirehoseData", Dataset = "ACC", Mutation = full)
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
    ## [BUG-13] no mutation data returns NULL silently
    expect_null(getMutationRate(new("FirehoseData", Dataset = "X")))
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
