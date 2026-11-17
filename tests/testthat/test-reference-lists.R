## Static reference lists of cohorts and Firehose run dates

test_that("getFirehoseDatasets lists the 38 Firehose cohorts", {
    ds <- getFirehoseDatasets()
    expect_type(ds, "character")
    expect_length(ds, 38L)
    expect_false(anyDuplicated(ds) > 0L)
    expect_true(all(c("ACC", "BRCA", "COADREAD", "GBMLGG", "KIPAN", "STES",
        "FPPP", "UVM") %in% ds))
})

test_that("getFirehoseRunningDates lists stddata dates, newest first", {
    d <- getFirehoseRunningDates()
    expect_length(d, 62L)
    expect_true(all(grepl("^20[0-9]{6}$", d)))
    expect_identical(d[1L], "20160128")
    expect_identical(d[length(d)], "20111026")
    expect_false(is.unsorted(rev(d)))
    expect_identical(getFirehoseRunningDates(last = 2), c("20160128", "20151101"))
    expect_identical(getFirehoseRunningDates(last = 0), character(0L))
})

test_that("getFirehoseAnalyzeDates lists analyses dates, newest first", {
    d <- getFirehoseAnalyzeDates()
    expect_length(d, 35L)
    expect_true(all(grepl("^20[0-9]{6}$", d)))
    expect_identical(d[1L], "20160128")
    expect_identical(d[length(d)], "20101223")
    expect_false(is.unsorted(rev(d)))
    expect_identical(getFirehoseAnalyzeDates(last = 2), c("20160128", "20150821"))
})
