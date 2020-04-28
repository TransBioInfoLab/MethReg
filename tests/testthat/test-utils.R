test_that("make_granges_from_names works", {
    library(GenomicRanges)
    gr <- make_granges_from_names("chr22:11-20")
    expect_true(is(gr, "GRanges"))
    expect_equal(width(gr), 10)
})
