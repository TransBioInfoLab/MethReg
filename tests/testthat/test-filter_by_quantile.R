test_that("filter_dnam_by_quant_diff remove correctly regions without variance", {
    dnam <- t(matrix(c(rep(0,15), rep(1,5)),ncol = 2))
    rownames(dnam) <- c("to_be_removed","to_keep")
    colnames(dnam) <- paste0("Samples",1:10)
    filtered <- filter_dnam_by_quant_diff(dnam, min.IQR.threshold = 0.2)
    expect_true(is(filtered, "matrix"))
    expect_true("to_keep" %in% rownames(filtered))
    expect_false("to_be_removed" %in% rownames(filtered))
})


test_that("filter_dnam_by_quant_diff handles NAs", {
    dnam <- t(matrix(c(NA,rep(0,13),NA, rep(1,4),NA,rep(NA,10)),ncol = 3))
    rownames(dnam) <- c("to_be_removed","to_keep","all_na")
    colnames(dnam) <- paste0("Samples",1:10)
    filtered <- filter_dnam_by_quant_diff(dnam, min.IQR.threshold = 0.2)
    expect_true(is(filtered, "matrix"))
    expect_true("to_keep" %in% rownames(filtered))
    expect_false("to_be_removed" %in% rownames(filtered))
    expect_false("all_na" %in% rownames(filtered))
})



test_that("filter_dnam_by_quant_diff handles NAs for a SE", {
    dnam <- t(matrix(c(NA,rep(0,13),NA, rep(1,4),NA,rep(NA,10)),ncol = 3))
    nrows <- 3; ncols <- 10
    rowRanges <- GRanges(
        rep(c("chr1"), c(3)),
        IRanges(floor(runif(3, 1e5, 1e6)), width = 100),
        strand = sample(c("+", "-"), 3, TRUE)
    )
    rownames(dnam) <-  c("to_be_removed","to_keep","all_na")

    colData <- DataFrame(
        Treatment = rep(c("ChIP"), 10),
        row.names = LETTERS[1:10]
    )

    rse <- SummarizedExperiment(
        assays = SimpleList(counts = dnam),
        rowRanges = rowRanges,
        colData = colData
    )

    filtered <- filter_dnam_by_quant_diff(rse, diff.mean.th = 0.2)
    expect_true(is(filtered, "SummarizedExperiment"))
    expect_true("to_keep" %in% rownames(filtered))
    expect_false("to_be_removed" %in% rownames(filtered))
    expect_false("all_na" %in% rownames(filtered))
})





test_that("filter_exp_by_quant_mean_FC remove correctly genes without variance", {
    exp <- t(matrix(c(rep(0,15), rep(1,5)),ncol = 2))
    rownames(exp) <- c("to_be_removed","to_keep")
    colnames(exp) <- paste0("Samples",1:10)
    filtered <- filter_exp_by_quant_mean_FC(exp, fold.change = 1.5)
    expect_true(is(filtered, "matrix"))
    expect_true("to_keep" %in% rownames(filtered))
    expect_false("to_be_removed" %in% rownames(filtered))
})


test_that("filter_genes_zero_expression remove correctly genes 0", {
    exp <- t(matrix(c(rep(0,15), rep(1,5)),ncol = 2))
    rownames(exp) <- c("one_hundred_percent_zero","fifty_percent_zero")
    colnames(exp) <- paste0("Samples",1:10)

    filtered <- filter_genes_zero_expression(exp, max.samples.percentage = 1)
    expect_true(is(filtered, "matrix"))
    expect_false("one_hundred_percent_zero" %in% rownames(filtered))
    expect_true("fifty_percent_zero" %in% rownames(filtered))

    filtered <- filter_genes_zero_expression(exp, max.samples.percentage = 0)
    expect_true(nrow(filtered) == 0)
})


