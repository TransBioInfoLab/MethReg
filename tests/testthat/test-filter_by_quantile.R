test_that("filter_regions_by_mean_quantile_difference remove correctly regions without variance", {
    dnam <- t(matrix(c(rep(0,15), rep(1,5)),ncol = 2))
    rownames(dnam) <- c("to_be_removed","to_keep")
    colnames(dnam) <- paste0("Samples",1:10)
    filtered <- filter_regions_by_mean_quantile_difference(dnam, diff.mean.th = 0.2)
    expect_true(is(filtered, "matrix"))
    expect_true("to_keep" %in% rownames(filtered))
    expect_false("to_be_removed" %in% rownames(filtered))
})


test_that("filter_genes_by_quantile_mean_fold_change remove correctly genes without variance", {
    exp <- t(matrix(c(rep(0,15), rep(1,5)),ncol = 2))
    rownames(exp) <- c("to_be_removed","to_keep")
    colnames(exp) <- paste0("Samples",1:10)
    filtered <- filter_genes_by_quantile_mean_fold_change(exp, fold.change = 1.5)
    expect_true(is(filtered, "matrix"))
    expect_true("to_keep" %in% rownames(filtered))
    expect_false("to_be_removed" %in% rownames(filtered))
})

