test_that("classification for non-signficant results", {

    # Not significant
    res <- get_tf_dnam_classification(
        low.estimate = 0.8, low.pval = 1,
        high.estimate = -0.2, high.pval = 1,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role, NA)
    expect_equal(res$DNAm.effect, NA)

    res <- get_tf_dnam_classification(
        low.estimate = 0.8, low.pval = NA,
        high.estimate = -0.2, high.pval = 1,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role, NA)
    expect_equal(res$DNAm.effect, NA)

    res <- get_tf_dnam_classification(
        low.estimate = -0.1, low.pval = NA,
        high.estimate = -0.8, high.pval = 0.0001,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role, "Repressor")
    expect_equal(res$DNAm.effect, "Enhancing")

    res <- get_tf_dnam_classification(
        low.estimate = 0.1, low.pval = NA,
        high.estimate = 0.8, high.pval = 0.0001,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role, "Activator")
    expect_equal(res$DNAm.effect, "Enhancing")

})

test_that("classification for signficant results in same direction", {

    # Same direction significants
    res <- get_tf_dnam_classification(
        low.estimate = 0.2, low.pval = 0.01,
        high.estimate = 0.8, high.pval = 0.01,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Activator")
    expect_equal(res$DNAm.effect,"Enhancing")

    res <- get_tf_dnam_classification(
        low.estimate = 0.8, low.pval = 0.01,
        high.estimate = 0.2, high.pval = 0.05,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Activator")
    expect_equal(res$DNAm.effect,"Attenuating")

    res <- get_tf_dnam_classification(
        low.estimate = -0.8, low.pval = 0.01,
        high.estimate = -0.2, high.pval = 0.05,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Repressor")
    expect_equal(res$DNAm.effect,"Attenuating")

    res <- get_tf_dnam_classification(
        low.estimate = -0.2, low.pval = 0.05,
        high.estimate = -0.8, high.pval = 0.01,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Repressor")
    expect_equal(res$DNAm.effect,"Enhancing")
})

test_that("classification just one signficant results in same direction", {

    res <- get_tf_dnam_classification(
        low.estimate = 0.2, low.pval = 1,
        high.estimate = 0.8, high.pval = 0.05,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Activator")
    expect_equal(res$DNAm.effect,"Enhancing")

    res <- get_tf_dnam_classification(
        low.estimate = 0.8, low.pval = 0.05,
        high.estimate = 0.2, high.pval = 1,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Activator")
    expect_equal(res$DNAm.effect,"Attenuating")

    res <- get_tf_dnam_classification(
        low.estimate = -0.8, low.pval = 0.05,
        high.estimate = -0.2, high.pval = 1,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Repressor")
    expect_equal(res$DNAm.effect,"Attenuating")

    res <- get_tf_dnam_classification(
        low.estimate = -0.2, low.pval = 1,
        high.estimate = -0.8, high.pval = 0.05,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Repressor")
    expect_equal(res$DNAm.effect,"Enhancing")

})

test_that("classification just one signficant results in different direction", {
    res <- get_tf_dnam_classification(
        low.estimate = -0.2, low.pval = 1,
        high.estimate = 0.8, high.pval = 0.05,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Activator")
    expect_equal(res$DNAm.effect,"Enhancing")

    res <- get_tf_dnam_classification(
        low.estimate = 0.8, low.pval = 0.05,
        high.estimate = -0.2, high.pval = 1,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Activator")
    expect_equal(res$DNAm.effect,"Attenuating")

    res <- get_tf_dnam_classification(
        low.estimate = -0.8, low.pval = 0.05,
        high.estimate = 0.2, high.pval = 1,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Repressor")
    expect_equal(res$DNAm.effect,"Attenuating")

    res <- get_tf_dnam_classification(
        low.estimate = 0.2, low.pval = 1,
        high.estimate = -0.8, high.pval = 0.05,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Repressor")
    expect_equal(res$DNAm.effect,"Enhancing")

})

test_that("classification both signficant results in different direction", {

    res <- get_tf_dnam_classification(
        low.estimate = -0.8, low.pval = 0.05,
        high.estimate = 0.2, high.pval = 0.05,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Dual")
    expect_equal(res$DNAm.effect,"Invert")

    res <- get_tf_dnam_classification(
        low.estimate = 0.8, low.pval = 0.05,
        high.estimate = -0.2, high.pval = 0.05,
        pvalue.threshold = 0.05
    )
    expect_equal(res$TF.role,"Dual")
    expect_equal(res$DNAm.effect,"Invert")

})
