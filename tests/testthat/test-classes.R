test_that("make_granges_from_names works", {
    res <- getClassification(0,1)
    expect_equal(res$TF,"Activator")
    expect_equal(res$DNAm,"M-plus")

    res <- getClassification(1,0.5)
    expect_equal(res$TF,"Activator")
    expect_equal(res$DNAm,"M-minus")

    res <- getClassification(0,-1)
    expect_equal(res$TF,"Repressor")
    expect_equal(res$DNAm,"M-plus")

    res <- getClassification(-0.5,0.4)
    expect_equal(res$TF,"Repressor")
    expect_equal(res$DNAm,"M-minus")

})