test_that("get_human_tfs return list of TFs", {
    skip_if_offline()
    skip_on_bioc()
    human.tfs <- get_human_tfs()
    expect_true("FOXA1" %in% human.tfs$external_gene_name)
    expect_true(!"COL6A1" %in% human.tfs$external_gene_name)
})


