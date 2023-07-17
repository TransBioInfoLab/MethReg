test_that("create_triplet_distance_based runs", {
    skip_on_bioc()
    triplet <- create_triplet_distance_based(
        region = "chr3:189631389-189631389",
        genome = "hg38",
        target.method = "genes.promoter.overlap",
        motif.search.window.size = 50,
        motif.search.p.cutoff = 1
    )
    expect_true("TP63" %in% triplet$target_symbol)
    expect_true(
        all(
            c("regionID",
              "target_symbol",
              "target",
              "TF_symbol",
              "TF") %in% colnames(triplet)
        )
    )
    expect_type(triplet$TF_symbol,"character")
    expect_type(triplet$target_symbol,"character")
})


