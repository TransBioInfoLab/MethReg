test_that("create_triplet_distance_based runs", {
    triplet <- create_triplet_distance_based(
        region = "chr3:189631389-189631389",
        genome = "hg38",
        target.method = "genes.promoter.overlap",
        motif.search.window.size = 50,
        motif.search.p.cutoff = 1
    )
    expect_true("TP63" %in% triplet$target_gene_name)
    expect_true(
        all(
            c("regionID",
              "target_gene_name",
              "target",
              "TF_external_gene_name",
              "TF") %in% colnames(triplet)
        )
    )
    expect_type(triplet$TF_external_gene_name,"character")
    expect_type(triplet$target_gene_name,"character")
})


