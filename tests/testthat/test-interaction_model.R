test_that("interaction_model works", {
    data("dna.met.chr21")
    dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
    data("gene.exp.chr21")
    triplet <- data.frame("regionID" = rownames(dna.met.chr21)[1:10],
                          "TF" = rownames(gene.exp.chr21)[11:20],
                          "target" = rownames(gene.exp.chr21)[1:10])
    results <- interaction_model(triplet, dna.met.chr21, gene.exp.chr21)
    expect_true("regionID" %in% colnames(results))
    expect_true("TF" %in% colnames(results))
    expect_true("target" %in% colnames(results))
    expect_true("pval_met" %in% colnames(results))
    expect_true("pval_rna.tf" %in% colnames(results))
    expect_true("pval_met:rna.tf" %in% colnames(results))
    expect_true("estimate_met" %in% colnames(results))
    expect_true("estimate_rna.tf" %in% colnames(results))
    expect_true("estimate_met:rna.tf" %in% colnames(results))
})

test_that("interaction_model handles 0 cases", {
    data("dna.met.chr21")
    dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
    data("gene.exp.chr21")
    triplet <- data.frame(
        "regionID" = rownames(dna.met.chr21)[1],
        "TF" = rownames(gene.exp.chr21)[11],
        "target" = rownames(gene.exp.chr21)[1]
    )
    gene.exp.chr21[1,-1] <- 0 # at least one of the two groups will have only 0 values
    results <- interaction_model(triplet, dna.met.chr21, gene.exp.chr21)
    expect_true(is.na(results$quant_estimate_metGrp))
})
