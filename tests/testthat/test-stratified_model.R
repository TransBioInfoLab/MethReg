test_that("stratified_model works", {
    data("dna.met.chr21")
    dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
    data("gene.exp.chr21.log2")
    triplet <- data.frame("regionID" = rownames(dna.met.chr21)[1:5],
                          "TF" = rownames(gene.exp.chr21.log2)[11:15],
                          "target" = rownames(gene.exp.chr21.log2)[1:5])
    results <- stratified_model(triplet, dna.met.chr21, gene.exp.chr21.log2)
    expect_true("regionID" %in% colnames(results))
    expect_true("TF" %in% colnames(results))
    expect_true("target" %in% colnames(results))
    expect_true("DNAmlow_pval_rna.tf" %in% colnames(results))
    expect_true("DNAmlow_estimate_rna.tf" %in% colnames(results))
    expect_true("DNAmhigh_pval_rna.tf" %in% colnames(results))
    expect_true("DNAmhigh_estimate_rna.tf" %in% colnames(results))
    expect_true("TF.affinity" %in% colnames(results))
    expect_true("TF.role" %in% colnames(results))
})

test_that("stratified_model handles 0 cases", {
    data("dna.met.chr21")
    dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
    data("gene.exp.chr21.log2")
    triplet <- data.frame(
        "regionID" = rownames(dna.met.chr21)[1],
        "TF" = rownames(gene.exp.chr21.log2)[11],
        "target" = rownames(gene.exp.chr21.log2)[1]
    )
    gene.exp.chr21.log2[1,-1] <- 0 # at least one of the two groups will have only 0 values
    results <- stratified_model(triplet, dna.met.chr21, gene.exp.chr21.log2)
    expect_true(is.na(results$TF.affinity))
    expect_true(is.na(results$TF.role))
})
