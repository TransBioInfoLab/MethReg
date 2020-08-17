test_that("plot_interaction_model return a ggplot object", {
    data("dna.met.chr21")
    dna.met.chr21 <- make_se_from_dnam_probes(dna.met.chr21)
    data("gene.exp.chr21.log2")
    triplet <- data.frame(
        "regionID" = rownames(dna.met.chr21)[1:5],
        "TF" = rownames(gene.exp.chr21.log2)[11:15],
        "target" = rownames(gene.exp.chr21.log2)[1:5])
    results <- stratified_model(triplet, dna.met.chr21, gene.exp.chr21.log2)
    plots <- plot_stratified_model(results[1,], dna.met.chr21, gene.exp.chr21.log2)
    # Adding color to samples
    metadata <- clinical[,"sample_type",drop = FALSE]
    plots <- plot_stratified_model(results[1,], dna.met.chr21, gene.exp.chr21.log2,metadata)
    expect_true(is(plots[[1]],"ggplot"))
})


