test_that("plot_interaction_model return a ggplot object", {
    data("dna.met.chr21")
    dna.met.chr21 <- make_se_from_dnam_probes(dna.met.chr21)
    data("gene.exp.chr21.log2")
    triplet <- data.frame(
        "regionID" = rownames(dna.met.chr21)[1:10],
        "TF" = rownames(gene.exp.chr21.log2)[11:20],
        "target" = rownames(gene.exp.chr21.log2)[1:10]
    )
    results <- interaction_model(triplet, dna.met.chr21, gene.exp.chr21.log2)
    plots <-  plots <- plot_interaction_model(
        triplet.results = results[1,],
        dnam = dna.met.chr21,
        exp = gene.exp.chr21.log2
    )
    # Adding color to samples
    metadata <- clinical[,"sample_type",drop = FALSE]
    plots <- plot_interaction_model(
        triplet.results = results[1,],
        dnam = dna.met.chr21,
        exp = gene.exp.chr21.log2,
        metadata =  metadata
    )
    expect_true(is(plots[[1]],"ggplot"))
})


