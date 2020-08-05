test_that("cor_region_dnam_target_gene works with SE and is filtering", {
    # dna.met.chr21 <- get(data("dna.met.chr21"))
    # dna.met.chr21 <- make_se_from_dnam_probes(dna.met.chr21)
    # gene.exp.chr21 <- get(data("gene.exp.chr21"))
    # gene.exp.chr21 <- make_se_from_gene_matrix(gene.exp.chr21)
    # # Map example region to closest gene
    # links <- get_region_target_gene(
    #     regions.gr = dna.met.chr21 %>% rowRanges(),
    #     genome = "hg19",
    #     method = "closest.gene"
    # )
    #
    # # Correalted DNAm and gene expression, display only significant associations
    # results <- cor_region_dnam_target_gene(
    #     links = links,
    #     dnam = dna.met.chr21,
    #     exp = gene.exp.chr21,
    #     filter.results = FALSE
    # )
    #
    # results.filtered <- cor_region_dnam_target_gene(
    #     links = links,
    #     dnam = dna.met.chr21,
    #     exp = gene.exp.chr21,
    #     filter.results = TRUE
    # )
    #
    # results.filtered.stringent <- cor_region_dnam_target_gene(
    #     links = links,
    #     dnam = dna.met.chr21,
    #     exp = gene.exp.chr21,
    #     min.cor.estimate = 0.8,
    #     filter.results = TRUE
    # )
    #
    # results.filtered.less.stringent <- cor_region_dnam_target_gene(
    #     links = links,
    #     dnam = dna.met.chr21,
    #     exp = gene.exp.chr21,
    #     min.cor.estimate = 0.2,
    #     min.cor.pval = 1,
    #     filter.results = TRUE
    # )
    #
    # expect_true(nrow(results) > nrow(results.filtered))
    # expect_true(nrow(results[results$met_exp_cor_fdr < 0.05 &  abs(results$met_exp_cor_estimate) >= 0,]) == nrow(results.filtered))
    # expect_true(nrow(results.filtered.stringent) == 0)
    # expect_true(min(abs(results.filtered.less.stringent$met_exp_cor_estimate)) > 0.2)
})

test_that("cor_region_dnam_target_gene correlation signal is correct", {
    dnam <- t(matrix(sort(c(runif(20))), ncol = 1))
    rownames(dnam) <- c("chr3:203727581-203728580")
    colnames(dnam) <- paste0("Samples",1:20)

    exp <- dnam
    rownames(exp) <- c("ENSG00000232886")
    colnames(exp) <- paste0("Samples",1:20)

    # Map example region to closest gene
    links <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886"
    )

    # Correalted DNAm and gene expression, display only significant associations
    results.cor.pos <- cor_region_dnam_target_gene(
        links = links,
        dnam = dnam,
        exp = exp,
        filter.results = FALSE
    )

    dnam <- t(matrix(sort(c(runif(20))), ncol = 1))
    rownames(dnam) <- c("chr3:203727581-203728580")
    colnames(dnam) <- paste0("Samples",1:20)

    exp <- t(matrix(sort(dnam,decreasing = TRUE)))
    rownames(exp) <- c("ENSG00000232886")
    colnames(exp) <- paste0("Samples",1:20)

    # Map example region to closest gene
    links <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886"
    )

    # Correalted DNAm and gene expression, display only significant associations
    results.cor.neg <- cor_region_dnam_target_gene(
        links = links,
        dnam = dnam,
        exp = exp,
        filter.results = FALSE
    )

    expect_true(results.cor.neg$met_exp_cor_estimate < 0)
    expect_true(results.cor.pos$met_exp_cor_estimate > 0)
    expect_true(results.cor.neg$met_exp_cor_fdr < 0.05)
    expect_true(results.cor.pos$met_exp_cor_fdr < 0.05)
})





