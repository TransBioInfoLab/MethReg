test_that("cor_dnam_target_gene correlation signal is correct", {
    dnam <- t(matrix(sort(c(runif(20))), ncol = 1))
    rownames(dnam) <- c("chr3:203727581-203728580")
    colnames(dnam) <- paste0("Samples",1:20)

    exp <- dnam
    rownames(exp) <- c("ENSG00000232886")
    colnames(exp) <- paste0("Samples",1:20)

    # Map example region to closest gene
    pair.dnam.target <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886"
    )

    # Correalted DNAm and gene expression, display only significant associations
    results.cor.pos <- cor_dnam_target_gene(
        pair.dnam.target = pair.dnam.target,
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
    pair.dnam.target <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886"
    )

    # Correalted DNAm and gene expression, display only significant associations
    results.cor.neg <- cor_dnam_target_gene(
        pair.dnam.target = pair.dnam.target,
        dnam = dnam,
        exp = exp,
        filter.results = FALSE
    )

    expect_true(results.cor.neg$met_exp_cor_estimate < 0)
    expect_true(results.cor.pos$met_exp_cor_estimate > 0)
    expect_true(results.cor.neg$met_exp_cor_fdr < 0.05)
    expect_true(results.cor.pos$met_exp_cor_fdr < 0.05)
})

test_that("cor_region_dnam_target_gene filter results", {
    dnam <- t(matrix(sort(c(runif(20))), ncol = 1))
    rownames(dnam) <- c("chr3:203727581-203728580")
    colnames(dnam) <- paste0("Samples",1:20)

    exp <- t(matrix(rep(0.5,20), ncol = 1))
    rownames(exp) <- c("ENSG00000232886")
    colnames(exp) <- paste0("Samples",1:20)

    # Map example region to closest gene
    pair.dnam.target <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886"
    )

    # Correalted DNAm and gene expression, display only significant associations
    results.cor.pos <- cor_dnam_target_gene(
        pair.dnam.target = pair.dnam.target,
        dnam = dnam,
        exp = exp,
        filter.results = TRUE
    )

    expect_true(nrow(results.cor.pos) == 0)
})

test_that("cor_dnam_target_gene checks input", {

    dnam <- t(matrix(sort(c(runif(20))), ncol = 1))
    rownames(dnam) <- c("chr3:203727581-203728580")
    colnames(dnam) <- paste0("Samples",1:20)

    exp <- t(matrix(rep(0.5,20), ncol = 1))
    rownames(exp) <- c("ENSG00000232886")
    colnames(exp) <- paste0("Samples",1:20)

    # Map example region to closest gene
    pair.dnam.target <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886"
    )
    # Correalted DNAm and gene expression, display only significant associations
    expect_error(cor_dnam_target_gene())
    expect_error(cor_dnam_target_gene(pair.dnam.target = pair.dnam.target))
    expect_error(cor_dnam_target_gene(pair.dnam.target = pair.dnam.target,dnam = dnam))
    expect_error(cor_dnam_target_gene(pair.dnam.target = pair.dnam.target,exp = exp))
    expect_error(cor_dnam_target_gene(pair.dnam.target = pair.dnam.target,dnam, "error", exp = exp))
})





