test_that("interaction_model handles 0 cases in more than 25% of samples with zeroinflated model", {
    dnam <- t(matrix(sort(c(runif(20))), ncol = 1))
    rownames(dnam) <- c("chr3:203727581-203728580")
    colnames(dnam) <- paste0("Samples",1:20)

    exp.target <-  c(runif(10,min = 0,max = 0),
                     runif(10,min = 0,max = 10)) %>%
        matrix(ncol = 1) %>%  t
    rownames(exp.target) <- c("ENSG00000232886")
    colnames(exp.target) <- paste0("Samples",1:20)

    exp.tf <-  t(matrix(sort(c(runif(20))), ncol = 1))
    rownames(exp.tf) <- c("ENSG00000232888")
    colnames(exp.tf) <- paste0("Samples",1:20)

    exp <- rbind(exp.tf, exp.target)
    # Map example region to closest gene
    triplet <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886",
        "TF" = "ENSG00000232888"
    )

    results <- interaction_model(triplet, dnam, exp)
    expect_equal(results$Model.quantile, "Zero-inflated Negative Binomial Model")
    expect_equal(results$Model.interaction, "Zero-inflated Negative Binomial Model")
    expect_true("regionID" %in% colnames(results))
    expect_true("TF" %in% colnames(results))
    expect_true("target" %in% colnames(results))
    expect_true("pval_met" %in% colnames(results))
    expect_true("pval_rna.tf" %in% colnames(results))
    expect_true("pval_met:rna.tf" %in% colnames(results))
    expect_true("estimate_met" %in% colnames(results))
    expect_true("estimate_rna.tf" %in% colnames(results))
    expect_true("estimate_met:rna.tf" %in% colnames(results))
    expect_true("quant_pval_metGrp" %in% colnames(results))
    expect_true("quant_pval_rna.tf" %in% colnames(results))
    expect_true("quant_pval_metGrp:rna.tf" %in% colnames(results))
    expect_true("quant_estimate_metGrp" %in% colnames(results))
    expect_true("quant_estimate_rna.tf" %in% colnames(results))
    expect_true("quant_estimate_metGrp:rna.tf" %in% colnames(results))
})


test_that("interaction_model performs rlm if no 0", {
    dnam <- runif(n = 20, min = 0, max = 1) %>%
        sort %>%
        matrix(ncol = 1) %>%  t
    rownames(dnam) <- c("chr3:203727581-203728580")
    colnames(dnam) <- paste0("Samples",1:20)

    exp.target <- runif(n = 20, min = 0, max = 1) %>%
        sort %>%
        matrix(ncol = 1) %>%  t
    rownames(exp.target) <- c("ENSG00000232886")
    colnames(exp.target) <- paste0("Samples",1:20)

    exp.tf <- runif(n = 20, min = 1, max = 1) %>%
        sort %>%
        matrix(ncol = 1) %>%  t
    rownames(exp.tf) <- c("ENSG00000232888")
    colnames(exp.tf) <- paste0("Samples",1:20)

    exp <- rbind(exp.tf, exp.target)
    # Map example region to closest gene
    triplet <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886",
        "TF" = "ENSG00000232888"
    )

    results <- interaction_model(triplet, dnam, exp)
    expect_equal(results$Model.interaction, "Robust Linear Model")
    expect_equal(results$Model.quantile, "Robust Linear Model")
    expect_true("regionID" %in% colnames(results))
    expect_true("TF" %in% colnames(results))
    expect_true("target" %in% colnames(results))
    expect_true("pval_met" %in% colnames(results))
    expect_true("pval_rna.tf" %in% colnames(results))
    expect_true("pval_met:rna.tf" %in% colnames(results))
    expect_true("estimate_met" %in% colnames(results))
    expect_true("estimate_rna.tf" %in% colnames(results))
    expect_true("estimate_met:rna.tf" %in% colnames(results))
    expect_true("quant_pval_metGrp" %in% colnames(results))
    expect_true("quant_pval_rna.tf" %in% colnames(results))
    expect_true("quant_pval_metGrp:rna.tf" %in% colnames(results))
    expect_true("quant_estimate_metGrp" %in% colnames(results))
    expect_true("quant_estimate_rna.tf" %in% colnames(results))
    expect_true("quant_estimate_metGrp:rna.tf" %in% colnames(results))
})


