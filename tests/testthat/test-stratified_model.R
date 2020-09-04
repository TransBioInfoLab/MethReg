test_that("stratified_model works", {
    library(dplyr)
    dnam <- runif(20,min = 0,max = 1) %>%
        matrix(ncol = 1) %>%  t
    rownames(dnam) <- c("chr3:203727581-203728580")
    colnames(dnam) <- paste0("Samples",1:20)

    exp.target <-  runif(20,min = 0,max = 10) %>%
        matrix(ncol = 1) %>%  t
    rownames(exp.target) <- c("ENSG00000232886")
    colnames(exp.target) <- paste0("Samples",1:20)

    exp.tf <- runif(20,min = 0,max = 10) %>%
        matrix(ncol = 1) %>%  t
    rownames(exp.tf) <- c("ENSG00000232888")
    colnames(exp.tf) <- paste0("Samples",1:20)

    exp <- rbind(exp.tf, exp.target)

    triplet <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886",
        "TF" = "ENSG00000232888"
    )

    results <- stratified_model(
        triplet = triplet,
        dnam = dnam,
        exp =  exp
    )

    expect_true("regionID" %in% colnames(results))
    expect_true("TF" %in% colnames(results))
    expect_true("target" %in% colnames(results))
    expect_true("DNAmlow_pval_rna.tf" %in% colnames(results))
    expect_true("DNAmlow_estimate_rna.tf" %in% colnames(results))
    expect_true("DNAmhigh_pval_rna.tf" %in% colnames(results))
    expect_true("DNAmhigh_estimate_rna.tf" %in% colnames(results))
    expect_true("DNAm.effect" %in% colnames(results))
    expect_true("TF.role" %in% colnames(results))
})

test_that("stratified_model handles 0 cases", {
    library(dplyr)
    dnam <- runif(20,min = 0,max = 1) %>%
        matrix(ncol = 1) %>%  t
    rownames(dnam) <- c("chr3:203727581-203728580")
    colnames(dnam) <- paste0("Samples",1:20)

    exp.target <-  c(1,runif(19,min = 0,max = 0)) %>%
        matrix(ncol = 1) %>%  t
    rownames(exp.target) <- c("ENSG00000232886")
    colnames(exp.target) <- paste0("Samples",1:20)

    exp.tf <- runif(20,min = 0,max = 10) %>%
        matrix(ncol = 1) %>%  t
    rownames(exp.tf) <- c("ENSG00000232888")
    colnames(exp.tf) <- paste0("Samples",1:20)

    exp <- rbind(exp.tf, exp.target)

    triplet <- data.frame(
        "regionID" =  c("chr3:203727581-203728580"),
        "target" = "ENSG00000232886",
        "TF" = "ENSG00000232888"
    )

    results <- stratified_model(
        triplet = triplet,
        dnam = dnam,
        exp =  exp
    )

    expect_true(is.na(results$DNAm.effect))
    expect_true(is.na(results$TF.role))
})
