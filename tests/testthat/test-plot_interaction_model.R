test_that("plot_interaction_model return a ggplot object", {

    skip_on_bioc()
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

    results <- interaction_model(
        triplet = triplet,
        dnam =  dnam,
        exp =  exp,
        stage.wise.analysis = FALSE,
        filter.correlated.tf.exp.dna = FALSE,
        sig.threshold = 1,
        filter.triplet.by.sig.term = FALSE
    )

    plots <- plot_interaction_model(
        triplet.results = results,
        dnam = dnam,
        exp = exp
    )
    expect_true(is(plots[[1]],"ggplot"))

    metadata <- data.frame("Gender" = c(rep("Male",10),rep("Female",10)))
    rownames(metadata) <- colnames(dnam)

    plots <- plot_interaction_model(
        triplet.results = results,
        dnam = dnam,
        exp = exp,
        metadata =  metadata
    )
    expect_true(is(plots[[1]],"ggplot"))
})


