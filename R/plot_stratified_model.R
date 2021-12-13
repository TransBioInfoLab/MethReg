#' @title Plot stratified model results
#' @description Create several plots to show interaction data
#' TF expression with target gene interaction using a linear model
#' \deqn{log2(RNA target) ~ log2(TF)}
#' to samples with highest DNAm values (top 25 percent) and lowest DNAm values (bottom 25 percent), separately.
#' @param triplet.results Output from function stratified_model
#' with Region ID, TF  (column name: TF),  and target gene  (column name: target),
#' p-values and estimates of interaction
#' @param dnam DNA methylation matrix  or SummarizedExperiment object
#'  (columns: samples same order as met, rows: regions/probes)
#' @param exp A gene expression matrix or SummarizedExperiment object
#'  (columns: samples same order as met, rows: genes)
#' @param metadata A data frame with samples as row names and one columns that will be used to
#' color the samples
#' @param tf.activity.es A matrix with normalized enrichment scores for each TF across all samples
#' to be used in linear models instead of TF gene expression.
#' @param label.dnam Used for label text. Option "beta-value" and "residuals"
#' @param label.exp Used for label text. Option "expression" and "residuals"
#' @return A ggplot object, includes a table with results from fitting stratified model,
#' and the following scatter plots: 1) TF vs DNAm, 2) Target vs DNAm,
#' 3) Target vs TF, 4) Target vs TF for samples in Q1 and Q4 for DNA methylation,
#' 5) Target vs DNAm for samples in Q1 and Q4 for the TF
#' @param dnam.group.threshold DNA methylation threshold percentage to define samples 
#' in the low methylated group and high methylated group. For example, 
#' setting the threshold to 0.3 (30\%) will assign samples with the lowest 30\% 
#' methylation in the low group and the highest 30\% methylation in the high group. 
#' Default is 0.25 (25\%), accepted threshold range (0.0,0.5].
#' @importFrom ggpubr ggscatter ggarrange ggtexttable ttheme
#' @importFrom ggplot2 xlab ylab geom_smooth
#' @importFrom tibble as_tibble
plot_stratified_model <-  function(
    triplet.results,
    dnam,
    exp,
    metadata,
    label.dnam = "beta-value",
    label.exp = "expression",
    tf.activity.es = NULL,
    dnam.group.threshold = 0.25
){
    
    label.dnam <- match.arg(label.dnam, choices = c("beta-value","residuals"))
    label.exp <- match.arg(label.exp, choices = c("expression","residuals"))
    
    #---------------------------------------------------------------------------
    # Input checking
    #---------------------------------------------------------------------------
    if (missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
    if (missing(exp)) stop("Please set exp argument with gene expression matrix")
    if (missing(triplet.results)) stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
    if (!all(c("regionID","TF","target") %in% colnames(triplet.results))) {
        stop("triplet must have the following columns names: regionID, TF, target")
    }
    if (is(dnam,"SummarizedExperiment")) dnam <- assay(dnam)
    if (is(exp,"SummarizedExperiment")) exp <- assay(exp)
    
    check_data(dnam, exp, metadata)
    
    out <- plyr::alply(
        .data = triplet.results,
        .margins = 1,
        .fun = function(row.triplet,metadata){
            
            df <- get_triplet_data(
                exp = exp,
                dnam = dnam,
                row.triplet = row.triplet,
                tf.es =  tf.activity.es
            )
            
            color <- NULL
            if(!missing(metadata)) {
                df <- cbind(df,metadata)
                color <- colnames(metadata)[1]
            }
            
            plots <- get_plot_results(
                df = df,
                row.triplet = row.triplet,
                color =  color,
                label.dnam =label.dnam,
                label.exp = label.exp,
                use_tf_enrichment_scores = is.null(tf.activity.es),
                dnam.group.threshold = dnam.group.threshold
            )
            save(plots,file = "test.rda")
            
            # Reformat p-values for better looking on the plots
            for (idx in grep("pval|fdr|value",colnames(row.triplet))) {
                row.triplet[,idx] <- format.pval(
                    row.triplet[,idx]  %>% as.data.frame(),
                    digits = 3)
            }
            for (idx in grep("estimate|median|minus",colnames(row.triplet))) {
                row.triplet[,idx] <- format(
                    row.triplet[,idx]  %>% as.data.frame(),
                    digits = 3)
            }
            
            table.plots <- get_table_stratified_plot(row.triplet)
            
            # Arrange the plots on the same page
            plot.table <- ggarrange(
                ggarrange(
                    table.plots$table.plot.metadata,
                    ggarrange(
                        table.plots$table.plot.lm.dna.low,
                        table.plots$table.plot.lm.dna.high,
                        nrow = 2
                    ),
                    ncol = 2
                ),
                ggarrange(
                    plots$tf.target,
                    plots$dnam.target,
                    #plots$dnam.tf,
                    ncol = 2
                ),
                plots$tf.target.quantile,
                plots$dnam.target.quantile,
                nrow = 4,
                heights = c(2,2,2.5,2.5))
            plot.table
        }, .progress = "time", metadata = metadata)
    attr(out,"split_type") <- NULL
    attr(out,"split_labels") <- NULL
    
    names(out) <- paste0(
        triplet.results$regionID,"_TF_",triplet.results$TF,
        "_target_",triplet.results$target
    )
    out
}

get_table_stratified_plot <- function(row.triplet){
    
    base_size <- 9
    tab <- row.triplet %>%
        dplyr::select(
            c("regionID",
              "target",
              "target_symbol",
              "TF",
              "TF_symbol",
              "TF.role",
              "DNAm.effect"
            )
        ) %>% t() %>% as_tibble(rownames = "Variable")
    
    tab$Variable <- c(
        "Region ID",
        "Target gene ID",
        "Target gene Symbol",
        "TF gene ID",
        "TF gene Symbol",
        "TF role",
        "DNAm effect"
    )
    table.plot.metadata <- ggtexttable(
        tab,
        rows = NULL,
        cols = NULL,
        theme = ttheme("mOrange", base_size = base_size)
    )
    
    # Get results for linear model with all samples
    table.plot.lm.dna.low <- get_table_plot_results(row.triplet, type = "DNAmlow")
    
    # Get results for linear model with DNAm high samples
    table.plot.lm.dna.high <- get_table_plot_results(row.triplet, type = "DNAmhigh")
    
    table.plot.list <- list(
        "table.plot.metadata" = table.plot.metadata,
        "table.plot.lm.dna.low" = table.plot.lm.dna.low,
        "table.plot.lm.dna.high" = table.plot.lm.dna.high
    )
    
    return(table.plot.list)
}

