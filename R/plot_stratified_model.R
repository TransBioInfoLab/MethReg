#' @title Plot stratifed model results
#' @description Create several plots to show interaction data
#' TF expression with target gene interaction using a linear model
#' \deqn{log2(RNA target) ~ log2(TF)}
#' to samples with higest DNAm values (top 25 percent) and lowest DNAm values (bottom 25 percent), separately.
#' @param triplet.results Output from function stratified_model
#' with Region ID, TF  (column name: TF),  and target gene  (column name: target),
#' p-values and estimates of interaction
#' @param dnam DNA methylation matrix  (columns: samples same order as met, rows: regions/probes)
#' @param exp gene expression matrix (columns: samples same order as met, rows: genes)
#' @param metadata A data frame with samples as rownames and one columns that will be used to
#' color the samples
#' @return A ggplot object with a table with the results and the
#' the following scatter plots: 1) TF vs DNAm, 2) Target vs DNAm,
#' 3) Target vs TF, 4) Target vs TF for samples in Q1 and Q4 for DNA methylation,
#' 5) Target vs DNAm for samples in Q1 and Q4 for the TF
#' @examples
#' library(dplyr)
#' dnam <- runif(20,min = 0,max = 1) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(dnam) <- c("chr3:203727581-203728580")
#' colnames(dnam) <- paste0("Samples",1:20)
#'
#' exp.target <-  runif(20,min = 0,max = 10) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(exp.target) <- c("ENSG00000232886")
#' colnames(exp.target) <- paste0("Samples",1:20)
#'
#' exp.tf <- runif(20,min = 0,max = 10) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(exp.tf) <- c("ENSG00000232888")
#' colnames(exp.tf) <- paste0("Samples",1:20)
#'
#' exp <- rbind(exp.tf, exp.target)
#'
#' triplet <- data.frame(
#'    "regionID" =  c("chr3:203727581-203728580"),
#'    "target" = "ENSG00000232886",
#'    "TF" = "ENSG00000232888"
#')
#'
#' results <- stratified_model(triplet = triplet,dnam = dnam, exp = exp)
#' plots <- plot_stratified_model(
#'     triplet.results = results,
#'     dnam = dnam,
#'     exp = exp
#' )
#' \dontrun{
#' data("dna.met.chr21")
#' dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
#' data("gene.exp.chr21")
#' triplet <- data.frame(
#'    "regionID" = rownames(dna.met.chr21)[1:5],
#'    "TF" = rownames(gene.exp.chr21)[11:15],
#'    "target" = rownames(gene.exp.chr21)[1:5]
#' )
#' results <- stratified_model(
#'    triplet = triplet,
#'    dnam = dna.met.chr21,
#'    exp = gene.exp.chr21
#' )
#' plots <- plot_stratified_model(
#'    triplet.results = results[1,],
#'    dnam = dna.met.chr21,
#'    exp = gene.exp.chr21
#' )
#' # Adding color to samples
#' metadata <- clinical[,"sample_type",drop = FALSE]
#' plots <- plot_stratified_model(
#'    triplet.results = results[1,],
#'    dnam = dna.met.chr21,
#'    exp = gene.exp.chr21,
#'    metadata = metadata
#' )
#' }
#' @export
#' @importFrom ggpubr ggscatter ggarrange ggtexttable ttheme
#' @importFrom ggplot2 xlab ylab geom_smooth
#' @importFrom tibble as_tibble
plot_stratified_model <-  function(
    triplet.results,
    dnam,
    exp,
    metadata
){

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

            rna.target <- exp[rownames(exp) == as.character(row.triplet$target), , drop = FALSE]
            rna.tf <- exp[rownames(exp) == as.character(row.triplet$TF), , drop = FALSE]
            met <- dnam[rownames(dnam) == as.character(row.triplet$regionID), , drop = FALSE]

            df <- data.frame(
                rna.target = rna.target %>% as.numeric,
                met = met %>% as.numeric,
                rna.tf = rna.tf %>% as.numeric,
                stringsAsFactors = FALSE
            )

            color <- NULL
            if(!missing(metadata)){
                df <- cbind(df,metadata)
                color <- colnames(metadata)[1]
            }

            plots <- get_plot_results(df, row.triplet, color)

            # Reformat p-values for better looking on the plots
            for(idx in grep("pval|fdr|value",colnames(row.triplet))) {
                row.triplet[,idx] <- format.pval(
                    row.triplet[,idx]  %>% as.data.frame(),
                    digits = 3)
            }
            for(idx in grep("estimate|median|minus",colnames(row.triplet))) {
                row.triplet[,idx] <- format(
                    row.triplet[,idx]  %>% as.data.frame(),
                    digits = 3)
            }

            table.plots <- get_table_stratified_plot(row.triplet)

            # Arrange the plots on the same page
            plot.table <- ggarrange(
                ggarrange(table.plots$table.plot.metadata,
                          ggarrange(
                              table.plots$table.plot.lm.dna.low,
                              table.plots$table.plot.lm.dna.high,
                              nrow = 2),
                          ncol = 2),
                ggarrange(plots$tf.target,
                          plots$dnam.target,
                          plots$dnam.tf,
                          ncol = 3),
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
              "TF.affinity")
        ) %>% t() %>% as_tibble(rownames = "Variable")

    tab$Variable <- c(
        "Region ID",
        "Target gene ID",
        "Target gene Symbol",
        "TF gene ID",
        "TF gene Symbol",
        "TF role",
        "TF affinity"
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

