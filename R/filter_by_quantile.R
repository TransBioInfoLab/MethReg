#' @title Filter regions with high DNA methylation variance
#' @description For each region, compares the mean DNA methylation of samples in high (Q4) vs. low (Q1)
#' and requires difference to be above a certain threshold.
#' @param dnam DNA methylation matrix
#' @param diff.mean.th Diff. mean Threshold between samples in Q4 and Q1
#' @param cores Number of CPU cores to be used in the analysis. Default: 1
#' @export
#' @examples
#' data("dna.met.chr21")
#' dna.met.chr21.filtered <- filter_regions_by_mean_quantile_difference(dna.met.chr21)
filter_regions_by_mean_quantile_difference <- function(dnam, diff.mean.th = 0.2, cores = 1){

    parallel <- register_cores(cores)
    diff_mean <- plyr::adply(dnam,.margins = 1,.fun = function(row){
        quant.met <-  quantile(row, na.rm = TRUE)

        low.cutoff <- quant.met[2]
        upper.cutoff <- quant.met[4]

        mean.q1 <- row %>% .[. <= low.cutoff] %>% mean(na.rm = TRUE)
        mean.q4 <- row %>% .[. >= upper.cutoff] %>% mean(na.rm = TRUE)
        data.frame("diff_mean" = mean.q4 - mean.q1, stringsAsFactors = FALSE)
    }, .progress = "time", .parallel = parallel)


    tab <- plyr::count(diff_mean$diff_mean > diff.mean.th)
    colnames(tab)[1] <- "Status"
    tab$Status[which(tab$Status == FALSE)] <- "Regions below threshold"
    tab$Status[which(tab$Status == TRUE)] <- "Regions above threshold"
    print(tab)

    diff.regions <- c(diff_mean %>% filter(diff_mean > 0.2) %>% pull(X1) %>% as.character())
    dnam[diff.regions,]
}


#' @title Filter regions with high DNA methylation variance
#' @description For each region, compares the mean DNA methylation of samples in high (Q4) vs. low (Q1)
#' and requires difference to be above a certain threshold.
#' @param dnam DNA methylation matrix
#' @param fold.change Mean fold change threshold between samples in Q4 and Q1. Default 1.5
#' @param cores Number of CPU cores to be used in the analysis. Default: 1
#' @export
#' @examples
#' data("gene.exp.chr21")
#' gene.exp.chr21.filtered <- filter_genes_by_quantile_mean_fold_change(gene.exp.chr21)
filter_genes_by_quantile_mean_fold_change <- function(exp, fold.change = 1.5, cores = 1){

    parallel <- register_cores(cores)
    diff.genes <- plyr::adply(exp,.margins = 1,.fun = function(row){
        quant <-  quantile(row, na.rm = TRUE)
        quant.fold.change <- data.frame("q4_div_q1" = quant[4] / quant[2])

        low.cutoff <- quant[2]
        upper.cutoff <- quant[4]

        mean.q1 <- row %>% .[. <= low.cutoff] %>% mean(na.rm = TRUE)
        mean.q4 <- row %>% .[. >= upper.cutoff] %>% mean(na.rm = TRUE)
        data.frame("diff_fold_change" = mean.q4 / mean.q1, stringsAsFactors = FALSE)
    }, .progress = "time",.parallel = parallel)


    tab <- plyr::count(diff.genes$diff_fold_change > fold.change)
    colnames(tab)[1] <- "Status"
    tab$Status[which(tab$Status == FALSE)] <- "Genes below threshold"
    tab$Status[which(tab$Status == TRUE)] <- "Genes above threshold"
    print(tab)

    diff.genes <- c(diff.genes %>% filter(diff_fold_change > fold.change) %>% pull(X1) %>% as.character())
    exp[diff.genes,]
}


#' @title Remove genes with gene expression 0 in more than 25% of the samples
#' @param exp Gene expression matrix
#' @param max.samples.percentage Max percentage of samples with gene expression as 0.
#' @noRd
filter_genes_zero_expression <- function(exp, max.samples.percentage = 0.25){
    genes.keep <- (rowSums(exp == 0) / ncol(exp) <= max.samples.percentage) %>% which %>% names
    message("Removing ", nrow(exp) - length(genes.keep), " out of ", nrow(exp), " genes")
    exp[genes.keep,]
}
