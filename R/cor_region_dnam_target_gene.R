#' @title Evaluate correlation of DNA methylation region and target gene expression
#' @description Evaluate correlation of the DNA methylation region and target gene expression
#' using spearman correlation test
#' @param links A dataframe with the following columns: regionID (DNA methylation) and target (target gene)
#' @param met DNA methylation matrix (rows are regions and columns are samples). Samples should be in the
#' same order as gene expression.
#' @param exp Gene expression matrix (rows are genes, columns are samples) log2-normalized (log2(exp + 1)).
#' Samples should be in the same order as the DNA methylation matrix.
#' @param min.cor.pval Filter of significant correlations (default: 0.05)
#' @param min.cor.estimate Filter of significant correlations (default: not applied)
#' @param file.out If provided, name of a csv file which will be used to save the results.
#' @param cores Number of CPU cores to be used. Default 1.
#' @importFrom plyr adply
#' @importFrom tibble tibble
#' @importFrom stats p.adjust cor.test
#' @importFrom dplyr filter
#' @export
#' @examples
#' library(dplyr)
#' data("dna.met.chr21")
#' dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
#' regions.gr <- make_granges_from_names(rownames(dna.met.chr21))
#' data("gene.exp.chr21")
#'
#' # Map example region to closest gene
#' links <- get_region_target_gene(regions.gr = regions.gr, genome = "hg19", method = "closest.gene")
#'
#' # Samples in met and exp datasets should be in the same order.
#' identical (colnames (dna.met.chr21), colnames(gene.exp.chr21))
#'
#' # Correalted DNAm and gene expression, display only significant associations
#' results <- cor_region_dnam_target_gene(links = links, met = dna.met.chr21, exp = gene.exp.chr21)
#'
#' # display all associations
#' results.all <- cor_region_dnam_target_gene(
#'    links = links,
#'    met = dna.met.chr21,
#'    exp = gene.exp.chr21,
#'    min.cor.pval = 1)
cor_region_dnam_target_gene <- function(
    links,
    met,
    exp,
    min.cor.pval = 0.05,
    min.cor.estimate = 0.0,
    file.out,
    cores = 1
){

    if(is.null(exp)) stop("Please set exp matrix")
    if(is.null(met)) stop("Please set met matrix")
    if(ncol(met) != ncol(exp)) stop("exp and met does not have the same size")
    if(!all(c("target","regionID") %in% colnames(links))) stop("links object must have target and regionID columns")

    # remove triplet with RNA expression equal to 0 for more than 25% of the samples
    genes.keep <- (rowSums(exp == 0) < 0.25) %>% which %>% names
    exp <- exp[genes.keep,]

    links <- links[links$target %in% rownames(exp),]
    links <- links[links$regionID %in% rownames(met),]
    if(nrow(links) == 0) stop("links not found in data. Please check rownames and links provided.")

    parallel <- register_cores(cores)

    correlation.df <- plyr::adply(
        .data = links,
        .margins = 1,
        .fun = function(link){
            exp <- exp[link$target,]
            met <- met[rownames(met) == link$regionID,]
            res <- cor.test(exp %>% as.numeric,
                            met %>% as.numeric,
                            method = "spearman",
                            exact = FALSE)
            return(tibble("met_exp_cor_pvalue" = res$p.value,
                          "met_exp_cor_estimate" = res$estimate))
        },.progress = "time",.parallel = parallel)


    correlation.df <- na.omit(correlation.df)
    correlation.df$met_exp_cor_fdr <- p.adjust(correlation.df$met_exp_cor_pvalue, method = "fdr")

    correlation.df <- correlation.df %>%
        dplyr::filter(.data$met_exp_cor_fdr <= min.cor.pval & abs(.data$met_exp_cor_estimate) >= min.cor.estimate)

    if(!missing(file.out)) readr::write_tsv(x = correlation.df, path = file.out)

    return(correlation.df)
}
