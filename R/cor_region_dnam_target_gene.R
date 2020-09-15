#' @title Evaluate correlation of DNA methylation region and
#' target gene expression
#' @description This function evaluate the correlation of the DNA methylation
#' and target gene expression using spearman rank correlation test.
#' Note that genes with RNA expression equal to 0 for all samples
#' will not be evaluated.
#' @return A data frame with the following information: regionID, target gene,
#' correlation pvalue and estimate between
#' DNA methylation and target gene expression, FDR corrected p-values.
#' @param links A dataframe with the following columns:
#' regionID (DNA methylation) and target (target gene)
#' @param dnam DNA methylation matrix or SummarizedExperiment object
#' with regions/cpgs in rows and samples in columns are samples.
#' Samples should be in the same order as gene expression matrix (exp).
#' @param exp Gene expression matrix  or SummarizedExperiment object
#' (rows are genes, columns are samples) log2-normalized (log2(exp + 1)).
#' Samples should be in the same order as the DNA methylation matrix.
#' @param filter.results
#' Filter results using min.cor.pval and min.cor.estimate thresholds
#' @param min.cor.pval
#' P-value threshold filter (default: 0.05)
#' @param min.cor.estimate
#' Correlation estimate threshold filter (default: not applied)
#' @param cores Number of CPU cores to be used. Default 1.
#' @importFrom plyr adply
#' @importFrom tibble tibble
#' @importFrom stats p.adjust cor.test
#' @export
#' @examples
#' dnam <- t(matrix(sort(c(runif(20))), ncol = 1))
#' rownames(dnam) <- c("chr3:203727581-203728580")
#' colnames(dnam) <- paste0("Samples",1:20)
#' exp <- dnam
#' rownames(exp) <- c("ENSG00000232886")
#' colnames(exp) <- paste0("Samples",1:20)
#'
#' links <- data.frame(
#'    "regionID" =  c("chr3:203727581-203728580"),
#'    "target" = "ENSG00000232886"
#' )
#'
#' # Correlated DNAm and gene expression, display only significant associations
#' results.cor.pos <- cor_dnam_target_gene(
#'    links = links,
#'    dnam = dnam,
#'    exp = exp,
#'    filter.results = TRUE,
#'    min.cor.pval = 0.05,
#'    min.cor.estimate = 0.0
#')
cor_dnam_target_gene <- function(
    links,
    dnam,
    exp,
    filter.results = TRUE,
    min.cor.pval = 0.05,
    min.cor.estimate = 0.0,
    cores = 1
){

    #-------------------------------------------------------------------------
    # Checking input
    #-------------------------------------------------------------------------
    if (missing(exp) || is.null(exp)) {
        stop("Please set an R matrix/SE to exp argument")
    }

    if (missing(dnam) || is.null(dnam)) {
        stop("Please set an R matrix/SE to dnam argument")
    }

    if (is(dnam,"SummarizedExperiment")) {
        dnam <- assay(dnam)
    }

    if (!is(dnam,"matrix")) {
        stop("dnam input is wrong")
    }

    if (is(exp,"SummarizedExperiment")) {
        exp <- assay(exp)
    }

    if (!is(exp,"matrix")) {
        stop("exp input is wrong")
    }

    if (ncol(dnam) != ncol(exp)) {
        stop("exp and dnam does not have the same size")
    }

    if (!all(c("target","regionID") %in% colnames(links))) {
        stop("links object must have target and regionID columns")
    }
    #-------------------------------------------------------------------------

    message("Removing genes with RNA expression equal to 0/NA for all samples")
    exp <- filter_genes_zero_expression(exp = exp, max.samples.percentage = 100)

    message("Removing regions with beta-values equal to NA for all samples")
    regions.keep <- (rowSums(is.na(dnam)) < ncol(dnam)) %>% which %>% names
    dnam <- dnam[regions.keep,, drop = FALSE]

    links <- links[links$target %in% rownames(exp),]
    links <- links[links$regionID %in% rownames(dnam),]
    if (nrow(links) == 0) stop("links not found in data. Please check rownames and links provided.")

    # reducing object sizes in case we will make it parallel
    exp <- exp[rownames(exp) %in% links$target,,drop = FALSE]
    dnam <- dnam[rownames(dnam) %in% links$regionID,, drop = FALSE]

    parallel <- register_cores(cores)

    correlation.df <- plyr::adply(
        .data = links,
        .margins = 1,
        .fun = function(link){
            tryCatch({
                exp <- exp[link$target,]
                dnam <- dnam[rownames(dnam) == link$regionID,]
                suppressWarnings({
                    res <- cor.test(
                        x = exp %>% as.numeric,
                        y = dnam %>% as.numeric,
                        method = "spearman",
                        exact = TRUE
                    )
                })
                return(
                    tibble(
                        "met_exp_cor_pvalue" = res$p.value,
                        "met_exp_cor_estimate" = res$estimate
                    )
                )
            }, error = function(e){
                return(
                    tibble(
                        "met_exp_cor_pvalue" = NA,
                        "met_exp_cor_estimate" = NA
                    )
                )
            })
        },
        .progress = "time",
        .parallel = parallel,
        .inform = TRUE
    )

    correlation.df <- na.omit(correlation.df)
    correlation.df$met_exp_cor_fdr <- p.adjust(
        correlation.df$met_exp_cor_pvalue,
        method = "fdr"
    )

    if (filter.results) {
        correlation.df <- correlation.df %>%
            dplyr::filter(
                .data$met_exp_cor_fdr <= min.cor.pval &
                    abs(.data$met_exp_cor_estimate) >= min.cor.estimate
            )
    }

    return(correlation.df)
}
