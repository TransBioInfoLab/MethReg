#' @title Evaluate correlation of TF expression and
#' target gene expression
#' @description This function evaluate the correlation of a TF
#' and target gene expression using spearman rank correlation test.
#' Note that genes with RNA expression equal to 0 for all samples
#' will not be evaluated.
#' @return A data frame with the following information: TF, target gene,
#' correlation p-value and estimate between
#' TF and target gene expression, FDR corrected p-values.
#' @param pair.tf.target A dataframe with the following columns:
#' TF and target (target gene)
#' @param tf.activity.es A matrix with normalized enrichment
#' scores for each TF across all samples to be used in linear models instead
#' of TF gene expression. See \code{\link{get_tf_ES}}.
#' @param exp Gene expression matrix  or SummarizedExperiment object
#' (rows are genes, columns are samples) log2-normalized (log2(exp + 1)).
#' Samples should be in the same order as the tf.activity.es matrix
#' @param cores Number of CPU cores to be used. Default 1.
#' @param verbose Show messages ?
#' @importFrom plyr adply
#' @importFrom tibble tibble
#' @importFrom stats p.adjust cor.test
#' @export
#' @examples
#' exp <- t(matrix(sort(c(runif(40))), ncol = 2))
#' rownames(exp) <- c("ENSG00000232886","ENSG00000232889")
#' colnames(exp) <- paste0("Samples",1:20)
#'
#' pair.tf.target <- data.frame(
#'    "TF" = "ENSG00000232889",
#'    "target" = "ENSG00000232886"
#' )
#'
#' # Correlated TF and gene expression
#' results.cor.pos <- cor_tf_target_gene(
#'    pair.tf.target = pair.tf.target,
#'    exp = exp,
#')
#' # Correlated TF and gene expression
#' results.cor.pos <- cor_tf_target_gene(
#'    pair.tf.target = pair.tf.target,
#'    exp = exp,
#'    tf.activity.es = exp
#')
cor_tf_target_gene <- function(
    pair.tf.target,
    exp,
    tf.activity.es,
    cores = 1,
    verbose = FALSE
){

    #-------------------------------------------------------------------------
    # Checking input
    #-------------------------------------------------------------------------
    if (missing(exp) || is.null(exp)) {
        stop("Please set an R matrix/SE to exp argument")
    }

    if (is(exp,"SummarizedExperiment")) {
        exp <- assay(exp)
    }

    if (!is(exp,"matrix")) {
        stop("exp input is wrong")
    }

    if(!missing(tf.activity.es)){
        if (ncol(tf.activity.es) != ncol(exp)) {
            stop("exp and tf.activity.es does not have the same size")
        }
    }

    if (!all(c("target","TF") %in% colnames(pair.tf.target))) {
        stop("pair.tf.target object must have target and tf columns")
    }
    #-------------------------------------------------------------------------

     if(verbose) {
         message("Removing genes with RNA expression equal to 0/NA for all samples")
     }

    exp <- filter_genes_zero_expression(exp = exp, max.samples.percentage = 100)

    pair.tf.target <- pair.tf.target %>% as.data.frame() %>%
        dplyr::filter(.data$target %in% rownames(exp))

    if (missing(tf.activity.es)){
        pair.tf.target <- pair.tf.target %>%
            dplyr::filter(.data$TF %in% rownames(exp))
    } else {
        pair.tf.target <- pair.tf.target %>%
            dplyr::filter(.data$TF %in% rownames(tf.activity.es))
    }

    if (nrow(pair.tf.target) == 0) {
        stop(
            "pair.tf.target not found in data.",
            "Please check rownames and pair.tf.target provided."
        )
    }
    # reducing object sizes in case we will make it parallel
    if (!missing(tf.activity.es)){
        idx <- rownames(tf.activity.es) %in% pair.tf.target$TF
        tf.activity.es <- tf.activity.es[idx,,drop = FALSE]
        idx <- rownames(exp) %in% pair.tf.target$target
        exp <- exp[idx,,drop = FALSE]
    } else {
        idx <- rownames(exp) %in% c(pair.tf.target$target,pair.tf.target$TF)
        exp <- exp[idx,,drop = FALSE]
    }

    parallel <- register_cores(cores)

    correlation.df <- plyr::adply(
        .data = pair.tf.target,
        .margins = 1,
        .fun = function(pair, exp, tf.activity.es){

            tryCatch({
                target <- exp[pair$target,]
                if (missing(tf.activity.es)) {
                    tf <- exp[pair$TF,]
                } else {
                    tf <- tf.activity.es[pair$TF,]
                }
                suppressWarnings({
                    res <- cor.test(
                        x = target %>% as.numeric,
                        y = tf %>% as.numeric,
                        method = "spearman",
                        exact = TRUE
                    )
                })
                return(
                    tibble(
                        "TF_vs_target_gene_spearman_cor_pvalue" = res$p.value,
                        "TF_vs_target_gene_spearman_cor_estimate" = res$estimate
                    )
                )
            }, error = function(e){
                print(e)
                return(
                    tibble(
                        "TF_vs_target_gene_spearman_cor_pvalue" = NA,
                        "TF_vs_target_gene_spearman_cor_estimate" = NA
                    )
                )
            })
        },
        .progress = "time",
        .parallel = parallel,
        .inform = TRUE,
        tf.activity.es = tf.activity.es,
        exp = exp
    )

    #correlation.df <- na.omit(correlation.df)
    correlation.df$TF_vs_target_gene_spearman_cor_fdr <- p.adjust(
        correlation.df$TF_vs_target_gene_spearman_cor_pvalue,
        method = "fdr"
    )

    return(correlation.df)
}
