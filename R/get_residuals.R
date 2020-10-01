#' @title Get residuals from regression model
#' @description Compute studentized residuals from fitting linear regression models to expression values
#' in a data matrix
#' @param data A matrix or SummarizedExperiment object
#' with samples as columns and features (gene, probes)
#' as rows. Note that expression values should typically be log2(expx + 1)
#' transformed before fitting linear regression models.
#' @param metadata.samples A data frame with samples as rows and columns the covariates.
#' No NA values are allowed, otherwise residual of the corresponding sample will be NA.
#' @param metadata.genes A data frame with genes (covariates) as rows and samples as columns.
#' For each evaluated gene, each column (e.g. CNA) that corresponds to the same gene
#' will be set as a single covariate variable. This can be used to correct copy number alterations for each gene.
#' @param cores Number of CPU cores to be used. Defaults to 1.
#' @return A residuals matrix with samples as columns and features (gene, probes) as rows
#' @details When only \code{metadata.samples} are provided, this function computes
#' residuals for expression values in a data matrix by fitting model
#'
#' \code{features ~ Sample_covariate1 + Sample_covariate2 ... + Sample_covariateN}
#' where \code{N} is the index of the columns in the metadata provided, \code{features} are
#' (typically log transformed) expression values.
#'
#' When the user additionally provide \code{metadata.genes},
#' that is, gene metadata (e.g. gene_covariate = copy number variations/alterations)
#' residuals are computed by fitting the following model:
#'
#' \code{features ~ Sample_covariate1 + Sample_covariate2 ... + Sample_covariateN + gene_covariate}
#'
#' @examples
#' data("gene.exp.chr21.log2")
#'
#' data("clinical")
#' metadata <- clinical[,c( "gender", "sample_type")]
#'
#' cnv <- matrix(
#'    sample(x = c(-2,-1,0,1,2),
#'    size = ncol(gene.exp.chr21.log2) * nrow(gene.exp.chr21.log2),replace = TRUE),
#'    nrow = nrow(gene.exp.chr21.log2),
#'    ncol = ncol(gene.exp.chr21.log2)
#' )
#' rownames(cnv) <- rownames(gene.exp.chr21.log2)
#' colnames(cnv) <- colnames(gene.exp.chr21.log2)
#'
#' gene.exp.residuals <- get_residuals(
#'    data = gene.exp.chr21.log2[1:3,],
#'    metadata.samples = metadata,
#'    metadata.genes = cnv
#' )
#' gene.exp.residuals <- get_residuals(
#'    data = gene.exp.chr21.log2[1:3,],
#'    metadata.samples = metadata,
#'    metadata.genes = cnv[1:2,]
#' )
#' gene.exp.residuals <- get_residuals(
#'    data = gene.exp.chr21.log2[1:3,],
#'    metadata.samples = metadata
#' )
#' @export
#' @importFrom stats rstudent na.exclude na.omit
get_residuals <- function(
    data,
    metadata.samples = NULL,
    metadata.genes = NULL,
    cores = 1
){

    if (missing(data) || is.null(data)) {
        stop("Please data argument with a matrix/SE")
    }

    if (is(data,"SummarizedExperiment")) {
        data <- assay(data)
    }

    if (is.null(metadata.samples) & is.null(metadata.genes)) {
        stop("Please set at least metadata.samples or metadata.genes argument with metadata information")
    }

    if (!all(colnames(data) == rownames(metadata.samples))) {
        stop("data columns names should be the same as metadata row names")
    }

    if (any(is.na(metadata.samples))) {
        message("There are NA's within the metadata, residuals for those samples will be NA.")
    }

    if (!is.null(metadata.genes)) {

        if (ncol(metadata.genes) != ncol(data)) {
            stop("metadata.genes and data should have the number of columns")
        }

        if (!all(colnames(metadata.genes) == colnames(data))) {
            stop("metadata.genes columns names should be the same as data columns names")
        }
        data <- data[rownames(data) %in% rownames(metadata.genes),]

    }
    parallel <- register_cores(cores)


    if (is.null(metadata.genes)) {
        cov_char <- stringr::str_c(colnames(metadata.samples), collapse = " + ")
        form <- stringr::str_c("exp ~ ", cov_char)
        message("Formula used: ", form)
    } else if (is.null(metadata.samples)) {
        cov_char <- ""
        form <- stringr::str_c("exp ~ ", cov_char)
        message("Formula used: ", form, "gene.covariate")
    } else {
        cov_char <- stringr::str_c(colnames(metadata.samples), collapse = " + ")
        form <- stringr::str_c("exp ~ ", cov_char)
        message("Formula used: ", form, " + gene.covariate")
    }

    resid <- plyr::adply(
        .data = data,
        .margins = 1,
        .fun = function(row, genes.names, metadata.genes){
            exp <- row %>% matrix
            colnames(exp) <- "exp"

            if(!is.null(metadata.samples)){
                dat <- cbind(exp, metadata.samples)
            } else {
                dat <- exp %>% as.data.frame()
            }
            dat$exp <- as.numeric(dat$exp)
            gene.name <- genes.names[parent.frame()$i[]]

            if (!is.null(metadata.genes)) {

                if (gene.name %in% rownames(metadata.genes)){
                    df <- data.frame(metadata.genes[gene.name,,drop = FALSE])
                    if (ncol(df) > 1) df <- df %>% t
                    colnames(df) <- gene.name
                    dat <- cbind(dat, df)

                    if(!is.null(metadata.samples)){
                        form <- paste0(form," + ",gene.name)
                    } else {
                        form <- paste0(form, gene.name)
                    }
                    #message("\nFormula used: ",form)
                }
            }
            #print(str(dat))
            fitE <- lm(form, data = dat, na.action = na.exclude)
            rstudent(fitE)
        },
        genes = rownames(data),
        metadata.genes = metadata.genes,
        .progress = "time",
        .inform = TRUE,
        .id = NULL,
        .parallel = parallel)

    rownames(resid) <- rownames(data)
    colnames(resid) <- colnames(data)
    return(resid %>% as.matrix())
}
