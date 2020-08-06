#' @title Get studentized residual from regression model
#' @description Get Residuals using a data matrix and a sample metadata data frame as
#' follows: features ~ Sample_covariate1 + Sample_covariate2 ... + Sample_covariateN
#' where n is the index of the columns of the metadata provided
#' In case the user also provide a gene metadata (i.e copy number variation/alteration)
#' get Residuals using:
#' features ~ Sample_covariate1 + Sample_covariate2 ... + Sample_covariateN + gene_covariate
#' @param data.matrix A matrix with samples as columns and features (gene, probes)
#' as rows.
#' @param metadata.samples A data frame with samples as rows and columns the covariates.
#' Observation: No NA values are allowed, otherwise residual of the sample will be NA.
#' @param metadata.genes A data frame with genes (covariates) as rows and samples on columns.
#' For each evaluated gene the columns will be set as a single covariates.
#' This can be used to correct each gene for copy number alterations.
#' @param cores Number of CPU cores to be used. Default 1.
#' @return A residuals matrix with samples as columns and features (gene, probes)
#' @examples
#' data("gene.exp.chr21")
#' data("clinical")
#' metadata <- clinical[,c( "gender", "sample_type")]
#' cnv <- matrix(
#'    sample(x = c(-2,-1,0,1,2),
#'    size = ncol(gene.exp.chr21) * nrow(gene.exp.chr21),replace = TRUE),
#'    nrow = nrow(gene.exp.chr21),
#'    ncol = ncol(gene.exp.chr21)
#' )
#' rownames(cnv) <- rownames(gene.exp.chr21)
#' colnames(cnv) <- colnames(gene.exp.chr21)
#' gene.exp.residuals <- get_residuals(
#'    data.matrix = gene.exp.chr21[1:3,],
#'    metadata.samples = metadata,
#'    metadata.genes = cnv
#' )
#' gene.exp.residuals <- get_residuals(
#'    data.matrix = gene.exp.chr21[1:3,],
#'    metadata.samples = metadata,
#'    metadata.genes = cnv[1:2,]
#' )
#' gene.exp.residuals <- get_residuals(
#'    data.matrix = gene.exp.chr21[1:3,],
#'    metadata.samples = metadata
#' )
#' @export
#' @importFrom stats rstudent na.exclude na.omit
get_residuals <- function(
    data.matrix,
    metadata.samples = NULL,
    metadata.genes = NULL,
    cores = 1
){

    if(is.null(data.matrix)) stop("Please data.matrix dnam argument with a matrix")
    if(is.null(metadata.samples)) stop("Please set metadata argument with metadata information")
    if(!all(colnames(data.matrix) == rownames(metadata.samples))) {
        stop("data.matrix columns names should be the same as metadata row names")
    }

    if(any(is.na(metadata.samples))){
        message("There are NA's within the metadata, residuals for those samples will be NA.")
    }

    if(!missing(metadata.genes)) {

        if(ncol(metadata.genes) != ncol(data.matrix)){
            stop("metadata.genes and data.matrix should have the number of columns")
        }

        if(!all(colnames(metadata.genes) == colnames(data.matrix))){
            stop("metadata.genes columns names should be the same as data.matrix columns names")
        }
    }
    parallel <- register_cores(cores)

    cov_char <- stringr::str_c(colnames(metadata.samples), collapse = " + ")
    form <- stringr::str_c("val ~ ", cov_char)

    if(missing(metadata.genes)) {
        message("Formula used: ", form)
    } else {
        message("Formula used: ", form, " + gene.covariate")
    }
    resid <- plyr::adply(
        .data = data.matrix,
        .margins = 1,
        .fun = function(row, genes.names,metadata.genes){
            val <- row %>% matrix
            colnames(val) <- "val"
            dat <- cbind(val, metadata.samples)
            dat$val <- as.numeric(dat$val)
            gene.name <- genes.names[parent.frame()$i[]]
            if(!is.null(metadata.genes)) {
                if(gene.name %in% rownames(metadata.genes)){
                    df <- data.frame(metadata.genes[gene.name,])
                    colnames(df) <- gene.name
                    dat <- cbind(dat, df)
                    form <- paste0(form," + ",gene.name)
                    #message("\nFormula used: ",form)
                }
            }
            fitE <- lm(form, data = dat, na.action = na.exclude)
            rstudent(fitE)
        },
        genes = rownames(data.matrix),
        metadata.genes = metadata.genes,
        .progress = "time",
        .inform = TRUE,
        .id = NULL,
        .parallel = parallel)

    rownames(resid) <- rownames(data.matrix)
    colnames(resid) <- colnames(data.matrix)
    return(resid %>% as.matrix())
}
