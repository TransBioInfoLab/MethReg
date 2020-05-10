#' @title Get Residuals
#' @description Get Residuals using a data matrix and a metadata data frame as
#' follows: features ~ covariate1 + covariate2 ... + covariateN
#' where n is the index of the columns of the metadata provided
#' @param data.matrix A matrix with samples as columns and features (gene, probes)
#' as rows.
#' @param metadata A data frame with samples as rows and columns the covariates.
#' Observation: No NA values are allowed, otherwise residual of the sample will be NA.
#' @return A residuals matrix with samples as columns and features (gene, probes)
#' @examples
#' data("gene.exp.chr21")
#' data("clinical")
#' metadata <- clinical[,c( "gender", "sample_type")]
#' gene.exp.residuals <- get_residuals(gene.exp.chr21[1:3,], metadata)
#' @export
#' @importFrom stats residuals na.exclude na.omit
get_residuals <- function(data.matrix,
                          metadata
){

    if(missing(data.matrix)) stop("Please data.matrix dnam argument with a matrix")
    if(missing(metadata)) stop("Please set metadata argument with metadata information")
    if(!all(colnames(data.matrix) == rownames(metadata))) {
        stop("data.matrix columns names should be the same as metadata row names")
    }

    cov_char <- stringr::str_c(colnames(metadata), collapse = " + ")
    form <- stringr::str_c("val ~ ", cov_char)
    message("Formula used: ", form)
    resid <- plyr::adply(
        .data = data.matrix,
        .margins = 1,
        .fun = function(row){
            val <- row %>% matrix
            colnames(val) <- "val"
            dat <- cbind(val, metadata)
            dat$val <- as.numeric(dat$val)
            fitE <- lm(form, data = dat, na.action = na.exclude)
            residuals(fitE)
        }, .progress = "time",.inform = TRUE,.id = NULL)
    rownames(resid) <- rownames(data.matrix)
    colnames(resid) <- colnames(data.matrix)
    return(resid)
}
