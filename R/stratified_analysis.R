#' @title Fits linear model to triplet data (Target, TF, DNAm) separately for
#' samples with DNAm high or low groups.
#' @description Should be used only for triplet data with significant
#'  \code{TF*DNAm} interaction from fitting models in \code{interaction_model}.
#' These models can be used to examine how TF activities differ in
#' samples with high DNAm or low DNAm values.
#' @param triplet Data frame with columns for DNA methylation region (regionID), TF  (TF), and target gene  (target)
#' @param dnam DNA methylation matrix  (columns: samples in the same order as \code{exp} matrix, rows: regions/probes)
#' @param exp A log2 (gene expression + 1) matrix (columns: samples in the same order as \code{dnam} matrix,
#' rows: genes represented by ensembl IDs (e.g. ENSG00000239415))
#' @return A dataframe with Region, TF, Estimates and P-value from linear model
#' @details This function fits linear model
#' \code{log2(RNA target) ~ log2(TF)}
#'
#' to samples with higest DNAm values (top 25 percent) and lowest DNAm values (bottom 25 percent), separately.
#'
#' To account for confounding effects from covariate variables, first use the \code{get_residuals} function to obtain
#' RNA residual values which have covariate effects removed, then fit interaction model. Note that no
#' log2 transformation is needed when \code{interaction_model} is applied to residuals data.
#'
#' @examples
#' data("dna.met.chr21")
#' dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
#' data("gene.exp.chr21")
#' triplet <- data.frame("regionID" = rownames(dna.met.chr21)[1:10],
#'                       "TF" = rownames(gene.exp.chr21)[11:20],
#'                       "target" = rownames(gene.exp.chr21)[1:10])
#' results <- stratified_model(triplet, dna.met.chr21, gene.exp.chr21)
#' @export
#' @importFrom rlang .data
stratified_model <- function(triplet,
                              dnam,
                              exp
){

    if(missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
    if(missing(exp)) stop("Please set exp argument with gene expression matrix")

    if(!all(grepl("ENSG", rownames(exp)))){
        stop("exp must have the following row names as ENSEMBL IDs (i.e. ENSG00000239415)")
    }

    if(missing(triplet)) stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
    if(!all(c("regionID","TF","target") %in% colnames(triplet))) {
        stop("triplet must have the following columns names: regionID, TF, target")
    }

    # remove triplet with RNA expression equal to 0 for more than 25% of the samples
    # remove triplet with RNA expression equal to 0 for more than 25% of the samples
    message("Removing triplet with RNA expression equal to 0 for more than 25% of the samples")
    genes.keep <- (rowSums(exp == 0)/ncol(exp) < 0.25) %>% which %>% names
    exp <- exp[genes.keep,]

    message("Removing triplet with no DNA methylation information for more than 25% of the samples")
    regions.keep <- (rowSums(is.na(dnam)) < (ncol(dnam) * 0.75)) %>% which %>% names
    dnam <- dnam[regions.keep,]

    triplet <- triplet %>% dplyr::filter(
        .data$target %in% rownames(exp) &
            .data$TF %in% rownames(exp) &
            .data$regionID %in% rownames(dnam))

    triplet$TF_symbol <- map_ensg_to_symbol(triplet$TF)
    triplet$target_symbol <- map_ensg_to_symbol(triplet$target)


    if(nrow(triplet) == 0){
        stop("We were not able to find the same rows from triple in the data, please check the input.")
    }

    out <- plyr::adply(
        .data = triplet,
        .margins = 1,
        .fun = function(row.triplet){
            rna.target <- exp[rownames(exp) == row.triplet$target, , drop = FALSE]
            met <- dnam[rownames(dnam) == as.character(row.triplet$regionID), ]
            rna.tf <- exp[rownames(exp) == row.triplet$TF, , drop = FALSE]

            data <- data.frame(
                rna.target = rna.target %>% as.numeric,
                met = met %>% as.numeric,
                rna.tf = rna.tf %>% as.numeric
            )

            low.cutoff <- quantile(data$met, na.rm = TRUE)[2]
            upper.cutoff <- quantile(data$met, na.rm = TRUE)[4]

            data.low <- data %>% dplyr::filter(met <= low.cutoff)
            data.high <- data %>% dplyr::filter(met >= upper.cutoff)

            results.low <- MASS::rlm(
                rna.target ~ rna.tf,
                data = data.low,
                psi = MASS::psi.bisquare,
                maxit = 100) %>% summary %>% coef %>% data.frame

            degrees.freedom.value <- nrow(data.low) - 2
            results.low$pval <- 2 * (1 - pt( abs(results.low$t.value), df = degrees.freedom.value) )

            results.low.pval <- results.low[-1,4,drop = F] %>% t %>% as.data.frame()
            colnames(results.low.pval) <- paste0("DNAmlow_pval_",colnames(results.low.pval))

            results.low.estimate <- results.low[-1,1,drop = F] %>% t %>% as.data.frame()
            colnames(results.low.estimate) <- paste0("DNAmlow_estimate_",colnames(results.low.estimate))

            results.high <- MASS::rlm(
                rna.target ~ rna.tf,
                data = data.high,
                psi = MASS::psi.bisquare,
                maxit = 100) %>% summary %>% coef %>% data.frame

            degrees.freedom.value <- nrow(data.high) - 2
            results.high$pval <- 2 * (1 - pt( abs(results.high$t.value), df = degrees.freedom.value) )

            results.high.pval <- results.high[-1,4,drop = F] %>% t %>% as.data.frame()
            colnames(results.high.pval) <- paste0("DNAmhigh_pval_",colnames(results.high.pval))

            results.high.estimate <- results.high[-1,1,drop = F] %>% t %>% as.data.frame()
            colnames(results.high.estimate) <- paste0("DNAmhigh_estimate_",colnames(results.high.estimate))

            out <- cbind(results.low.pval,
                         results.low.estimate,
                         results.high.pval,
                         results.high.estimate
            ) %>% data.frame()

        }, .progress = "time")

    return(out)
}

