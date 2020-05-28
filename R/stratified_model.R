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
#' @param cores Number of CPU cores to be used. Default 1.
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
stratified_model <- function(
    triplet,
    dnam,
    exp,
    cores = 1
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
    # message("Removing triplet with RNA expression equal to 0 for more than 25% of the samples")
    exp <- filter_genes_zero_expression(exp,max.samples.percentage = 0)

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

    parallel <- register_cores(cores)

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


            results.low <- stratified_model_aux(data.low,"DNAmlow")
            results.low.pval <- results.low$pval
            results.low.estimate <- results.low$estimate

            results.high <- stratified_model_aux(data.high,"DNAmhigh")
            results.high.pval <- results.high$pval
            results.high.estimate <- results.high$pval

            class <- getClassification(results.low.estimate, results.high.estimate)

            out <- cbind(
                results.low.pval,
                results.low.estimate,
                results.high.pval,
                results.high.estimate,
                class$TF,
                class$DNAm
            ) %>% data.frame()

        }, .progress = "time", .parallel = parallel)

    return(out)
}

stratified_model_aux <- function(data, prefix = ""){
    pct.zeros.samples <- sum(data$rna.target == 0) / nrow(data)

    if (pct.zeros.samples > 0.25) {
        print(data)
        results <- pscl::zeroinfl(
            trunc(rna.target) ~ rna.tf | 1,
            data = data,
            dist = "negbin",
            EM = FALSE) %>% summary %>% coef
        results <- results$count %>% data.frame
        results.pval <- results[c(-1,-5),4,drop = F] %>%
            t %>%
            as.data.frame()
        colnames(results.pval) <- paste0(prefix,"_pval_",colnames(results.pval))

        results.estimate <- results[c(-1,-5),1,drop = F] %>%
            t %>%
            as.data.frame()
        colnames(results.estimate) <- paste0(prefix,"_estimate_",colnames(results.estimate))
    } else {
        results <- MASS::rlm(
            rna.target ~ rna.tf,
            data = data,
            psi = MASS::psi.bisquare,
            maxit = 100) %>% summary %>% coef %>% data.frame

        degrees.freedom.value <- nrow(data) - 2
        results$pval <- 2 * (1 - pt( abs(results$t.value), df = degrees.freedom.value) )

        results.pval <- results[-1,4,drop = F] %>% t %>% as.data.frame()
        colnames(results.pval) <- paste0(prefix,"_pval_",colnames(results.pval))

        results.estimate <- results[-1,1,drop = F] %>% t %>% as.data.frame()
        colnames(results.estimate) <- paste0(prefix,"_estimate_",colnames(results.estimate))
    }

    return(
        list("estimate" = results.estimate,
             "pval" = results.pval)
    )
}

getClassification <- function(low.estimate, high.estimate){

    estimate.vector <- c(low.estimate %>% as.numeric, high.estimate %>% as.numeric)
    slope_estimate <- estimate.vector[which.max(abs(estimate.vector))]
    TFclass <- ifelse(slope_estimate > 0, "Activator", "Repressor")

    if(TFclass == "Repressor") {
        if(low.estimate < high.estimate) {
            DNAmClass <- "M-minus"
        } else {
            DNAmClass <- "M-plus"
        }
    } else {
        if(low.estimate < high.estimate) {
            DNAmClass <- "M-plus"
        } else {
            DNAmClass <- "M-minus"
        }
    }
    return(list("DNAm" = DNAmClass,"TF" = TFclass))
}

