#' @title Fits linear models to triplet data (Target, TF, DNAm) for
#' samples with high DNAm or low DNAm separately, and annotates TF (activator/repressor,
#' methyl-plus/methyl-minus).
#' @description Should be used after fitting \code{interaction_model}, and only
#' for triplet data with significant \code{TF*DNAm} interaction. This analysis
#' examines in more details on how TF activities differ in samples with high DNAm or low DNAm values.
#' @param triplet Data frame with columns for DNA methylation region (regionID), TF  (TF), and target gene  (target)
#' @param dnam DNA methylation matrix  (columns: samples in the same order as \code{exp} matrix, rows: regions/probes)
#' @param exp A log2 (gene expression + 1) matrix (columns: samples in the same order as \code{dnam} matrix,
#' rows: genes represented by ensembl IDs (e.g. ENSG00000239415))
#' @param cores Number of CPU cores to be used. Default 1.
#' @return A dataframe with \code{Region, TF, target, TF_symbol target_symbol}, results for
#' fitting linear models to samples with low methylation (\code{DNAmlow_pval_rna.tf},
#' \code{DNAmlow_estimate_rna.tf}), or samples with high methylation (\code{DNAmhigh_pval_rna.tf},
#' \code{DNAmhigh_pval_rna.tf.1}), annotations for TF (\code{class.TF}) and (\code{class.TF.DNAm}).
#'
#' @details This function fits linear model
#' \code{log2(RNA target) = log2(TF)}
#'
#' to samples with higest DNAm values (top 25 percent) or lowest DNAm values (bottom 25 percent), separately.
#'
#' There are two implementations of these models, depending on whether there are an excessive
#' amount (i.e. more than 25 percent) of samples with zero counts in RNAseq data:
#'
#' \itemize{
#' \item When percent of zeros in RNAseq data is less than
#' 25 percent, robust linear models are implemented using \code{rlm} function from \code{MASS} package. This
#' gives outlier gene expression values reduced weight. We used \code{"psi.bisqure"}
#' option in function \code{rlm} (bisquare weighting,
#' https://stats.idre.ucla.edu/r/dae/robust-regression/).
#'
#' \item When percent of zeros in RNAseq data is more than 25 percent, zero inflated negative binomial models
#' are implemented using \code{zeroinfl} function from \code{pscl} package. This assumes there are
#' two processes that generated zeros (1) one where the counts are always zero
#' (2) another where the count follows a negative binomial distribution.
#'}
#'
#' To account for confounding effects from covariate variables, first use the \code{get_residuals} function to obtain
#' RNA residual values which have covariate effects removed, then fit interaction model. Note that no
#' log2 transformation is needed when \code{interaction_model} is applied to residuals data.
#'
#' This function also provides annotations for TFs. A TF is annotated as \code{activator} if
#' increasing amount of TF (higher TF gene expression) corresponds to increased target gene expression. A TF
#' is annotated as \code{repressor} if increasing amount of TF (higher TF gene expression) corresponds to
#' decrease in target gene expression.
#'
#' In addition, a TF is annotated as \code{M-plus} (methyl-plus) if more TF regulation on gene transcription
#' is observed in samples with high DNAm. That is, DNAm enhances TF regulation on target gene expression.
#' On the other hand, a TF is annotated as \code{M-minus} (methyl-minus) if more TF regulation on gene
#' transcription is observed in samples with low DNAm. That is, DNAm reduces TF regulation
#' on target gene expression.
#'
#'
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
#' \dontrun{
#' data("dna.met.chr21")
#' dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
#' data("gene.exp.chr21")
#' triplet <- data.frame("regionID" = rownames(dna.met.chr21)[1:10],
#'                       "TF" = rownames(gene.exp.chr21)[11:20],
#'                       "target" = rownames(gene.exp.chr21)[1:10])
#' results <- stratified_model(triplet, dna.met.chr21, gene.exp.chr21)
#' }
#' @export
#' @importFrom tibble tibble
#' @importFrom rlang .data
stratified_model <- function(
    triplet,
    dnam,
    exp,
    cores = 1
){

    if(missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
    if(missing(exp)) stop("Please set exp argument with gene expression matrix")

    if(is(dnam,"SummarizedExperiment")) dnam <- assay(dnam)
    if(is(exp,"SummarizedExperiment")) exp <- assay(exp)

    if(!all(grepl("ENSG", rownames(exp)))){
        stop("exp must have the following row names as ENSEMBL IDs (i.e. ENSG00000239415)")
    }

    if(missing(triplet)) stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
    if(!all(c("regionID","TF","target") %in% colnames(triplet))) {
        stop("triplet must have the following columns names: regionID, TF, target")
    }

    message("Removing genes with RNA expression equal to 0 for all samples from triplets")
    exp <- filter_genes_zero_expression_all_samples(exp)

    message("Removing triplet with no DNA methylation information for more than 25% of the samples")
    regions.keep <- (rowSums(is.na(dnam)) < (ncol(dnam) * 0.75)) %>% which %>% names
    dnam <- dnam[regions.keep,,drop = FALSE]

    triplet <- triplet %>% dplyr::filter(
        .data$target %in% rownames(exp) &
            .data$TF %in% rownames(exp) &
            .data$regionID %in% rownames(dnam)
    )

    triplet$TF_symbol <- map_ensg_to_symbol(triplet$TF)
    triplet$target_symbol <- map_ensg_to_symbol(triplet$target)

    # Remove cases where target is also the TF if it exists
    triplet <- triplet %>% dplyr::filter(
            .data$TF != .data$target
    )

    if(nrow(triplet) == 0){
        stop("We were not able to find the same rows from triple in the data, please check the input.")
    }

    parallel <- register_cores(cores)

    out <- plyr::adply(
        .data = triplet,
        .margins = 1,
        .fun = function(row.triplet){

            data <- make_df_from_triple(exp, dnam, row.triplet)

            low.cutoff <- quantile(data$met, na.rm = TRUE)[2]
            upper.cutoff <- quantile(data$met, na.rm = TRUE)[4]

            data.low <- data %>% dplyr::filter(.data$met <= low.cutoff)
            data.high <- data %>% dplyr::filter(.data$met >= upper.cutoff)

            results.low <- stratified_model_aux(data.low,"DNAmlow")
            results.low.pval <- results.low$pval
            results.low.estimate <- results.low$estimate

            results.high <- stratified_model_aux(data.high,"DNAmhigh")
            results.high.pval <- results.high$pval
            results.high.estimate <- results.high$estimate

            classification <- getClassification(results.low.estimate, results.high.estimate)

            tibble::tibble(
                "DNAmlow_pval_rna.tf" = results.low.pval %>% as.numeric(),
                "DNAmlow_estimate_rna.tf" = results.low.estimate %>% as.numeric(),
                "DNAmhigh_pval_rna.tf" = results.high.pval %>% as.numeric(),
                "DNAmhigh_estimate_rna.tf" = results.high.estimate %>% as.numeric(),
                "TF.affinity" = classification$TF.affinity,
                "TF.role" = classification$TF.role
            )
        }, .progress = "time", .parallel = parallel, .inform = TRUE)

    return(out)
}


#' @importFrom MASS rlm psi.bisquare
#' @importFrom stats coef pt
stratified_model_aux <- function(data, prefix = ""){
    pct.zeros.samples <- sum(data$rna.target == 0, na.rm = TRUE) / nrow(data)

    if (pct.zeros.samples > 0.25) {
        results <-  tryCatch({
            pscl::zeroinfl(
            trunc(rna.target) ~ rna.tf | 1,
            data = data,
            dist = "negbin",
            EM = FALSE) %>% summary %>% coef
        }, error = function(e){
            # message("Binary model: ", e)
            return(NULL)
        })

        if(is.null(results)) return(stratified_model_aux_no_results(pct.zeros.samples))

        results <- results$count %>% data.frame

        results.pval <- results["rna.tf","Pr...z..",drop = F] %>%
            t %>%
            as.data.frame()
        colnames(results.pval) <- paste0(prefix,"_pval_",colnames(results.pval))

        results.estimate <- results["rna.tf","Estimate",drop = F] %>%
            t %>%
            as.data.frame()
        colnames(results.estimate) <- paste0(prefix,"_estimate_",colnames(results.estimate))

    } else {

        results <- tryCatch({

            rlm(rna.target ~ rna.tf,
                data = data,
                psi = psi.bisquare,
                maxit = 100) %>% summary %>% coef %>% data.frame

        }, error = function(e){
            # message("Binary model: ", e)
            return(NULL)
        })

        if(is.null(results)) return(stratified_model_aux_no_results(pct.zeros.samples))

        degrees.freedom.value <- nrow(data) - 2
        results$pval <- 2 * (1 - pt( abs(results$t.value), df = degrees.freedom.value) )

        results.pval <- results[-1,4,drop = F] %>% t %>% as.data.frame()
        colnames(results.pval) <- paste0(prefix,"_pval_",colnames(results.pval))

        results.estimate <- results[-1,1,drop = F] %>% t %>% as.data.frame()
        colnames(results.estimate) <- paste0(prefix,"_estimate_",colnames(results.estimate))
    }

    return(
        list("estimate" = results.estimate,
             "pval" = results.pval,
             "Model" = ifelse(pct.zeros.samples > 0.25,
                              "Zero-inflated Negative Binomial Model",
                              "Robust Linear Model"),
             "percet_zero_target_genes" = paste0(round(pct.zeros.samples * 100, digits = 2)," %")
        )
    )
}
stratified_model_aux_no_results <- function(pct.zeros.samples){
    list("estimate" = NA,
         "pval" = NA,
         "Model" = "Robust Linear Model",
         "percet_zero_target_genes" = paste0(round(pct.zeros.samples * 100, digits = 2)," %")
    )

}

getClassification <- function(low.estimate, high.estimate){

    estimate.vector <- c(low.estimate %>% as.numeric, high.estimate %>% as.numeric)

    if(any(is.na(estimate.vector))){
        return(list("TF.affinity" = NA,"TF.role" = NA))
    }

    slope_estimate <- estimate.vector[which.max(abs(estimate.vector))]
    TF.role <- ifelse(slope_estimate > 0, "Activator", "Repressor")

    if(TF.role == "Repressor") {

        if(low.estimate < high.estimate) {
            TF.affinity <- "M-minus"
        } else {
            TF.affinity <- "M-plus"
        }

    } else {

        if(low.estimate < high.estimate) {
            TF.affinity <- "M-plus"
        } else {
            TF.affinity <- "M-minus"
        }
    }
    return(list("TF.affinity" = TF.affinity,"TF.role" = TF.role))
}

