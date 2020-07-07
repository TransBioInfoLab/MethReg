#' @title Fits linear models with interaction to triplet data (Target, TF, DNAm)
#' @description Identify DNA methylation (DNAm) changes that work synergistically with TFs in regulating
#' target gene expression, by fitting robust linear model or zero inflated negative binomial model to triplet data.
#' @param triplet Data frame with columns for DNA methylation region (regionID), TF  (TF), and target gene  (target)
#' @param dnam DNA methylation matrix  (columns: samples in the same order as \code{exp} matrix, rows: regions/probes)
#' @param exp A log2 (gene expression count + 1) matrix (columns: samples in the same order as \code{dnam} matrix,
#' rows: genes represented by ensembl IDs (e.g. ENSG00000239415))
#' @param cores Number of CPU cores to be used. Default 1.
#' @return A dataframe with \code{Region, TF, target, TF_symbo, target_symbol, estimates and P-values},
#' after fitting robust linear models or zero-inflated negative binomial models (see Details above).
#'
#' Model 1 (considering DNAm values as a continuous variable) generates \code{pval_met}, \code{pval_rna.tf},
#' \code{pval_met.rna.tf} and \code{estimates_met}, \code{estimates_rna.tf}, \code{estimates_met.rna.tf}.
#' Model 2 (considering DNAm values as a binary variable)
#' generates \code{quant_pval_metGrp}, \code{quant_pval_rna.tf}, \code{quant_pval_metGrp.rna.tf},
#' \code{quant_estimates_metGrp}, \code{quant_estimates_rna.tf}, \code{quant_estimates_metGrp.rna.tf}
#'
#' \code{Model.interaction} indicates which model (robust linear model or zero inflated model)
#' was used to fit Model 1, and \code{Model.quantile} indicates which model(robust linear model or zero
#' inflated model) was used to fit Model 2.
#'
#'@details This function fits the linear model
#'
#' \code{log2(RNA target) ~ log2(TF) + DNAm + log2(TF) * DNAm}
#'
#' to triplet data in two ways:
#'
#' Model 1 : by considering \code{DNAm} as a continuous variable
#'
#' Model 2:  by considering \code{DNAm} as a binary variable - we defined a binary group for
#' DNA methylation values (high = 1, low = 0). That is, samples with the highest
#' DNAm levels (top 25 percent) has high = 1, samples with lowest
#' DNAm levels (bottom 25 pecent) has high = 0. Note that in this
#' implementation, only samples wih DNAm values in the first and last quartiles
#' are considered.
#'
#' We can then identify significant DNAm by TF interactions by selecting those that are significant in both
#' Model 1 and Model 2.
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
#' }
#'
#' To account for confounding effects from covariate variables, first use the \code{get_residuals} function to obtain
#' RNA or DNAm residual values which have covariate effects removed, then fit interaction model. Note that no
#' log2 transformation is needed when \code{interaction_model} is applied to residuals data.
#'
#' @examples
#' data("dna.met.chr21")
#' dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
#' data("gene.exp.chr21")
#' triplet <- data.frame("regionID" = rownames(dna.met.chr21)[1:10],
#'                       "TF" = rownames(gene.exp.chr21)[11:20],
#'                       "target" = rownames(gene.exp.chr21)[1:10])
#' results <- interaction_model(triplet, dna.met.chr21, gene.exp.chr21)
#'
#' # select those that are significant in both models
#' # results <- results[results$`pval_met:rna.tf`< 0.05 & results$`quant_pval_metGrp:rna.tf`< 0.05 ,]
#'
#' @export
#' @importFrom rlang .data
#' @importFrom MASS rlm psi.bisquare
interaction_model <- function(
    triplet,
    dnam,
    exp,
    cores = 1
){

    if(missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
    if(missing(exp)) stop("Please set exp argument with gene expression matrix")

    check_data(dnam, exp)

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
    dnam <- dnam[regions.keep,]

    triplet <- triplet %>% dplyr::filter(
        .data$target %in% rownames(exp) &
            .data$TF %in% rownames(exp) &
            .data$regionID %in% rownames(dnam))

    if(nrow(triplet) == 0){
        stop("We were not able to find the same rows from triple in the data, please check the input.")
    }

    triplet$TF_symbol <- map_ensg_to_symbol(triplet$TF)
    triplet$target_symbol <- map_ensg_to_symbol(triplet$target)

    parallel <- register_cores(cores)

    plyr::adply(
        .data = triplet,
        .margins = 1,
        .fun = function(row.triplet){

            data <- make_df_from_triple(exp, dnam, row.triplet)

            quant.met <-  quantile(data$met,na.rm = TRUE)
            quant.diff <- data.frame("met.q4_minus_q1" = quant.met[4] - quant.met[2])

            low.cutoff <- quant.met[2]
            upper.cutoff <- quant.met[4]

            data.high.low <- data %>% dplyr::filter(.data$met <= low.cutoff | .data$met >= upper.cutoff)
            data.high.low$metGrp <- ifelse(data.high.low$met <= low.cutoff, 0, 1)

            pct.zeros.in.samples <- sum(data$rna.target == 0, na.rm = TRUE) / nrow(data)

            if(pct.zeros.in.samples > 0.25){
                itx.all <- interaction_model_zeroinfl(data)
            } else {
                itx.all <- interaction_model_rlm(data)
            }

            # Quantile model: we will use data.high.low (Q4 and Q1 only)
            pct.zeros.in.quant.samples <- sum(data.high.low$rna.target == 0,na.rm = TRUE) / nrow(data.high.low)
            if(pct.zeros.in.quant.samples > 0.25){
                itx.quant <- interaction_model_quant_zeroinfl(data.high.low)
            } else {
                itx.quant <- interaction_model_quant_rlm(data.high.low)
            }

            # Create output
            interaction_model_output(
                itx.all,
                pct.zeros.in.samples,
                quant.diff,
                itx.quant,
                pct.zeros.in.quant.samples
            )
        },
        .progress = "time",
        .parallel = parallel,
        .inform = TRUE,
        .paropts = list(.errorhandling = 'pass'))
}


interaction_model_no_results <- function(){
    cbind(
        #"Model.interaction" = NA,
        #"met.q4_minus_q1" = NA,
        "quant_pval_metGrp" = NA,
        "quant_pval_rna.tf" = NA,
        "quant_pval_metGrp:rna.tf" = NA,
        "quant_estimate_metGrp" = NA,
        "quant_estimate_rna.tf" = NA,
        "quant_estimate_metGrp:rna.tf" = NA)
    #"Model.quantile" = NA,
    #"% 0 target genes (All samples)" = NA,
    #"% of 0 target genes (Q1 and Q4)" = NA) %>% as.data.frame
}

make_df_from_triple <- function(exp, dnam, row.triplet){
    rna.target <- exp[rownames(exp) == row.triplet$target, , drop = FALSE]
    met <- dnam[rownames(dnam) == as.character(row.triplet$regionID), ]
    rna.tf <- exp[rownames(exp) == row.triplet$TF, , drop = FALSE]

    data <- data.frame(
        rna.target = rna.target %>% as.numeric,
        met = met %>% as.numeric,
        rna.tf = rna.tf %>% as.numeric
    )
    data
}

interaction_model_output <- function(itx.all,
                                     pct.zeros.in.samples,
                                     quant.diff,
                                     itx.quant,
                                     pct.zeros.in.quant.samples){
    if(is.null(itx.quant)) itx.quant <- interaction_model_no_results()
    if(is.null(itx.all)) itx.all <- data.frame(rep(NA,6) %>% t)

    cbind(
        itx.all,
        data.frame(
            "Model interaction" =
                ifelse(pct.zeros.in.samples > 0.25,
                       "Zero-inflated Negative Binomial Model",
                       "Robust Linear Model")
        ),
        quant.diff,
        itx.quant,
        data.frame(
            "Model quantile" =
                ifelse(pct.zeros.in.quant.samples > 0.25,
                       "Zero-inflated Negative Binomial Model",
                       "Robust Linear Model")
        ),
        "% 0 target genes (All samples)" = paste0(round(pct.zeros.in.samples * 100,digits = 2)," %"),
        "% of 0 target genes (Q1 and Q4)" = paste0(round(pct.zeros.in.quant.samples * 100,digits = 2)," %")
    )
}


interaction_model_rlm <- function(data){
    rlm.bisquare <- tryCatch({
        # 2) fit linear model: target RNA ~ DNAm + RNA TF
        rlm (
            rna.target ~ met + rna.tf + rna.tf * met,
            data = data,
            psi = MASS::psi.bisquare,
            maxit = 100) %>% summary %>% coef %>% data.frame
    }, error = function(e){
        # message("Continuous model: ", e)
        return(NULL)
    })

    # if(is.null(rlm.bisquare)) return(interaction_model_no_results())
    if(is.null(rlm.bisquare)) return(NULL)

    degrees.freedom.value <- nrow(data) - 4
    rlm.bisquare$pval <- 2 * (1 - pt( abs(rlm.bisquare$t.value), df = degrees.freedom.value) )

    #mod1 <- lm("rna ~ met + tf + met * tf", data = df)
    all.pval <- rlm.bisquare[-1,4,drop = F] %>% t %>% as.data.frame()
    colnames(all.pval) <- paste0("pval_",colnames(all.pval))

    all.estimate <- rlm.bisquare[-1,1,drop = F] %>% t %>% as.data.frame()
    colnames(all.estimate) <- paste0("estimate_",colnames(all.estimate))
    return(cbind(all.pval, all.estimate))
}

#' @importFrom pscl zeroinfl
interaction_model_zeroinfl <- function(data){
    zinb <- tryCatch({
        pscl::zeroinfl(
            trunc(rna.target) ~ met + rna.tf + rna.tf * met | 1,
            data = data,
            dist = "negbin",
            EM = FALSE) %>% summary %>% coef
    }, error = function(e){
        # message("Continuous model: ", e)
        return(NULL)
    })
    if(is.null(zinb)) return(interaction_model_no_results())

    zinb <- zinb$count %>% data.frame

    all.pval <- zinb[c(-1,-5),4,drop = F] %>% t %>% as.data.frame()
    colnames(all.pval) <- paste0("pval_",colnames(all.pval))

    all.estimate <- zinb[c(-1,-5),1,drop = F] %>% t %>% as.data.frame()
    colnames(all.estimate) <- paste0("estimate_",colnames(all.estimate))
    return(cbind(all.pval, all.estimate))
}

interaction_model_quant_zeroinfl <- function(data){
    zinb.quant <- tryCatch({
        pscl::zeroinfl(
            trunc(rna.target) ~ metGrp + rna.tf + metGrp * rna.tf | 1,
            data = data,
            dist = "negbin",
            EM = FALSE) %>% summary %>% coef
    }, error = function(e){
        # message("Continuous model: ", e)
        return(NULL)
    })
    if(is.null(zinb.quant)) return(interaction_model_no_results())

    zinb.quant <- zinb.quant$count %>% data.frame
    quant.pval <- zinb.quant[c(-1,-5),4,drop = F] %>%
        t %>%
        as.data.frame()
    colnames(quant.pval) <- paste0("quant_pval_",colnames(quant.pval))

    quant.estimate <- zinb.quant[c(-1,-5),1,drop = F] %>%
        t %>%
        as.data.frame()
    colnames(quant.estimate) <- paste0("quant_estimate_",colnames(quant.estimate))

    return(cbind(quant.pval, quant.estimate))
}



interaction_model_quant_rlm <- function(data){
    rlm.bisquare.quant <- tryCatch({
        rlm (
            rna.target ~ metGrp + rna.tf + metGrp * rna.tf,
            data = data,
            psi = MASS::psi.bisquare,
            maxit = 100) %>% summary %>% coef %>% data.frame
    }, error = function(e){
        #message("Binary model: ", e)
        return(NULL)
    })

    if(is.null(rlm.bisquare.quant)) return(NULL)

    # if(is.null(rlm.bisquare.quant)) return(interaction_model_no_results())

    degrees.freedom.value <- nrow(data) - 4
    rlm.bisquare.quant$pval <- 2 * (1 - pt( abs(rlm.bisquare.quant$t.value),
                                            df = degrees.freedom.value) )

    quant.pval <- rlm.bisquare.quant[-1,4,drop = F] %>%
        t %>%
        as.data.frame()
    colnames(quant.pval) <- paste0("quant_pval_",colnames(quant.pval))

    quant.estimate <- rlm.bisquare.quant[-1,1,drop = F] %>%
        t %>%
        as.data.frame()
    colnames(quant.estimate) <- paste0("quant_estimate_",colnames(quant.estimate))
    return(cbind(quant.pval, quant.estimate))
}

