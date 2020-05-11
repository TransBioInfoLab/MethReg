#' @title Fits robust linear models with interaction to triplet data (Target, TF, DNAm)
#' @description Identify DNA methylation (DNAm) changes that work synergistically with TFs in regulating
#' target gene expression.
#' @param triplet Data frame with columns for DNA methylation region (regionID), TF  (TF), and target gene  (target)
#' @param dnam DNA methylation matrix  (columns: samples in the same order as \code{exp} matrix, rows: regions/probes)
#' @param exp A log2 (gene expression + 1) matrix (columns: samples in the same order as \code{dnam} matrix,
#' rows: genes represented by ensembl IDs (e.g. ENSG00000239415))
#' @return A dataframe with Region, TF, Estimates and P-values, after fitting robust linear
#' models using two approaches(see Details above).
#'
#' Approch 1 (considering DNAm values as a continuous variable) generates \code{pval_met, pval_rna.tf, pval_met.rna.tf
#' and estimates_met, estimates_rna.tf, estimates_met.rna.tf}. Approach 2 (considering DNAm values as a binary variable)
#' generates \code{quant_pval_metGrp, quant_pval_rna.tf, quant_pval_metGrp.rna.tf,
#' quant_estimates_metGrp, quant_estimates_rna.tf, quant_estimates_metGrp.rna.tf}
#'
#'
#'
#'
#' @details This function fits robust linear model
#'
#' \code{log2(RNA target) ~ log2(TF) + DNAm + log2(TF) * DNAm}
#'
#' The robust linear model, implemented using \code{rlm} function from \code{MASS},
#' gives outlier gene expression values reduced weight. We used \code{"psi.bisqure"}
#' option in function \code{rlm} (bisquare weighting,
#' https://stats.idre.ucla.edu/r/dae/robust-regression/).
#'
#' The above linear model were fit to triplet data frame in two ways:
#'
#' Model 1 : by considering DNAm as a continuous variable
#'
#' Model 2:  by considering DNAm as a binary variable - we defined a binary group for
#' DNA methylation values (high = 1, low = 0). That is, samples with the highest
#' DNAm levels (top 25 percent) has high = 1, samples with lowest
#' DNAm levels (bottom 25 pecent) has high = 0. Note that in this
#' implementation, only samples wih DNAm values in the first and last quartiles
#' are considered.
#'
#' We can identify significant DNAm by TF interactions by selecting those that are significant in both models.
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
interaction_model <- function(triplet,
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
    genes.keep <- (rowSums(exp == 0) < 0.25) %>% which %>% names
    exp <- exp[genes.keep,]

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

            # 2) fit linear model: target RNA ~ DNAm + RNA TF
            results <- MASS::rlm(
                rna.target ~ met + rna.tf + rna.tf * met,
                data = data
            )

            rlm.bisquare <- rlm (
                rna.target ~ met + rna.tf + rna.tf * met,
                data = data,
                psi = MASS::psi.bisquare,
                maxit = 100) %>% summary %>% coef %>% data.frame


            degrees.freedom.value <- nrow(data) - 4
            rlm.bisquare$pval <- 2 * (1 - pt( abs(rlm.bisquare$t.value), df = degrees.freedom.value) )

            #mod1 <- lm("rna ~ met + tf + met * tf", data = df)
            all.pval <- rlm.bisquare[-1,4,drop = F] %>% t %>% as.data.frame()
            colnames(all.pval) <- paste0("pval_",colnames(all.pval))

            all.estimate <- rlm.bisquare[-1,1,drop = F] %>% t %>% as.data.frame()
            colnames(all.estimate) <- paste0("estimates_",colnames(all.estimate))


            quant.met <-  quantile(data$met,na.rm = TRUE)
            quant.diff <- data.frame("met.q4_minus_q1" = quant.met[4] - quant.met[2])

            low.cutoff <- quant.met[2]
            upper.cutoff <- quant.met[4]

            data.high.low <- data %>% dplyr::filter(met <= low.cutoff | met >= upper.cutoff)
            data.high.low$metGrp <- ifelse(data.high.low$met <= low.cutoff,0,1)

            rlm.bisquare.quant <-   rlm (
                rna.target ~ metGrp + rna.tf + metGrp * rna.tf,
                data = data.high.low,
                psi = MASS::psi.bisquare,
                maxit = 100) %>% summary %>% coef %>% data.frame

            degrees.freedom.value <- nrow(data.high.low) - 4
            rlm.bisquare.quant$pval <- 2 * (1 - pt( abs(rlm.bisquare.quant$t.value), df = degrees.freedom.value) )

            quant.pval <- rlm.bisquare.quant[-1,4,drop = F] %>%
                t %>%
                as.data.frame()
            colnames(quant.pval) <- paste0("quant_pval_",colnames(quant.pval))

            quant.estimate <- rlm.bisquare.quant[-1,1,drop = F] %>%
                t %>%
                as.data.frame()
            colnames(quant.estimate) <- paste0("quant_estimates_",colnames(quant.estimate))

            out <- cbind(all.pval,
                         all.estimate,
                         quant.diff,
                         quant.pval, quant.estimate
            )
            out
        }, .progress = "time")

    return(out)
}

