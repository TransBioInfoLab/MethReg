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
#' (1) by considering DNAm as a continuous variable
#'
#' (2) by considering DNAm as a binary variable - we defined a binary group for
#' DNA methylation values (high = 1, low = 0). That is, samples with the highest
#' DNAm levels (top 25 percent) has high = 1, samples with lowest
#' DNAm levels (bottom 25 pecent) has high = 0. Note that in this
#' implementation, only samples wih DNAm values in the first and last quartiles
#' are considered.
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
#' @export
#' @importFrom rlang .data
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

            df <- data.frame(
                rna.target = rna.target %>% as.numeric,
                met = met %>% as.numeric,
                rna.tf = rna.tf %>% as.numeric
            )

            # filter rows with excessive zeros
            prop.zeros <- sum(rna.target == 0)/length(rna.target)

            # fit linear models: target RNA ~ met + TF + met*TF
            # use met as continous variable
            if (prop.zeros < 0.25){

                rlm.bisquare <- data.frame (
                    coef (
                        summary(
                            MASS::rlm ( rna.target ~ met + rna.tf + met * rna.tf,
                                  data = df,
                                  psi = psi.bisquare, maxit = 100)
                        )
                    )
                )
                df.value <- nrow(df) - 4
                rlm.bisquare$pval <- 2 * (1 - pt( abs(rlm.bisquare$t.value), df = df.value) )

                #mod1 <- lm("rna ~ met + tf + met * tf", data = df)
                results <- rlm.bisquare[-1,4,drop = F] %>% t %>% as.data.frame()
                colnames(results) <- paste0("pval_",colnames(results))
                estimates <- rlm.bisquare[-1,1,drop = F] %>% t %>% as.data.frame()
                colnames(estimates) <- paste0("estimates_",colnames(estimates))

                # fit linear models: target RNA ~ met + TF + met*TF
                # use met as binary variable

                quant.met <-  quantile(df$met, na.rm = TRUE)
                quant.tf <-  quantile(df$rna.tf,na.rm = TRUE)
                quant.diff <- data.frame("met.q4_minus_q1" = quant.met[4] - quant.met[2],
                                         "tf.q4_minus_q1" = quant.tf[4] - quant.tf[2])


                # Keep samples in low andf high quintile
                df2 <- df[df$met < quant.met[2] | df$met > quant.met[4],]
                # Low quintile group is 0, high is 1
                # df2$metGrp <- NA
                df2$metGrp <- ifelse(df2$met <= quant.met[2],0,1)

                rlm.bisquare.quant <- data.frame (
                    coef (
                        summary(
                            MASS::rlm ( rna.target ~ metGrp + rna.tf + metGrp * rna.tf,
                                  data = df2,
                                  psi = psi.bisquare, maxit = 100)
                        )
                    )
                )
                df.value <- nrow(df2) - 4
                rlm.bisquare.quant$pval <- 2 * (1 - pt( abs(rlm.bisquare.quant$t.value), df = df.value) )
                results.quant <- rlm.bisquare.quant[-1,4,drop = F] %>%
                    t %>%
                    as.data.frame()
                colnames(results.quant) <- paste0("quant_pval_",colnames(results.quant))
                estimates.quant <- rlm.bisquare.quant[-1,1,drop = F] %>%
                    t %>%
                    as.data.frame()
                colnames(estimates.quant) <- paste0("quant_estimates_",colnames(estimates.quant))

                out <- data.frame(cbind(results,
                                 estimates,
                                 quant.diff,
                                 results.quant,
                                 estimates.quant
                                 ),
                           row.names = NULL,stringsAsFactors = FALSE)

                out
        }
            }, .progress = "time")

    return(out)
}

