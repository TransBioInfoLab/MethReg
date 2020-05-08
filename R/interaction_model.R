#' @title Calculate interaction linear model Target ~ TF + DNAm
#' @description Evaluate interaction of DNA methylation,
#' TF expression with target gene interaction using a linear model
#' \deqn{log2(RNA target) ~ log2(TF) + DNAm + log2(TF) * DNAm}
#'
#' To consider covariates, RNA can also be the residuals.
#' \deqn{log2(RNA target residuals) ~ log2(TF residual) + DNAm + log2(TF residual) * DNAm}
#'
#' @param triplet Dataframe with region (column name: regionID),
#' TF  (column name: TF),  and target gene  (column name: target),
#' @param dnam DNA methylation matrix  (columns: samples same order as met, rows: regions/probes)
#' @param exp A log2 (gene expression + 1) matrix with samples as columns in the same order as met
#' and genes as rows represented by ensembl IDs ENSG00000239415)
#' @return A dataframe with Region, TF, Estimates and P-value from linear model
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

            data <- data.frame(
                rna.target = rna.target %>% as.numeric,
                met = met %>% as.numeric,
                rna.tf = rna.tf %>% as.numeric
            )

            # 2) fit linear model: target RNA ~ DNAm + RNA TF
            results <- lm (
                rna.target ~ met + rna.tf + rna.tf * met,
                data = data
            )

            results.pval <- summary(results)$coefficients[-1, 4, drop = F] %>% t %>% as.data.frame()
            colnames(results.pval) <- paste0("pval_", colnames(results.pval))

            results.estimate <- summary(results)$coefficients[-1, 1, drop = F] %>% t %>% as.data.frame()
            colnames(results.estimate) <- paste0("estimate_", colnames(results.estimate))

            low.cutoff <- quantile(data$met, na.rm = TRUE)[2]
            upper.cutoff <- quantile(data$met, na.rm = TRUE)[4]

            data.low <- data %>% dplyr::filter(met <= low.cutoff)
            data.high <- data %>% dplyr::filter(met >= upper.cutoff)

            results.low <- lm (
                rna.target ~ met + rna.tf + rna.tf * met,
                data = data.low
            )
            results.high <- lm (
                rna.target ~ met + rna.tf + rna.tf * met,
                data = data.high
            )

            results.low.pval <- summary(results.low)$coefficients[-1,4,drop = F] %>% t %>% as.data.frame()
            colnames(results.low.pval) <- stringr::str_c("DNAmlow_pval_", colnames(results.low.pval))

            results.low.estimate <- summary(results.low)$coefficients[-1,1,drop = F] %>% t %>% as.data.frame()
            colnames(results.low.estimate) <- stringr::str_c("DNAmlow_estimate_", colnames(results.low.estimate))

            results.high.pval <- summary(results.high)$coefficients[-1,4,drop = F] %>% t %>% as.data.frame()
            colnames(results.high.pval) <- stringr::str_c("DNAmhigh_pval_", colnames(results.high.pval))

            results.high.estimate <- summary(results.high)$coefficients[-1,1,drop = F] %>% t %>% as.data.frame()
            colnames(results.high.estimate) <- stringr::str_c("DNAmhigh_estimate_", colnames(results.high.estimate))

            out <- cbind(results.pval, results.estimate,
                         results.low.pval, results.low.estimate,
                         results.high.pval, results.high.estimate
            ) %>% data.frame()

        }, .progress = "time")

    return(out)
}

