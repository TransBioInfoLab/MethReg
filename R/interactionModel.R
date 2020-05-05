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
#' @param exp gene expression matrix (columns: samples same order as met, rows: genes)
#' @return A dataframe with Region, TF, Estimates and P-value from linear model
#' @examples
#' triplet <- data.frame("regionID" = paste0("region_",1:9),
#'                       "TF" = paste0("ESNG",10),
#'                       "target" = paste0("ESNG",1:9))
#' dnam <- runif(200) %>% matrix(10)
#' colnames(dnam) <- paste0("Sample_",1:20)
#' rownames(dnam) <- paste0("region_",1:10)
#' exp <- rexp(200) %>% matrix(10)
#' colnames(exp) <- paste0("Sample_",1:20)
#' rownames(exp) <- paste0("ESNG",1:10)
#' results <- interaction_model(triplet, dnam, exp)
#' \dontrun{
#' # TODO create a more real example
#'  human.tfs <- get_human_tfs()
#'  TF.target <- get_tf_targets_cistrome("COAD_READ*")
#'  DNAm.target <- get_region_target_gene("hg38", method = "closest.gene")
#' }
#' @export
interaction_model <- function(triplet,
                              dnam,
                              exp
){

    if(missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
    if(missing(exp)) stop("Please set exp argument with gene expression matrix")
    if(missing(triplet)) stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
    if(!all(c("regionID","TF","target") %in% colnames(triplet))) {
        stop("triplet must have the following columns names: regionID, TF, target")
    }

    triplet <- triplet %>% dplyr::filter(target %in% rownames(exp) &
                                             TF %in% rownames(exp) &
                                             regionID %in% rownames(dnam))
    if(nrow(triplet) == 0){
        stop("We were not able to find the same rows from triple in the data, please check the input.")
    }
    out <- plyr::adply(.data = triplet,.margins = 1,.fun = function(row.triplet){

        rna.target <- exp[rownames(exp) == row.triplet$target, , drop = FALSE]
        met <- dnam[rownames(dnam) == as.character(row.triplet$regionID), ]
        rna.tf <- exp[rownames(exp) == row.triplet$TF, , drop = FALSE]

        data <- data.frame(
            rna.target = log2(rna.target + 1) %>% as.numeric,
            met = met %>% as.numeric,
            rna.tf = log2(rna.tf + 1) %>% as.numeric
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

        out <- cbind(results.pval, results.estimate) %>% data.frame()

    }, .progress = "time")

    return(out)
}

