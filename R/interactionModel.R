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
#' @return A dataframe with Region, TF, Estimates and Pvalue from linear model
#' @examples
#' \dontrun{
#'  human.tfs <- get_human_tfs()
#'  TF.target <- get_tf_targets_cistrome("COAD_READ*")
#'  DNAm.target <- get_dnam_target_gene("hg38", method = "closest.gene")
#' }
#' @export
interaction_model <- function(triplet,
                              dnam,
                              exp
){

    if(missing(met)) stop("Please set met argument with DNA methylation matrix")
    if(missing(exp)) stop("Please set exp argument with gene expression matrix")
    if(missing(triplet)) stop("Please set triplet argument with interactors (region,TF, target gene) data frame")

    rna.target <- exp[rownames(exp) == triplet$target, , drop = FALSE]
    met <- dnam[rownames(dnam) == as.character(triplet$regionID), ]
    rna.tf <- exp[rownames(exp) == triplet$TF, , drop = FALSE]

    data <- data.frame(
        rna.target = log2(rna.target + 1) %>% as.numeric,
        met = met %>% as.numeric,
        rna.tf = log2(rna.tf + 1) %>% as.numeric
    )

    # 2) fit linear model: target RNA ~ DNAm + RNA TF
    results <- lm (
        rna.target ~ met + rna.target + rna.target * met,
        data = data
    )

    results.pval <- summary(results.cases)$coefficients[-1, 4, drop = F] %>% t %>% as.data.frame()
    colnames(results.pval) <- paste0("cases_pval_", colnames(results.pval))

    results.estimate <- summary(results.cases)$coefficients[-1, 1, drop = F] %>% t %>% as.data.frame()
    colnames(results.estimate) <- paste0("cases_estimate_", colnames(results.estimate))

    out <- cbind(results.pval, results.estimate) %>% data.frame()

    return(out)
}

