#' @title Calculate interaction linear model Target ~ TF + DNAm
#' @description Evaluate interaction of DNA methylation,
#' TF expression with target gene interaction using a linear model
#' log2(RNA target ) ~ log2(TF) + DNAm + log2(TF) * DNAm
#' @param met DNA methylation matrix
#' @param exp gene EXpression matrix
#' @param TF.target Dataframe with TF and target gene
#' @param DNAm.target Dataframe with DNAm and target gene
#' @return A dataframe with Region, TF, Estimates and Pvalue from linear model
#' @examples
#' \dontrun{
#'  human.tfs <- get_human_tfs()
#'  TF.target <- get_tf_targets_cistrome("COAD_READ*")
#'  DNAm.target <- get_dnam_target_gene("hg38", method = "closest.gene")
#' }
#' @export
interaction_model <- function(met,
                              exp,
                              TF.target,
                              DNAm.target){

    # changing RNA of target to log2
    data$rna <- log2(data[[target]] + 1)

    # transforming TF to log2 TF + 1
    data[,tfs] <- log2(data[,tfs] + 1)


    fmla.mod2 <- as.formula(paste("rna ~ ",  as.character(pairs[idx,2])))
    mod2 <- lm(fmla.mod2, data = data.high)
    return()
}

