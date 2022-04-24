#' @title Wrapper for MethReg functions
#' @description Wrapper for the following MethReg functions:
#'  1) DNAm vs Target gene spearman correlation
#'  2) TF vs Target gene spearman correlation
#'  3) interaction_model
#'  4) stratified model
#' @param triplet Data frame with columns for DNA methylation region (regionID), 
#' TF  (TF), and target gene  (target)
#' @param dnam DNA methylation matrix or SummarizedExperiment object
#' (columns: samples in the same order as \code{exp} matrix, rows: regions/probes)
#' @param exp A matrix or SummarizedExperiment object object
#'  (columns: samples in the same order as \code{dnam},
#' rows: genes represented by ensembl IDs (e.g. ENSG00000239415))
#' @param tf.activity.es A matrix with normalized enrichment scores for each TF across all samples
#' to be used in linear models instead of TF gene expression. See \code{\link{get_tf_ES}}.
#' @param dnam.group.percent.threshold DNA methylation threshold percentage 
#' to define samples  in the low methylated group and high methylated group. 
#' For example, 
#' setting the threshold to 0.3 (30\%) will assign samples with the lowest 30\% 
#' methylation in the low group and the highest 30\% methylation in the high group. 
#' Default is 0.25 (25\%), accepted threshold range (0.0,0.5].
#' @param perform.correlation.analaysis Perform correlation analysis ?
#' @param remove.sig.correlated.tf.exp.dnam  
#' If wilcoxon test of TF expression Q1 and Q4 is significant (pvalue < 0.05),
#' triplet will be removed.
#' @param remove.nonsig.correlated.dnam.target.gene
#' If spearman correlation of target expression and DNAm for all samples 
#' is not significant (pvalue > 0.05), triplet will be removed
#' If wilcoxon test of target expression Q1 and Q4 is not significant (pvalue > 0.05),
#' triplet will be removed.
#' @param remove.nonsig.correlated.dnam.target.gene.threshold.pvalue
#' Cut-off for remove.nonsig.correlated.dnam.target.gene in the spearman test
#' @param remove.nonsig.correlated.dnam.target.gene.threshold.estimate
#' Cut-off for remove.nonsig.correlated.dnam.target.gene in the spearman test
#' @param filter.triplet.by.sig.term Filter significant triplets ?
#' Select triplets if any term is significant
#' 1) interaction (TF x DNAm) p-value < 0.05 or 
#' 2) DNAm p-value < 0.05 or 
#' 3) TF p-value < 0.05 in binary model
#' @param filter.triplet.by.sig.term.pvalue.threshold 
#' P-values/FDR Threshold to filter significant triplets.
#' @param filter.triplet.by.sig.term.using.fdr
#'  Uses FRD instead of p-value when using filter.triplet.by.sig.term.
#' @param multiple.correction.by.stage.wise.analysis 
#' A boolean indicating if stagewise analysis should be performed
#' to correct for multiple comparisons. If set to FALSE then FDR analysis is performed.
#' @param cores Number of CPU cores to be used. Default 1.
#' @param verbose A logical argument indicating if
#' messages output should be provided.
#' @return A dataframe with \code{Region, TF, target, TF_symbo, target_symbol, estimates and P-values},
#' after fitting robust linear models or zero-inflated negative binomial models (see Details above).
#'
#' Model considering DNAm values as a binary variable generates \code{quant_pval_metGrp},
#' \code{quant_pval_rna.tf}, \code{quant_pval_metGrp.rna.tf},
#' \code{quant_estimates_metGrp}, \code{quant_estimates_rna.tf}, \code{quant_estimates_metGrp.rna.tf}.
#'
#' \code{Model.interaction} indicates which model (robust linear model or zero inflated model)
#' was used to fit Model 1, and \code{Model.quantile} indicates which model(robust linear model or zero
#' inflated model) was used to fit Model 2.
#'
#'@details This function fits the linear model
#'
#' \code{log2(RNA target) ~ log2(TF) + DNAm + log2(TF) * DNAm}
#'
#' to triplet data as follow:
#'
#' Model by considering \code{DNAm} as a binary variable - we defined a binary group for
#' DNA methylation values (high = 1, low = 0). That is, samples with the highest
#' DNAm levels (top 25 percent) has high = 1, samples with lowest
#' DNAm levels (bottom 25 percent) has high = 0. Note that in this
#' implementation, only samples with DNAm values in the first and last quartiles
#' are considered.
#'
#' In these models, the term \code{log2(TF)} evaluates direct effect of TF on
#' target gene expression, \code{DNAm} evaluates direct effect of DNAm on target
#' gene expression, and \code{log2(TF)*DNAm} evaluates synergistic effect of DNAm
#' and TF, that is, if TF regulatory activity is modified by DNAm.
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
#' Note that only triplets with TF expression not significantly different in high vs. low
#' methylation groups will be evaluated (Wilcoxon test, p > 0.05).
#'
#' @export
methReg_analysis <- function(
    triplet,
    dnam,
    exp,
    tf.activity.es = NULL,
    dnam.group.percent.threshold = 0.25,
    perform.correlation.analaysis = TRUE,
    remove.nonsig.correlated.dnam.target.gene = FALSE,
    remove.nonsig.correlated.dnam.target.gene.threshold.pvalue = 0.01,
    remove.nonsig.correlated.dnam.target.gene.threshold.estimate = 0.2,
    remove.sig.correlated.tf.exp.dnam = TRUE,
    filter.triplet.by.sig.term = TRUE,
    filter.triplet.by.sig.term.using.fdr = TRUE,
    filter.triplet.by.sig.term.pvalue.threshold = 0.05,
    multiple.correction.by.stage.wise.analysis = TRUE,
    tf.dnam.classifier.pval.threshold = 0.001,
    verbose = FALSE,
    cores = 1
){
  if(perform.correlation.analaysis){
    
    message("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo")
    message("o Performing correlation analysis of Target expression and TF")
    message("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo")
    triplet <- triplet %>% 
      cor_tf_target_gene(
        exp = exp,
        tf.activity.es = tf.activity.es,
        cores = cores
      ) 
    
    message("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo")
    message("o Performing correlation analysis of Target expression and DNAm")
    message("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo")
    triplet <- triplet %>% cor_dnam_target_gene(
        dnam = dnam,
        exp = exp,
        cores = cores,
        filter.results = remove.nonsig.correlated.dnam.target.gene, 
        min.cor.estimate = remove.nonsig.correlated.dnam.target.gene.threshold.estimate,
        min.cor.pval = remove.nonsig.correlated.dnam.target.gene.threshold.pvalue
      ) 
    
  }
  
  message("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo")
  message("o Performing interaction analysis of Target expression, DNAm and TF")
  message("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo")
  
  results <- triplet %>%  interaction_model(
    dnam = dnam,
    exp = exp,
    tf.activity.es = tf.activity.es,
    cores = cores,
    filter.correlated.tf.exp.dnam = remove.sig.correlated.tf.exp.dnam,
    filter.triplet.by.sig.term = filter.triplet.by.sig.term,
    sig.threshold = filter.triplet.by.sig.term.pvalue.threshold,
    fdr = filter.triplet.by.sig.term.using.fdr,
    stage.wise.analysis = multiple.correction.by.stage.wise.analysis,
    dnam.group.threshold = dnam.group.percent.threshold,
    filter.correlated.target.exp.dnam = remove.nonsig.correlated.dnam.target.gene
  )  
  
  if(nrow(results) == 0){
    stop("No significant results from the interaction model analysis after filtering.")
  }
  message("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo")
  message("o Classification of significant interaction of Target expression, DNAm and TF")
  message("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo")
  results %>%  stratified_model(
    dnam = dnam,
    exp = exp,
    tf.activity.es = tf.activity.es,
    cores = cores,
    tf.dnam.classifier.pval.thld = tf.dnam.classifier.pval.threshold,
    dnam.group.threshold = dnam.group.percent.threshold
  )
}