#' To provide functional annotations for differentially methylated regions (DMRs)
#' and differentially methylated CpG sites (DMS),`coMethTF` performs
#' integrative analyses using matched DNA methylation and gene expression
#' along with Transcription Factor Binding Sites (TFBS) data.
#' coMethTF evaluates, prioritizes and annotates DNA methylation regions
#' (or sites) with high regulatory potential that works synergistically with
#' TFs to regulate target gene expressions, without any additional ChIP-seq data.
#' @docType package
#' @name coMethTF
NULL


#' A DNA methylation matrix for 50 samples (only chrmossome 21)
#' @docType data
#' @name dna.met.chr21
#' @format A matrix: 50 samples
NULL

#' A gene expression matrix for 50 samples (only chrmossome 21)
#' @docType data
#' @name gene.exp.chr21
#' @format A matrix: 50 samples
NULL
