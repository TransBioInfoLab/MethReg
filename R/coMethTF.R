#' @title coMethTF: functional annotation of DMRs identified in epigenome-wide association studies
#' @description To provide functional annotations for differentially methylated regions (DMRs)
#' and differentially methylated CpG sites (DMS),`coMethTF` performs
#' integrative analyses using matched DNA methylation and gene expression
#' along with Transcription Factor Binding Sites (TFBS) data.
#' coMethTF evaluates, prioritizes and annotates DNA methylation regions
#' (or sites) with high regulatory potential that works synergistically with
#' TFs to regulate target gene expressions, without any additional ChIP-seq data.
#' @docType package
#' @name coMethTF
NULL

#' TCGA-COAD clinical matrix for 38 samples
#' Retrieved from GDC using TCGAbiolinks
#' @docType data
#' @name clinical
#' @format A matrix: 38 samples
NULL

#' TCGA-COAD DNA methylation matrix for 38 samples (only chrmossome 21)
#' Retrieved from GDC using TCGAbiolinks
#' @docType data
#' @name dna.met.chr21
#' @format A matrix: 38 samples
NULL

#' TCGA-COAD gene expression matrix for 38 samples (only chrmossome 21)
#' Retrieved from GDC using TCGAbiolinks
#' @docType data
#' @name gene.exp.chr21
#' @format A matrix: 38 samples
NULL
