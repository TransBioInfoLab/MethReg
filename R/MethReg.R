#' @title MethReg: functional annotation of DMRs identified in epigenome-wide
#' association studies
#' @description To provide functional annotations for
#' differentially methylated regions (DMRs)
#' and differentially methylated CpG sites (DMS), \code{MethReg} performs
#' integrative analyses using matched DNA methylation and gene expression
#' along with Transcription Factor Binding Sites (TFBS) data.
#' MethReg evaluates, prioritizes and annotates DNA methylation regions
#' (or sites) with high regulatory potential that works synergistically with
#' TFs to regulate target gene expressions, without any
#' additional ChIP-seq data.
#' @docType package
#' @name MethReg
NULL

#' TCGA-COAD clinical matrix for 38 samples
#' retrieved from GDC using TCGAbiolinks
#' @docType data
#' @name clinical
#' @format A matrix: 38 samples (rows) and variables (columns)
#' patient, sample, gender and sample_type
NULL

#' TCGA-COAD DNA methylation matrix (beta-values) for 38 samples (only chr21)
#' retrieved from GDC using TCGAbiolinks.
#' @docType data
#' @name dna.met.chr21
#' @format A beta-value matrix with 38 samples, includes CpG IDs in the rows and
#' TCGA sample identifiers in the columns
NULL

#' TCGA-COAD gene expression matrix (log2 (FPKM-UQ + 1))
#' for 38 samples (only chromosome 21) retrieved from GDC using TCGAbiolinks
#' @docType data
#' @name gene.exp.chr21.log2
#' @format A log2 (FPKM-UQ + 1) gene expression matrix with 38 samples,
#' includes Ensembl IDs in the rows and TCGA sample identifiers in the columns
NULL
