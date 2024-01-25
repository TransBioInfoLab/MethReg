#' TCGA-COAD clinical matrix for 38 samples
#' retrieved from GDC using \link[TCGAbiolinks]{TCGAbiolinks}
#' @docType data
#' @name clinical
#' @format A matrix: 38 samples (rows) and variables (columns)
#' patient, sample, gender and sample_type
"clinical"

#' TCGA-COAD DNA methylation matrix (beta-values) for 38 samples (only chr21)
#' retrieved from GDC using \link[TCGAbiolinks]{TCGAbiolinks}
#' @docType data
#' @name dna.met.chr21
#' @format A beta-value matrix with 38 samples, includes CpG IDs in the rows and
#' TCGA sample identifiers in the columns
"dna.met.chr21"

#' TCGA-COAD gene expression matrix (log2 (FPKM-UQ + 1))
#' for 38 samples (only chromosome 21) retrieved from GDC using TCGAbiolinks
#' @docType data
#' @name gene.exp.chr21.log2
#' @format A log2 (FPKM-UQ + 1) gene expression matrix with 38 samples,
#' includes Ensembl IDs in the rows and TCGA sample identifiers in the columns
"gene.exp.chr21.log2"
