#' For a given dnam region maps candidate target gene and scan for TF binding
#' @description This function wraps two other functions:
#' get_region_target_gene and get_tf_in_region from the package.
#' This function will map a region to a target gene using two methods
#' one mapping to the closest gene, and another one mapping to any gene within
#' a given window.
#' To find TF binding to the region JASPAR2020 is used to look into the region.
#'
#' @importFrom tidyr separate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @param dnam DNA methylation matrix  (columns: samples in the same order as \code{exp} matrix, rows: regions/probes)
#' @param motif.search.window.size Integer value to extend the regions. For example, a value of 50 will
#' extend 25 bp upstream and 25 downstream the region. Default is no increase
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param motif.search.p.cutoff Motifmatcher pvalue cut-off. Default 1e-8.
#' @param cores Number of CPU cores to be used. Default 1.
#' @param method How genes are mapped to regions: closest gene promoter to the region ("closest.gene"); or
#' genes within a window around the region ("window").
#' @param window.size When \code{method = "window"}, number of base pairs to extend the region (+- window.size/2).
#' Default is 500kbp (or +/- 250kbp, i.e. 250k bp from start or end of the region)
#' @examples
#' data("dna.met.chr21")
#' triplet <- get_triplet(dna.met.chr21[1:100,], motif.search.window.size = 50)
#' @export
get_triplet <- function(
    dnam,
    genome = c("hg38","hg19"),
    target.method = c("closest.gene","window"),
    target.window.size = 500 * 10^3,
    motif.search.window.size = 0,
    motif.search.p.cutoff = 1e-8,
    cores = 1
){

    target.method <- match.arg(target.method)
    genome <- match.arg(genome)

    if(any(grepl("cg",rownames(dnam)))){
        dnam.regions <- map_probes_to_regions(dnam)
    }

    dnam.regions.gr <- make_granges_from_names(rownames(dnam.regions))

    message("Finding target genes")
    region.target <- get_region_target_gene(
        regions.gr = dnam.regions.gr,
        genome = genome,
        method = target.method
    )

    message("Looking for TFBS")
    region.tf <- get_tf_in_region(
        region = dnam.regions.gr,
        genome = genome,
        window.size = motif.search.window.size,
        p.cutoff = motif.search.p.cutoff,
        cores = cores
    )
    triplet <- dplyr::inner_join(region.target,region.tf)
    return(triplet)
}
