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
#' @param region A Granges or a named vector with regions (i.e "chr21:100002-1004000")
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param target.method How genes are mapped to regions: closest gene promoter to the region ("closest.gene"); or
#' genes within a window around the region ("window").
#' @param target.window.size When \code{method = "window"}, number of base pairs to extend the region (+- window.size/2).
#' Default is 500kbp (or +/- 250kbp, i.e. 250k bp from start or end of the region)
#' @param motif.search.window.size Integer value to extend the regions. For example, a value of 50 will
#' extend 25 bp upstream and 25 downstream the region. Default is no increase
#' @param motif.search.p.cutoff Motifmatcher pvalue cut-off. Default 1e-8.
#' @param cores Number of CPU cores to be used. Default 1.
#' @examples
#' \dontrun{
#' data("dna.met.chr21")
#' dnam.regions <- map_probes_to_regions(dna.met.chr21)
#' triplet <- get_triplet(
#'    region = rownames(dnam.regions)[1:100],
#'    motif.search.window.size = 50
#' )
#' }
#' @export
get_triplet <- function(
    region,
    genome = c("hg38","hg19"),
    target.method = c("closest.gene","window"),
    target.window.size = 500 * 10^3,
    motif.search.window.size = 0,
    motif.search.p.cutoff = 1e-8,
    cores = 1
){

    target.method <- match.arg(target.method)
    genome <- match.arg(genome)

    if(is(region,"character") | is(region,"factor")){
        region.gr <- make_granges_from_names(region)
        region.names <- region
    } else if(is(region,"GenomicRanges")){
        region.gr <- region
        region.names <- make_names_from_granges(region)
    }

    message("Finding target genes")
    region.target <- get_region_target_gene(
        regions.gr = region.gr,
        genome = genome,
        method = target.method
    )

    message("Looking for TFBS")
    region.tf <- get_tf_in_region(
        region = region.gr,
        genome = genome,
        window.size = motif.search.window.size,
        p.cutoff = motif.search.p.cutoff,
        cores = cores
    )
    triplet <- dplyr::inner_join(region.target, region.tf)
    return(triplet)
}
