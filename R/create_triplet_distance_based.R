#' @title Map DNAm to target genes using distance approaches, and
#' TF to the DNAm region using JASPAR2020 TFBS.
#' @description This function wraps two other functions
#' \code{get_region_target_gene} and \code{get_tf_in_region} from the package.
#' This function will map a region to a target gene using three methods
#' (mapping to the closest gene,
#' mapping to any gene within a given window of distance, or mapping to a
#' fixed number of nearby genes upstream or downstream).
#' To find TFs binding to the region, JASPAR2020 is used.
#'
#' @importFrom tidyr separate
#' @import GenomicRanges
#' @param region A Granges or a named vector with
#' regions (i.e "chr21:100002-1004000")
#' @param genome Human genome reference "hg38" or "hg19"
#' @param target.method How genes are mapped to regions: regions
#' overlapping gene promoter ("genes.promoter.overlap");
#' genes within a window around the region ("window"); or fixed number of
#' nearby genes upstream and downstream from the region
#' @param target.window.size When \code{method = "window"}, number of base
#' pairs to extend the region (+- window.size/2).
#' Default is 500kbp (or +/- 250kbp, i.e. 250k bp from start or end of the region)
#' @param target.num.flanking.genes Number of flanking genes upstream and
#' downstream to search.
#' For example, if \code{target.num.flanking.genes = 5}, it will return the
#' 5 genes upstream and 5 genes downstream
#' @param motif.search.window.size Integer value to extend the regions.
#' For example, a value of 50 will
#' extend 25 bp upstream and 25 downstream the region. Default is no increase
#' @param motif.search.p.cutoff motifmatchr pvalue cut-off. Default 1e-8.
#' @param cores Number of CPU cores to be used. Default 1.
#' @return A data frame with TF, target and RegionID information.
#' @examples
#' regions.names <- c("chr3:189631389-189632889","chr4:43162098-43163498")
#' triplet <- create_triplet_distance_based(
#'    region = regions.names,
#'    motif.search.window.size = 500,
#'    target.method = "closest.gene"
#' )
#' \dontrun{
#' data("dna.met.chr21")
#' dnam.regions <- make_se_from_dnam_probes(dna.met.chr21)
#' triplet <- create_triplet_distance_based(
#'    region = rownames(dnam.regions)[1:100],
#'    motif.search.window.size = 50
#' )
#' }
#' @export
create_triplet_distance_based <- function(
    region,
    genome = c("hg38","hg19"),
    target.method = c(
        "genes.promoter.overlap",
        "window",
        "nearby.genes",
        "closest.gene"
    ),
    target.window.size = 500 * 10^3,
    target.num.flanking.genes = 5,
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
    } else if(is(region,"SummarizedExperiment")){
        region.gr <- rowRanges(region)
        region.names <- make_names_from_granges(region.gr)
    }

    message("Finding target genes")
    region.target <- get_region_target_gene(
        regions.gr = region.gr,
        genome = genome,
        window.size = target.window.size,
        method = target.method,
        num.flanking.genes = target.num.flanking.genes
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
