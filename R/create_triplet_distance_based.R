#' @title Map DNAm to target genes using distance approaches, and
#' TF to the DNAm region using JASPAR2024 TFBS.
#' @description This function wraps two other functions
#' \code{get_region_target_gene} and \code{get_tf_in_region} from the package.
#' This function will map a region to a target gene using three methods
#' (mapping to the closest gene,
#' mapping to any gene within a given window of distance, or mapping to a
#' fixed number of nearby genes upstream or downstream).
#' To find TFs binding to the region, JASPAR2024 is used.
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
#' @param target.rm.promoter.regions.from.distal.linking When performing distal linking
#' with method = "windows" or method = "nearby.genes", or "closest.gene.tss",
#' if set to TRUE (default), probes in promoter regions will be removed
#' from the input.
#' @param target.promoter.upstream.dist.tss Number of base pairs (bp) upstream of
#' TSS to consider as promoter regions. Defaults to 2000 bp.
#' @param target.promoter.downstream.dist.tss Number of base pairs (bp) downstream of
#' TSS to consider as promoter regions. Defaults to 2000 bp.
#' @param motif.search.window.size Integer value to extend the regions.
#' For example, a value of 50 will
#' extend 25 bp upstream and 25 downstream the region. Default is no increase
#' @param motif.search.p.cutoff motifmatchr pvalue cut-off. Default 1e-8.
#' @param TF.peaks.gr  A granges with TF peaks to be overlaped with input region
#' Metadata column expected "id" with TF name. Default NULL. Note that Remap catalog
#' can be used as shown in the examples.
#' @param max.distance.region.target Max distance between region and target gene. Default 1Mbp.
#' @param cores Number of CPU cores to be used. Default 1.
#' @return A data frame with TF, target and RegionID information.
#' @examples
#' regions.names <- c("chr3:189631389-189632889","chr4:43162098-43163498")
#' triplet <- create_triplet_distance_based(
#'    region = regions.names,
#'    motif.search.window.size = 500,
#'    target.method = "closest.gene"
#' )
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
    target.promoter.upstream.dist.tss = 2000,
    target.promoter.downstream.dist.tss = 2000,
    target.rm.promoter.regions.from.distal.linking = TRUE,
    motif.search.window.size = 0,
    motif.search.p.cutoff = 1e-8,
    TF.peaks.gr = NULL,
    max.distance.region.target =  10^6,
    cores = 1
){

    target.method <- match.arg(target.method)
    genome <- match.arg(genome)
    probeID <- NULL
    if (is(region,"character") | is(region,"factor")) {
        region.gr <- make_granges_from_names(region)
    } else if (is(region,"GenomicRanges")) {
        region.gr <- region

        if ("probeID" %in% colnames(values(region))) {
            probeID <- values(region)$probeID
        }

    } else if (is(region,"SummarizedExperiment")) {
        region.gr <- rowRanges(region)
        if ("probeID" %in% colnames(values(region))) {
            probeID <- values(region)$probeID
        }

    } else {
        stop(
            "region input is a ", class(region),
            ". Expecting a charcter, GRanges or SE"
        )
    }
    region.names <- make_names_from_granges(region.gr)

    message("Finding target genes")
    region.target <- get_region_target_gene(
        regions.gr = region.gr,
        genome = genome,
        window.size = target.window.size,
        method = target.method,
        num.flanking.genes = target.num.flanking.genes,
        promoter.upstream.dist.tss = target.promoter.upstream.dist.tss,
        promoter.downstream.dist.tss = target.promoter.downstream.dist.tss,
        rm.promoter.regions.from.distal.linking = target.rm.promoter.regions.from.distal.linking
    )

    message("Looking for TFBS")
    region.tf <- get_tf_in_region(
        region = region.gr,
        genome = genome,
        window.size = motif.search.window.size,
        p.cutoff = motif.search.p.cutoff,
        cores = cores,
        TF.peaks.gr = TF.peaks.gr
    )
    triplet <- dplyr::inner_join(region.target, region.tf)
    triplet <- triplet %>% dplyr::filter(!is.na(.data$TF))
    
    message("Removing regions and target genes from different chromosomes")
    triplet <- triplet %>% dplyr::filter(!is.na(.data$distance_region_target_tss))

    message("Removing regions and target genes with ditance higher than ", max.distance.region.target, " bp")
    triplet <- triplet %>% dplyr::filter(abs(.data$distance_region_target_tss) < max.distance.region.target)


    triplet <- triplet %>%
        dplyr::relocate(.data$distance_region_target_tss, .after = dplyr::last_col()) %>%
        dplyr::relocate(contains("pos"), .after = dplyr::last_col())

    if (!is.null(probeID)) {
        triplet$probeID <- probeID[match(triplet$regionID, region.names)]
        triplet <- triplet %>%
            dplyr::relocate(.data$probeID, .after = "regionID")
    }

    return(triplet)
}
