#' @title Mapping regions to gene
#' @description To map a region to genes there are two options: 1) closest gene
#' 2) map to all genes within a window around the region (default window.width = 500kbp
#' (i.e. +/- 250kbp from start or end of the region)).
#' @param regions.gr A Genomic Ranges objec GRanges
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param method How regions are mapped to genes: closest gene ("closest.gene"); or
#' genes within a window around the region ("window").
#' @param window.width When \code{method = "window"}, number of base pairs to extend the region (+- window.width/2).
#' Default is 500kbp (or +/- 250kbp, i.e. 250k bp from start or end of the region)
#' @importFrom GenomicRanges distanceToNearest nearest ranges makeGRangesFromDataFrame values seqnames distance
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom tidyr unite
#' @importFrom ELMER getTSS
#' @examples
#' library(GenomicRanges)
#' library(dplyr)
#'
#' # Create example region
#' regions.gr <- data.frame(
#'               chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'               start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'               end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'               stringsAsFactors = FALSE)  %>%
#'      makeGRangesFromDataFrame
#'
#'  # map to closest gene
#'  region.closest <- get_region_target_gene(
#'                       regions.gr = regions.gr,
#'                       genome = "hg19",
#'                       method = "closest.gene")
#'
#'  # map to all gene within region +- 250kbp
#'  region.window <- get_region_target_gene(
#'                       regions.gr = regions.gr,
#'                       genome = "hg19",
#'                       method = "window")
#' @export
get_region_target_gene <- function(
    regions.gr,
    genome = c("hg38","hg19"),
    method = c("closest.gene","window"),
    window.width = 500000
){

    method <- match.arg(method)
    genome <- match.arg(genome)

    if(!is(regions.gr,"GRanges")) stop("regions.gr must be a GRanges")

    if(method == "closest.gene"){
        tssAnnot <- ELMER::getTSS(genome = genome)

        neargenes <- tssAnnot[nearest(regions.gr,tssAnnot)] %>% as.data.frame()
        distance.region.tss <- values(distanceToNearest(regions.gr,tssAnnot))$distance

        neargenes <- cbind(neargenes[,c("seqnames","start","end","external_gene_name","ensembl_gene_id")],
                           distance.region.tss)

        colnames(neargenes)[1:3] <- c("gene_chrom","gene_start","gene_end")
        colnames(neargenes)[4] <- "target_gene_name"
        colnames(neargenes)[5] <- "target"

        regionID <- paste0(
            regions.gr %>% seqnames %>% as.character(),
            ":",
            regions.gr %>% start,
            "-",
            regions.gr %>% end)
        out <- dplyr::bind_cols(
            data.frame("regionID" = regionID, stringsAsFactors = FALSE),
            neargenes
        ) %>% tibble::as_tibble()

    } else {
        geneAnnot <- get_gene_information(genome = genome,as.granges = TRUE)
        geneAnnot$entrezgene <- NULL
        geneAnnot <- unique(geneAnnot)
        regions.gr.extend <- regions.gr + (window.width/2)

        overlap <- findOverlaps(regions.gr.extend,geneAnnot)

        regionID <- paste0(
            regions.gr[queryHits(overlap)] %>% seqnames %>% as.character(),
            ":",
            regions.gr[queryHits(overlap)] %>% start,
            "-",
            regions.gr[queryHits(overlap)] %>% end)


        regionID.extended <- paste0(
            regions.gr.extend[queryHits(overlap)] %>% seqnames %>% as.character(),
            ":",
            regions.gr.extend[queryHits(overlap)] %>% start,
            "-",
            regions.gr.extend[queryHits(overlap)] %>% end
        )

        genes.overlapping <- geneAnnot[subjectHits(overlap)] %>% as.data.frame()
        colnames(genes.overlapping) <- paste0("gene_",colnames(genes.overlapping))

        colnames(genes.overlapping)[grep("ensembl_gene_id",colnames(genes.overlapping))] <- "target"

        out <- dplyr::bind_cols(
            data.frame("regionID" = regionID,stringsAsFactors = FALSE),
            data.frame("regionID.extended" = regionID.extended %>% as.character(),
                       "window.extended.width" = window.width,
                       "Distance region-gene" = distance(regions.gr[queryHits(overlap)],
                                                         geneAnnot[subjectHits(overlap)]),
                       stringsAsFactors = FALSE),
            genes.overlapping)  %>% tibble::as_tibble()
    }
    return(out)
}


