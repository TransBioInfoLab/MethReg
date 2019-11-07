#' @title Mapping regions to closes gene
#' @description Maps a given region to th closest gene
#' @param regions A dataframe with chrom, start, end or a GRanges
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param method For the moment mapping region to closest gene
#' @importFrom coMethDMR AnnotateResults
#' @importFrom GenomicRanges distanceToNearest nearest ranges
#' @importFrom tidyr unite
#' @importFrom ELMER getTSS
#' @examples
#' regions.gr <- data.frame(
#'  chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'   start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'   end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'                          stringsAsFactors = FALSE)  %>%
#'                          makeGRangesFromDataFrame
#'  getDNAm.target(regions.gr = regions.gr ,genome = "hg19")
getDNAm.target <- function(
    regions.gr,
    genome = c("hg38","hg19"),
    method = c("closest.gene")
){

    method <- match.arg(method)
    genome <- match.arg(genome)

    if(!is(regions.gr,"GRanges")) stop("regions.gr must be a GRanges")

    tssAnnot <- ELMER::getTSS(genome = genome)
    neargenes <- tssAnnot[nearest(regions.gr,tssAnnot)] %>% as.data.frame()
    distance.region.tss <- values(distanceToNearest(regions.gr,tssAnnot))$distance

    neargenes <- cbind(neargenes[,c("seqnames","start","end","external_gene_name","ensembl_gene_id")],
                       distance.region.tss)

    colnames(neargenes)[1:3] <- c("gene_chrom","gene_start","gene_end")
    regionID <- ranges(regions.gr) %>%
        as.data.frame() %>%
        tidyr::unite(col = "regionID")

    region.df <- regions.gr %>% as.data.frame()
    colnames(region.df)[1:3] <- c("region_chrom","region_start","region_end")

    out <- cbind(region.df, neargenes)
    return(out)
}
