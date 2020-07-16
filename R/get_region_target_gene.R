#' @title Obtain target genes that are close to input regions
#' @description To map genes to a region there are two options: 1) closest gene
#' 2) map to all genes within a window around the region (default window.size = 500kbp
#' (i.e. +/- 250kbp from start or end of the region)).
#' @param regions.gr A Genomic Ranges object (GRanges)
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param method How genes are mapped to regions: closest gene promoter to the region ("closest.gene"); or
#' genes within a window around the region ("window"); or a fixed number genes upstream
#' and downstream of the region "nearest.genes"
#' @param window.size When \code{method = "window"}, number of base pairs to extend the region (+- window.size/2).
#' Default is 500kbp (or +/- 250kbp, i.e. 250k bp from start or end of the region)
#' @param num.flanking.genes Number of flanking genes upstream and downstream to search.
#' For example, if num.flanking.genes = 5, it will return the 5 genes upstream
#' and 5 genes dowstream of the given region.
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom tidyr unite
#' @importFrom dplyr select
#' @importFrom methods is
#' @examples
#' library(GenomicRanges)
#' library(dplyr)
#'
#' # Create example region
#' regions.gr <- data.frame(
#'        chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'        start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'        end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'        stringsAsFactors = FALSE)  %>%
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
#'                       method = "window",
#'                       window.size = 500 * 10^3)
#'
#'  # map regions to n upstream and n dowstream genes
#'  region.nearest <- get_region_target_gene(
#'                       regions.gr = regions.gr,
#'                       genome = "hg19",
#'                       method = "nearest.genes",
#'                       num.flanking.genes = 5)
#' @export
get_region_target_gene <- function(
    regions.gr,
    genome = c("hg38","hg19"),
    method = c("closest.gene","window","nearest.genes"),
    window.size = 500 * 10^3,
    num.flanking.genes = 5
){

    method <- match.arg(method)
    genome <- match.arg(genome)

    if(!is(regions.gr,"GRanges")) stop("regions.gr must be a GRanges")

    if(method == "closest.gene"){
        out <- get_region_target_gene_closest(regions.gr, genome)
    } else if(method == "window.size"){
        out <- get_region_target_gene_window(regions.gr, genome, window.size)
    } else {
        out <- get_region_target_gene_nearest.genes(
            regions.gr = regions.gr,
            genome = genome,
            num.flanking.genes = num.flanking.genes
        )
    }
    return(out)
}

get_region_target_gene_closest <- function(
    regions.gr,
    genome
){

    # Get gene information
    gene.info <- get_gene_information(genome = genome, as.granges = TRUE)

    # get gene promoter
    gene.promoters <-  gene.info %>%
        promoters(upstream = 2000, downstream = 2000)

    hits <- findOverlaps(regions.gr, gene.promoters, ignore.strand = FALSE, select = "all")

    # overlap region and promoter
    neargenes <- gene.info[subjectHits(hits)] %>% as.data.frame(row.names = NULL)

    regions.gr <- regions.gr[queryHits(hits)]
    neargenes <- cbind(
        neargenes[,c("seqnames",
                     "start",
                     "end",
                     "external_gene_name",
                     "ensembl_gene_id")])

    colnames(neargenes)[1:3] <- c("target_gene_chrom","target_gene_start","target_gene_end")
    colnames(neargenes)[4] <- "target_gene_name"
    colnames(neargenes)[5] <- "target"

    regionID <- regions.gr %>% as.data.frame %>% dplyr::select(1:3)
    regionID <- paste0(regionID[[1]],":",regionID[[2]],"-",regionID[[3]])
    out <- dplyr::bind_cols(
        data.frame("regionID" = regionID, stringsAsFactors = FALSE),
        neargenes[,4:5]
    ) %>% tibble::as_tibble()
    out
}

get_region_target_gene_window <- function(
    regions.gr,
    genome,
    window.size
){
    geneAnnot <- get_gene_information(genome = genome, as.granges = TRUE)
    geneAnnot$entrezgene <- NULL
    geneAnnot <- unique(geneAnnot)
    regions.gr.extend <- regions.gr + (window.size/2)

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

    genes.overlapping <- geneAnnot[subjectHits(overlap)] %>% as.data.frame(row.names = NULL)
    colnames(genes.overlapping) <- paste0("gene_",colnames(genes.overlapping))

    colnames(genes.overlapping)[grep("ensembl_gene_id",colnames(genes.overlapping))] <- "target"
    colnames(genes.overlapping)[grep("external_gen",colnames(genes.overlapping))] <- "target_gene_name"
    genes.overlapping <- genes.overlapping[,grep("target",colnames(genes.overlapping))]

    out <- dplyr::bind_cols(
        data.frame("regionID" = regionID,stringsAsFactors = FALSE),
        #data.frame("regionID.extended" = regionID.extended %>% as.character(),
        #           "window.extended.width" = window.width,
        #           "Distance region-gene" = distance(regions.gr[queryHits(overlap)],
        #                                             geneAnnot[subjectHits(overlap)]),
        #           stringsAsFactors = FALSE),
        genes.overlapping)  %>% tibble::as_tibble()
    out
}

#' @examples
#' regions.gr <- data.frame(
#'        chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'        start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'        end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'        stringsAsFactors = FALSE)  %>%
#'      makeGRangesFromDataFrame
#' genome <- "hg38"
#' num.flanking.genes <- 10
#' @noRd
get_region_target_gene_nearest.genes <- function(
    regions.gr,
    genome,
    num.flanking.genes
){

    names(regions.gr) <- make_names_from_granges(regions.gr)
    print( names(regions.gr))
    genes.gr <- get_gene_information(genome = genome, as.granges = TRUE)
    genes.gr$entrezgene <- NULL
    genes.gr <- unique(genes.gr)

    # Optimized version
    # Idea: vectorize search
    # 1) For all regions, get nearest gene
    # 2) check follow and overlapping genes recursively
    # 3) check precede and overlapping genes recursively
    # 4) map the positions based on min distance (L1)
    # The input data has to be at gene level and not transcript which would broke
    # some of the optimizations for which we remove the genes already evaluated

    # 1) For all regions, get nearest gene
    message("Identifying genes close to the region")
    nearest.idx <- nearest(
        regions.gr,
        genes.gr,
        select = "all",
        ignore.strand = TRUE
    )

    closest.genes <- tibble::tibble(
        genes.gr[nearest.idx %>% subjectHits] %>% as.data.frame(), # nearest gene to the region
        "ID" = names(regions.gr)[nearest.idx %>% queryHits], # region being evaluated
        "Distance" = get_region_gene_distance(nearest.idx, regions.gr, genes.gr)
    )

    message("Identifying ",num.flanking.genes, " genes downstream to the region")
    precede.genes <- get_region_target_gene_nearest.genes_aux(
        direction.fun = GenomicRanges::precede,
        nearest.idx = nearest.idx,
        genes.gr = genes.gr,
        regions.gr = regions.gr,
        num.flanking.genes = num.flanking.genes
    )

    # 2) check follow and overlapping genes recursively
    message("Identifying ",num.flanking.genes, " genes upstream of the region")
    follow.genes <- get_region_target_gene_nearest.genes_aux(
        direction.fun = GenomicRanges::follow,
        nearest.idx = nearest.idx,
        genes.gr = genes.gr,
        regions.gr = regions.gr,
        num.flanking.genes = num.flanking.genes
    )
    ret <- rbind(closest.genes,precede.genes,follow.genes) %>% unique
    ret <- ret[order(ret$Distance),]

    ret <- ret[,
               c("ID",
                 "ensembl_gene_id",
                 grep("external_gene_", colnames(ret), value = TRUE),
                 "Distance")
               ]

    message("Identifying gene position for each region")
    print(ret[order(ret$ID),])
    ret <- get_region_target_gene_nearest.genes_addPos(ret, num.flanking.genes)
    colnames(ret)[1:3] <- c("regionID", "target", "target_gene_name")
    return(ret)
}

get_region_gene_distance <- function(idx, regions.gr,geneAnnot){
    signal <- ifelse(
        start(regions.gr[idx %>% queryHits]) < start(geneAnnot[idx  %>% subjectHits]), 1, -1
    )
    distance(
        regions.gr[idx %>% queryHits],
        geneAnnot[idx %>% subjectHits],
        select = "all",
        ignore.strand = TRUE) * signal

}

#' @importFrom progress progress_bar
#' @importFrom GenomicRanges findOverlaps distance nearest
get_region_target_gene_nearest.genes_aux <- function(
    direction.fun,
    nearest.idx,
    genes.gr,
    regions.gr = regions.gr,
    num.flanking.genes
){

    pb <- progress::progress_bar$new(total = num.flanking.genes)
    evaluating <- nearest.idx %>% queryHits
    # we will start our search from the nearest gene from the region
    idx <- nearest.idx %>% as.data.frame()

    ret <- NULL
    # 2) check precede and overlapping genes recursively
    for (i in 1:(num.flanking.genes)) {
        idx <- unique(
            rbind(
                tibble::as_tibble(
                    findOverlaps(
                        genes.gr[idx$subjectHits],
                        genes.gr,
                        ignore.strand = TRUE,
                        type = "any",
                        select = "all"
                    )
                ),
                tibble::as_tibble(
                    direction.fun(
                        genes.gr[idx$subjectHits],
                        genes.gr,
                        select = "all",
                        ignore.strand = TRUE
                    )
                )
            )
        )
        idx$evaluating <-  evaluating[idx$queryHits]

        # remove same target gene and region if counted twice
        idx <- idx[!duplicated(idx[, c("subjectHits","evaluating")]), ]

        # todo remove already evaluated previously (we don't wanna do it again)
        idx <- idx[
            !paste0(genes.gr[idx$subjectHits]$ensembl_gene_id, names(regions.gr)[idx$evaluating]) %in%
                paste0(ret$ensembl_gene_id, ret$ID),
            ]
        evaluating <- evaluating[idx$queryHits]

        ret <- rbind(ret, # keep old results
                     tibble::tibble(
                         genes.gr[idx$subjectHits] %>% as.data.frame(),
                         "ID" = names(regions.gr)[evaluating],
                         "Distance" = ifelse(start(regions.gr[evaluating]) < start(genes.gr[idx$subjectHits]), 1,-1) *
                             distance(regions.gr[evaluating],
                                      genes.gr[idx$subjectHits],
                                      select = "all",
                                      ignore.strand = TRUE)
                     )
        )
        pb$tick()
    }
    ret <- ret[!duplicated(ret[,c("ensembl_gene_id","ID")]),]
    pb$terminate()
    return(ret)
}

#' @importFrom dplyr group_by
get_region_target_gene_nearest.genes_addPos <- function(ret, num.flanking.genes){
    f <- function(pairs) {
        center <- which(abs(pairs$Distance) == min(abs(pairs$Distance)))[1]
        pos <- setdiff(-center:(nrow(pairs) - center), 0)
        pairs$Side <- ifelse(pos > 0, paste0("R", abs(pos)), paste0("L", abs(pos)))

        out <- pairs %>%
            dplyr::filter(
                pairs$Side %in% c(paste0("R", 1:(num.flanking.genes)),
                                  paste0("L", 1:(num.flanking.genes))
                )
            )

        if (nrow(out) < num.flanking.genes) {
            if (paste0("R", floor(num.flanking.genes)) %in% out$Side) {
                cts <- length(grep("L", sort(pairs$Side), value = TRUE))
                out <- pairs %>%
                    dplyr::filter(
                        Side %in% c(paste0("R", 1:(num.flanking.genes - cts)),
                                    grep("L", sort(out$Side), value = TRUE))
                    )
            } else {
                cts <- length(grep("R", sort(pairs$Side), value = TRUE))
                out <- pairs %>%
                    dplyr::filter(
                        pairs$Side %in%
                            c(paste0("L", 1:(num.flanking.genes - cts)),
                              grep("R", sort(out$Side), value = TRUE))
                    )
            }
        }
        out <- out[order(out$Distance), ]
        return(out)
    }

    ret %>% group_by(ID) %>% do(f(.))
}
