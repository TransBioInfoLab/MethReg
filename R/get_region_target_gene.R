#' @title Obtain target genes of input regions based on distance
#' @description To map an input region to genes there are three options:
#' 1) map region to closest gene
#' 2) map region to all genes within a window around the region
#' (default window.size = 500kbp
#' (i.e. +/- 250kbp from start or end of the region)).
#' 3) map region to a fixed number of nearby genes (upstream/downstream)
#' @param regions.gr A Genomic Ranges object (GRanges) or a
#' SummarizedExperiment object (rowRanges will be used)
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param method How genes are mapped to regions: region overlapping gene
#' promoter ("genes.promoter.overlap"); or
#' genes within a window around the region ("window"); or a fixed number genes upstream
#' and downstream of the region ("nearby.genes");
#' or closest gene to the region ("closest.gene")
#' @param window.size When \code{method = "window"}, number of base pairs
#' to extend the region (+- window.size/2).
#' Default is 500kbp (or +/- 250kbp, i.e. 250k bp from start or end of the region)
#' @param num.flanking.genes When \code{method = "nearby.genes"}, set the number
#' of flanking genes upstream and downstream to search.Defaults to 5.
#' For example, if num.flanking.genes = 5, it will return the 5 genes upstream
#' and 5 genes downstream of the given region.
#' @param rm.promoter.regions.from.distal.linking When performing distal linking
#' with method = "windows" or method = "nearby.genes", if set to TRUE (default),
#' probes in promoter regions will be removed from the input.
#' @param promoter.upstream.dist.tss Number of base pairs (bp) upstream of
#' TSS to consider as promoter regions. Defaults to 2000 bp.
#' @param promoter.downstream.dist.tss Number of base pairs (bp) downstream of
#' TSS to consider as promoter regions. Defaults to 2000 bp.
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
#'  region.genes.promoter.overlaps <- get_region_target_gene(
#'                       regions.gr = regions.gr,
#'                       genome = "hg19",
#'                       method = "genes.promoter.overlap")
#'
#'  # map to all gene within region +- 250kbp
#'  region.window.genes <- get_region_target_gene(
#'                       regions.gr = regions.gr,
#'                       genome = "hg19",
#'                       method = "window",
#'                       window.size = 500 * 10^3
#'  )
#'
#'  # map regions to n upstream and n downstream genes
#'  region.nearby.genes <- get_region_target_gene(
#'                       regions.gr = regions.gr,
#'                       genome = "hg19",
#'                       method = "nearby.genes",
#'                       num.flanking.genes = 5
#'  )
#' @export
#' @return A data frame with the following information:
#' regionID, Target symbol, Target ensembl ID
#' @details For the analysis of probes in promoter regions (promoter analysis),
#'  we recommend setting
#'  \code{method = "genes.promoter.overlap"}.
#'
#'  For the analysis of probes in distal regions (distal analysis),
#'  we recommend setting either \code{method = "window"} or \code{method = "nearby.genes"}.
#'
#'  Note that because \code{method = "window"} or \code{method = "nearby.genes"} are
#'  mainly used for analyzing distal probes,
#'  by default \code{rm.promoter.regions.from.distal.linking = TRUE} to
#'  remove probes in promoter regions.
get_region_target_gene <- function(
    regions.gr,
    genome = c("hg38","hg19"),
    method = c(
        "genes.promoter.overlap",
        "window",
        "nearby.genes",
        "closest.gene"
    ),
    promoter.upstream.dist.tss = 2000,
    promoter.downstream.dist.tss = 2000,
    window.size = 500 * 10^3,
    num.flanking.genes = 5,
    rm.promoter.regions.from.distal.linking = TRUE
){

    method <- match.arg(method)
    genome <- match.arg(genome)

    if(is(regions.gr,"SummarizedExperiment")) {
        regions.gr <- SummarizedExperiment::rowRanges(regions.gr)
    }

    if(!is(regions.gr,"GRanges")) stop("regions.gr must be a GRanges")
    if(method != "genes.promoter.overlap" & rm.promoter.regions.from.distal.linking){
        message("Removing regions overlapping promoter regions")
        regions.gr <- subset_by_non_promoter_regions(
            regions.gr = regions.gr,
            genome =  genome,
            upstream = promoter.upstream.dist.tss,
            downstream = promoter.downstream.dist.tss
        )
    }

    if(method == "genes.promoter.overlap"){
        message("Mapping regions to the closest gene")
        out <- get_region_target_gene_by_promoter_overlap(
            regions.gr = regions.gr,
            genome = genome,
            upstream = promoter.upstream.dist.tss,
            downstream = promoter.downstream.dist.tss
        )
    } else if(method == "window"){
        message("Mapping regions to genes within a window of size: ", window.size, " bp")
        out <- get_region_target_gene_window(regions.gr, genome, window.size)
    } else if(method == "nearby.genes"){
        out <- get_region_target_gene_nearby.genes(
            regions.gr = regions.gr,
            genome = genome,
            num.flanking.genes = num.flanking.genes
        )
    } else if(method == "closest.gene"){
        out <- get_region_target_gene_closest(
            regions.gr = regions.gr,
            genome = genome
        )
    }
    return(out)
}

#' @noRd
#' @examples
#' regions.gr <- data.frame(
#'        chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'        start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'        end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'        stringsAsFactors = FALSE)  %>%
#'      makeGRangesFromDataFrame
#' @importFrom GenomicRanges distanceToNearest
#' @importFrom S4Vectors values
get_region_target_gene_closest <- function(
    regions.gr,
    genome
){
    gene.info <- get_gene_information(genome = genome, as.granges = TRUE)
    hits <- nearest(regions.gr,gene.info)
    nearest.genes <- gene.info[hits] %>% as.data.frame(row.names = NULL)
    distance <- distanceToNearest(regions.gr,gene.info)


    nearest.genes <- dplyr::bind_cols(
        nearest.genes[,c(
            "seqnames",
            "start",
            "end",
            "external_gene_name",
            "ensembl_gene_id"
        )],
        data.frame("distance_region_to_target_gene" = values(distance)[["distance"]])
    )

    colnames(nearest.genes)[1:3] <- c(
        "target_gene_chrom",
        "target_gene_start",
        "target_gene_end"
    )
    colnames(nearest.genes)[4] <- "target_gene_name"
    colnames(nearest.genes)[5] <- "target"

    regionID <- regions.gr %>% as.data.frame %>% dplyr::select(1:3)
    regionID <- paste0(regionID[[1]],":",regionID[[2]],"-",regionID[[3]])
    out <- dplyr::bind_cols(
        data.frame("regionID" = regionID, stringsAsFactors = FALSE),
        nearest.genes[,c("target_gene_name", "target","distance_region_to_target_gene")]
    ) %>% tibble::as_tibble()
    out
}

get_region_target_gene_by_promoter_overlap <- function(
    regions.gr,
    genome,
    upstream,
    downstream
){

    # Get gene information
    gene.promoters <- get_promoter_regions(
        genome = genome,
        upstream = upstream,
        downstream = downstream
    )

    hits <- findOverlaps(
        query = regions.gr,
        subject = gene.promoters,
        ignore.strand = TRUE,
        select = "all"
    )

    # overlap region and promoter
    neargenes <- gene.promoters[subjectHits(hits)] %>% as.data.frame(row.names = NULL)

    regions.gr <- regions.gr[queryHits(hits)]
    neargenes <- cbind(
        neargenes[,c("seqnames",
                     "start",
                     "end",
                     "external_gene_name",
                     "ensembl_gene_id")]
    )

    colnames(neargenes)[1:3] <- c(
        "target_gene_chrom",
        "target_gene_start",
        "target_gene_end"
    )
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
get_region_target_gene_nearby.genes <- function(
    regions.gr,
    genome,
    num.flanking.genes
){

    names(regions.gr) <- make_names_from_granges(regions.gr)
    genes.gr <- get_gene_information(genome = genome, as.granges = TRUE)
    genes.gr$entrezgene <- NULL
    genes.gr <- unique(genes.gr)

    # Optimized version
    # Idea: vectorize search
    # 1) For all regions, get nearby gene
    # 2) check follow and overlapping genes recursively
    # 3) check precede and overlapping genes recursively
    # 4) map the positions based on min distance (L1)
    # The input data has to be at gene level and not transcript which would broke
    # some of the optimizations for which we remove the genes already evaluated

    # 1) For all regions, get nearby gene
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
    precede.genes <- get_region_target_gene_nearby.genes_aux(
        direction.fun = GenomicRanges::precede,
        nearest.idx = nearest.idx,
        genes.gr = genes.gr,
        regions.gr = regions.gr,
        num.flanking.genes = num.flanking.genes
    )

    # 2) check follow and overlapping genes recursively
    message("Identifying ",num.flanking.genes, " genes upstream of the region")
    follow.genes <- get_region_target_gene_nearby.genes_aux(
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
    ret <- get_region_target_gene_nearby.genes_addPos(ret, num.flanking.genes)
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
get_region_target_gene_nearby.genes_aux <- function(
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
    for (i in seq_len(num.flanking.genes)) {
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

        ret <- rbind(
            ret, # keep old results
            tibble::tibble(
                genes.gr[idx$subjectHits] %>% as.data.frame(),
                "ID" = names(regions.gr)[evaluating],
                "Distance" = ifelse(
                    start(regions.gr[evaluating]) < start(genes.gr[idx$subjectHits]), 1,-1) *
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
#' @importFrom plyr .
get_region_target_gene_nearby.genes_addPos <- function(
    ret,
    num.flanking.genes
){

    f <- function(pairs) {
        center <- which(abs(pairs$Distance) == min(abs(pairs$Distance)))[1]
        pos <- setdiff(-center:(nrow(pairs) - center), 0)
        pairs$Side <- ifelse(pos > 0, paste0("R", abs(pos)), paste0("L", abs(pos)))

        out <- pairs %>%
            dplyr::filter(
                pairs$Side %in% c(paste0("R", seq_len(num.flanking.genes)),
                                  paste0("L", seq_len(num.flanking.genes))
                )
            )

        if (nrow(out) < num.flanking.genes) {
            if (paste0("R", floor(num.flanking.genes)) %in% out$Side) {
                cts <- length(grep("L", sort(pairs$Side), value = TRUE))
                out <- pairs %>%
                    dplyr::filter(
                        .data$Side %in% c(
                            paste0("R", seq_len(num.flanking.genes - cts)),
                            grep("L", sort(out$Side), value = TRUE))
                    )
            } else {
                cts <- length(grep("R", sort(pairs$Side), value = TRUE))
                out <- pairs %>%
                    dplyr::filter(
                        pairs$Side %in%
                            c(paste0("L", seq_len(num.flanking.genes - cts)),
                              grep("R", sort(out$Side), value = TRUE)
                            )
                    )
            }
        }
        out <- out[order(out$Distance), ]
        return(out)
    }

    ret %>% group_by(.data$ID) %>% do(f(.))
}
