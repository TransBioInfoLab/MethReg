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
#' genes within a window around the region ("window");
#' or a fixed number genes upstream
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
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom tidyr unite
#' @importFrom dplyr select filter bind_rows
#' @importFrom methods is
#' @importFrom utils head
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
#'                       method = "genes.promoter.overlap"
#'  )
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
#'  we recommend setting either \code{method = "window"}
#'  or \code{method = "nearby.genes"}.
#'
#'  Note that because \code{method = "window"} or
#'  \code{method = "nearby.genes"} are
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

    if (is(regions.gr,"SummarizedExperiment")) {
        regions.gr <- SummarizedExperiment::rowRanges(regions.gr)
    }

    if (!is(regions.gr,"GRanges")) stop("regions.gr must be a GRanges")

    if (method != "genes.promoter.overlap" & rm.promoter.regions.from.distal.linking) {

        message("Removing regions overlapping promoter regions")
        regions.gr <- subset_by_non_promoter_regions(
            regions.gr = regions.gr,
            genome =  genome,
            upstream = promoter.upstream.dist.tss,
            downstream = promoter.downstream.dist.tss
        )

        if (length(regions.gr) == 0) {
            stop("After removing promoter regions, regions.gr is empty")
        }
    }

    if (method == "genes.promoter.overlap") {

        message("Mapping regions to the closest gene")
        out <- get_region_target_gene_by_promoter_overlap(
            regions.gr = regions.gr,
            genome = genome,
            upstream = promoter.upstream.dist.tss,
            downstream = promoter.downstream.dist.tss
        )

    } else if (method == "window") {

        message("Mapping regions to genes within a window of size: ", window.size, " bp")
        out <- get_region_target_gene_window(regions.gr, genome, window.size)

    } else if (method == "nearby.genes") {

        out <- get_region_target_gene_nearby.genes(
            regions.gr = regions.gr,
            genome = genome,
            num.flanking.genes = num.flanking.genes
        )

    } else if (method == "closest.gene") {
        out <- get_region_target_gene_closest(
            regions.gr = regions.gr,
            genome = genome
        )
    }
    out <- get_distance_region_target(out, genome = genome)
    out$target_tss_pos_in_relation_to_region <- NULL
    out$region_pos_in_relation_to_gene_tss <- NULL

    out <- out %>% dplyr::rename(target_symbol = .data$target_gene_name)

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

    nearest.genes <- dplyr::bind_cols(
        nearest.genes[,c(
            "seqnames",
            "start",
            "end",
            "external_gene_name",
            "ensembl_gene_id"
        )],
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
        nearest.genes[,c("target_gene_name", "target")]
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
        neargenes[,c(
            "seqnames",
            "start",
            "end",
            "external_gene_name",
            "ensembl_gene_id"
        )]
    )

    colnames(neargenes)[1:3] <- c(
        "target_gene_chrom",
        "target_gene_start",
        "target_gene_end"
    )
    colnames(neargenes)[4] <- "target_gene_name"
    colnames(neargenes)[5] <- "target"

    regionID <- regions.gr %>% data.frame %>% dplyr::select(1:3)
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
        regions.gr[queryHits(overlap)] %>% end
    )


    genes.overlapping <- geneAnnot[subjectHits(overlap)] %>% as.data.frame(row.names = NULL)
    colnames(genes.overlapping) <- paste0("gene_",colnames(genes.overlapping))
    colnames(genes.overlapping)[grep("ensembl_gene_id",colnames(genes.overlapping))] <- "target"
    colnames(genes.overlapping)[grep("external_gen",colnames(genes.overlapping))] <- "target_gene_name"
    genes.overlapping <- genes.overlapping[,grep("target",colnames(genes.overlapping))]

    out <- dplyr::bind_cols(
        data.frame(
            "regionID" = regionID,stringsAsFactors = FALSE),
        genes.overlapping
    )  %>% tibble::as_tibble()
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
        # nearest gene to the region
        genes.gr[nearest.idx %>% subjectHits] %>% data.frame(),
        # region being evaluated
        "ID" = names(regions.gr)[nearest.idx %>% queryHits]
    )

    message("Identifying ",num.flanking.genes, " genes downstream of the region")
    precede.genes <- get_region_target_gene_nearby.genes_aux(
        direction.fun = GenomicRanges::precede,
        nearest.idx = nearest.idx,
        genes.gr = genes.gr,
        regions.gr = regions.gr,
        num.flanking.genes = num.flanking.genes
    )

    # 2) check follow and overlapping genes recursively
    message("Identifying ", num.flanking.genes, " genes upstream of the region")
    follow.genes <- get_region_target_gene_nearby.genes_aux(
        direction.fun = GenomicRanges::follow,
        nearest.idx = nearest.idx,
        genes.gr = genes.gr,
        regions.gr = regions.gr,
        num.flanking.genes = num.flanking.genes
    )

    ret <- rbind(closest.genes, precede.genes, follow.genes) %>% unique

    ret <- ret[,c(
        "ID",
        "ensembl_gene_id",
        grep("external_gene_", colnames(ret), value = TRUE)
    )]

    message("Identifying gene position for each region")
    colnames(ret)[1:3] <- c("regionID", "target", "target_gene_name")

    ret <- get_distance_region_target(ret, genome = genome)

    ret <- ret %>% dplyr::group_by(.data$regionID,.data$target_tss_pos_in_relation_to_region) %>%
        filter(
            abs(.data$distance_region_target_tss) <=
                (abs(.data$distance_region_target_tss) %>% sort %>% head(num.flanking.genes) %>% max)
        ) %>% dplyr::ungroup()

    return(ret)
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
    idx <- nearest.idx %>% data.frame()

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
                genes.gr[idx$subjectHits] %>% data.frame(),
                "ID" = names(regions.gr)[evaluating],
            )
        )
        pb$tick()
    }
    ret <- ret[!duplicated(ret[,c("ensembl_gene_id","ID")]),]
    pb$terminate()
    return(ret)
}



#' @title Calcule distance (in bp) between DNAm region and target gene TSS
#' @description Given a dataframe with a region ("regionID") and a
#' target gene ("target"), returns a data frame with the distance included.
#' @examples
#' region.target <- data.frame(
#'   "regionID" = "chr21:17511749-17511750",
#'   "target" =  "ENSG00000215326"
#' )
#' region.target.with.distance <- get_distance_region_target(
#'     region.target = region.target
#')
#' @noRd
get_distance_region_target <- function(
    region.target,
    genome = c("hg38","hg19")
){
    genome <- match.arg(genome)

    if ( !all( c("target","regionID") %in%  colnames(region.target))) {
        stop("Input requires columns regionID (chrx:start:end) and target (ENSG)")
    }
    region.target.only <- region.target %>%
        dplyr::filter(!is.na(.data$regionID)) %>%
        dplyr::filter(!is.na(.data$target)) %>%
        dplyr::select(c("target","regionID")) %>%
        unique()


    # resize is used to keep only the TSS to calculate the distance to TSS
    genes.gr <- get_gene_information(
        genome = genome,
        as.granges =  TRUE
    ) %>% resize(1)

    # We only need to calculate the distance of genes in the input
    genes.gr <- genes.gr[genes.gr$ensembl_gene_id %in% region.target.only$target]

    # Adding new information
    region.target.only <- region.target.only %>%
        cbind(
            data.frame(
                "distance_region_target_tss" = NA,
                "target_tss_pos_in_relation_to_region" = NA,
                "region_pos_in_relation_to_gene_tss" = NA
            )
        )

    # If the gene has no information the distance is NA
    region.target.no.info <- region.target.only %>%
        dplyr::filter(!.data$target %in% genes.gr$ensembl_gene_id)

    # If the gene has information the distance will be calculated
    region.target.info <- region.target.only %>%
        dplyr::filter(.data$target %in% genes.gr$ensembl_gene_id)

    regions.gr <- make_granges_from_names(
        names = region.target.info$regionID
    )


    idx <- match(region.target.info$target,genes.gr$ensembl_gene_id)
    dist <- distance(
        regions.gr,
        genes.gr[idx]
    )

    region.target.info$distance_region_target_tss <- dist
    region.target.info$region_pos_in_relation_to_gene_tss <- ifelse(
        as.logical(strand(genes.gr[idx]) != "-"),
        ifelse(start(regions.gr) < start(genes.gr[idx]), "upstream", "downstream"),
        ifelse(start(regions.gr) < start(genes.gr[idx]), "downstream", "upstream")
    )

    region.target.info$target_tss_pos_in_relation_to_region <- ifelse(
        start(regions.gr) < start(genes.gr[idx]), "right", "left"
    )


    region.target.info$distance_region_target_tss <-
        region.target.info$distance_region_target_tss * ifelse(region.target.info$region_pos_in_relation_to_gene_tss == "downstream",1, -1)

    # output both results together
    region.target.only <- plyr::rbind.fill(region.target.info, region.target.no.info)

    # using target and region keys, map the distance to the original input
    region.target$distance_region_target_tss <-
        region.target.only$distance_region_target_tss[
            match(
                paste0(
                    region.target$regionID,
                    region.target$target
                ),
                paste0(
                    region.target.only$regionID,
                    region.target.only$target
                )
            )
            ]

    region.target$target_tss_pos_in_relation_to_region <-
        region.target.only$target_tss_pos_in_relation_to_region[
            match(
                paste0(
                    region.target$regionID,
                    region.target$target
                ),
                paste0(
                    region.target.only$regionID,
                    region.target.only$target
                )
            )
            ]

    region.target$region_pos_in_relation_to_gene_tss <-
        region.target.only$region_pos_in_relation_to_gene_tss[
            match(
                paste0(
                    region.target$regionID,
                    region.target$target
                ),
                paste0(
                    region.target.only$regionID,
                    region.target.only$target
                )
            )
            ]

    return(region.target)
}




