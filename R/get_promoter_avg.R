#' @title Access human TF from Lambert et al 2018
#' @description For each probe within the same promoter region,
#' take the mean or median beta-values
#' @return A list with a data frame of the new regions and merged ones,
#' and the calculated mean/median beta-value matrix
#' @export
#' @examples
#' data("dna.met.chr21")
#' promoter.avg <- get_promoter_avg(dna.met.chr21, genome = "hg19", arrayType = "450k")
get_promoter_avg <- function(
    dnam,
    genome,
    arrayType,
    cores = 1)
{

    # We will start by defining the promoter regions
    message("o Get promoter regions for ", genome)
    promoter.gr <- get_promoter_regions(genome)
    message("oo Number of promoter regions: ", length(promoter.gr))

    # For each promoter region we will then
    # take the mean DNA methylation beta-values of all
    # probes within it

    # Get probes regions for mapping the motifs
    message("o Get DNA methylation regions overlapping promoter regions")

    # If inout are probes we need to map to regions
    if(any(grepl("cg", rownames(dnam)))){
        dnam <- map_probes_to_regions(dnam, genome = genome, arrayType = arrayType)
    }

    probes.gr <- make_granges_from_names(rownames(dnam))


    # Find which probes overlap with the regions
    hits <- findOverlaps(promoter.gr, probes.gr, ignore.strand = TRUE) %>% as.data.frame()
    if(nrow(hits) == 0) stop("No overlap found between promoter regions and DNA methylation array found")

    region.with.more.than.one.probe <- unique(hits$queryHits[duplicated(hits$queryHits)])
    unique.hits <- hits[!hits$queryHits %in% region.with.more.than.one.probe,]

    promoter.matrix <- NULL
    # Do we have probes mapped to unique promoter regions, if so copy probes and rename
    # probes to regions
    if(nrow(unique.hits) > 0){
        promoter.matrix <- dnam[unique.hits$subjectHits,, drop = FALSE] %>% as.matrix()
        rownames(promoter.matrix) <- make_names_from_granges(promoter.gr[unique.hits$queryHits])
    }

    message("o Get mean DNA methylation of regions overlapping each promoter region")
    parallel <- register_cores(cores)

    # Do we have regions overlapping with multiple probes ?
    if(length(region.with.more.than.one.probe) > 0){
        non.unique.hits <- hits[hits$queryHits %in% region.with.more.than.one.probe,]

        non.unique.promoter <- plyr::adply(
            unique(non.unique.hits$queryHits),
            .margins = 1,
            function(x){
                idx <- hits %>% filter(queryHits == x) %>% pull(subjectHits)
                rows <- make_names_from_granges(probes.gr[idx])
                Matrix::colMeans(dnam[rows,,drop = FALSE],na.rm = TRUE)
            }, .id = NULL, .parallel = parallel , .progress = "time", .inform = TRUE)

        rownames(non.unique.promoter) <- make_names_from_granges(promoter.gr[unique(non.unique.hits$queryHits)])

        if(is.null(promoter.matrix)) {
            promoter.matrix <- non.unique.promoter
        } else {
            promoter.matrix <-  rbind(promoter.matrix, non.unique.promoter)
        }
    }

    promoter.matrix %>% as.matrix
}


#' @importFrom GenomicRanges promoters strand strand<-
get_promoter_regions <- function(
    genome,
    upstream = 2000,
    downstream = 2000){

    genes <- get_gene_information(genome = genome, as.granges = TRUE)
    promoters.gr <- promoters(genes, upstream = upstream, downstream = downstream)
    strand(promoters.gr) <- "*"
    return(promoters.gr %>% unique)
}
