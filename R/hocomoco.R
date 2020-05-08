#' @title Get human TF list for a region using HOCOMOCO prediction
#' @description This function uses a pre-computed dataset for EPIC and HM450 Array that was created as follow:
#' For each HOCOMOCO human TF, the motif was search around the probe (+-250bp), and a binary matrix was created, with 1
#' if the motif was found, 0 if not. This function uses this pre-computed dataset to extend the probes to region, using the overlap probes
#' overlapping the regions as follows: for each region, get the probes within it and a motif will be one if at least one of the
#' overlapping probes has the motif (value 1 in the original dataset).
#' Then for each TF motifs found within the region, we select the TFs within the same TF family/subfamily since they
#' have similar binding motifs.
#' @importFrom dplyr pull filter %>% tbl
#' @importFrom plyr adply
#' @importFrom sesameData sesameDataGet
#' @importFrom GenomicRanges findOverlaps
#' @return A dataframe with region, a TF name and TF gene ensembl ID
#' @examples
#' \dontrun{
#'  regions.names <- c("chr1:60591-79592","chr4:40162197-43162198")
#'  region.tf <- get_tf_in_region(regions.names,
#'                  genome = "hg19",
#'                  arrayType = "450k",
#'                  classification = "subfamily")
#'
#'  regions.gr <- make_granges_from_names(regions.names)
#'  region.tf <- get_tf_in_region(regions.gr,
#'                  genome = "hg38",
#'                  arrayType = "450k",
#'                  classification = "family")
#' }
#' @export
#' @param region Region to map. Either a Granges or a named vertor
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param arrayType DNA methylation array "450k" or "EPIC"
#' @param classification TF classification to be used.
#' Either "subfamily" or "family".
get_tf_in_region <- function(region,
                            genome = c("hg19","hg38"),
                            arrayType = c("450k","EPIC"),
                            classification = c("subfamily","family")) {

    if (!requireNamespace("ELMER.data", quietly = TRUE)) {
        stop("ELMER.data is needed. Please install it.",
             call. = FALSE)
    }

    arrayType <- toupper(match.arg(arrayType))
    genome <- match.arg(genome)
    classification <- match.arg(classification)

    if(is(region,"character") | is(region,"factor")){
        regions.gr <- make_granges_from_names(region)
        region.names <- region
    } else if(is(region,"GenomicRanges")){
        regions.gr <- region
        region.names <- make_names_from_granges(regions.gr)
    }

    # get pre-computed data
    motifs.probes <- get(data(list = paste0("Probes.motif.",genome,".",arrayType),package = "ELMER.data"))
    motifs.tfs <- get(data(list = paste0("TF.",classification),package = "ELMER.data"))

    # Get probes regions for mapping the motifs
    probes.gr <- sesameDataGet(paste0(gsub("450K","HM450",arrayType),".",genome,".manifest"))
    probes.gr <- probes.gr[names(probes.gr) %in% rownames(motifs.probes)]

    # Find which probes overlap with the regions
    hits <- findOverlaps(regions.gr,probes.gr) %>% as.data.frame()

    if(nrow(hits) == 0) stop("No overlap found between regions and DNA methylation array found")

    # For each region overlapping get all probes overallping and their motifs
    df <- plyr::adply(unique(hits$queryHits),1, function(x){
        idx <- hits %>% filter(queryHits == x) %>% pull(subjectHits)
        probes <- names(probes.gr)[idx]
        motifs <- colnames(motifs.probes)[Matrix::colSums(motifs.probes[probes,,drop = FALSE]) > 0]
        tfs <- unlist(motifs.tfs[motifs]) %>% as.character() %>% unique()
        tibble::tibble(region.names[x],tfs)
    },.id = NULL,.progress = "time",.inform = TRUE)

    df$TF_id <- map_symbol_to_ensg(df$tfs)
    colnames(df) <- c("regionID","TF_external_gene_name","TF_ensembl_gene_id")
    return(df)
}
