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
#' @import dbplyr
#' @return A dataframe with region, a TF name and TF gene ensembl ID
#' @examples
#' \dontrun{
#'  regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#'  get_region_tf("chr22:18267969-18268249",
#'                  genome = "hg19",
#'                  arrayType = "450k",
#'                  classification = "subfamily")
#' }
#' @export
#' @importFrom sesameData sesameDataGet
#' @importFrom GenomicRanges findOverlaps
get_region_tf <- function(region,
                            genome = c("hg19","hg38"),
                            arrayType = c("450k","EPIC"),
                            classification = c("subfamily","family")) {

    if (!requireNamespace("ELMER.data", quietly = TRUE)) {
        stop("ELMER.data is needed. Please install it.",
             call. = FALSE)
    }
    if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
        stop("TCGAbiolinks is needed. Please install it.",
             call. = FALSE)
    }

    arrayType <- toupper(match.arg(arrayType))
    genome <- match.arg(genome)
    classification <- match.arg(classification)

    if(is(region,"character") | is(region,"factor")){
        regions.gr <- make_granges_from_names(region)
    } else  if(is(region,"GenomicRanges")){
        stop("Region must be a list of character")
    }

    motifs.probes <- get(data(list = paste0("Probes.motif.",genome,".",arrayType),package = "ELMER.data"))
    motifs.tfs <- get(data(list = paste0("TF.",classification),package = "ELMER.data"))

    probes.gr <- sesameDataGet(paste0(gsub("450K","HM450",arrayType),".",genome,".manifest"))
    probes.gr <- probes.gr[names(probes.gr) %in% rownames(motifs.probes)]
    hits <- findOverlaps(regions.gr,probes.gr) %>% as.data.frame()
    df <- plyr::adply(unique(hits$queryHits),1, function(x){
        idx <- hits %>% filter(queryHits == x) %>% pull(subjectHits)
        probes <- names(probes.gr)[idx]
        motifs <- colnames(motifs.probes)[Matrix::colSums(motifs.probes[probes,,drop = FALSE]) > 0]
        tfs <- unlist(motifs.tfs[motifs]) %>% as.character() %>% unique()
        tibble::tibble(region[x],tfs)
    },.id = NULL,.progress = "time",.inform = TRUE)
    genome.info <- TCGAbiolinks::get.GRCh.bioMart(genome)
    df$TF_id <- genome.info$ensembl_gene_id[match(df$tfs,genome.info$external_gene_name)]
    colnames(df) <- c("regionID","TF_external_gene_name","TF_ensembl_gene_id")
    return(df)
}
