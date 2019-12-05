#' @title Get human TF list for region using hocomoco prediction
#' @description Get list of TF mapping to regions using
#' TF motifs to regions
#' @importFrom dplyr pull filter %>% tbl
#' @importFrom plyr adply
#' @import dbplyr
#' @return A dataframe with TF, target and correlation
#' @examples
#' \dontrun{
#'  regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#'  mapTFBSHocomoco("chr22:18267969-18268249",
#'                  genome = "hg19",
#'                  arrayType = "450k",
#'                  classification = "subfamily")
#' }
#' @export
#' @importFrom sesameData sesameDataGet
#' @importFrom GenomicRanges findOverlaps
mapTFBSHocomoco <- function(region,
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
        regions.gr <- makeGrangesFromNames(region)
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

#' Create granges from name
#' @importFrom tidyr separate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @examples
#' regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#' regions.gr <- makeGrangesFromNames(regions.names)
makeGrangesFromNames <- function(names){
    names %>%
        data.frame %>%
        separate(col = ".",into = c("chr","start","end")) %>%
        makeGRangesFromDataFrame()
}
