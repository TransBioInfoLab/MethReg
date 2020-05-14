#' @title Get human TF list for a region using HOCOMOCO prediction
#' @description Given a genomic region, this function obtains TFs that bind close to it (+-250bp).
#' To this end, we used a pre-computed dataset for EPIC and HM450 Array that was created as follows:
#' for each HOCOMOCO human TF, the motif was searched around the probe (+-250bp), and a binary matrix was created,
#' with 1 if the motif was found, 0 if not.
#' The function \code{get_tf_in_region} uses this pre-computed dataset to link regions to TFs:
#' for each region, obtain the probes within it and a motif will be selected
#' if at least one of the probes within the region has the motif (i.e. value 1 in the original dataset).
#' Then for each TF motifs found within the region, we select the TFs within the same TF family/subfamily since they
#' have similar binding motifs.
#' @importFrom dplyr pull filter %>% tbl
#' @importFrom plyr adply
#' @importFrom sesameData sesameDataGet
#' @importFrom GenomicRanges findOverlaps
#' @importFrom methods is
#' @importFrom utils data
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
#' @param cores Number of cores to be used. Default 1.
#' @importFrom Matrix colSums
get_tf_in_region <- function(region,
                              genome = c("hg19","hg38"),
                              arrayType = c("450k","EPIC"),
                              classification = c("subfamily","family"),
                              cores = 1) {

    check_package("ELMER.data")

    arrayType <- match.arg(arrayType)
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
    message("Loading pre-computed probes TFBS data")
    motifs.probes <- get(data(list = paste0("Probes.motif.",genome,".",toupper(arrayType)),package = "ELMER.data",envir = environment()))
    message("Transforming probes TFBS data to regions TFBS data")
    motifs.probes <- map_motif_probes_to_regions(motifs.probes = motifs.probes,
                                                 genome = genome,
                                                 arrayType = arrayType,
                                                 regions.gr = regions.gr)

    message("Mapping regions to HOCOMOCO TFBS. This may take a while...")
    parallel <- register_cores(cores)

    message("Mapping regions to TFs")
    motifs.tfs <- get(data(list = paste0("TF.",classification),package = "ELMER.data",envir = environment()))

    # For each region/probe get the TFBS in region as a data.frame
    motifs.probes.df <- plyr::alply(
        colnames(motifs.probes),
        .margins = 1,
        function(colum.name){
            colum <- motifs.probes[,colum.name,drop = FALSE]
            regions <- rownames(colum)[which(colum %>% pull > 0)];
            tfs <- unlist(motifs.tfs[colum.name]) %>% as.character() %>% unique()
            expand.grid(regions,tfs)
        }, .progress = "time",.parallel = parallel)
    motifs.probes.df <- dplyr::bind_rows(motifs.probes.df)
    colnames(motifs.probes.df) <- c("motif","TF_external_gene_name")

    message("Preparing output")
    # Merge results of region and motif to motif TF

    motifs.probes.df$TF_ensembl_gene_id <- map_symbol_to_ensg(motifs.probes.df$TF_external_gene_name)
    motifs.probes.df <- motifs.probes.df %>% na.omit
    motifs.probes.df %>% unique
}

# Since motifs.probes is not region based, we need to
# combine the rows of the sparce matrix to the overlapping regions
map_motif_probes_to_regions <- function(motifs.probes,
                                        genome,
                                        arrayType,
                                        regions.gr){

    # Get probes regions for mapping the motifs
    probes.gr <- get_met_probes_info(genome = genome,arrayType = arrayType)
    probes.gr <- probes.gr[match(rownames(motifs.probes),names(probes.gr))]
    regions.gr <- unique(regions.gr) # we cannot have duplicated regions

    # Find which probes overlap with the regions
    hits <- findOverlaps(regions.gr, probes.gr, ignore.strand = TRUE) %>% as.data.frame()
    if(nrow(hits) == 0) stop("No overlap found between regions and DNA methylation array found")

    region.with.more.than.one.probe <- unique(hits$queryHits[duplicated(hits$queryHits)])
    unique.hits <- hits[!hits$queryHits %in% region.with.more.than.one.probe,]

    motif.matrix <- NULL
    # Do we have probes mapped to unique regions, if so copy probes and rename
    # probes to regions
    if(nrow(unique.hits) > 0){
        motif.matrix <- motifs.probes[unique.hits$subjectHits,, drop = FALSE] %>% as.matrix()
        rownames(motif.matrix) <- make_names_from_granges(regions.gr[unique.hits$queryHits])
    }

    # Do we have regions overlapping with multiple probes ?
    if(length(region.with.more.than.one.probe) > 0){
        non.unique.hits <- hits[hits$queryHits %in% region.with.more.than.one.probe,]
        non.unique.motifs <- plyr::adply(unique(non.unique.hits$queryHits),1, function(x){
            idx <- hits %>% filter(queryHits == x) %>% pull(subjectHits)
            probes <- names(probes.gr)[idx]
            Matrix::colSums(motifs.probes[probes,,drop = FALSE])
        },.id = NULL,.progress = "time",.inform = TRUE)

        rownames(non.unique.motifs) <- make_names_from_granges(regions.gr[unique(non.unique.hits$queryHits)])

        if(is.null(motif.matrix)) {
            motif.matrix <- non.unique.motifs
        } else {
            motif.matrix <-  rbind(motif.matrix, non.unique.motifs)
        }
    }

    if(is(motifs.probes,"ngCMatrix")){
        motif.matrix <-  motif.matrix %>% as.matrix() %>% as.data.frame()
    }

    # remove motifs not found in any regions
    motif.matrix <- motif.matrix[,colSums(motif.matrix) > 0]
    motif.matrix
}
