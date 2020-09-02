#' @title Map TF and target genes using regulon databases or
#' any user provided target-tf, to
#' TF to the DNAm region using JASPAS2020 TFBS.
#' @description This function wraps two other functions
#' \code{get_region_target_gene} and \code{get_tf_in_region} from the package.
#' This function will map a region to a target gene using three methods
#' (mapping to the closest gene,
#' mapping to any gene within a given window of distance, or mapping to a
#' fixed number of nearby genes upstream or downstream).
#' To find TFs binding to the region, JASPAR2020 is used.
#'
#' @importFrom tidyr separate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @param region A Granges or a named vector with
#' regions (i.e "chr21:100002-1004000")
#' @param genome Human genome reference "hg38" or "hg19"
#' @param regulons.min.confidence Minimun confidence score  ("A", "B","C","D", "E")
#' classifying regulons based on their quality from Human DoRothEA database
#'  \link[dorothea]{dorothea_hs}.
#' @param motif.search.window.size Integer value to extend the regions.
#' For example, a value of 50 will
#' extend 25 bp upstream and 25 downstream the region. Default is no increase
#' @param motif.search.p.cutoff motifmatchr pvalue cut-off. Default 1e-8.
#' @param cores Number of CPU cores to be used. Default 1.
#' @param tf.target A dataframe with tf and target columns. If not provided,
#' \link[dorothea]{dorothea_hs} will be used.
#' @return A data frame with TF, target and RegionID information.
#' @examples
#' \dontrun{
#' data("dna.met.chr21")
#' dnam.regions <- make_se_from_dnam_probes(dna.met.chr21)
#' triplet <- create_triplet_regulon_based(
#'    region = rownames(dnam.regions)[1:100],
#'    motif.search.window.size = 50,
#'    min.confidence = "B"
#' )
#' }
#' @noRd
#' @importFrom SummarizedExperiment rowRanges
create_triplet_regulon_based <- function(
    region,
    genome = c("hg38","hg19"),
    min.confidence = c("A", "B","C","D", "E"),
    motif.search.window.size = 0,
    motif.search.p.cutoff = 1e-8,
    cores = 1,
    tf.target
){

    min.confidence <- match.arg(min.confidence)
    genome <- match.arg(genome)

    if (is(region, "character") | is(region, "factor")) {
        region.gr <- make_granges_from_names(region)
        region.names <- region
    } else if (is(region, "GenomicRanges")) {
        region.gr <- region
        region.names <- make_names_from_granges(region.gr)
    } else if (is(region, "SummarizedExperiment")) {
        region.gr <- rowRanges(region)
        region.names <- make_names_from_granges(region.gr)
    }

    message("Mapping target and TF genes")
    if(missing(tf.target)){
        tf.target <- get_regulon_dorothea(min.confidence = min.confidence)
    } else {
        # check regulons input data
        cols <- c("tf", "target")
        if(!all(cols %in% colnames(tf.target))){
            stop("regulons must have columns tf and target")
        }
        tf.target$tf_ensg <- map_symbol_to_ensg(tf.target$tf)
        tf.target$target_ensg <- map_symbol_to_ensg(tf.target$target)
    }

    tf.target$target_name <- tf.target$target
    tf.target$target <- tf.target$target_ensg
    tf.target$TF <- tf.target$tf_ensg

    message("Looking for TFBS")
    region.tf <- get_tf_in_region(
        region = region.gr,
        genome = genome,
        window.size = motif.search.window.size,
        p.cutoff = motif.search.p.cutoff,
        cores = cores
    )
    triplet <- dplyr::inner_join(tf.target, region.tf)

    triplet <- get_distance_region_target(triplet)

    message("Removing regions and target genes from different chromosomes")
    triplet <- triplet %>% dplyr::filter(!is.na(.data$distance_region_target))

    return(triplet)
}
