#' @title Map TF and target genes using regulon databases or
#' any user provided target-tf. Maps
#' TF to the DNAm region with TFBS using JASPAR2020 TFBS.
#' @description This function wraps two other functions
#' \code{get_region_target_gene} and \code{get_tf_in_region} from the package.
#'
#' @importFrom tidyr separate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @param region A Granges or a named vector with
#' regions (i.e "chr21:100002-1004000")
#' @param genome Human genome reference "hg38" or "hg19"
#' @param regulons.min.confidence Minimun confidence score  ("A", "B","C","D", "E")
#' classifying regulons based on their quality from Human DoRothEA database
#'  \link[dorothea]{dorothea_hs}. The default minimun confidence score is "B".
#' @param motif.search.window.size Integer value to extend the regions.
#' For example, a value of 50 will
#' extend 25 bp upstream and 25 downstream the region. Default is no increase
#' @param motif.search.p.cutoff motifmatchr pvalue cut-off. Default 1e-8.
#' @param cores Number of CPU cores to be used. Default 1.
#' @param tf.target A dataframe with tf and target columns. If not provided,
#' \link[dorothea]{dorothea_hs} will be used.
#' @return A data frame with TF, target and RegionID information.
#' @examples
#' triplet <- create_triplet_regulon_based(
#'    region = c("chr1:69591-69592", "chr1:898803-898804"),
#'    motif.search.window.size = 50,
#'    regulons.min.confidence = "B",
#'      motif.search.p.cutoff = 0.05
#' )
#' @importFrom SummarizedExperiment rowRanges
#' @export
create_triplet_regulon_based <- function(
    region,
    genome = c("hg38","hg19"),
    regulons.min.confidence = "B",
    motif.search.window.size = 0,
    motif.search.p.cutoff = 1e-8,
    cores = 1,
    tf.target
){

    regulons.min.confidence <- match.arg(
        arg = regulons.min.confidence,
        choices =   c("A", "B", "C", "D", "E")
    )

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
        tf.target <- get_regulon_dorothea(min.confidence = regulons.min.confidence)
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

    if(nrow(triplet) == 0){
        stop("No triplets found")
    }

    triplet <- get_distance_region_target(triplet)

    message("Removing regions and target genes from different chromosomes")
    triplet <- triplet %>% dplyr::filter(!is.na(.data$distance_region_target))

    return(triplet)
}
