#' @title Get human TFs for regions by scanning it with motifmatchr using
#' JASPAR 2020 database
#' @description Given a genomic region, this function uses motifmatchr
#' and JASPAR2020
#' to scan the region for 554 human transcription factors binding sites. There is also
#' an option (argument \code{window.size}) to extend the scanning region
#' before performing the search, which
#' by default is 0 (do not extend)
#' @return A data frame with the following information: regionID, TF symbol, TF ensembl ID
#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedArray colSums
#' @importFrom IRanges width
#' @param region A vector of region names or GRanges object with the DNA
#' methylation regions to be scanned for the motifs
#' @param window.size Integer value to extend the regions.
#' For example, a value of 50 will
#' extend 25 bp upstream and 25 bp downstream the region.
#' The default is not to increase the scanned region.
#' @param genome Human genome of reference "hg38" or "hg19".
#' @param p.cutoff motifmatchr p.cutoff. Default 1e-8.
#' @param cores Number of CPU cores to be used. Default 1.
#' @param verbose A logical argument indicating if
#' messages output should be provided.
#' @examples
#'  regions.names <- c("chr3:189631389-189632889","chr4:43162098-43163498")
#'  region.tf <- get_tf_in_region(
#'                  region = regions.names,
#'                  genome = "hg38"
#'  )
#' @export
get_tf_in_region <- function(
    region,
    window.size = 0,
    genome = c("hg19","hg38"),
    p.cutoff = 1e-8,
    cores = 1,
    verbose = FALSE
) {

    check_package("JASPAR2020")
    check_package("TFBSTools")

    parallel <- register_cores(cores)

    if (is(region,"character") | is(region,"factor")) {
        region.gr <- make_granges_from_names(region)
        region.names <- region
    } else if (is(region,"GenomicRanges")) {
        region.gr <- region
        region.names <- make_names_from_granges(region)
    }

    region.gr <- region.gr + (window.size/2)
    # region <- resize(region,width = 50,fix = "center")

    if (min(IRanges::width(region.gr)) < 8) {
        stop("Minimun region size is 8, please set window.size argument")
    }

    genome <- match.arg(genome)

    opts <- list()
    opts[["species"]] <- 9606 # homo sapies
    # opts[["all_versions"]] <- TRUE
    PFMatrixList <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
    motifs.names <- lapply(PFMatrixList, function(x)(TFBSTools::name(x)))
    names(PFMatrixList) <- motifs.names
    PFMatrixList <- PFMatrixList[grep("::|var",motifs.names,invert = TRUE)]

     if(verbose)  message("Evaluating ", length(PFMatrixList), " JASPAR Human TF motifs")
     if(verbose)  message("This may take a while...")
    suppressWarnings({
        motif.matrix <- motifmatchr::matchMotifs(
            pwms = PFMatrixList,
            subject = region.gr,
            genome = genome,
            p.cutoff = p.cutoff
        ) %>% SummarizedExperiment::assay()
    })
    rownames(motif.matrix) <- region.names

    # remove motifs not found in any regions
    motif.matrix <- motif.matrix[,DelayedArray::colSums(motif.matrix) > 0, drop = FALSE]

    if (ncol(motif.matrix) == 0) {
        message("No motifs found")
        return(NULL)
    }

    if (is(motif.matrix, "lgCMatrix")) {
        motif.matrix <- motif.matrix[!duplicated(rownames(motif.matrix)),]
        motif.matrix <- motif.matrix %>% as.matrix() %>% as.data.frame()
    }

     if(verbose)  message("Preparing output")
    motifs.probes.df <- plyr::alply(
        colnames(motif.matrix),
        .margins = 1,
        function(colum.name){
            colum <- motif.matrix[,colum.name, drop = FALSE]
            regions <- rownames(colum)[which(colum %>% pull > 0)];
            tfs <- colum.name
            expand.grid(regions,tfs,stringsAsFactors = FALSE)
        }, .progress = "time",.parallel = parallel)
    motifs.probes.df <- dplyr::bind_rows(motifs.probes.df)
    colnames(motifs.probes.df) <- c("regionID","TF_external_gene_name")

    motifs.probes.df$TF <- map_symbol_to_ensg(
        motifs.probes.df$TF_external_gene_name
    )

    motifs.probes.df <- motifs.probes.df %>% na.omit
    return(motifs.probes.df %>% unique)

}
