#' @title Get human TF list for a region using JASPAR 2020 database and motifmatchr
#' @description Given a genomic region, this function obtains TFs within it using a motif search.
#' To this end, we use  a pre-computed dataset for EPIC and HM450 Array that was created as follows:
#' each JASPAR 2020 human TF motif is searched within region and a binary matrix is created,
#' with 1 if the motif was found, 0 if not.
#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedArray colSums
#' @param region A vector of region names or GRanges object with the DNA methylation regions to be scanned for the motifs
#' @param window.size Integer value to extend the regions. For example, a value of 50 will
#' extend 25 bp upstream and 25 downstream the region. Default is no increase
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param p.cutoff Motifmatcher p.cutoff. Default 1e-8.
#' @param cores Number of CPU cores to be used. Default 1.
#' @examples
#' \dontrun{
#'  regions.names <- c("chr1:79502-79592","chr4:43162098-43162198")
#'  region.tf <- get_tf_in_region(
#'                  region = regions.names,
#'                  genome = "hg38"
#'  )
#'  regions.names <- c("chr1:79592-79592","chr4:43162198-43162198")
#'  regions.gr <- make_granges_from_names(regions.names)
#'  region.tf <- get_tf_in_region(
#'                  region = regions.gr,
#'                  window.size = 25,
#'                  genome = "hg38"
#'  )
#' }
#' @export
get_tf_in_region <- function(
    region,
    window.size = 0,
    genome = c("hg19","hg38"),
    p.cutoff = 1e-8,
    cores = 1)
{

    check_package("JASPAR2020")
    check_package("TFBSTools")

    parallel <- register_cores(cores)

    if(is(region,"character") | is(region,"factor")){
        region.gr <- make_granges_from_names(region)
        region.names <- region
    } else if(is(region,"GenomicRanges")){
        region.gr <- region
        region.names <- make_names_from_granges(region)
    }

    region.gr <- region.gr + (window.size/2)
    # region <- resize(region,width = 50,fix = "center")

    if (min(IRanges::width(region.gr)) < 2)
        stop("Minimun region size is 2, please set window.size argument")

    genome <- match.arg(genome)

    opts <- list()
    opts[["species"]] <- 9606 # homo sapies
    # opts[["all_versions"]] <- TRUE
    PFMatrixList <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
    motifs.names <- lapply(PFMatrixList, function(x)(TFBSTools::name(x)))
    names(PFMatrixList) <- motifs.names
    PFMatrixList <- PFMatrixList[grep("::|var",motifs.names,invert = TRUE)]

    message("Evaluating ", length(PFMatrixList), " JASPAR Human TF motifs")
    message("This may take a while...")
    motif.matrix <- motifmatchr::matchMotifs(
        pwms = PFMatrixList,
        subject = region.gr,
        genome = genome,
        p.cutoff = p.cutoff
    ) %>% SummarizedExperiment::assay()
    rownames(motif.matrix) <- region.names

    # remove motifs not found in any regions
    motif.matrix <- motif.matrix[,DelayedArray::colSums(motif.matrix) > 0, drop = FALSE]

    if(is(motif.matrix, "lgCMatrix")){
        motif.matrix <-  motif.matrix %>% as.matrix() %>% as.data.frame()
    }

    message("Preparing output")
    motifs.probes.df <- plyr::alply(
        colnames(motif.matrix),
        .margins = 1,
        function(colum.name){
            colum <- motif.matrix[,colum.name, drop = FALSE]
            regions <- rownames(colum)[which(colum %>% pull > 0)];
            tfs <- colum.name
            expand.grid(regions,tfs)
        }, .progress = "time",.parallel = parallel)
    motifs.probes.df <- dplyr::bind_rows(motifs.probes.df)
    colnames(motifs.probes.df) <- c("regionID","TF_external_gene_name")

    motifs.probes.df$TF <- map_symbol_to_ensg(motifs.probes.df$TF_external_gene_name)
    motifs.probes.df <- motifs.probes.df %>% na.omit
    return(motifs.probes.df %>% unique)

}



