#' @title Get human TFs for regions by either scanning it with motifmatchr using
#' JASPAR 2020 database or overlapping with TF chip-seq from user input
#' @description Given a genomic region, this function maps TF in regions
#' using two methods: 1) using motifmatchr nd JASPAR2022 to scan the
#' region for 554 human transcription factors
#' binding sites. There is also  an option (argument \code{window.size})
#' to extend the scanning region before performing the search, which
#' by default is 0 (do not extend).
#' 2) Using user input  TF chip-seq to check for overlaps between region
#' and TF peaks.
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
#' @param TF.peaks.gr  A granges with TF peaks to be overlaped with input region
#' Metadata column expected "id" with TF name. Default NULL. Note that Remap catalog
#' can be used as shown in the examples.
#' @examples
#'  regions.names <- c("chr3:189631389-189632889","chr4:43162098-43163498")
#'  region.tf <- get_tf_in_region(
#'                  region = regions.names,
#'                  genome = "hg38"
#'  )
#'
#' \dontrun{
#'    library(ReMapEnrich)
#'    demo.dir <- "~/ReMapEnrich_demo"
#'    dir.create(demo.dir, showWarnings = FALSE, recursive = TRUE)
#'    # Use the function DowloadRemapCatalog
#'    remapCatalog2018hg38 <- downloadRemapCatalog(demo.dir, assembly = "hg38")
#'    # Load the ReMap catalogue and convert it to Genomic Ranges
#'    remapCatalog <- bedToGranges(remapCatalog2018hg38)
#'    regions.names <- c("chr3:189631389-189632889","chr4:43162098-43163498")
#'    region.tf.remap <- get_tf_in_region(
#'                    region = regions.names,
#'                    genome = "hg38",
#'                    TF.peaks.gr = remapCatalog
#'    )
#'  }
#' @export
get_tf_in_region <- function(
    region,
    window.size = 0,
    genome = c("hg19","hg38"),
    p.cutoff = 1e-8,
    cores = 1,
    TF.peaks.gr = NULL,
    verbose = FALSE
) {

    check_package("JASPAR2022")
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
    # region <- region.gr %>% resize(width(.) + window.size, fix = "center")

    genome <- match.arg(genome)

    if (is.null(TF.peaks.gr)) {

        if (min(IRanges::width(region.gr)) < 8) {
            stop("Minimun region size is 8, please set window.size argument")
        }

        opts <- list()
        opts[["species"]] <- 9606 # homo sapiens
        # opts[["all_versions"]] <- TRUE
        PFMatrixList <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
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
        colnames(motifs.probes.df) <- c("regionID","TF_symbol")

        motifs.probes.df$TF <- map_symbol_to_ensg(
            motifs.probes.df$TF_symbol
        )

        motifs.probes.df <- motifs.probes.df %>% na.omit

    } else {

       hits <- findOverlaps(TF.peaks.gr, region.gr, ignore.strand = TRUE)
       motifs.probes.df <- data.frame(
           "regionID" = region.names[subjectHits(hits)],
           "TF_symbol" = TF.peaks.gr$id[queryHits(hits)]
       )
       motifs.probes.df$TF <- map_symbol_to_ensg(motifs.probes.df$TF_symbol)
    }
    return(motifs.probes.df %>% unique)
}
