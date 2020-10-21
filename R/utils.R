#' Create a Granges object from a genmic region string
#' @description Given a region name such as chr22:18267969-18268249, we will create a Granges
#' object
#' @importFrom tidyr separate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @param names A region name as "chr22:18267969-18268249" or a vector of region names.
#' @examples
#' regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#' regions.gr <- make_granges_from_names(regions.names)
#' @export
make_granges_from_names <- function(names){
    names %>%
        data.frame %>%
        separate(col = ".",into = c("chr","start","end")) %>%
        makeGRangesFromDataFrame()
}


#' Create region name from Granges
#' @description Given a GRanges returns region name such as chr22:18267969-18268249
#' @importFrom stringr str_c
#' @importFrom dplyr %>%
#' @importFrom GenomicRanges start end seqnames
#' @param region A GenomicRanges object
#' @examples
#' regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#' regions.gr <- make_granges_from_names(regions.names)
#' make_names_from_granges(regions.gr)
#' @export
make_names_from_granges <- function(region){
    str_c(
        region %>% seqnames %>% as.character,":",
        region %>% start %>% as.character,"-",
        region %>% end %>% as.character)
}


#' Change row names of a matrix from DNAm probes (CpGs) to genomic region.
#' @description Given a DNA methylation matrix with probes as row names,
#' map probes to genomic regions
#' @param dnam A DNA methylation matrix
#' @param genome Human genome of reference hg38 or hg19
#' @param arrayType DNA methylation array type (450k or EPIC)
#' @param rm.masked.probes Remove masked probes ? Default: TRUE
#' @examples
#' data(dna.met.chr21)
#' dna.met.chr21.with.region.name <- map_probes_to_regions(dna.met.chr21)
#' @importFrom sesameData sesameDataCacheAll sesameDataGet
#' @noRd
map_probes_to_regions <- function(
    dnam,
    genome = c("hg38","hg19"),
    arrayType = c("450k","EPIC"),
    rm.masked.probes = TRUE
){
    genome <- match.arg(genome)
    arrayType <- match.arg(arrayType)

    probe.info <- get_met_probes_info(genome, arrayType)

    if (rm.masked.probes) {
        # Remove probes that should be masked
        probe.info <- probe.info[!probe.info$MASK_general,]
        # Keep non-masked probes
        dnam <- dnam[rownames(dnam) %in% names(probe.info),]
    }

    rownames(dnam) <- make_names_from_granges(probe.info[rownames(dnam)])
    return(dnam)
}

#' Get HM450/EPIC manifest files from Sesame package
#' @description
#' Returns a data frame with HM450/EPIC manifest information
#' files from Sesame package
#' @examples
#' regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#' regions.gr <- make_granges_from_names(regions.names)
#' make_names_from_granges(regions.gr)
#' @export
get_met_probes_info <- function(
    genome = c("hg38","hg19"),
    arrayType = c("450k","EPIC")
){
    genome <- match.arg(genome)
    arrayType <- match.arg(arrayType)

    check_package("sesameData")
    check_package("sesame")

    sesameDataCacheAll()
    sesameDataGet(
        str_c(
            ifelse(arrayType == "450k","HM450","EPIC"),
            ".",
            genome,
            ".manifest"
        )
    )
}


#' @param genome Human genome of reference. Options: hg38, hg19.
#' @param ensembl.gene.id Gene ensembl ID. A character vectors
#' @description Given a GRanges returns region name such as chr22:18267969-18268249
#' @examples
#' data(gene.exp.chr21.log2)
#' gene.symbols <- map_ensg_to_symbol(rownames(gene.exp.chr21.log2))
#' @noRd
map_ensg_to_symbol <- function(
    ensembl.gene.id,
    genome = "hg38"
){
    gene.location <- get_gene_information(genome)
    symbols <- gene.location[match(ensembl.gene.id,gene.location$ensembl_gene_id),]$external_gene_name
    return(symbols)
}

#' @param genome Human genome of reference. Options: hg38, hg19.
#' @param ensembl.gene.id Gene ensembl ID. A character vectors
#' @description Given a GRanges returns region name such as chr22:18267969-18268249
#' @examples
#' gene.symbols <- map_symbol_to_ensg("TP63"s)
#' @noRd
map_symbol_to_ensg <- function(
    gene.symbol,
    genome = "hg38"
){
    gene.location <- get_gene_information(genome)
    ensembl_gene_id <- gene.location[match(gene.symbol,gene.location$external_gene_name),]$ensembl_gene_id
    return(ensembl_gene_id)
}


#' @title Get human genome information from biomaRt
#' @param genome Human genome of reference. Options: hg38, hg19.
#' @param TSS add TSS information
#' @noRd
get_gene_information_biomart <- function(
    genome = c("hg38","hg19"),
    TSS = FALSE
){
    check_package("biomaRt")
    genome <- match.arg(genome)
    tries <- 0L
    msg <- character()
    while (tries < 3L) {
        gene.location <- tryCatch({
            host <- ifelse(
                genome == "hg19",
                "grch37.ensembl.org",
                "www.ensembl.org"
            )
            mirror <- list(NULL, "useast", "uswest", "asia")[[tries + 1]]
            ensembl <- tryCatch({
                message(
                    ifelse(
                        is.null(mirror),
                        paste0("Accessing ", host, " to get gene information"),
                        paste0("Accessing ", host, " (mirror ", mirror, ")")
                    )
                )
                biomaRt::useEnsembl(
                    "ensembl",
                    dataset = "hsapiens_gene_ensembl",
                    host = host,
                    mirror = mirror
                )
            }, error = function(e) {
                message(e)
                return(NULL)
            })

            # Column values we will recover from the database
            attributes <- c(
                "ensembl_gene_id",
                "external_gene_name",
                "chromosome_name",
                "strand",
                "end_position",
                "start_position",
                "gene_biotype"
            )

            if (TSS)  attributes <- c(attributes, "transcription_start_site")

            db.datasets <- biomaRt::listDatasets(ensembl)
            description <- db.datasets[db.datasets$dataset == "hsapiens_gene_ensembl", ]$description
            message(paste0("Downloading genome information (try:", tries, ") Using: ", description))
            gene.location <- biomaRt::getBM(
                attributes = attributes,
                filters = "chromosome_name",
                values = c(seq_len(22),"X","Y"),
                mart = ensembl
            )
            gene.location
        }, error = function(e) {
            msg <<- conditionMessage(e)
            tries <<- tries + 1L
            NULL
        })
        if (!is.null(gene.location)) break
        if (tries == 3L) stop("failed to get URL after 3 tries:", "\n  error: ", msg)
    }
}

#' @noRd
get_gene_information <- function(
    genome = "hg38",
    as.granges = FALSE
){

    if (genome == "hg19") {
        gene.location <- gene.location.hg19
    } else {
        gene.location <- gene.location.hg38
    }

    if (as.granges) {
        gene.location$strand[gene.location$strand == 1] <- "+"
        gene.location$strand[gene.location$strand == -1] <- "-"
        gene.location$chromosome_name <- paste0("chr", gene.location$chromosome_name)
        gene.location <-  gene.location %>%
            makeGRangesFromDataFrame(
                seqnames.field = "chromosome_name",
                start.field = "start_position",
                end.field = "end_position", keep.extra.columns = TRUE
            )
    }

    return(gene.location)
}


check_data <- function(dnam, exp, metadata){

    if (!is(dnam,"matrix")) {
        stop("DNA methylation should be a matrix object")
    }

    if (!is(exp,"matrix")) {
        stop("Gene expression data should be a matrix object")
    }

    if (ncol(dnam) != ncol(exp)) {
        stop("DNA methylation and gene expression don't have the same number of samples")
    }

    if (!all(colnames(dnam) == colnames(exp))) {
        stop("DNA methylation and gene expression don't have the column names")
    }

    if (!missing(metadata)) {

        if (nrow(metadata) != ncol(exp)) {
            stop("Metadata and data don't have the same number of samples")
        }

        if (all(rownames(metadata) != colnames(exp))) {
            stop("Metadata rownames and data columns don't have the same names of samples")
        }
    }

}

#' @title check if package is avaliable
#' @param package Package name
#' @noRd
check_package <- function(package){
    suppressMessages({
        if (!requireNamespace(package, quietly = TRUE)) {
            stop(
                package,
                " package is needed for this function to work. Please install it.",
                call. = FALSE
            )
        }
    })
}

#' @title register cores
#' @param cores A interger which defines the number of cores to be used in parallel
#' @noRd
register_cores <- function(cores){

    check_package("parallel")
    check_package("doParallel")

    parallel <- FALSE
    if (cores > 1) {
        if (cores > parallel::detectCores()) cores <- parallel::detectCores()
        doParallel::registerDoParallel(cores)
        parallel = TRUE
    }
    return(parallel)
}

#' @title Subset regions to those not overlapping promoter regions
#' @description Subset a granges object to those not overlapping promoter regions
#' (default +- 2kb away from TSS)
#' @importFrom IRanges subsetByOverlaps
#' @noRd
subset_by_non_promoter_regions <- function(
    regions.gr,
    genome,
    upstream = 2000,
    downstream = 2000
){
    message("o Get promoter regions for ", genome)
    promoter.gr <- get_promoter_regions(
        genome = genome,
        upstream = upstream,
        downstream = downstream
    )
    promoter.regions <- IRanges::subsetByOverlaps(
        regions.gr,
        promoter.gr,
        ignore.strand = TRUE
    )

    message("o Remove promoter regions")
    GenomicRanges::setdiff(regions.gr, promoter.regions)
}

#' @title Subset regions to those  overlapping promoter regions
#' @description Subset a granges object to those overlapping promoter regions
#' (default +- 2kb away from TSS)
#' @importFrom IRanges subsetByOverlaps
#' @noRd
subset_by_promoter_regions <- function(
    regions.gr,
    genome,
    upstream = 2000,
    downstream = 2000
){
    message("o Get promoter regions for ", genome)
    promoter.gr <- get_promoter_regions(
        genome = genome,
        upstream = upstream,
        downstream = downstream
    )
    promoter.regions <- IRanges::subsetByOverlaps(regions.gr, promoter.gr)
    return(promoter.regions)
}

#' @title Get promoter genes using biomart
#' @description Subset a granges object to those overlapping promoter regions
#' (default +- 2kb away from TSS)
#' @noRd
#' @importFrom GenomicRanges promoters strand strand<-
get_promoter_regions <- function(
    genome,
    upstream = 2000,
    downstream = 2000
){

    genes <- get_gene_information(genome = genome, as.granges = TRUE)
    promoters.gr <- promoters(genes, upstream = upstream, downstream = downstream)
    strand(promoters.gr) <- "*"
    return(promoters.gr %>% unique)
}


#' @title Transform DNA methylation array into a summarized Experiment object
#' @param dnam DNA methylation matrix with beta-values or m-values as data,
#' row as cpgs "cg07946458" or regions ("chr1:232:245") and column as samples
#' @param genome Human genome of reference: hg38 or hg19
#' @param arrayType DNA methylation array type (450k or EPIC)
#' @param betaToM indicates if converting methylation beta values to mvalues
#' @param verbose A logical argument indicating if
#' messages output should be provided.
#' @export
#' @examples
#' library(dplyr)
#' dnam <- runif(20, min = 0,max = 1) %>% sort %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(dnam) <- c("chr3:203727581-203728580")
#' colnames(dnam) <- paste0("Samples",1:20)
#'  se <- make_dnam_se(dnam)
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @return A summarized Experiment object with DNA methylation probes mapped to
#' genomic regions
make_dnam_se <- function(
    dnam,
    genome = c("hg38","hg19"),
    arrayType = c("450k","EPIC"),
    betaToM = FALSE,
    verbose = FALSE
) {
    genome <- match.arg(genome)
    arrayType <- match.arg(arrayType)

    check_package("SummarizedExperiment")
    check_package("S4Vectors")

    if(verbose)  message("o Creating a SummarizedExperiment from DNA methylation input")

    # Get probes annotation
    # Prepare all data matrices
    if (any(grepl("chr", dnam %>% rownames()))) {
        rowRanges <- dnam %>% rownames() %>% make_granges_from_names()
        colData <- S4Vectors::DataFrame(samples = colnames(dnam))
    } else {
        if(verbose)  message("oo Fetching probes metadata")
        annotation <- get_met_probes_info(genome = genome, arrayType = arrayType)
        strand(annotation) <- "*" # don't add strand info for CpGs

        # Keep only annotation with information in the methylation array
        rowRanges <- annotation[names(annotation) %in% rownames(dnam),, drop = FALSE]
        if (length(rowRanges) == 0) {
            message("We were not able to map the rownames to cpgs probes identifiers. Please, check your input.")
            return(NULL)
        }

        # remove masked probes
        if(verbose)  message("oo Removing masked probes")
        rowRanges <- rowRanges[!rowRanges$MASK_general]
        # rowRanges <- rowRanges[grep("cg",names(rowRanges))] # remove rs probes
        dnam <- dnam[rownames(dnam) %in% names(rowRanges), , drop = FALSE]
        dnam <- dnam[names(rowRanges), , drop = FALSE]
    }
    # Prepare all data matrices
    colData <- S4Vectors::DataFrame(samples = colnames(dnam))
    assay <- data.matrix(dnam)

    if (betaToM) {
        ### Compute M values
        assay <- log2(assay / (1 - assay))
    }

    rowRanges$probeID <- names(rowRanges)
    regions.name <-  make_names_from_granges(rowRanges)
    names(rowRanges) <- regions.name
    rownames(dnam) <- regions.name

    # Create SummarizedExperiment
    if(verbose)  message("oo Preparing SummarizedExperiment object")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = assay,
        rowRanges = rowRanges,
        colData = colData,
        metadata = list(
            "genome" = genome,
            "arrayType" = arrayType
        )
    )
    return(se)
}

#' @title Transform gene expression matrix into a Summarized Experiment object
#' @param exp Gene expression matrix with gene expression counts,
#' row as ENSG gene IDS and column as samples
#' @param genome Human genome of reference: hg38 or hg19
#' @param verbose A logical argument indicating if
#' messages output should be provided.
#' @export
#' @examples
#' gene.exp.chr21.log2 <- get(data("gene.exp.chr21.log2"))
#' gene.exp.chr21.log2.se <- make_exp_se(gene.exp.chr21.log2)
#' @return A summarized Experiment object
make_exp_se <- function(
    exp,
    genome = c("hg38","hg19"),
    verbose = FALSE
) {
    # Data checking
    genome <- match.arg(genome)

    if (!all(grepl("ENSG", rownames(exp)))) {
        stop("Please the gene expression matrix should receive ENSEMBLE IDs (ENSG)")
    }

    if(verbose)  message("o Creating a SummarizedExperiment from gene expression input")
    gene.info <- get_gene_information(genome = genome, as.granges = TRUE)
    rowRanges <- gene.info[match(exp %>% rownames(),gene.info$ensembl_gene_id),]
    names(rowRanges) <- rowRanges$ensembl_gene_id
    colData <- S4Vectors::DataFrame(samples = colnames(exp))

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = exp %>% data.matrix(),
        rowRanges = rowRanges,
        colData = colData,
        metadata = list("genome" = genome)
    )
    return(se)
}
