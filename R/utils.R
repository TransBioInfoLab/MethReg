#' Create Granges from name
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
#' @examples
#' regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#' regions.gr <- make_granges_from_names(regions.names)
#' make_names_from_granges(regions.gr)
#' @noRd
make_names_from_granges <- function(region){
    str_c(
        region %>% seqnames %>% as.character,":",
        region %>% start %>% as.character,"-",
        region %>% end %>% as.character)
}


#' Change probes names to region names
#' @description Given a DNA methylation matrix with probes as row names,
#' map probes to genomic regions
#' @param dnam A DNA methylation matrix
#' @param genome Human genome of reference hg38 or hg19
#' @param arrayType DNA methylation array type (450k or EPIC)
#' @param rm.masked.probes Remove masked probes ? Default: TRUE
#' @examples
#' data(dna.met.chr21)
#' dna.met.chr21.with.region.name <- map_probes_to_regions(dna.met.chr21)
#' @export
#' @importFrom sesameData sesameDataCacheAll sesameDataGet
map_probes_to_regions <- function(
    dnam,
    genome = c("hg38","hg19"),
    arrayType = c("450k","EPIC"),
    rm.masked.probes = TRUE
){
    genome <- match.arg(genome)
    arrayType <- match.arg(arrayType)

    probe.info <- get_met_probes_info(genome, arrayType)

    if(rm.masked.probes){
        # Remove probes that should be masked
        probe.info <- probe.info[!probe.info$MASK_general,]
        # Keep non-masked probes
        dnam <- dnam[rownames(dnam) %in% names(probe.info),]
    }

    rownames(dnam) <- make_names_from_granges(probe.info[rownames(dnam)])
    return(dnam)
}

get_met_probes_info <- function(
    genome = c("hg38","hg19"),
    arrayType = c("450k","EPIC")
){
    genome <- match.arg(genome)
    arrayType <- match.arg(arrayType)

    sesameDataCacheAll()
    sesameDataGet(
        str_c(ifelse(arrayType == "450k","HM450","EPIC"),".",
              genome,".manifest")
    )
}


#' @param genome Human genome of reference. Options: hg38, hg19.
#' @param ensembl.gene.id Gene ensembl ID. A character vectors
#' @description Given a GRanges returns region name such as chr22:18267969-18268249
#' @examples
#' data(gene.exp.chr21)
#' gene.symbols <- map_ensg_to_symbol(rownames(gene.exp.chr21))
#' @noRd
#' @importFrom biomaRt useEnsembl listDatasets getBM
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
#' @importFrom biomaRt useEnsembl listDatasets getBM
map_symbol_to_ensg <- function(
    gene.symbol,
    genome = "hg38"
){
    gene.location <- get_gene_information(genome)
    ensembl_gene_id <- gene.location[match(gene.symbol,gene.location$external_gene_name),]$ensembl_gene_id
    return(ensembl_gene_id)
}


get_gene_information_biomart <- function(genome = "hg38"){
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
                message(ifelse(is.null(mirror),
                               paste0("Accessing ",
                                      host, " to get gene information"),
                               paste0("Accessing ",
                                      host, " (mirror ", mirror, ")")))
                useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                           host = host, mirror = mirror)
            }, error = function(e) {
                message(e)
                return(NULL)
            })
            attributes <- c(
                "ensembl_gene_id",
                "external_gene_name",
                "chromosome_name",
                "strand",
                "end_position",
                "start_position"
            )
            db.datasets <- listDatasets(ensembl)
            description <- db.datasets[db.datasets$dataset == "hsapiens_gene_ensembl", ]$description
            message(paste0("Downloading genome information (try:", tries, ") Using: ", description))
            gene.location <- getBM(
                attributes = attributes,
                filters = "chromosome_name",
                values = c(1:22,"X","Y"),
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
get_gene_information <- function(genome = "hg38", as.granges = FALSE){
    if(genome == "hg19"){
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
                end.field = "end_position", keep.extra.columns = TRUE)
    }

    return(gene.location)
}


check_data <- function(dnam, exp, metadata){

    if(!is(dnam,"matrix")) {
        stop("DNA methylation should be a matrix object")
    }

    if(!is(exp,"matrix")) {
        stop("Gene expression data should be a matrix object")
    }

    if(ncol(dnam) != ncol(exp)){
        stop("DNA methylation and gene expression do not has the same number of samples")
    }

    if(!all(colnames(dnam) == colnames(exp))){
        stop("DNA methylation and gene expression do not has the column names")
    }

    if(!missing(metadata)){

        if(nrow(metadata) != ncol(exp)){
            stop("Metadata and data do not has the same number of samples")
        }

        if(all(rownames(metadata) != colnames(exp))){
            stop("Metadata and data do not has the same number of samples")
        }
    }

}

#' @title register cores
#' @param package Package name
#' @noRd
check_package <- function(package){
    if (!requireNamespace(package, quietly = TRUE)) {
        stop(package, " package is needed for this function to work. Please install it.",
             call. = FALSE)
    }
}

#' @title register cores
#' @param cores A interger which defines the number of cores to be used in parallel
#' @noRd
register_cores <- function(cores){

    check_package("parallel")
    check_package("doParallel")

    parallel <- FALSE
    if (cores > 1){
        if (cores > parallel::detectCores()) cores <- parallel::detectCores()
        doParallel::registerDoParallel(cores)
        parallel = TRUE
    }
    return(parallel)
}

#' @importFrom IRanges subsetByOverlaps
#' @export
get_non_promoter_regions <- function(regions.gr, genome){
    message("o Get promoter regions for ", genome)
    promoter.gr <- get_promoter_regions(genome) %>% reduce
    promoter.regions <- IRanges::subsetByOverlaps(regions.gr,promoter.gr)
    message("o Remove promoter regions")
    GenomicRanges::setdiff(regions.gr, promoter.regions)
}

#' @title Transform DNA methylation array to a summarized Experiment object
#' @param met DNA methylation matrix with beta-values or m-values as data,
#' row as regions or regions and column as samples
#' @export
#' @examples
#' dna.met.chr21 <- get(data("dna.met.chr21"))
#' se <- make_se_from_dnam_probes(dna.met.chr21)
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
make_se_from_dnam_probes <- function (
    met,
    genome = c("hg38","hg19"),
    arrayType = c("450k","EPIC")
) {
    genome <- match.arg(genome)
    arrayType <- match.arg(arrayType)

    check_package("SummarizedExperiment")
    check_package("S4Vectors")

    message("o Creating a SummarizedExperiment from DNA methylation input")

    # Get probes annotation
    message("oo Fetching probes metadata")
    annotation <- get_met_probes_info(genome = genome, arrayType =  arrayType)

    # Keep only annotation with information in the methylation array
    rowRanges <- annotation[names(annotation) %in% rownames(met),, drop = FALSE]
    if(length(rowRanges) == 0){
        message("We were not able to map the rownames to cpgs probes identifiers. Please, check your input.")
        return(NULL)
    }

    # remove masked probes
    message("oo Removing masked probes")
    rowRanges <- rowRanges[!rowRanges$MASK_general]

    # Prepare all data matrices
    colData <- S4Vectors::DataFrame(samples = colnames(met))
    met <- met[rownames(met) %in% names(rowRanges), , drop = FALSE]
    met <- met[names(rowRanges), , drop = FALSE]
    assay <- data.matrix(met)

    # Create SummarizedExperiment
    message("oo Preparing SummarizedExperiment object")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = assay,
        rowRanges = rowRanges,
        colData = colData
    )
    return(se)
}


#' @title Transform DNA methylation array to a summarized Experiment object
#' @param met DNA methylation matrix with beta-values or m-values as data,
#' row as cpgs and column as samples
#' @export
#' @examples
#' dna.met.chr21 <- get(data("dna.met.chr21"))
#' dna.met.chr21.regions <- map_probes_to_regions(dna.met.chr21)
#' dnam.se <- make_se_from_dnam_regions(dna.met.chr21.regions)
make_se_from_dnam_regions <- function (
    met
) {

    check_package("SummarizedExperiment")
    check_package("S4Vectors")

    message("o Creating a SummarizedExperiment from DNA methylation input")

    # Prepare all data matrices
    rowRanges <- met %>% rownames() %>% make_granges_from_names()
    colData <- S4Vectors::DataFrame(samples = colnames(met))
    assay <- data.matrix(met)

    # Create SummarizedExperiment
    message("oo Preparing SummarizedExperiment object")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = assay,
        rowRanges = rowRanges,
        colData = colData
    )
    return(se)
}

#' @title Transform gene expression matrix to a Summarized Experiment object
#' @param exp Gene expressio  matrix with gene expression counts,
#' row as ENSG gene IDS and column as samples
#' @export
#' @examples
#' gene.exp.chr21 <- get(data("gene.exp.chr21"))
#' gene.exp.chr21.se <- make_se_from_gene_matrix(gene.exp.chr21)
make_se_from_gene_matrix <- function (
    exp,
    genome = c("hg38","hg19")
) {
    # Data checking
    genome <- match.arg(genome)

    if (!all(grepl("ENSG", rownames(exp)))) {
        stop("Please the gene expression matrix should receive ENSEMBLE IDs (ENSG)")
    }

    message("o Creating a SummarizedExperiment from gene expression input")
    gene.info <- get_gene_information(genome = genome, as.granges = TRUE)
    rowRanges <- gene.info[match(exp %>% rownames(),gene.info$ensembl_gene_id),]
    colData <- S4Vectors::DataFrame(samples = colnames(exp))

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = exp %>% data.matrix(),
        rowRanges = rowRanges,
        colData = colData
    )
    return(se)
}
