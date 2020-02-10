#' @title Mapping regions to closes gene
#' @description Maps a given region to th closest gene
#' @param regions A dataframe with chrom, start, end or a GRanges
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param method How to map regions to genes: closest gene, gens overlapping a window
#' around the genomic region input.
#' @param window.width Number of base pairs to extend the region (+-window.width/2)
#' @importFrom coMethDMR AnnotateResults
#' @importFrom GenomicRanges distanceToNearest nearest ranges makeGRangesFromDataFrame values
#' @importFrom tidyr unite
#' @importFrom ELMER getTSS
#' @examples
#' library(GenomicRanges)
#' library(dplyr)
#' regions.gr <- data.frame(
#'  chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'   start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'   end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'                          stringsAsFactors = FALSE)  %>%
#'                          makeGRangesFromDataFrame
#'  getDNAm.target(regions.gr = regions.gr ,genome = "hg19",)
#'  getDNAm.target(regions.gr = regions.gr ,genome = "hg19",method = "window")
getDNAm.target <- function(
    regions.gr,
    genome = c("hg38","hg19"),
    method = c("closest.gene","window"),
    window.width = 500000
){

    method <- match.arg(method)
    genome <- match.arg(genome)

    if(!is(regions.gr,"GRanges")) stop("regions.gr must be a GRanges")


    if(method == "closest.gene"){
        tssAnnot <- ELMER::getTSS(genome = genome)

        neargenes <- tssAnnot[nearest(regions.gr,tssAnnot)] %>% as.data.frame()
        distance.region.tss <- values(distanceToNearest(regions.gr,tssAnnot))$distance

        neargenes <- cbind(neargenes[,c("seqnames","start","end","external_gene_name","ensembl_gene_id")],
                           distance.region.tss)

        colnames(neargenes)[1:3] <- c("gene_chrom","gene_start","gene_end")
        regionID <- paste0( seqnames(regions.gr) %>% as.character(),
                            ":",
                            start(regions.gr),
                            "-",
                            end(regions.gr))
        out <- cbind(regionID, neargenes) %>% tibble::as_tibble()
    } else {
        geneAnnot <- ELMER:::get.GRCh(genome = genome,as.granges = TRUE)
        regions.gr.extend <- regions.gr + (window.width/2)

        overlap <- findOverlaps(regions.gr.extend,geneAnnot)

        regionID <- paste0(
            seqnames(  regions.gr[queryHits(overlap)]) %>% as.character(),
            ":",
            start(  regions.gr[queryHits(overlap)]),
            "-",
            end(  regions.gr[queryHits(overlap)]))


        regionID.extended <- paste0(
            seqnames(  regions.gr.extend[queryHits(overlap)]) %>% as.character(),
            ":",
            start(  regions.gr.extend[queryHits(overlap)]),
            "-",
            end(  regions.gr.extend[queryHits(overlap)]))

        genes.overlapping <- tssAnnot[subjectHits(overlap)] %>% as.data.frame()
        colnames(genes.overlapping) <- paste0("gene_",colnames(genes.overlapping))

        out <- dplyr::bind_cols(data.frame("regionID" = regionID),
                                data.frame("regionID.extended" = regionID.extended,
                                           "window.extended.width" = window.width),
                                genes.overlapping)
    }
    return(out)
}


#' @title Evaluate correlation of region and gene using spearman test
#' @description Evaluate correlation of region and gene using spearman test
#' @param links.df A dataframe with the following collumns
#' region and gene
#' @param met DNA methylation matrix (rows are regions, columns samples). Samples should be in the
#' same order as gene expression.
#' @param exp Gene expression matrix (rows are genes, columns samples)
#' Samples should be in the same order as gene expression.
#' @param file.out A csv file name to save the results.
#' @param min.cor.pval Filter of significant correlations (default: 0.05)
#' @param min.cor.estimate Filter of significant correlations (default: not applied)
#' @importFrom plyr adply
#' @export
#' @examples
#' regions.gr <- data.frame(
#'  chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'   start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'   end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'                          stringsAsFactors = FALSE)  %>%
#'                          makeGRangesFromDataFrame
#' map <- getDNAm.target(regions.gr = regions.gr ,genome = "hg19")
#' map <- unite(map,col = "regionID",c("region_chrom", "region_start", "region_end" ))
#' links <- tibble::tibble(regionID = map$regionID, geneID = map$ensembl_gene_id)
#' met <- matrix(rep(0,length(links$regionID) * 4),
#'               nrow = length(links$regionID),
#'               dimnames = list(c(links$regionID),c(paste0("S",c(1:4)))))
#'
#' exp <- matrix(rep(0,length(links$geneID) * 4),
#'               nrow = length(links$geneID),
#'               dimnames = list(c(links$geneID),c(paste0("S",c(1:4)))))
#' corMetGene(links, met, exp)
corMetGene <- function(links,
                       met,
                       exp,
                       file.out,
                       min.cor.pval = 0.05,
                       min.cor.estimate = 0.0){

    if(is.null(exp)) stop("Please set exp matrix")
    if(is.null(met)) stop("Please set met matrix")
    if(ncol(met) != ncol(exp)) stop("exp and met does not have the same size")
    if(!all(c("geneID","regionID") %in% colnames(links))) stop("links object must have geneID and regionID columns")

    links <- links[links$geneID %in% rownames(exp),]
    links <- links[links$regionID %in% rownames(met),]
    if(nrow(links) == 0) stop("links not found in data. Please check rownames and links provided.")

    correlation.df <- plyr::adply(.data = links,
                                  .margins = 1,
                                  .fun = function(link){
                                      exp <- log2(exp[link$geneID,] + 1)
                                      met <- met[rownames(met) == link$regionID,]
                                      res <- cor.test(exp %>% as.numeric,
                                                      met %>% as.numeric,
                                                      method = "spearman",
                                                      exact = FALSE)
                                      return(tibble(met_exp_cor_pvalue = res$p.value,
                                                    met_exp_cor_estimate = res$estimate))
                                  },.progress = "time")


    correlation.df <- na.omit(correlation.df)
    correlation.df$met_exp_cor_fdr <- p.adjust(correlation.df$met_exp_cor_pvalue, method = "fdr")

    correlation.df <- correlation.df %>%
        dplyr::filter(met_exp_cor_fdr <= min.cor.pval & abs(met_exp_cor_estimate) >= min.cor.estimate)


    if(!missing(file.out)) readr::write_tsv(x = correlation.df, path = file.out)

    return(correlation.df)
}
