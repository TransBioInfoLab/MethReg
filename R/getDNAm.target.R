#' @title Mapping regions to closes gene
#' @description Maps a given region to th closest gene
#' @param regions A dataframe with chrom, start, end or a GRanges
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param method For the moment mapping region to closest gene
#' @importFrom coMethDMR AnnotateResults
#' @importFrom GenomicRanges distanceToNearest nearest ranges makeGRangesFromDataFrame values
#' @importFrom tidyr unite
#' @importFrom ELMER getTSS
#' @examples
#' regions.gr <- data.frame(
#'  chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'   start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'   end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'                          stringsAsFactors = FALSE)  %>%
#'                          makeGRangesFromDataFrame
#'  getDNAm.target(regions.gr = regions.gr ,genome = "hg19")
getDNAm.target <- function(
    regions.gr,
    genome = c("hg38","hg19"),
    method = c("closest.gene")
){

    method <- match.arg(method)
    genome <- match.arg(genome)

    if(!is(regions.gr,"GRanges")) stop("regions.gr must be a GRanges")

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
    out <- cbind(regionID, neargenes)
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
corMetGene <- function(links, met, exp, file.out){

    if(is.null(exp)) stop("Please set exp matrix")
    if(is.null(met)) stop("Please set met matrix")
    if(ncol(met) != ncol(exp)) stop("exp and met does not have the same size")

    correlation.df <- plyr::adply(.data = links,
                                  .margins = 1,
                                  .fun = function(link){

                                      idx <- which(rownames(exp) %in% link$geneID)

                                      if(length(idx) < 1) {
                                          return(tibble(p.value = NA, estimate_rho = NA))
                                      }

                                      exp <- log2(exp[idx,] + 1)
                                      met <- met[link$regionID,]
                                      res <- cor.test(exp,
                                                      met,
                                                      method = "spearman",
                                                      exact = FALSE)
                                      return(tibble(p.value = res$p.value,
                                                    estimate_rho = res$estimate))
                                  },.progress = "time")


    correlation.df <- na.omit(correlation.df)
    correlation.df$fdr <- p.adjust(correlation.df$p.value,method = "fdr")
    if(!missing(file.out)) readr::write_tsv(x = correlation.df, path = file.out)

    return(correlation.df)
}
