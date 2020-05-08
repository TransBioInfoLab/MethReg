#' @title Mapping regions to gene
#' @description To map a region to genes there are two options: 1) closest gene
#' 2) map to all genes within a window around the region (default window.width = 500kbp
#' (i.e. +/- 250kbp from start or end of the region)).
#' @param regions.gr A Genomic Ranges objec GRanges
#' @param genome Human genome of reference "hg38" or "hg19"
#' @param method How regions are mapped to genes: closest gene ("closest.gene); or
#' genes within a window around the region ("window").
#' @param window.width When \code{method = "window"}, number of base pairs to extend the region (+- window.width/2).
#' Default is 500kbp (+- 250kbp, i.e. 250k bp from start or end of the region)
#' @importFrom GenomicRanges distanceToNearest nearest ranges makeGRangesFromDataFrame values seqnames distance
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom tidyr unite
#' @importFrom ELMER getTSS
#' @examples
#' library(GenomicRanges)
#' library(dplyr)
#'
#' # Create example region
#' regions.gr <- data.frame(
#'               chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'               start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'               end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'               stringsAsFactors = FALSE)  %>%
#'      makeGRangesFromDataFrame
#'
#'  # map to closest gene
#'  region.closest <- get_region_target_gene(regions.gr = regions.gr, genome = "hg19", method = "closest.gene")
#'
#'  # map to all gene within region +- 250kbp
#'  region.window <- get_region_target_gene(regions.gr = regions.gr, genome = "hg19", method = "window")
#' @export
get_region_target_gene <- function(
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
        regionID <- paste0(
            regions.gr %>% seqnames %>% as.character(),
            ":",
            regions.gr %>% start,
            "-",
            regions.gr %>% end)
        out <- cbind(regionID, neargenes) %>% tibble::as_tibble()
    } else {
        geneAnnot <- get_gene_information(genome = genome,as.granges = TRUE)
        geneAnnot$entrezgene <- NULL
        geneAnnot <- unique(geneAnnot)
        regions.gr.extend <- regions.gr + (window.width/2)

        overlap <- findOverlaps(regions.gr.extend,geneAnnot)

        regionID <- paste0(
            regions.gr[queryHits(overlap)] %>% seqnames %>% as.character(),
            ":",
            regions.gr[queryHits(overlap)] %>% start,
            "-",
            regions.gr[queryHits(overlap)] %>% end)


        regionID.extended <- paste0(
            regions.gr.extend[queryHits(overlap)] %>% seqnames %>% as.character(),
            ":",
            regions.gr.extend[queryHits(overlap)] %>% start,
            "-",
            regions.gr.extend[queryHits(overlap)] %>% end
        )

        genes.overlapping <- geneAnnot[subjectHits(overlap)] %>% as.data.frame()
        colnames(genes.overlapping) <- paste0("gene_",colnames(genes.overlapping))

        out <- dplyr::bind_cols(data.frame("regionID" = regionID),
                                data.frame("regionID.extended" = regionID.extended,
                                           "window.extended.width" = window.width,
                                           "Distance region-gene" = distance(regions.gr[queryHits(overlap)],
                                                                             geneAnnot[subjectHits(overlap)])),
                                genes.overlapping)  %>% tibble::as_tibble()
    }
    return(out)
}


#' @title Evaluate correlation of DNA methylation region and target gene expressions
#' @description Evaluate correlation of the DNA methylation region and target gene expression
#' using spearman correlation test
#' @param links A dataframe with the following columns: regionID (DNA methylation) and  geneID (target gene)
#' @param met DNA methylation matrix (rows are regions and columns are samples). Samples should be in the
#' same order as gene expression.
#' @param exp Gene expression matrix (rows are genes, columns are samples)
#' Samples should be in the same order as the DNA methylation matrix.
#' @param min.cor.pval Filter of significant correlations (default: 0.05)
#' @param min.cor.estimate Filter of significant correlations (default: not applied)
#' @param file.out If provided, name of a csv file which will be used to save the results.
#' @importFrom plyr adply
#' @importFrom tibble tibble
#' @importFrom stats p.adjust cor.test
#' @export
#' @examples
#' library(dplyr)
#' # Create example region
#' regions.gr <- data.frame(
#'                 chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'                 start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'                 end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'                 stringsAsFactors = FALSE)  %>%
#'               GenomicRanges::makeGRangesFromDataFrame()
#'
#' # Map example region to closest gene
#' map <- get_region_target_gene(regions.gr = regions.gr, genome = "hg19", method = "closest.gene")
#' map <- tidyr::unite(map,col = "regionID",c("gene_chrom", "gene_start", "gene_end"))
#' links <- tibble::tibble(regionID = map$regionID, geneID = map$ensembl_gene_id)
#'
#' # Create data example
#' met <- matrix(runif(length(links$regionID) * 4, 0, 1),
#'               nrow = length(links$regionID),
#'               dimnames = list(c(links$regionID),c(paste0("S",c(1:4)))))
#'
#' exp <- matrix(runif(length(links$regionID) * 4, 0, 10),
#'               nrow = length(links$geneID),
#'               dimnames = list(c(links$geneID),c(paste0("S",c(1:4)))))
#'
#' # Samples in met and exp datasets should be in the same order.
#' identical (colnames (met), colnames(exp))
#'
#' # Correalted DNAm and gene expression, display only significant associations
#' cor_region_dnam_target_gene(links = links, met = met, exp = exp)
#'
#' # display all associations
#' cor_region_dnam_target_gene(links = links, met = met, exp = exp, min.cor.pval = 1)
cor_region_dnam_target_gene <- function(
    links,
    met,
    exp,
    min.cor.pval = 0.05,
    min.cor.estimate = 0.0,
    file.out
){

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
                                      return(tibble("met_exp_cor_pvalue" = res$p.value,
                                                    "met_exp_cor_estimate" = res$estimate))
                                  },.progress = "time")


    correlation.df <- na.omit(correlation.df)
    correlation.df$met_exp_cor_fdr <- p.adjust(correlation.df$met_exp_cor_pvalue, method = "fdr")

    correlation.df <- correlation.df %>%
        dplyr::filter(met_exp_cor_fdr <= min.cor.pval & abs(met_exp_cor_estimate) >= min.cor.estimate)


    if(!missing(file.out)) readr::write_tsv(x = correlation.df, path = file.out)

    return(correlation.df)
}
