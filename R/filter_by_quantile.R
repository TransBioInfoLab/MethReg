#' @title Select regions with variations in DNA methylation
#' levels above a threshold
#' @description
#' For each region, compares the mean DNA methylation (DNAm) levels
#' in samples with high DNAm  (Q4) vs. low DNAm (Q1) and requires
#' the difference to be above a threshold.
#' @param dnam DNA methylation matrix or SumarizedExperiment object
#' @param diff.mean.th
#' Threshold for difference in mean DNAm levels for samples in Q4 and Q1
#' @param cores
#' Number of CPU cores to be used in the analysis. Default: 1
#' @export
#' @examples
#' data("dna.met.chr21")
#' dna.met.chr21.filtered <- filter_dnam_by_quant_diff(
#'   dna.met.chr21
#' )
#' @return
#' A subset of the original matrix only with the
#' rows passing the filter threshold.
filter_dnam_by_quant_diff <- function(
    dnam,
    diff.mean.th = 0.2,
    cores = 1
){
    if(is(dnam,"SummarizedExperiment")){
        matrix <- assay(dnam)
    } else {
        matrix <- dnam
    }

    # remove all NA rows
    keep.rows <- which(rowSums(is.na(matrix)) != ncol(matrix))
    if(length(keep.rows) < nrow(matrix)){
        message("Removing rows with NAs for all samples")
        matrix <- matrix[keep.rows,]
    }

    IQR <- calculate_IQR(matrix)
    tab <- plyr::count(IQR$IQR > diff.mean.th)
    colnames(tab)[1] <- "Status"
    tab$Status[which(tab$Status == FALSE)] <- "Regions below threshold"
    tab$Status[which(tab$Status == TRUE)] <- "Regions above threshold"
    print(tab)

    diff.regions <- c(
        IQR %>%
            filter(IQR > diff.mean.th) %>%
            pull(.data$ID) %>%
            as.character()
    )
    dnam[diff.regions,,drop = FALSE]
}

#' @examples
#' library(dplyr)
#' 10 %>% rnorm %>%
#'   matrix(nrow = 1,dimnames = list(c("row1"), LETTERS[1:10])) %>%
#'   calculate_IQR
#' @noRd
calculate_IQR <- function(matrix){
    check_package("matrixStats")
    tibble::tibble(
        "ID" = rownames(matrix),
        "IQR" = matrixStats::rowIQRs(matrix, na.rm = TRUE)
    )
}

#' @examples
#' library(dplyr)
#' 10 %>% rnorm %>%
#'   matrix(nrow = 1,dimnames = list(c("row1"), LETTERS[1:10])) %>%
#'   calculate_mean_q4_minus_mean_q1
#' @noRd
calculate_mean_q4_minus_mean_q1 <- function(matrix, cores = 1){
    parallel <- register_cores(cores)

    plyr::adply(.data = matrix,.margins = 1,.fun = function(row){
        qs <- quantile(row, probs =  c(0.25,0.75), na.rm = TRUE)
        qs.mean <- tapply(row, findInterval(row, qs), mean, rm.na = TRUE)
        tibble::tibble("diff.mean" = max(qs.mean) - min(qs.mean))
    },.parallel = parallel, .progress = "time",.id = "ID")
}


#' @title Select genes with variations above a threshold
#' @description For each gene, compares the mean gene expression
#' levels in samples in high expression (Q4)
#' vs. samples with low gene expression (Q1),
#' and requires the fold change to be above a certain threshold.
#' @param exp Gene expression matrix or SumarizedExperiment object
#' @param fold.change
#' Threshold for fold change of mean gene
#' expression levels in samples with high
#' (Q4) and low (Q1) gene expression levels. Defaults to 1.5.
#' @param cores
#' Number of CPU cores to be used in the analysis. Default: 1
#' @export
#' @examples
#' data("gene.exp.chr21.log2")
#' gene.exp.chr21.log2.filtered <- filter_exp_by_quant_mean_FC(
#'   gene.exp.chr21.log2
#' )
#' @return
#' A subset of the original matrix only with the rows passing
#' the filter threshold.
filter_exp_by_quant_mean_FC <- function(
    exp,
    fold.change = 1.5,
    cores = 1
){

    if(is(exp,"SummarizedExperiment")){
        matrix <- assay(exp)
    } else {
        matrix <- exp
    }


    parallel <- register_cores(cores)
    diff.genes <- plyr::adply(matrix,.margins = 1,.fun = function(row){
        quant <-  quantile(row, na.rm = TRUE)
        quant.fold.change <- data.frame("q4_div_q1" = quant[4] / quant[2])

        low.cutoff <- quant[2]
        upper.cutoff <- quant[4]

        mean.q1 <- row[row <= low.cutoff] %>% mean(na.rm = TRUE)
        mean.q4 <- row[row >= upper.cutoff] %>% mean(na.rm = TRUE)
        data.frame(
            "diff_fold_change" = mean.q4 / mean.q1,
            stringsAsFactors = FALSE
        )
    }, .progress = "time",.parallel = parallel)


    tab <- plyr::count(diff.genes$diff_fold_change > fold.change)
    colnames(tab)[1] <- "Status"
    tab$Status[which(tab$Status == FALSE)] <- "Genes below threshold"
    tab$Status[which(tab$Status == TRUE)] <- "Genes above threshold"
    print(tab)

    diff.genes <- c(
        diff.genes %>%
            filter(.data$diff_fold_change > fold.change) %>%
            pull(.data$X1) %>%
            as.character()
    )
    exp[diff.genes,,drop = FALSE]
}


#' @title Remove genes with gene expression level equal to 0 in a
#' substantial percentage of the samples
#' @param exp Gene expression matrix or SumarizedExperiment object
#' @param max.samples.percentage Max percentage of samples with gene
#' expression as 0, for genes to be selected.
#' If max.samples.percentage 100, remove genes with 0 for 100% samples.
#' If max.samples.percentage 25, remove genes with 0 for more
#' than 25% of the samples.
#'
#' @noRd
#' @return A subset of the original matrix only with the rows
#' passing the filter threshold.
filter_genes_zero_expression <- function(
    exp,
    max.samples.percentage = 0.25
){

    if (is(exp,"SummarizedExperiment")) {
        matrix <- assay(exp)
    } else {
        matrix <- exp
    }

    na.or.zeros <- matrix == 0 | is.na(matrix)
    percent.na.or.zeros <- rowSums(na.or.zeros) / ncol(matrix)

    genes.keep <- (percent.na.or.zeros < max.samples.percentage) %>% which %>% names
    message("Removing ", nrow(matrix) - length(genes.keep), " out of ", nrow(matrix), " genes")
    exp[genes.keep,, drop = FALSE]
}

#' @title
#' Remove genes with gene expression level equal to 0 or NA in a all samples
#' @param exp Gene expression matrix or a Summarized Experiment object
#' @noRd
#' @examples
#' data("gene.exp.chr21.log2")
#' gene.exp.chr21.log2.filtered <- filter_genes_zero_expression_all_samples(
#'   gene.exp.chr21.log2
#' )
#' @return
#' A subset of the original matrix only with the rows
#' passing the filter threshold.
filter_genes_zero_expression_all_samples <- function(
    exp
){
    if(is(exp,"SummarizedExperiment")){
        exp <- assay(exp)
    }
    idx.all.zero <- rowSums(exp == 0, na.rm = TRUE) == ncol(exp)
    idx.all.na <- rowSums(is.na(exp)) == ncol(exp)

    # do not keep if it is all zero or all NA
    genes.keep <- rownames(exp)[!(idx.all.zero | idx.all.na)] %>% na.omit()
    if(length(genes.keep) < nrow(exp) & length(genes.keep) > 0){
        message(
            "Removing ", nrow(exp) - length(genes.keep),
            " out of ", nrow(exp), " genes"
        )
    }
    exp <- exp[genes.keep,,drop = FALSE]
}



#' @description  Give dnam and expression outputs the % of samples with
#' 0 in Q1 and Q4. This function is useful if the analysis is performed using
#' residuals. This function will add the following column to the input data
#' "% of 0 target genes (Q1 and Q4)".
#' @noRd
add_percent_zero_q1_q4 <- function(
    region.target,
    dnam,
    exp
){

    if (is(exp,"SummarizedExperiment")) {
        exp <- assay(exp)
    }

    if (is(dnam,"SummarizedExperiment")) {
        dnam <- assay(dnam)
    }

    aux <- plyr::adply(
        unique(region.target[,c("regionID","target")]),
        .margins = 1,
        .fun = function(row) {

            rna.target <- exp[rownames(exp) == row$target, , drop = FALSE]
            met <- dnam[rownames(dnam) == as.character(row$regionID), ]

            data <- data.frame(
                rna.target = rna.target %>% as.numeric,
                met = met %>% as.numeric
            )

            quant.met <-  quantile(data$met,na.rm = TRUE)
            low.cutoff <- quant.met[2]
            upper.cutoff <- quant.met[4]

            data.high.low <- data %>%
                filter(.data$met <= low.cutoff | .data$met >= upper.cutoff)
            data.high.low$metGrp <- ifelse(data.high.low$met <= low.cutoff, 0, 1)

            pct.zeros.in.quant.samples <- sum(
                data.high.low$rna.target == 0,
                na.rm = TRUE) / nrow(data.high.low)

            data.frame("% of 0 target genes (Q1 and Q4)" = paste0(
                round(pct.zeros.in.quant.samples * 100,digits = 2),
                " %")
            )
        }
    )
    region.target$`% of 0 target genes (Q1 and Q4)` <- aux$X..of.0.target.genes..Q1.and.Q4.[
        match(paste0(region.target$regionID,region.target$target),paste0(aux$regionID,aux$target))
        ]
    return(region.target)
}
