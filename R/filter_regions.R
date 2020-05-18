#' @title Filter regions with high DNA methylation variance
#' @description For each region, compares the mean DNA methylation of samples in high (Q4) vs. low (Q1)
#' and requires difference to be above a certain threshold.
#' @param dnam DNA methylation matrix
#' @param diff.mean.th Diff. mean Threshold between samples in Q4 and Q1
#' @export
#' @examples
#' data("dna.met.chr21")
#' dna.met.chr21.filtered <- filter_regions_by_mean_quantile_difference(dna.met.chr21)
filter_regions_by_mean_quantile_difference <- function(dnam, diff.mean.th = 0.2){

    diff_mean <- plyr::adply(dnam,.margins = 1,.fun = function(row){
        quant.met <-  quantile(row, na.rm = TRUE)
        quant.diff <- data.frame("met.q4_minus_q1" = quant.met[4] - quant.met[2])

        low.cutoff <- quant.met[2]
        upper.cutoff <- quant.met[4]

        data.low <- row %>% .[. <= low.cutoff] %>% mean(na.rm = TRUE)
        data.high <- row %>% .[. >= upper.cutoff] %>% mean(na.rm = TRUE)
        data.frame("diff_mean" = data.high - data.low, stringsAsFactors = FALSE)
    }, .progress = "time")


    tab <- plyr::count(diff_mean$diff_mean > diff.mean.th)
    colnames(tab)[1] <- "Status"
    tab$Status[which(tab$Status == FALSE)] <- "Regions below threshold"
    tab$Status[which(tab$Status == TRUE)] <- "Regions above threshold"
    print(tab)

    diff.regions <- c(diff_mean %>% filter(diff_mean > 0.2) %>% pull(X1) %>% as.character())
    dnam[diff.regions,]
}
