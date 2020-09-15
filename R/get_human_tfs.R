#' @title Access human TF from Lambert et al 2018
#' @description Access human TF from Lambert et al 2018 (PMID: 29425488)
#' @importFrom readr read_csv
#' @return A dataframe with Human TF
#' @examples
#'  human.tfs <- get_human_tfs()
#' @export
get_human_tfs <- function() {
    ans <- get_lambert_tf()
    human.TF <- suppressWarnings(
        readr::read_csv(ans, col_types = readr::cols())
    )
    human.TF <- human.TF %>% filter(.data$`Is TF?` == "Yes")
    colnames(human.TF)[2:3] <- c("ensembl_gene_id","external_gene_name")
    human.TF$X1 <- NULL
    return(human.TF)
}

#' @title Download supplemental file from Lambert et al 2018
#' @description Download supplemental file from Lambert et al 2018 (PMID: 29425488)
#' @param bfc A BiocFileCache instance
#' @return A named vector with the url for the local file
#' @examples
#' if (interactive()) file.url <- get_lambert_tf()
#' @noRd
get_lambert_tf <- function(
    bfc = BiocFileCache::BiocFileCache(ask = FALSE)
) {

    check_package("BiocFileCache")
    url <-  "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"

    # check if url is being tracked
    res <- BiocFileCache::bfcquery(bfc, url)

    if (BiocFileCache::bfccount(res) == 0L) {
        # if it is not in cache, add
        ans <- BiocFileCache::bfcadd(bfc, rname = "Lambert_TF", fpath = url)
    } else {
        # if it is in cache, get path to load
        ans <- BiocFileCache::bfcrpath(bfc, rnames = "Lambert_TF")
    }
    ans
}
