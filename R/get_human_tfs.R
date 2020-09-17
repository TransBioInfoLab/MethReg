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



#' @title Access to TF target gene from DMTDB
#' @description Access DMTD for all cancer from
#' http://bio-bigdata.hrbmu.edu.cn/DMTDB/download.jsp
#' and filter TF and target gene for a given cancer.
#' @param cancer A TCGA cancer identifier
#' @noRd
get_DMTD_target_tf <- function(
    cancer = c(
        "BLCA",
        "BRCA",
        "CESC",
        "COAD",
        "ESCA",
        "HNSC",
        "KIRC",
        "KIRP",
        "LGG",
        "LIHC",
        "LUAD",
        "LUSC",
        "OV",
        "PAAD",
        "PCPG",
        "PRAD",
        "SARC",
        "SKCM",
        "STAD",
        "TGCT",
        "THCA",
        "UCEC"
    )
){

    cancer <- match.arg(cancer)
    database <- get_DMTD() %>% readr::read_tsv() %>% filter(.data$Cancer == cancer)
    return(database)
}


#' @title Download supplemental file from Lambert et al 2018
#' @description Download supplemental file from Lambert et al 2018 (PMID: 29425488)
#' @param bfc A BiocFileCache instance
#' @return A named vector with the url for the local file
#' @examples
#' if (interactive()) file.url <- get_lambert_tf()
#' @noRd
get_DMTD <- function(
    bfc = BiocFileCache::BiocFileCache(ask = FALSE)
) {

    check_package("BiocFileCache")
    url.root <- "http://bio-bigdata.hrbmu.edu.cn/DMTDB/download_loading.jsp"
    url.options <- "?path=download/DMTD_V2/all.txt&name=All_DMTD_V2.txt"
    url <- paste0(url.root, url.options)

    # check if url is being tracked
    res <- BiocFileCache::bfcquery(bfc, url)

    if (BiocFileCache::bfccount(res) == 0L) {
        # if it is not in cache, add
        ans <- BiocFileCache::bfcadd(bfc, rname = "DMTD_V2", fpath = url)
    } else {
        # if it is in cache, get path to load
        ans <- BiocFileCache::bfcrpath(bfc, rnames = "DMTD_V2")
    }
    ans
}




