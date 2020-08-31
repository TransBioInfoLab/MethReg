#' @title Access human TF from Lambert et al 2018
#' @description Access human TF from Lambert et al 2018 (PMID: 29425488)
#' @importFrom readr read_csv
#' @return A dataframe with Human TF
#' @examples
#'  human.tfs <- get_human_tfs()
#' @export
get_human_tfs <- function() {
    human.TF <- readr::read_csv(
        "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"
    )
    human.TF <- human.TF[human.TF$`Is TF?` == "Yes",]
    colnames(human.TF)[2:3] <- c("ensembl_gene_id","external_gene_name")
    human.TF$X1 <- NULL
    return(human.TF)
}
