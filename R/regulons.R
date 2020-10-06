#' @title Calculate enrichment scores for each TF across all samples using
#' \code{\link[dorothea]{dorothea}} and \code{\link[viper]{viper}}.
#' @param exp Gene expression matrix with gene expression counts,
#' row as ENSG gene IDS and column as samples
#' @param min.confidence Minimun confidence score  ("A", "B","C","D", "E")
#' classifying regulons based on their quality from Human DoRothEA database.
#' @examples
#' regulons <- get_regulon_dorothea(min.confidence = "E")
#' @return A dataframe with tf, confidence and target gene from dorothea package.
#' @noRd
get_regulon_dorothea <- function(
    min.confidence = c("A", "B","C","D", "E")
){
    check_package("dorothea")
    min.confidence <- match.arg(min.confidence)

    dorothea_hs <- get(utils::data(dorothea_hs, package = "dorothea"))
    confidence.set <- LETTERS[1:5][1:which(LETTERS[1:5] == min.confidence)]
    regulons = dorothea_hs %>%
        filter(.data$confidence %in% confidence.set)

    regulons$tf_ensg <- map_symbol_to_ensg(regulons$tf)
    regulons$target_ensg <- map_symbol_to_ensg(regulons$target)

    return(regulons)
}

#' @title Calculate enrichment scores for each TF across all samples using
#' dorothea and viper.
#' @param exp Gene expression matrix with gene expression counts,
#' row as ENSG gene IDS and column as samples
#' @param min.confidence Minimun confidence score ("A", "B","C","D", "E")
#' classifying regulons based on their quality from Human DoRothEA database.
#' The default minimun confidence score is "B"
#' @param regulons DoRothEA regulons in table format. Same as \link[dorothea]{run_viper}.
#' If not specified Bioconductor (human) dorothea regulons besed on GTEx will be.
#' used \link[dorothea]{dorothea_hs}.
#' @examples
#' gene.exp.chr21.log2 <- get(data("gene.exp.chr21.log2"))
#' tf_es <- get_tf_ES(gene.exp.chr21.log2)
#' @return A matrix of normalized enrichment scores for each TF across all samples
#' @export
get_tf_ES <- function(
    exp,
    min.confidence = "B",
    regulons
){
    check_package("dorothea")
    check_package("viper")

    min.confidence <- match.arg(
        arg = min.confidence,
        choices =   c("A", "B", "C", "D", "E")
    )

    if(missing(regulons)){
        regulons <- get_regulon_dorothea(min.confidence = min.confidence)
    } else {
        cols <- c("tf", "target", "mor")
        if(!all(cols %in% colnames(regulons))){
            stop("regulons must have columns tf, target, and mor")
        }
    }

    if(all(grepl("ENSG",rownames(exp)))){
        rownames(exp) <- map_ensg_to_symbol(rownames(exp))
    }

    tf_activities <- tryCatch({
        tf_activities <- dorothea::run_viper(
            input = exp,
            regulons = regulons,
            options =  list(
                method = "scale",
                minsize = 4,
                eset.filter = FALSE,
                cores = 1,
                verbose = FALSE
            )
        )
        rownames(tf_activities) <- map_symbol_to_ensg(rownames(tf_activities))
        tf_activities <- tf_activities[!is.na(rownames(tf_activities)),]

        tf_activities
    }, error = function(e) {
        message(e)
        return(NULL)
    })
    tf_activities

}
