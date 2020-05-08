#' @title Access TF targets from Cistrome cancer database
#' @description Access cistrome TF cancer database
#' which is a comprehensive resource for predicted transcription factor (TF)
#' targets and enhancer profiles in cancers.
#' For more information: http://cistrome.org/CistromeCancer/
#' @importFrom dplyr pull filter %>% tbl
#' @importFrom plyr adply
#' @param tcga.study A TCGA study from the database
#' (i.e. "ACC","BLCA", "BRCA_1", "BRCA_2", "CESC", "COAD_READ*")
#' Please check function for a complete list: get_cistrome_studies()
#' @param minCor Cistrome minimum correlation between TF and target gene
#' @importFrom rlang .data
#' @return A dataframe with TF, target and correlation
#' @examples
#' \dontrun{
#'  get_tf_targets_cistrome("ACC")
#'  get_tf_targets_cistrome("COAD_READ*")
#' }
#' @export
get_tf_targets_cistrome <- function(tcga.study, minCor = 0.2) {

    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop("reshape2 package is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    con <- get_cistrome_dbconn()
    check_cistrome_study(con, tcga.study)

    db <- dplyr::tbl(con, "targets_tf")
    study.tf <- db %>%
        dplyr::filter(.data$study == tcga.study) %>%
        dplyr::pull(.data$tf) %>% unique

    # find the TF within database table
    results <- plyr::adply(study.tf,
                           .margins = 1,
                           .fun = function(x){
                               # get tf targets
                               tars <- dplyr::tbl(con, "targets_tf")  %>%
                                   dplyr::filter(.data$tf == x & .data$study == tcga.study) %>%
                                   dplyr::pull(2)

                               dplyr::tbl(con, "cor_tf") %>%
                                   dplyr::select(c("tf", "feature", tcga.study)) %>%
                                   dplyr::filter(.data$tf == x & .data$feature %in% tars) %>%
                                   as.data.frame
                           },.id = NULL, .progress = "time",.parallel = FALSE)
    results <- results %>% reshape2::melt()
    colnames(results) <- c("TF", "Target", "study","cor")
    results$cor <- results$cor / 100
    results <- na.omit(results)
    message("Cistrome TF target # links: ", nrow(results))
    message("Applying minCor filter (minCor: ", minCor)
    results <- results %>% dplyr::filter(.data$cor > minCor)

    message("Cistrome TF target with minCor # links: ", nrow(results))
    return(results)
}

get_cistrome_dbconn <- function(){
    if (!requireNamespace("RSQLite", quietly = TRUE)) {
        stop("RSQLite package is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("DBI", quietly = TRUE)) {
        stop("DBI package is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("downloader", quietly = TRUE)) {
        stop("download package is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("R.utils", quietly = TRUE)) {
        stop("R.utils package is needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (!requireNamespace("gunzip", quietly = TRUE)) {
        stop("gunzip package is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    file <- "cRegulome.db"
    if(!file.exists(file)){
        downloader::download("https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/9537385/cRegulome.db.gz","cRegulome.db.gz")
        gunzip::gunzip(paste0(file,".gz"), remove = FALSE)
    }
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = file)
    return(con)
}

#' @importFrom knitr kable
check_cistrome_study <- function(conn, study){
    studies <- get_cistrome_studies()
    if(!study %in% studies){
        message("Available studies: ")
        print(knitr::kable(data.frame("Studies" = studies)))
        stop("Study not found. Please check list above.")
    }
}

get_cistrome_studies <- function(conn = NULL){


    if (!requireNamespace("DBI", quietly = TRUE)) {
        stop("DBI package is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if(!is.null(conn)) {
        studies <- DBI::dbListFields(conn, "cor_tf")
    } else {
        studies <- c("ACC",
                     "BLCA",
                     "BRCA_1",
                     "BRCA_2",
                     "CESC",
                     "COAD_READ*",
                     "GBM*",
                     "HNSC",
                     "KICH",
                     "KIRC",
                     "KIRP",
                     "LAML",
                     "LGG",
                     "LIHC",
                     "LUAD",
                     "LUSC",
                     "MESO",
                     "OV",
                     "PAAD",
                     "PCPG",
                     "PRAD",
                     "SARC",
                     "SKCM",
                     "STES*",
                     "TGCT",
                     "THCA",
                     "THYM",
                     "UCEC",
                     "UVM"
        )}
    studies
}

#' @title register cores
#' @param cores A interger which defines the number of cores to be used in parallel
#' @noRd
register_cores <- function(cores){

    if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("parallel package is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("doParallel package is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    parallel <- FALSE
    if (cores > 1){
        if (cores > parallel::detectCores()) cores <- parallel::detectCores()
        doParallel::registerDoParallel(cores)
        parallel = TRUE
    }
    return(parallel)
}
