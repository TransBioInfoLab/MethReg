#' Create Granges from name
#' @description Given a region name such as chr22:18267969-18268249, we will create a Granges
#' object
#' @importFrom tidyr separate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @examples
#' regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#' regions.gr <- make_granges_from_names(regions.names)
#' @noRd
make_granges_from_names <- function(names){
    names %>%
        data.frame %>%
        separate(col = ".",into = c("chr","start","end")) %>%
        makeGRangesFromDataFrame()
}
