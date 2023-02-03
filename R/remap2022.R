#' @title Access REMAP2022 non-redundant peaks
#' @description  Access REMAP2022 non-redundant peaks
#' @param cell_line filter peaks using cell line description field
#' @param species which species to access REMAP2022 non-redundant peaks.
#' Options are: homo_sapiens or mus_musculus
#' @export
readRemap2022 <- function(
    cell_line,
    species = c("homo_sapiens","mus_musculus")
  ){
  
  species <- match.arg(species)
  
  url <- dplyr::case_when(
    species == "homo_sapiens" ~  "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_nr_macs2_hg38_v1_0.bed.gz",
    species == "mus_musculus" ~  "https://remap.univ-amu.fr/storage/remap2022/mm39/MACS2/remap2022_nr_macs2_mm39_v1_0.bed.gz"
  )
  
  message("Downloading: ", basename(url))
  file <- paste0("readRemap2022/",basename(url))
  
  if(!file.exists(file)){
    dir.create(dirname(file),showWarnings = FALSE,recursive = TRUE)
    downloader::download(url, file)
  }
  
  remapCatalog <- bedToGranges(file)
  remapCatalog$id <- gsub("-|,","",remapCatalog$id)
  remapCatalog$cell_lines <- gsub("^[[:alnum:]]*:","",remapCatalog$id)
  remapCatalog$id <- gsub(":[[:alnum:]]*$","",remapCatalog$id)
  
  if(!missing(cell_line)){
    metadata <- openxlsx::read.xlsx(
      "https://remap.univ-amu.fr/storage/remap2022/biotypes/remap2022_hsap_biotypes.xlsx"
    )
    cell.lines <- metadata %>% 
      dplyr::filter(grepl(cell_line,.data$description,ignore.case = T)) %>% 
      dplyr::pull(biotype)
    remapCatalog <- remapCatalog[remapCatalog$cell_lines %in% gsub("-|,","",cell.lines),]
  }
  
  remapCatalog
}

bedToGranges <- function (path) {
  bedData <- bedImport(path)
  if (ncol(bedData) > 6) 
    bedData <- bedData[, -c(7:ncol(bedData))]
  if (ncol(bedData) < 3) 
    stop("File has less than 3 columns")
  if ("strand" %in% colnames(bedData)) 
    bedData$strand <- gsub(
      pattern = "[^+-]+", replacement = "*", 
      x = bedData$strand
    )
  if (ncol(bedData) == 3) {
    granges <- with(bedData, GenomicRanges::GRanges(
      chrom, 
      IRanges::IRanges(chromStart, chromEnd))
    )
  }
  else if (ncol(bedData) == 4) {
    granges = with(bedData, GenomicRanges::GRanges(
      chrom, 
      IRanges::IRanges(chromStart, chromEnd), id = name)
    )
  }
  else if (ncol(bedData) == 5) {
    granges <- with(bedData, GenomicRanges::GRanges(
      chrom, 
      IRanges::IRanges(chromStart, chromEnd), id = name, 
      score = score)
    )
  }
  else if (ncol(bedData) == 6) {
    granges <- with(bedData, GenomicRanges::GRanges(
      chrom, 
      IRanges::IRanges(chromStart, chromEnd), id = name, 
      score = score, strand = strand)
    )
  }
  else {
    stop("Error while constructing the GRanges object. \n             No number of columns have been matched.")
  }
  return(granges)
}

bedImport <- function (path) {
  regions <- as.data.frame(
    data.table::fread(
      path, 
      header = FALSE, 
      sep = "\t", 
      stringsAsFactors = FALSE, 
      quote = ""
    )
  )
  col.names <- c(
    "chrom", "chromStart", "chromEnd", "name", 
    "score", "strand", "thickStart", "thickEnd", "itemRgb", 
    "blockCount", "blockSizes", "blockStarts"
  )
  colnames(regions) <- col.names[1:ncol(regions)]
  return(regions)
}