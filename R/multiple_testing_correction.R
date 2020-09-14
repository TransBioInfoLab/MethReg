#' @examples
#' \dontrun{
#' data("dna.met.chr21")
#' dna.met.chr21 <- make_se_from_dnam_probes(dna.met.chr21)
#' data("gene.exp.chr21.log2")
#' triplet <- data.frame(
#'     "regionID" = rownames(dna.met.chr21)[1:10],
#'     "TF" = rownames(gene.exp.chr21.log2)[11:20],
#'     "target" = rownames(gene.exp.chr21.log2)[1:10]
#' )
#' results <- interaction_model(triplet, dna.met.chr21, gene.exp.chr21.log2)
#' calculate_stage_wise_adjustment(results)
#' }
#' @noRd
calculate_stage_wise_adjustment <- function(results){

    check_package("stageR")

    if(!"tripletID" %in% colnames(results)){
        results$tripletID <- create_triplet_ID(results)
    }

    interactiol.col <- grep("quant_pval_metGrp:",colnames(results),value = TRUE)
    results <- stage_wise_adjustment(results, interactiol.col)

    dnam.col <- grep("quant_pval_metGrp$",colnames(results),value = TRUE)
    results <- stage_wise_adjustment(results, dnam.col)

    tf.col <- grep("quant_pval_rna.tf$|quant_pval_es.tf$",colnames(results),value = TRUE)
    results <- stage_wise_adjustment(results, tf.col)

    return(results)
}

stage_wise_adjustment <- function(
    results,
    col
){

    min.pval <- results %>%
        dplyr::group_by(.data$regionID) %>%
        dplyr::summarise(min(.data[[col]]))

    results$tripletID <- paste0(
        gsub("[[:punct:]]", "_", results$regionID),
        "_TF_",results$TF_symbol,
        "_target_",results$target
    )

    # Preparing StageR input
    pScreen.pval <- min.pval[[2]]
    names(pScreen.pval) <- gsub("[[:punct:]]", "_", min.pval$regionID)

    pConfirmation <- results[,col] %>% as.matrix  ## LW: change to p-values
    rownames(pConfirmation) <-  results$tripletID
    colnames(pConfirmation) <- "transcript"

    triplet2region <- data.frame(
        row.names = results$tripletID,
        "transcript" = results$tripletID,
        "gene" = gsub("[[:punct:]]", "_", results$regionID)
    )

    pScreen.pval.stageRObj <- stageR::stageRTx(
        pScreen = pScreen.pval,
        pConfirmation = pConfirmation,
        pScreenAdjusted = FALSE,
        tx2gene = triplet2region
    )

    pScreen.pval.stageRObj <- stageR::stageWiseAdjustment(
        object = pScreen.pval.stageRObj,
        method = "dte",
        alpha = 0.05
    )

    padj <- stageR::getAdjustedPValues(
        pScreen.pval.stageRObj,
        onlySignificantGenes = FALSE,
        order = FALSE
    )

    colnames(padj) <- c(
        "regionID",
        "tripletID",
        gsub("pval","region_stage_wise_adj_pval",col),
        gsub("pval","triplet_stage_wise_adj_pval",col)
    )
    padj$regionID <- NULL

    results <- dplyr::left_join(results, padj, by = c("tripletID"))
    results <- results %>%
        relocate(
            c(
                gsub("pval","region_stage_wise_adj_pval",col),
                gsub("pval","triplet_stage_wise_adj_pval",col)
            ),
            .after = col
        )
    results$tripletID <- NULL

    return(results)
}


#' @examples
#' \dontrun{
#' data("dna.met.chr21")
#' dna.met.chr21 <- make_se_from_dnam_probes(dna.met.chr21)
#' data("gene.exp.chr21.log2")
#' triplet <- data.frame(
#'     "regionID" = rownames(dna.met.chr21)[1:10],
#'     "TF" = rownames(gene.exp.chr21.log2)[11:20],
#'     "target" = rownames(gene.exp.chr21.log2)[1:10]
#' )
#' results <- interaction_model(triplet, dna.met.chr21, gene.exp.chr21.log2)
#' results <- calculate_fdr_per_region_adjustment(results)
#' }
#' @noRd
calculate_fdr_per_region_adjustment <- function(results){

    if(!"tripletID" %in% colnames(results)){
        results$tripletID <- create_triplet_ID(results)
    }

    for(pval.col in grep("quant_pval_",colnames(results),value = TRUE)){
        fdr.col <- gsub("pval","fdr",pval.col)
        fdr.by.region <- results %>%
            group_by(.data$regionID) %>%
            summarise(
                "fdr.by.region" = p.adjust(.data[[pval.col]], method = "fdr"),
                "tripletID" = .data$tripletID
            )
        results[[fdr.col]] <-  fdr.by.region$fdr.by.region[match(results$tripletID,fdr.by.region$tripletID)]

        results <- results %>% relocate(fdr.col, .after = pval.col)

    }

    results$tripletID <- NULL
    results
}

#' @examples
#' \dontrun{
#' data("dna.met.chr21")
#' dna.met.chr21 <- make_se_from_dnam_probes(dna.met.chr21)
#' data("gene.exp.chr21.log2")
#' triplet <- data.frame(
#'     "regionID" = rownames(dna.met.chr21)[1:10],
#'     "TF" = rownames(gene.exp.chr21.log2)[11:20],
#'     "target" = rownames(gene.exp.chr21.log2)[1:10]
#' )
#' results <- interaction_model(triplet, dna.met.chr21, gene.exp.chr21.log2)
#' results <- calculate_fdr_per_region_adjustment(results)
#' }
#' @noRd
calculate_fdr_adjustment <- function(results){

    for(pval.col in grep("quant_pval_",colnames(results),value = TRUE)){
        fdr.col <- gsub("pval","fdr",pval.col)
        results[[fdr.col]] <- p.adjust(results[[pval.col]], method = "fdr")
        results <- results %>% relocate(fdr.col, .after = pval.col)
    }
    results
}

#' @examples
#' \dontrun{
#' data("dna.met.chr21")
#' dna.met.chr21 <- make_se_from_dnam_probes(dna.met.chr21)
#' data("gene.exp.chr21.log2")
#' triplet <- data.frame(
#'     "regionID" = rownames(dna.met.chr21)[1:10],
#'     "TF" = rownames(gene.exp.chr21.log2)[11:20],
#'     "target" = rownames(gene.exp.chr21.log2)[1:10]
#' )
#' triplet$tripletID <- create_triplet_ID(triplet)
#' }
#' @noRd
create_triplet_ID <- function(df){
    paste0(gsub("[[:punct:]]", "_", df$regionID),"_TF_",df$TF_symbol,"_target_",df$target)
}
