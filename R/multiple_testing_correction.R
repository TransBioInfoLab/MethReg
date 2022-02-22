
#' @examples
#' results <- data.frame(
#'   "regionID" = c("chr3:203727581-203728580","chr4:203727581-203728580"),
#'   "TF" = c("ENSG00000232886","ENSG00000232887"),
#'   "target" = c("ENSG00000232889","ENSG00000242887"),
#'   "quant_pval_metGrp" = c(0.96,0.96),
#'   "quant_pval_rna.tf" = c(0.2, 0.2),
#'   "quant_pval_metGrp:rna.tf" = c(0.2, 0.2)
#' )
#' results <- calculate_stage_wise_adjustment(results)
#' @noRd
calculate_stage_wise_adjustment <- function(results){

    check_package("stageR")

    interactiol.col <- grep("metGrp:.*_pvalue|metGrp\\.",colnames(results),value = TRUE)
    results <- stage_wise_adjustment(results, interactiol.col)

    dnam.col <- grep("metGrp_pvalue",colnames(results),value = TRUE)
    results <- stage_wise_adjustment(results, dnam.col)

    tf.col <- grep("tf_pvalue$",colnames(results),value = TRUE)
    tf.col <- grep("metGrp",tf.col,invert = TRUE,value = TRUE)
    results <- stage_wise_adjustment(results, tf.col)

    return(results)
}


#' @examples
#' results <- data.frame(
#'   "regionID" = c("chr3:203727581-203728580","chr4:203727581-203728580"),
#'   "TF" = c("ENSG00000232886","ENSG00000232887"),
#'   "target" = c("ENSG00000232889","ENSG00000242887"),
#'   "quant_pval_metGrp" = c(0.96,0.96),
#'   "quant_pval_rna.tf" = c(0.2, 0.2),
#'   "quant_pval_metGrp:rna.tf" = c(0.2, 0.2)
#' )
#' results <- stage_wise_adjustment(results,"quant_pval_metGrp.rna.tf")
#' @noRd
stage_wise_adjustment <- function(
    results,
    col
){
    check_package("stageR")

    if(!"tripletID" %in% colnames(results)){
        results$tripletID <- create_triplet_ID(results)
    }


    min.pval <- results %>%
        dplyr::group_by(.data$regionID) %>%
        dplyr::summarise(min(.data[[col]]))

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
        alpha = 0.05,
        allowNA = TRUE
    )

    padj <- stageR::getAdjustedPValues(
        pScreen.pval.stageRObj,
        onlySignificantGenes = FALSE,
        order = FALSE
    )
    # padj: geneID", "txID" ,"gene" ,"transcript"
    # equivalento to: region, triplet, region pval, triplet.pval

    results[[gsub("pval", "region_stage_wise_adj_pval", col)]] <- padj[[3]][match(results$tripletID, padj[[2]])]
    results[[gsub("pval", "triplet_stage_wise_adj_pval", col)]] <- padj[[4]][match(results$tripletID, padj[[2]])]
    results$tripletID <- NULL
    results <- results %>%
        relocate(
            c(
                gsub("pval","region_stage_wise_adj_pval",col),
                gsub("pval","triplet_stage_wise_adj_pval",col)
            ),
            .after = col
        )

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
#' @importFrom dplyr relocate
calculate_fdr_per_region_adjustment <- function(results){

    if(!"tripletID" %in% colnames(results)){
        results$tripletID <- create_triplet_ID(results)
    }

    for(pval.col in grep("RLM_.*pvalue",colnames(results),value = TRUE)){
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
    paste0(gsub("[[:punct:]]", "_", df$regionID),"_TF_",df$TF,"_target_",df$target)
}
