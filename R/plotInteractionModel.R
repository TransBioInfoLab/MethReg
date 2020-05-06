
#' @title Plot interaction data
#' @description Create several plots to show interaction data
#' TF expression with target gene interaction using a linear model
#' \deqn{log2(RNA target) ~ log2(TF) + DNAm + log2(TF) * DNAm}
#'
#' To consider covariates, RNA can also be the residuals.
#' \deqn{log2(RNA target residuals) ~ log2(TF residual) + DNAm + log2(TF residual) * DNAm}
#'
#' @param triplet.results Output from function interaction_model
#' with Region ID, TF  (column name: TF),  and target gene  (column name: target),
#' p-vallues and estimates of interaction
#' @param dnam DNA methylation matrix  (columns: samples same order as met, rows: regions/probes)
#' @param exp gene expression matrix (columns: samples same order as met, rows: genes)
#' @return A dataframe with Region, TF, Estimates and P-value from linear model
#' @examples
#' data("dna.met.chr21")
#' dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
#' data("gene.exp.chr21")
#' triplet <- data.frame("regionID" = rownames(dna.met.chr21)[1:10],
#'                       "TF" = rownames(gene.exp.chr21)[11:20],
#'                       "target" = rownames(gene.exp.chr21)[1:10])
#' results <- interaction_model(triplet, dna.met.chr21, gene.exp.chr21)
#' plots <- plot_interaction_model(results[1,], dna.met.chr21, gene.exp.chr21)
#' @export
#' @importFrom ggpubr ggscatter ggarrange ggtexttable ttheme
#' @importFrom ggplot2 xlab ylab geom_smooth
#' @importFrom tibble as_tibble
plot_interaction_model <-  function(triplet.results,
                                    dnam,
                                    exp
){

    if(missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
    if(missing(exp)) stop("Please set exp argument with gene expression matrix")
    if(missing(triplet.results)) stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
    if(!all(c("regionID","TF","target") %in% colnames(triplet))) {
        stop("triplet must have the following columns names: regionID, TF, target")
    }

    out <- plyr::alply(.data = triplet.results, .margins = 1, .fun = function(row.triplet){

        rna.target <- exp[rownames(exp) == row.triplet$target, , drop = FALSE]
        met <- dnam[rownames(dnam) == as.character(row.triplet$regionID), ]
        rna.tf <- exp[rownames(exp) == row.triplet$TF, , drop = FALSE]

        TCGAbiolinks::get.GRCh.bioMart()
        df <- data.frame(
            rna.target = rna.target %>% as.numeric,
            met = met %>% as.numeric,
            rna.tf = rna.tf %>% as.numeric
        )

        target.lab <- bquote(atop("Target" ~.(row.triplet$target_symbol %>% as.character())))
        region.lab <- "DNA methylation"
        tf.lab <- bquote(atop("TF" ~.(row.triplet$TF_symbol %>% as.character())))

        tf.target.plot <- ggscatter(df,
                                    x = "rna.tf",
                                    y = "rna.target",
                                    size = 1
        ) + xlab(tf.lab) +
            ylab(target.lab) +
            geom_smooth(method = MASS::rlm, se = FALSE)

        dnam.target.plot <- ggscatter(df,
                                      x = "met",
                                      y = "rna.target",
                                      size = 1
        ) + ylab(target.lab) +
            xlab(region.lab) +
            geom_smooth(method = MASS::rlm, se = FALSE)

        dnam.tf.plot <- ggscatter(df,
                                  x = "met",
                                  y = "rna.tf",
                                  size = 1
        ) + ylab(tf.lab)  +
            xlab(region.lab) +
            geom_smooth(method = MASS::rlm, se = FALSE)

        # quintile plots met
        quant <-  quantile(df$met,na.rm = TRUE)
        quantile_lower_cutoff <- quant[2]
        quantile_upper_cutoff <- quant[4]
        range1 <- paste0("[",paste(round(quant[1:2],digits = 3),collapse = ","),"]")
        range2 <- paste0("[",paste(round(quant[4:5],digits = 3),collapse = ","),"]")

        df$group <- NA
        df$group[df$met > quantile_upper_cutoff] <- paste0("DNAm high quartile ", range2)
        df$group[df$met < quantile_lower_cutoff] <- paste0("DNAm low quartile " , range1)

        df$group <- factor(df$group,
                           levels = c(paste0("DNAm low quartile " , range1),
                                      paste0("DNAm high quartile ", range2)
                           )
        )

        tf.target.quantile.plot <- ggscatter(df[!is.na(df$group),],
                                             x = "rna.tf",
                                             y = "rna.target",
                                             facet.by = "group",
                                             size = 1
        ) + xlab(tf.lab)  +
            ylab(target.lab) +
            geom_smooth(method = MASS::rlm, se = FALSE)


        # quintile plots TF
        quant <- quantile(df$rna.tf, na.rm = TRUE)
        quantile_lower_cutoff <- quant[2]
        quantile_upper_cutoff <- quant[4]
        range1 <- paste0("[",paste(round(quant[1:2],digits = 3),collapse = ","),"]")
        range2 <- paste0("[",paste(round(quant[4:5],digits = 3),collapse = ","),"]")

        df$groupTF <- NA
        df$groupTF[df$TF > quantile_upper_cutoff] <- paste0("TF high quartile ", range2)
        df$groupTF[df$TF < quantile_lower_cutoff] <- paste0("TF low quartile " , range1)
        df$groupTF <- factor(df$groupTF,
                             levels = c(
                                 paste0("TF low quartile " , range1),
                                 paste0("TF high quartile ", range2)
                             )
        )

        tf.target.quantile.plot <- ggscatter(
            df[!is.na(df$group),],
            x = "rna.tf",
            y = "rna.target",
            facet.by = "group",
            size = 1
        ) + xlab(tf.lab)  +
            ylab(target.lab) +
            geom_smooth(method = MASS::rlm, se = FALSE)


        # Reformat p-values for better looking on the plots
        for(idx in grep("pval|fdr|value",colnames(row.triplet))) {
            row.triplet[,idx] <- format.pval(row.triplet[,idx],digits = 3)
        }
        for(idx in grep("estimate|median|minus",colnames(row.triplet))) {
            row.triplet[,idx] <- format(row.triplet[,idx],digits = 3)
        }

        base_size <- 9
        table.plot1 <- ggtexttable(
            row.triplet[,c("regionID",
                           "target",
                           "target_symbol",
                           "TF",
                           "TF_symbol")] %>%
                t() %>%
                as_tibble(rownames = "Variable"),
            rows = NULL,
            cols = NULL,
            theme = ttheme("mOrange", base_size = base_size)
        )

        table.plot2a <- row.triplet[,grep("estimate",colnames(row.triplet),value = T),drop  = FALSE] %>%
            t() %>%
            as_tibble(rownames = "Variable")
        colnames(table.plot2a)[2] <- "Estimate"
        table.plot2a$Variable <- gsub("estimate_","", table.plot2a$Variable)
        table.plot2b <- row.triplet[,grep("pval", colnames(row.triplet), value = T), drop  = FALSE] %>%
            t() %>%
            as_tibble(rownames = "Variable")
        table.plot2b$Variable <- gsub("pval_","",table.plot2b$Variable)
        colnames(table.plot2b)[2] <- "P-value"
        table.plot2 <- merge(table.plot2a,table.plot2b, by = "Variable",sort = FALSE)

        table.plot2 <- ggtexttable(table.plot2,
                                   rows = NULL,
                                   cols = c("DNAm vs Target","Estimate","P-Values"),
                                   theme = ttheme("mOrange", base_size = base_size)
        )

        # Arrange the plots on the same page
        plot.table <- ggarrange(
            ggarrange(table.plot1,
                      table.plot2,
                      ncol = 2),
            ggarrange(tf.target.plot,
                      dnam.target.plot,
                      dnam.tf.plot,
                      ncol = 3),
            tf.target.quantile.plot,
            nrow = 3,
            heights = c(1,2,2))
        plot.table

    }, .progress = "time")
    names(out) <- paste0(triplet.results$regionID,"_TF_",triplet.results$TF,"_target_",triplet.results$target)
    out
}
