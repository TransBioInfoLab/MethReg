#' @title Plot interaction model results
#' @description Create several plots to show interaction data
#' TF expression with target gene interaction using a linear model
#' \deqn{log2(RNA target) = log2(TF) + DNAm + log2(TF) * DNAm}
#'
#' To consider covariates, RNA can also be the residuals.
#' \deqn{log2(RNA target residuals) = log2(TF residual) + DNAm + log2(TF residual) * DNAm}
#'
#' @param triplet.results Output from function interaction_model
#' with Region ID, TF  (column name: TF),  and target gene  (column name: target),
#' p-values and estimates of interaction
#' @param dnam DNA methylation matrix or SummarizedExperiment object
#'  (columns: samples same order as met, rows: regions/probes)
#' @param exp gene expression matrix or a SummarizedExperiment object
#' (columns: samples same order as met, rows: genes)
#' @param metadata A data frame with samples as rownames and one columns that will be used to
#' color the samples
#' @param tf.activity.es A matrix with normalized enrichment scores for each TF across all samples
#' to be used in linear models instead of TF gene expression.
#' @return A ggplot object, includes a table with results from fitting interaction model,
#' and the the following scatter plots: 1) TF vs DNAm, 2) Target vs DNAm,
#' 3) Target vs TF, 4) Target vs TF for samples in Q1 and Q4 for DNA methylation,
#' 5) Target vs DNAm for samples in Q1 and Q4 for the TF
#' @examples
#' library(dplyr)
#' dnam <- runif(20, min = 0,max = 1) %>% sort %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(dnam) <- c("chr3:203727581-203728580")
#' colnames(dnam) <- paste0("Samples",1:20)
#'
#' exp.target <-  runif(20,min = 0,max = 10) %>% sort %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(exp.target) <- c("ENSG00000232886")
#' colnames(exp.target) <- paste0("Samples",1:20)
#'
#' exp.tf <- runif(20,min = 0,max = 10) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(exp.tf) <- c("ENSG00000101412")
#' colnames(exp.tf) <- paste0("Samples",1:20)
#'
#' exp <- rbind(exp.tf, exp.target)
#'
#' triplet <- data.frame(
#'    "regionID" =  c("chr3:203727581-203728580"),
#'    "target" = "ENSG00000232886",
#'    "TF" = "ENSG00000101412"
#')
#'
#' results <- interaction_model(
#'   triplet = triplet,
#'   dnam = dnam,
#'   exp = exp,
#'   fdr = FALSE,
#'   filter.correlated.tf.exp.dna = FALSE
#' )
#' plots <- plot_interaction_model(
#'     triplet.results = results,
#'     dnam = dnam,
#'     exp = exp
#' )
#' @export
#' @importFrom ggpubr ggscatter ggarrange ggtexttable ttheme stat_cor
#' @importFrom ggplot2 xlab ylab geom_smooth
#' @importFrom tibble as_tibble
plot_interaction_model <-  function(
    triplet.results,
    dnam,
    exp,
    metadata,
    tf.activity.es = NULL
){

    if (missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")

    if (missing(exp)) stop("Please set exp argument with gene expression matrix")

    if (missing(triplet.results)) {
        stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
    }

    if (!all(c("regionID","TF","target") %in% colnames(triplet.results))) {
        stop("triplet must have the following columns names: regionID, TF, target")
    }

    if (is(dnam,"SummarizedExperiment")) dnam <- assay(dnam)
    if (is(exp,"SummarizedExperiment")) exp <- assay(exp)

    check_data(dnam, exp, metadata)

    out <- plyr::alply(
        .data = triplet.results,
        .margins = 1,
        .fun = function(row.triplet,metadata){


            df <- get_triplet_data(
                exp = exp,
                dnam = dnam,
                row.triplet = row.triplet,
                tf.es =  tf.activity.es
            )

            row.triplet <- cbind(
                row.triplet,
                stratified_model_results(
                    df
                )
            )

            color <- NULL
            if (!missing(metadata)) {
                df <- cbind(df,metadata)
                color <- colnames(metadata)[1]
            }

            # Reformat p-values for better looking on the plots
            idx <- grep("pval|fdr|value",colnames(row.triplet))
            row.triplet[,idx] <- format.pval(row.triplet[,idx],digits = 3)
            idx <- grep("estimate|median|minus",colnames(row.triplet))
            row.triplet[,idx] <- format(row.triplet[,idx],digits = 3)

            plots <- get_plot_results(
                df = df,
                row.triplet = row.triplet,
                color = color,
                use_tf_enrichment_scores = is.null(tf.activity.es)
            )

            table.plots <- get_table_plot(row.triplet)

            # Arrange the plots on the same page
            suppressWarnings({
                suppressMessages({

                    plot.table <-
                        ggarrange(
                            ggarrange(
                                table.plots$table.plot.metadata,
                                #table.plots$table.plot.wilcoxon,
                                #table.plots$table.plot.lm.all,
                                table.plots$table.plot.lm.quant,
                                table.plots$table.plot.lm.dna.low,
                                table.plots$table.plot.lm.dna.high,
                                heights = c(0.8,0.5,0.25,0.25),
                                ncol = 1),
                            ggarrange(
                                ggarrange(
                                    plots$tf.target,
                                    plots$dnam.target,
                                    #plots$dnam.tf,
                                    ncol = 2),
                                plots$tf.target.quantile,
                                plots$dnam.target.quantile,
                                nrow = 3
                            ), ncol = 2,widths = c(1,2)
                        )
                })
            })
            plot.table
        }, .progress = "time", metadata = metadata)
    attr(out,"split_type") <- NULL
    attr(out,"split_labels") <- NULL

    if (nrow(triplet.results) > 0) {
        names(out) <- paste0(
            triplet.results$regionID,
            "_TF_",
            triplet.results$TF,
            "_target_",
            triplet.results$target
        )
    }
    out
}

get_table_plot <- function(row.triplet){

    base_size <- 9
    tab <- row.triplet %>%
        dplyr::select(
            c("regionID",
              "target",
              "target_symbol",
              "TF",
              "TF_symbol",
              "met.IQR",
              "TF.role",
              "DNAm.effect")
        ) %>% t() %>% as_tibble(rownames = "Variable")

    tab$Variable <- c(
        "Region ID",
        "Target gene ID",
        "Target gene Symbol",
        "TF gene ID",
        "TF gene Symbol",
        "Diff. DNAm (Q4 - Q1)",
        "TF role",
        "DNAm effect"
    )

    table.plot.metadata <- ggtexttable(
        tab,
        rows = NULL,
        cols = NULL,
        theme = ttheme("mGreen", base_size = base_size)
    )

    tab <- row.triplet %>%
        dplyr::select(
            c("Wilcoxon_pval_tf_q4met_vs_q1met")
        ) %>% t() %>% as_tibble(rownames = "Variable")

    tab$Variable <- c(
        "TF Q4 vs TF Q1"
    )
    table.plot.wilcoxon <- ggtexttable(
        tab,
        rows = NULL,
        cols = c("Wilcoxon","P-Values"),
        theme = ttheme("mGreen", base_size = base_size)
    )

    # Get results for linear model with all samples
    #table.plot.lm.all <- get_table_plot_results(row.triplet, type = "all")

    # Get results for linear model with DNAm high samples
    table.plot.lm.quant <- get_table_plot_results(row.triplet, type = "quantile")

    # Get results for linear model with all samples
    table.plot.lm.dna.low <- get_table_plot_results(row.triplet, type = "DNAmlow")

    # Get results for linear model with DNAm high samples
    table.plot.lm.dna.high <- get_table_plot_results(row.triplet, type = "DNAmhigh")

    table.plot.list <- list(
        "table.plot.metadata" = table.plot.metadata,
        #"table.plot.lm.all" = table.plot.lm.all,
        "table.plot.lm.quantile" = table.plot.lm.quant,
        "table.plot.wilcoxon" = table.plot.wilcoxon,
        "table.plot.lm.dna.low" = table.plot.lm.dna.low,
        "table.plot.lm.dna.high" = table.plot.lm.dna.high
    )

    return(table.plot.list)
}

#' @importFrom stats quantile lm
#' @importFrom MASS rlm
get_plot_results <- function(
    df,
    row.triplet,
    color,
    use_tf_enrichment_scores = FALSE
){

    target.lab <- bquote(atop("Target" ~.(row.triplet$target_symbol %>% as.character())))
    region.lab <- "DNA methylation"

    if (use_tf_enrichment_scores) {
        tf.lab <- bquote(atop("Enrichment scores TF" ~.(row.triplet$TF_symbol %>% as.character())))
    } else {
        tf.lab <- bquote(atop("TF" ~.(row.triplet$TF_symbol %>% as.character())))
    }

    # quintile plots met
    quant <- quantile(df$met,na.rm = TRUE)
    quantile_lower_cutoff <- quant[2]
    quantile_upper_cutoff <- quant[4]
    range1 <- paste0("[",paste(round(quant[1:2], digits = 3), collapse = ","),"]")
    range2 <- paste0("[",paste(round(quant[4:5], digits = 3), collapse = ","),"]")

    df$DNAm.group <- NA
    df$DNAm.group[df$met >= quantile_upper_cutoff] <- paste0("DNAm high quartile ", range2)
    df$DNAm.group[df$met <= quantile_lower_cutoff] <- paste0("DNAm low quartile " , range1)

    df$DNAm.group <- factor(
        df$DNAm.group,
        levels = c(
            paste0("DNAm low quartile " , range1),
            paste0("DNAm high quartile ", range2)
        )
    )

    # quintile plots TF
    quant <- quantile(df$rna.tf, na.rm = TRUE)
    quantile_lower_cutoff <- quant[2]
    quantile_upper_cutoff <- quant[4]
    range1 <- paste0("[",paste(round(quant[1:2], digits = 3),collapse = ","),"]")
    range2 <- paste0("[",paste(round(quant[4:5], digits = 3),collapse = ","),"]")

    df$TF.group <- NA
    df$TF.group[df$rna.tf >= quantile_upper_cutoff] <- paste0("TF high quartile ", range2)
    df$TF.group[df$rna.tf <= quantile_lower_cutoff] <- paste0("TF low quartile ", range1)
    df$TF.group <- factor(
        df$TF.group,
        levels = c(
            paste0("TF low quartile " , range1),
            paste0("TF high quartile ", range2)
        )
    )


    tf.target.plot <- get_scatter_plot_results(
        df,
        x = "rna.tf",
        y = "rna.target",
        color = color,
        xlab = tf.lab,
        ylab = target.lab
    )

    #dnam.target.plot <- get_scatter_plot_results(
    #    df,
    #    x = "met",
    #    y = "rna.target",
    #    color = color,
    #    ylab = target.lab,
    #    xlab = region.lab
    #)

    dnam.target.plot <- get_box_plot_results(
        df,
        y = "rna.target",
        facet.by = "DNAm.group",
        ylab = target.lab
    )

    # dnam.target.plot <- get_histogram_plot_results(
    #     df,
    #    x = "rna.target",
    #    facet.by = "DNAm.group",
    #    xlab = target.lab
    # )

    dnam.tf.plot <- get_scatter_plot_results(
        df,
        x = "met",
        y = "rna.tf",
        color = color,
        ylab = tf.lab,
        xlab = region.lab
    )

    tf.target.quantile.plot <- get_scatter_plot_results(
        df[!is.na(df$DNAm.group),],
        x = "rna.tf",
        xlab = tf.lab,
        y = "rna.target",
        ylab = target.lab,
        facet.by = "DNAm.group",
        color = color # ifelse(is.null(color),"met",color)
    )

    dnam.target.quantile.plot <- get_scatter_plot_results(
        df[!is.na(df$TF.group),],
        x = "met",
        y = "rna.target",
        facet.by = "TF.group",
        color = color,
        ylab = target.lab,
        xlab = region.lab
    )

    return(
        list(
            "dnam.target.quantile" = dnam.target.quantile.plot,
            "tf.target.quantile" = tf.target.quantile.plot,
            "dnam.tf" = dnam.tf.plot,
            "dnam.target" = dnam.target.plot,
            "tf.target" = tf.target.plot
        )
    )
}

#' @noRd
#' @examples
#' df <- data.frame(
#'    x = runif(20),
#'    group = c(rep("low",10),rep("high",10))
#' )
#' get_histogram_plot_results(df = df, x =  "x",facet.by = "group", xlab = "expr")
get_histogram_plot_results <- function(
    df,
    x,
    facet.by,
    xlab
){

    df <- df[!is.na(df[[facet.by]]),]
    df[[facet.by]] <-  ifelse(grepl("low",df[[facet.by]]),"DNAm.low","DNAm.high")

    suppressWarnings({
        p <- ggpubr::gghistogram(
            df,
            x = x,
            add = "mean",
            rug = TRUE,
            fill = facet.by,
            palette = c("#00AFBB", "#E7B800"),
            add_density = FALSE
        ) + xlab(xlab) + ggplot2::theme(legend.title = ggplot2::element_blank())
    })
    p
}

#' @noRd
#' @examples
#' df <- data.frame(
#'    exp = runif(20),
#'    group = c(rep("high",10),rep("low",10))
#' )
#' get_box_plot_results(df = df, y =  "exp",facet.by = "group", ylab = "expr")
get_box_plot_results <- function(
    df,
    y,
    facet.by,
    ylab
){

    df <- df[!is.na(df[[facet.by]]),]
    df[[facet.by]] <- factor(
        ifelse(grepl("low",df[[facet.by]]),"DNAm.low","DNAm.high"),
        levels = c("DNAm.low","DNAm.high")
    )

    suppressWarnings({
        p <- ggpubr::ggboxplot(
            data = df,
            x = facet.by,
            y = y,
            add = "jitter",
            color = facet.by,
            palette = c("#477dbf", "#d04042"),
            add_density = FALSE
        ) + ylab(ylab) +
            ggplot2::theme(legend.title = ggplot2::element_blank()) +
            xlab("") +
            ggpubr::stat_compare_means(label.y = max(df[[y]]) * 1.1) +
            ggplot2::guides( color = FALSE) +
            ggplot2::ylim(min(df[[y]]),max(df[[y]]) * 1.2)
    })
    p
}

#' @importFrom sfsmisc f.robftest
#' @importFrom stats as.formula
#' @noRd
#' @examples
#' df <- data.frame(x = runif(20),y = runif(20))
#' get_scatter_plot_results(df, "x","y",NULL, xlab = "x", ylab = "y")
get_scatter_plot_results <- function(
    df,
    x,
    y,
    color,
    xlab,
    ylab,
    facet.by
){
    if (missing(facet.by)) {
        if (!is.null(color)) {
            p <- ggscatter(
                df,
                x = x,
                y = y,
                color = color,
                size = 1
            )
        } else {
            p <- ggscatter(
                df,
                x = x,
                y = y,
                size = 1
            )
        }
    } else {
        if (!is.null(color)) {
            p <- ggscatter(
                df,
                x = x,
                y = y,
                facet.by = facet.by,
                color = color,
                size = 1
            ) #+ ggplot2::theme(legend.position = "right")
        } else {
            p <- ggscatter(
                df,
                x = x,
                y = y,
                facet.by = facet.by,
                size = 1
            )
        }
    }

    p <- p + xlab(xlab) + ylab(ylab)

    suppressWarnings({
        suppressMessages({
            p <- p + geom_smooth(method = MASS::rlm, se = FALSE)
        })
    })

    if (missing(facet.by)) {
        rlm.res <- get_rlm_val_pval(df, x, y)

        p <- p + ggplot2::annotate(
            geom = "text",
            x = min(df[[x]], na.rm = TRUE),
            y = max(df[[y]] * 1.4, na.rm = TRUE),
            hjust = 0,
            vjust = 1,
            color = 'blue',
            label = paste0(
                #gsub("rna\\.","",y), " ~ ", gsub("rna\\.","",x),
                "RLM estimate = ",formatC(rlm.res$rlm.val, digits = 2, format = "e"),
                "\nRLM p-value = ",  formatC(rlm.res$rlm.p.value, digits = 2, format = "e"))
        )
    } else {
        # lower Annotation
        rlm.res.low <- get_rlm_val_pval(df %>% dplyr::filter(grepl("low",df[[facet.by]])),x , y)

        ann_text.low <- data.frame(
            x = min(df[[x]], na.rm = TRUE),
            y = max(df[[y]] * 1.4, na.rm = TRUE),
            facet.by = factor(grep("low",df[[facet.by]],value = TRUE),levels = unique(df[[facet.by]]))
        )
        colnames(ann_text.low) <- c(x,y,facet.by)

        # higher Annotation
        p <- p + ggplot2::geom_text(
            data = ann_text.low,
            hjust = 0,
            vjust = 1,
            color = 'blue',
            label = paste0(
                # gsub("rna\\.","",y), " ~ ", gsub("rna\\.","",x),
                "RLM estimate = ",formatC(rlm.res.low$rlm.val, digits = 2, format = "e"),
                "\nRLM p-value  = ",  formatC(rlm.res.low$rlm.p.value, digits = 2, format = "e"))
        )

        rlm.res.high <- get_rlm_val_pval(df %>% dplyr::filter(grepl("high",df[[facet.by]])), x , y)
        ann_text.high <- data.frame(
            x = min(df[[x]], na.rm = TRUE),
            y = max(df[[y]] * 1.4, na.rm = TRUE),
            facet.by = factor(grep("high",df[[facet.by]],value = TRUE),levels = unique(df[[facet.by]]))
        )
        colnames(ann_text.high) <- c(x,y,facet.by)

        # higher Annotation
        p <- p + ggplot2::geom_text(
            data = ann_text.high,
            hjust = 0,
            vjust = 1,
            color = 'blue',
            label = paste0(
                #gsub("rna\\.","",y), " ~ ", gsub("rna\\.","",x),
                "RLM estimate = ",formatC(rlm.res.high$rlm.val, digits = 2, format = "e"),
                "\nRLM p-value = ",  formatC(rlm.res.high$rlm.p.value, digits = 2, format = "e"))
        )
    }
    return(p)
    # stat_cor(method = "spearman",color = "blue")
}

get_rlm_val_pval <- function(df, x, y){
    suppressMessages({
        suppressWarnings({
            rls <- MASS::rlm(
                formula = as.formula(paste0(y, "~",x)),
                data = df,
                psi = psi.bisquare,
                maxit = 100
            )
            rlm <- rls %>% summary %>% coef %>% data.frame
            rlm.val <- rlm[-1,1]
        })
    })

    rlm.p.value <- tryCatch({
        degrees.freedom.value <- nrow(df) - 2
        pval <- 2 * (1 - pt( abs(rlm$t.value[-1]), df = degrees.freedom.value))
        pval
        #ftest <- sfsmisc::f.robftest(rls)
        #ftest$p.value
    }, error = function(e){
        #message(e);
        return("NA")
    })

    return(
        list(
            rlm.p.value = rlm.p.value,
            rlm.val = rlm.val
        )
    )
}

get_table_plot_results <- function(row.triplet, type){

    base.size <- 9

    if (type == "all") {
        pattern.estimate <- "^estimate"
        pattern.pval <- "^pval"
        title <- "Target ~ TF + DNAm +\n TF * DNAm"
        theme.color <- "mOrange"
    } else if (type == "quantile") {
        pattern.estimate <- "^quant_estimate"
        pattern.pval <- "^quant_pval"
        title <- "Target ~ TF + \nDNAm Quant. Group +\n TF * DNAm Quant. Group"
        theme.color <- "mGreen"
    } else if (type == "DNAmlow") {
        pattern.estimate <- "^DNAmlow_estimate"
        pattern.pval <- "^DNAmlow_pval"
        title <- "Target ~ TF\nDNAm low samples"
        theme.color <- "mBlue"
    } else if (type == "DNAmhigh") {
        pattern.estimate <- "^DNAmhigh_estimate"
        pattern.pval <- "^DNAmhigh_pval"
        title <- "Target ~ TF\nDNAm high samples"
        theme.color <- "mRed"
    }

    col.idx <- grep(pattern.estimate,colnames(row.triplet),value = TRUE)
    table.plot.estimate <- row.triplet[,col.idx,drop  = FALSE] %>%
        t() %>%
        as_tibble(rownames = "Variable")
    colnames(table.plot.estimate)[2] <- "Estimate"
    table.plot.estimate$Variable <- gsub(
        paste0(pattern.estimate,"|_"),"", table.plot.estimate$Variable
    )

    table.plot.estimate$Variable[grep("^rna.tf$|^es.tf$",table.plot.estimate$Variable)] <- "Direct effect of TF"
    table.plot.estimate$Variable[grep(":rna.tf$|:es.tf$",table.plot.estimate$Variable)] <- "Synergistic effect of DNA and TF"
    table.plot.estimate$Variable[grep("^met",table.plot.estimate$Variable)] <- "Direct effect of DNAm"


    col.idx <- grep(pattern.pval, colnames(row.triplet), value = TRUE)
    table.plot.pval <- row.triplet[,col.idx, drop  = FALSE] %>%
        t() %>%
        as_tibble(rownames = "Variable")

    table.plot.pval$Variable <- gsub(
        pattern = paste0(pattern.pval,"|_"),
        replacement = "",
        table.plot.pval$Variable
    )

    table.plot.pval$Variable[grep("^rna.tf$|^es.tf$",table.plot.pval$Variable)] <- "Direct effect of TF"
    table.plot.pval$Variable[grep(":rna.tf$|:es.tf$",table.plot.pval$Variable)] <- "Synergistic effect of DNA and TF"
    table.plot.pval$Variable[grep("^met",table.plot.pval$Variable)] <- "Direct effect of DNAm"

    colnames(table.plot.pval)[2] <- "P-value"

    table.plot <- merge(
        table.plot.estimate,
        table.plot.pval,
        by = "Variable",
        sort = FALSE
    )

    table.plot.lm.all <- ggtexttable(
        table.plot,
        rows = NULL,
        cols = c(title,"Estimate","P-Values"),
        theme = ttheme(theme.color, base_size = base.size)
    )

    table.plot.lm.all
}
