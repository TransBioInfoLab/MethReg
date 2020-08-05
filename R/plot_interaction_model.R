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
#' p-values and estimates of interaction
#' @param dnam DNA methylation matrix  (columns: samples same order as met, rows: regions/probes)
#' @param exp gene expression matrix (columns: samples same order as met, rows: genes)
#' @param metadata A data frame with samples as rownames and one columns that will be used to
#' color the samples
#' @return A dataframe with Region, TF, Estimates and P-value from linear model
#' @examples
#' dnam <- t(matrix(sort(c(runif(20))), ncol = 1))
#' rownames(dnam) <- c("chr3:203727581-203728580")
#' colnames(dnam) <- paste0("Samples",1:20)
#'
#' exp.target <-  c(runif(10,min = 0,max = 0),
#'                 runif(10,min = 0,max = 10)) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(exp.target) <- c("ENSG00000232886")
#' colnames(exp.target) <- paste0("Samples",1:20)
#'
#' exp.tf <-  t(matrix(sort(c(runif(20))), ncol = 1))
#' rownames(exp.tf) <- c("ENSG00000232888")
#' colnames(exp.tf) <- paste0("Samples",1:20)
#'
#' exp <- rbind(exp.tf, exp.target)
#' # Map example region to closest gene
#' triplet <- data.frame(
#'    "regionID" =  c("chr3:203727581-203728580"),
#'    "target" = "ENSG00000232886",
#'    "TF" = "ENSG00000232888"
#')
#'
#' results <- interaction_model(triplet, dnam, exp)
#' plots <- plot_interaction_model(
#'     triplet.results = results,
#'     dnam = dnam,
#'     exp = exp
#' )
#' \dontrun{
#' data("dna.met.chr21")
#' dna.met.chr21 <- map_probes_to_regions(dna.met.chr21)
#' data("gene.exp.chr21")
#' triplet <- data.frame("regionID" = rownames(dna.met.chr21)[1:10],
#'                       "TF" = rownames(gene.exp.chr21)[11:20],
#'                       "target" = rownames(gene.exp.chr21)[1:10])
#' results <- interaction_model(triplet, dna.met.chr21, gene.exp.chr21)
#' plots <- plot_interaction_model(results[1,], dna.met.chr21, gene.exp.chr21)
#' # Adding color to samples
#' metadata <- clinical[,"sample_type",drop = FALSE]
#' plots <- plot_interaction_model(results[1,], dna.met.chr21, gene.exp.chr21,metadata)
#' }
#' @export
#' @importFrom ggpubr ggscatter ggarrange ggtexttable ttheme stat_cor
#' @importFrom ggplot2 xlab ylab geom_smooth
#' @importFrom tibble as_tibble
plot_interaction_model <-  function(triplet.results,
                                    dnam,
                                    exp,
                                    metadata
){

    if(missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
    if(missing(exp)) stop("Please set exp argument with gene expression matrix")
    if(missing(triplet.results)) stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
    if(!all(c("regionID","TF","target") %in% colnames(triplet.results))) {
        stop("triplet must have the following columns names: regionID, TF, target")
    }

    if(is(dnam,"SummarizedExperiment")) dnam <- assay(dnam)
    if(is(exp,"SummarizedExperiment")) exp <- assay(exp)

    check_data(dnam, exp, metadata)

    out <- plyr::alply(
        .data = triplet.results,
        .margins = 1,
        .fun = function(row.triplet,metadata){

            rna.target <- exp[rownames(exp) == row.triplet$target, , drop = FALSE]
            met <- dnam[rownames(dnam) == as.character(row.triplet$regionID), , drop = FALSE]
            rna.tf <- exp[rownames(exp) == row.triplet$TF, , drop = FALSE]

            df <- data.frame(
                rna.target = rna.target %>% as.numeric,
                met = met %>% as.numeric,
                rna.tf = rna.tf %>% as.numeric
            )

            color <- NULL
            if(!missing(metadata)){
                df <- cbind(df,metadata)
                color <- colnames(metadata)[1]
            }

            plots <- get_plot_results(df, row.triplet, color)

            # Reformat p-values for better looking on the plots
            for(idx in grep("pval|fdr|value",colnames(row.triplet))) {
                row.triplet[,idx] <- format.pval(row.triplet[,idx],digits = 3)
            }
            for(idx in grep("estimate|median|minus",colnames(row.triplet))) {
                row.triplet[,idx] <- format(row.triplet[,idx],digits = 3)
            }

            table.plots <- get_table_plot(row.triplet)

            # Arrange the plots on the same page
            plot.table <- ggarrange(
                ggarrange(table.plots$table.plot.metadata,
                          table.plots$table.plot.lm.all,
                          table.plots$table.plot.lm.quant,
                          ncol = 3),
                ggarrange(plots$tf.target,
                          plots$dnam.target,
                          plots$dnam.tf,
                          ncol = 3),
                plots$tf.target.quantile,
                plots$dnam.target.quantile,
                nrow = 4,
                heights = c(2,2,2.5,2.5))
            plot.table
        }, .progress = "time", metadata = metadata)
    attr(out,"split_type") <- NULL
    attr(out,"split_labels") <- NULL
    names(out) <- paste0(triplet.results$regionID,"_TF_",triplet.results$TF,"_target_",triplet.results$target)
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
              "met.q4_minus_q1")
        ) %>% t() %>% as_tibble(rownames = "Variable")

    tab$Variable <- c(
        "Region ID",
        "Target gene ID",
        "Target gene Symbol",
        "TF gene ID",
        "TF gene Symbol",
        "Diff. DNAm (q4 - q1)"
    )
    table.plot.metadata <- ggtexttable(
        tab,
        rows = NULL,
        cols = NULL,
        theme = ttheme("mOrange", base_size = base_size)
    )

    # Get results for linear model with all samples
    table.plot.lm.all <- get_table_plot_results(row.triplet, type = "all")

    # Get results for linear model with DNAm high samples
    table.plot.lm.quant <- get_table_plot_results(row.triplet, type = "quantile")

    table.plot.list <- list(
        "table.plot.metadata" = table.plot.metadata,
        "table.plot.lm.all" = table.plot.lm.all,
        "table.plot.lm.quantile" = table.plot.lm.quant
    )

    return(table.plot.list)
}

#' @importFrom  stats quantile lm
#' @importFrom MASS rlm
get_plot_results <- function(df, row.triplet, color){

    target.lab <- bquote(atop("Target" ~.(row.triplet$target_symbol %>% as.character())))
    region.lab <- "DNA methylation"
    tf.lab <- bquote(atop("TF" ~.(row.triplet$TF_symbol %>% as.character())))

    # quintile plots met
    quant <-  quantile(df$met,na.rm = TRUE)
    quantile_lower_cutoff <- quant[2]
    quantile_upper_cutoff <- quant[4]
    range1 <- paste0("[",paste(round(quant[1:2],digits = 3),collapse = ","),"]")
    range2 <- paste0("[",paste(round(quant[4:5],digits = 3),collapse = ","),"]")

    df$DNAm.group <- NA
    df$DNAm.group[df$met >= quantile_upper_cutoff] <- paste0("DNAm high quartile ", range2)
    df$DNAm.group[df$met <= quantile_lower_cutoff] <- paste0("DNAm low quartile " , range1)

    df$DNAm.group <- factor(
        df$DNAm.group,
        levels = c(paste0("DNAm low quartile " , range1),
                   paste0("DNAm high quartile ", range2)
        )
    )

    # quintile plots TF
    quant <- quantile(df$rna.tf, na.rm = TRUE)
    quantile_lower_cutoff <- quant[2]
    quantile_upper_cutoff <- quant[4]
    range1 <- paste0("[",paste(round(quant[1:2],digits = 3),collapse = ","),"]")
    range2 <- paste0("[",paste(round(quant[4:5],digits = 3),collapse = ","),"]")

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


    tf.target.plot <- get_plot_results_aux(
        df,
        x = "rna.tf",
        y = "rna.target",
        color = color,
        xlab = tf.lab,
        ylab = target.lab)

    dnam.target.plot <- get_plot_results_aux(
        df,
        x = "met",
        y = "rna.target",
        color = color,
        ylab = target.lab,
        xlab = region.lab)

    dnam.tf.plot <- get_plot_results_aux(
        df,
        x = "met",
        y = "rna.tf",
        color = color,
        ylab = tf.lab,
        xlab = region.lab)

    tf.target.quantile.plot <- get_plot_results_aux(
        df[!is.na(df$DNAm.group),],
        x = "rna.tf",
        xlab = tf.lab,
        y = "rna.target",
        ylab = target.lab,
        facet.by = "DNAm.group",
        color = color
    )

    dnam.target.quantile.plot <- get_plot_results_aux(
        df[!is.na(df$TF.group),],
        x = "met",
        y = "rna.target",
        facet.by = "TF.group",
        color = color,
        ylab = target.lab,
        xlab = region.lab)

    return(list("dnam.target.quantile" = dnam.target.quantile.plot,
                "tf.target.quantile" = tf.target.quantile.plot,
                "dnam.tf" = dnam.tf.plot,
                "dnam.target" = dnam.target.plot,
                "tf.target" = tf.target.plot))
}


#' @importFrom sfsmisc f.robftest
#' @importFrom stats as.formula
#' @noRd
#' @examples
#' df <- data.frame(x = runif(20),y = runif(20))
#' get_plot_results_aux(df, "x","y",NULL, xlab = "x", ylab = "y")
get_plot_results_aux <- function(
    df,
    x,
    y,
    color,
    xlab,
    ylab,
    facet.by
){
    if(missing(facet.by)){
        if(!is.null(color)){
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
    } else{
        if(!is.null(color)){
            p <- ggscatter(
                df,
                x = x,
                y = y,
                facet.by = facet.by,
                color = color,
                size = 1
            )
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
        p <- p + geom_smooth(method = MASS::rlm, se = FALSE)

        rls <- MASS::rlm(
            as.formula(paste0(y, "~",x)),
            data = df,
            psi = psi.bisquare,
            maxit = 100)
    })
    rlm.val <- rls %>% summary %>% coef %>% data.frame
    rlm.val <- rlm.val[-1,1]
    rlm.p.value <- tryCatch({
        ftest <- sfsmisc::f.robftest(rls)
        ftest$p.value
    }, error = function(e){
        # message(e);
        return("NA")
    })

    p <- p + ggplot2::annotate(
        geom = "text",
        x = min(df[[x]], na.rm = TRUE),
        y = max(df[[y]], na.rm = TRUE),
        hjust = 0,
        vjust = 1,
        color = 'blue',
        label = paste0(x, ".", y, ".",
                       "rlm = ",formatC(rlm.val, digits = 2, format = "e"),
                       " pval.rlm = ",  formatC(rlm.p.value, digits = 2, format = "e"))
    )
    return(p)
    # stat_cor(method = "spearman",color = "blue")
}

get_table_plot_results <- function(row.triplet, type){

    base_size <- 9
    if(type == "all"){
        pattern.estimate <- "^estimate"
        pattern.pval <- "^pval"
        title <- "Target ~ TF + DNAm +\n TF * DNAm"
    } else if(type == "quantile"){
        pattern.estimate <- "^quant_estimate"
        pattern.pval <- "^quant_pval"
        title <- "Target ~ TF + \nDNAm Quant. Group +\n TF * DNAm Quant. Group"
    } else if(type == "DNAmlow"){
        pattern.estimate <- "^DNAmlow_estimate"
        pattern.pval <- "^DNAmlow_pval"
        title <- "Target ~ TF\nDNAm low samples"
    } else if(type == "DNAmhigh"){
        pattern.estimate <- "^DNAmhigh_estimate"
        pattern.pval <- "^DNAmhigh_pval"
        title <- "Target ~ TF\nDNAm high samples"
    }

    table.plot.estimate <- row.triplet[,grep(pattern.estimate,colnames(row.triplet),value = T),drop  = FALSE] %>%
        t() %>%
        as_tibble(rownames = "Variable")
    colnames(table.plot.estimate)[2] <- "Estimate"
    table.plot.estimate$Variable <- gsub(paste0(pattern.estimate,"|_"),"", table.plot.estimate$Variable)

    table.plot.pval <- row.triplet[,grep(pattern.pval, colnames(row.triplet), value = T), drop  = FALSE] %>%
        t() %>%
        as_tibble(rownames = "Variable")
    table.plot.pval$Variable <- gsub(paste0(pattern.pval,"|_"),"",table.plot.pval$Variable)
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
        theme = ttheme("mOrange", base_size = base_size)
    )

    table.plot.lm.all
}
