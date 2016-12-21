#' Make a MA plot of some limma output with pretty colors and shapes
#'
#' Yay pretty colors and shapes!
#'
#' @param pairwise  The result from all_pairwise(), which should be changed to handle other invocations too.
#' @param type  Type of table to use: deseq, edger, limma, basic.
#' @param table  Result from edger to use, left alone it chooses the first.
#' @param fc  Cutoff for log2(fold-change) significant.
#' @param pval_cutoff  Cutoff to define 'significant' by p-value.
#' @param ...  Extra arguments are passed to arglist.
#' @return a plot!
#' @seealso \link{plot_ma_de}
#' @examples
#'  \dontrun{
#'   prettyplot <- edger_ma(all_aprwise) ## [sic, I'm witty! and can speel]
#' }
#' @export
extract_de_ma <- function(pairwise, type="edger", table=NULL, fc=1, pval_cutoff=0.05, ...) {
    arglist <- list(...)

    expr_col <- NULL
    p_col <- NULL
    fc_col <- "logFC"
    if (type == "edger") {
        expr_col <- "logCPM"
        p_col <- "FDR"
    } else if (type == "deseq") {
        expr_col <- "baseMean"  ## This column will need to be changed from base 10 to log scale.
        p_col <- "adj.P.Val"
    } else if (type == "limma") {
        expr_col <- "AveExpr"
        p_col <- "adj.P.Val"
    } else if (type == "basic") {
        expr_col <- "numerator_median"
        p_col <- "p"
    } else {
        stop("Unable to determine the type of plot to perform: deseq, limma, edger, or basic.")
    }

    ## Step 1,  figure out of this is the result of: pairwise_edger(),
    ## all_pairwise(),  or combine_de_tables().
    ## Check if this is the result of all_pairwise()
    ## If it is the result of all_pairwise(),  then it will have
    ## the result of de_edger() sitting in the edger slot.
    ## Therefore we are really only testing between combine() or all()

    ## This is an explicit test against all_pairwise() and reduces it to the edger result.
    if (!is.null(pairwise[[type]])) {
        pairwise <- pairwise[[type]]
    }
    ## Unfortunately, the structure of a combine_de_tables() is different from that
    ## of de_edger().  Ergo,  a little work will be required to get the data out
    ## properly.

    ## combined_de_tables provides: slots including 'data' and 'de_summary'
    ## while de_edger provides: 'all_tables',  'contrasts_performed' and such.
    possible_tables <- NULL
    all_tables <- NULL
    the_table <- NULL
    if (is.null(table)) {
        table <- 1
    }
    if (!is.null(pairwise[["all_tables"]])) {
        ## If there is a 'all_tables', slot,  then this is de_edger output.
        possible_tables <- names(pairwise[["all_tables"]])
        all_tables <- pairwise[["all_tables"]]
    } else {
        ## Then this is combined_tables() output.
        possible_tables <- names(pairwise[["data"]])
        all_tables <- pairwise[["data"]]
        fc_col <- paste0(type, "_logfc")
        p_col <- paste0(type, "_adjp")
        expr_col <- NULL
        if (type == "deseq") {
            expr_col <- "deseq_basemean"
        } else if (type == "edger") {
            expr_col <- "edger_logcpm"
        } else if (type == "limma") {
            expr_col <- "limma_ave"
        } else if (type == "basic") {
            expr_col <- "basic_nummed"
        }
    }

    if (is.numeric(table)) {
        the_table <- possible_tables[[table]]
    } else {
        the_table <- table
        revname <- strsplit(x=table, split="_vs_")
        revname <- paste0(revname[[1]][2], "_vs_", revname[[1]][1])
        if (!(table %in% possible_tables) & revname %in% possible_tables) {
            message("Trey you doofus, you reversed the name of the table.")
            the_table <- revname
        } else if (!(table %in% possible_tables) & !(revname %in% possible_tables)) {
            stop("Unable to find the table in the set of possible tables.")
        }
    }
    the_table <- all_tables[[the_table]]
    ## DESeq2 returns the median values as base 10, but we are using log2 (or log10?)
    if (type == "deseq") {
        the_table[["log_basemean"]] <- log2(x=the_table[[expr_col]] + 1.0)
        expr_col <- "log_basemean"
    }

    ma_material <- plot_ma_de(table=the_table, expr_col=expr_col, fc_col=fc_col,
                              p_col=p_col, logfc_cutoff=fc, pval_cutoff=pval_cutoff) ##, ...)
    return(ma_material)
}

#' Perform a coefficient scatter plot of a limma/deseq/edger/basic table.
#'
#' Plot the gene abundances for two coefficients in a differential expression comparison.
#' By default, genes past 1.5 z scores from the mean are colored red/green.
#'
#' @param output  Result from the de_ family of functions, all_pairwise, or combine_de_tables().
#' @param toptable  Chosen table to query for abundances.
#' @param type  Query limma, deseq, edger, or basic outputs.
#' @param x  The x-axis column to use, either a number of name.
#' @param y  The y-axis column to use.
#' @param z  Define the range of genes to color (FIXME: extend this to p-value and fold-change).
#' @param p  Set a p-value cutoff for coloring the scatter plot (currently not supported).
#' @param fc  Set a fold-change cutoff for coloring points in the scatter plot (currently not supported.)
#' @param n  Set a top-n fold-change for coloring the points in the scatter plot (this should work, actually).
#' @param color_low  Color for the genes less than the mean.
#' @param color_high  Color for the genes greater than the mean.
#' @param ...  More arguments are passed to arglist.
#' @export
extract_coefficient_scatter <- function(output, toptable=NULL, type="limma", x=1, y=2, z=1.5,
                                        p=NULL, fc=NULL, n=NULL,
                                        color_low="#DD0000", color_high="#7B9F35", ...) {
    arglist <- list(...)
    ## This is an explicit test against all_pairwise() and reduces it to result from type.
    if (!is.null(output[[type]])) {
        output <- output[[type]]
    }

    qlimit <- 0.1
    if (!is.null(arglist[["qlimit"]])) {
        qlimit <- arglist[["qlimit"]]
    }
    fc_column <- paste0(type, "_logfc")
    if (!is.null(arglist[["fc_column"]])) {
        fc_column <- arglist[["fc_column"]]
    }
    p_column <- paste0(type, "_adjp")
    if (!is.null(arglist[["p_column"]])) {
        p_column <- arglist[["p_column"]]
    }
    gvis_filename <- NULL
    gvis_trendline <- TRUE
    tooltip_data <- NULL
    base_url <- NULL
    if (!is.null(arglist[["gvis_filename"]])) {
        gvis_filename <- arglist[["gvis_filename"]]
    }
    if (!is.null(arglist[["gvis_trendline"]])) {
        gvis_trendline <- arglist[["gvis_trendline"]]
    }
    if (!is.null(arglist[["tooltip_data"]])) {
        tooltip_data <- arglist[["tooltip_data"]]
    }
    if (!is.null(arglist[["base_url"]])) {
        base_url <- arglist[["base_url"]]
    }

    ## Extract the set of availble names -- FIXME this should be standardized!!
    thenames <- NULL
    if (type == "edger") {
        thenames <- names(output[["contrasts"]][["identities"]])
    } else if (type == "limma") {
        coefficients <- output[["pairwise_comparisons"]][["coefficients"]]
        thenames <- colnames(coefficients)
    } else if (type == "deseq") {
        thenames <- names(output[["coefficients"]])
    } else if (type == "basic") {
        thenames <- names(output[["conditions_table"]])
    } else {
        stop("I do not know what type you wish to query.")
    }

    message("This can do comparisons among the following columns in the pairwise result:")
    message(toString(thenames))
    xname <- ""
    yname <- ""
    if (is.numeric(x)) {
        xname <- thenames[[x]]
    } else {
        xname <- x
    }
    if (is.numeric(y)) {
        yname <- thenames[[y]]
    } else {
        yname <- y
    }
    message(paste0("Actually comparing ", xname, " and ", yname, "."))

    ## Now extract the coefficent df
    if (type == "edger") {
        coefficient_df <- output[["lrt"]][[1]][["coefficients"]]
        coefficient_df <- coefficient_df[, c(xname, yname)]
        if (max(coefficient_df) < 0) {
            coefficient_df <- coefficient_df * -1.0
        }
    } else if (type == "limma") {
        coefficient_df <- output[["pairwise_comparisons"]][["coefficients"]]
        coefficient_df <- coefficients[, c(x, y)]
    } else if (type == "deseq") {
        first_df <- output[["coefficients"]][[xname]]
        first_df[["delta"]] <- log2(as.numeric(first_df[["baseMean"]])) + first_df[["log2FoldChange"]]
        second_df <- output[["coefficients"]][[yname]]
        second_df[["delta"]] <- log2(as.numeric(second_df[["baseMean"]])) + second_df[["log2FoldChange"]]
        first_col <- first_df[, c("baseMean", "log2FoldChange", "delta")]
        colnames(first_col) <- c("mean.1", "fc.1", xname)
        second_col <- second_df[, c("baseMean", "log2FoldChange", "delta")]
        colnames(second_col) <- c("mean.2", "fc.2", yname)
        coefficient_df <- merge(first_col, second_col, by="row.names")
        rownames(coefficient_df) <- coefficient_df[["Row.names"]]
        coefficient_df <- coefficient_df[-1]
        coefficient_df <- coefficient_df[, c(xname, yname)]
        coefficient_df[is.na(coefficient_df)] <- 0
    } else if (type == "basic") {
        coefficient_df <- output[["medians"]]
        coefficient_df <- coefficient_df[, c(xname, yname)]
    }

    maxvalue <- max(coefficient_df) + 1.0
    plot <- sm(plot_linear_scatter(df=coefficient_df, loess=TRUE, gvis_filename=gvis_filename,
                                   gvis_trendline=gvis_trendline, first=xname, second=yname,
                                   tooltip_data=tooltip_data, base_url=base_url,
                                   pretty_colors=FALSE, color_low=color_low, color_high=color_high,
                                   p=p, fc=fc, n=n, z=z))
    plot[["scatter"]] <- plot[["scatter"]] +
        ggplot2::scale_x_continuous(limits=c(0, maxvalue)) +
        ggplot2::scale_y_continuous(limits=c(0, maxvalue))
    plot[["df"]] <- coefficient_df
    return(plot)
}

#' Create venn diagrams describing how well deseq/limma/edger agree.
#'
#' The sets of genes provided by limma and friends would ideally always agree,  but they do not.
#' Use this to see out how much the (dis)agree.
#'
#' @param table Which table to query?
#' @param adjp  Use adjusted p-values
#' @param ... More arguments are passed to arglist.
#' @export
de_venn <- function(table, adjp=FALSE, euler=FALSE, p=0.05, ...) {
    arglist <- list(...)
    combine_tables <- function(d, e, l) {
        ddf <- as.data.frame(l[, "limma_logfc"])
        rownames(ddf) <- rownames(l)
        colnames(ddf) <- c("limma_logfc")
        ddf <- merge(ddf, e, by="row.names", all=TRUE)
        rownames(ddf) <- ddf[["Row.names"]]
        ddf <- ddf[, -1]
        ddf <- ddf[, c("limma_logfc.x", "edger_logfc")]
        ddf <- merge(ddf, d, by="row.names", all=TRUE)
        rownames(ddf) <- ddf[["Row.names"]]
        ddf <- ddf[, -1]
        ddf <- ddf[, c("limma_logfc.x", "edger_logfc.x", "deseq_logfc")]
        colnames(ddf) <- c("limma", "edger", "deseq")
        return(ddf)
    }

    limma_p <- "limma_p"
    deseq_p <- "deseq_p"
    edger_p <- "edger_p"
    if (isTRUE(adjp)) {
        limma_p <- "limma_adjp"
        deseq_p <- "deseq_adjp"
        edger_p <- "edger_adjp"
    }

    limma_sig <- sm(get_sig_genes(table, column="limma_logfc", p_column=limma_p, p=p))
    edger_sig <- sm(get_sig_genes(table, column="edger_logfc", p_column=edger_p, p=p))
    deseq_sig <- sm(get_sig_genes(table, column="deseq_logfc", p_column=deseq_p, p=p))
    comp_up <- combine_tables(deseq_sig[["up_genes"]],
                              edger_sig[["up_genes"]],
                              limma_sig[["up_genes"]])
    comp_down <- combine_tables(deseq_sig[["down_genes"]],
                                edger_sig[["down_genes"]],
                                limma_sig[["down_genes"]])

    up_d <- sum(!is.na(comp_up[["deseq"]]) & is.na(comp_up[["edger"]]) & is.na(comp_up[["limma"]]))
    up_e <- sum(is.na(comp_up[["deseq"]]) & !is.na(comp_up[["edger"]]) & is.na(comp_up[["limma"]]))
    up_l <- sum(is.na(comp_up[["deseq"]]) & is.na(comp_up[["edger"]]) & !is.na(comp_up[["limma"]]))
    up_de <- sum(!is.na(comp_up[["deseq"]]) & !is.na(comp_up[["edger"]]) & is.na(comp_up[["limma"]]))
    up_dl <- sum(!is.na(comp_up[["deseq"]]) & is.na(comp_up[["edger"]]) & !is.na(comp_up[["limma"]]))
    up_el <- sum(is.na(comp_up[["deseq"]]) & !is.na(comp_up[["edger"]]) & !is.na(comp_up[["limma"]]))
    up_del <- sum(!is.na(comp_up[["deseq"]]) & !is.na(comp_up[["edger"]]) & !is.na(comp_up[["limma"]]))

    down_d <- sum(!is.na(comp_down[["deseq"]]) & is.na(comp_down[["edger"]]) & is.na(comp_down[["limma"]]))
    down_e <- sum(is.na(comp_down[["deseq"]]) & !is.na(comp_down[["edger"]]) & is.na(comp_down[["limma"]]))
    down_l <- sum(is.na(comp_down[["deseq"]]) & is.na(comp_down[["edger"]]) & !is.na(comp_down[["limma"]]))
    down_de <- sum(!is.na(comp_down[["deseq"]]) & !is.na(comp_down[["edger"]]) & is.na(comp_down[["limma"]]))
    down_dl <- sum(!is.na(comp_down[["deseq"]]) & is.na(comp_down[["edger"]]) & !is.na(comp_down[["limma"]]))
    down_el <- sum(is.na(comp_down[["deseq"]]) & !is.na(comp_down[["edger"]]) & !is.na(comp_down[["limma"]]))
    down_del <- sum(!is.na(comp_down[["deseq"]]) & !is.na(comp_down[["edger"]]) & !is.na(comp_down[["limma"]]))

    up_ones <- c("d" = up_d, "e" = up_e, "l" = up_l)
    up_twos <- c("d&e" = up_de, "d&l" = up_dl, "e&l" = up_el)
    up_threes <- c("d&e&l" = up_del)
    up_fun <- up_venneuler <- up_venn_data <- NULL
    if (isTRUE(euler)) {
        tt <- sm(require.auto("venneuler"))
        up_fun <- plot_fun_venn(ones=up_ones, twos=up_twos, threes=up_threes)
        up_venneuler <- up_fun[["plot"]]
        up_venn_data <- up_fun[["data"]]
    }
    tt <- sm(require.auto("js229/Vennerable"))
    up_venn <- Vennerable::Venn(SetNames = c("d", "e", "l"),
                                Weight = c(0, up_d, up_e, up_de,
                                           up_l, up_dl, up_el,
                                           up_del))
    up_res <- Vennerable::plot(up_venn, doWeights=FALSE)
    ##up_res <- plot(up_venn, doWeights=FALSE)
    up_venn_noweight <- grDevices::recordPlot()

    down_ones <- c("d" = down_d, "e" = down_e, "l" = down_l)
    down_twos <- c("d&e" = down_de, "d&l" = down_dl, "e&l" = down_el)
    down_threes <- c("d&e&l" = down_del)
    down_fun <- down_venneuler <- down_venn_data <- NULL
    if (isTRUE(euler)) {
        down_fun <- plot_fun_venn(ones=down_ones, twos=down_twos, threes=down_threes)
        down_venneuler <- down_fun[["plot"]]
        down_venn_data <- down_fun[["data"]]
    }
    down_venn <- Vennerable::Venn(SetNames = c("d", "e", "l"),
                                  Weight = c(0, down_d, down_e, down_de,
                                             down_l, down_dl, down_el,
                                             down_del))
    down_res <- Vennerable::plot(down_venn, doWeights=FALSE)
    ##down_res <- plot(down_venn, doWeights=FALSE)
    down_venn_noweight <- grDevices::recordPlot()

    retlist <- list(
        "up_venn" = up_venn,
        "up_noweight" = up_venn_noweight,
        "up_data" = comp_up,
        "down_venn" = down_venn,
        "down_noweight" = down_venn_noweight,
        "down_data" = comp_down)
    return(retlist)
}

#' Given a DE table with fold changes and p-values, show how 'significant' changes with changing cutoffs.
#'
#' Sometimes one might want to know how many genes are deemed significant while shifting the bars
#' which define significant.  This provides that metrics as a set of tables of numbers of
#' significant up/down genes when p-value is held constant, as well as number when fold-change is
#' held constant.
#'
#' @param table DE table to examine.
#' @param p_column Column in the DE table defining the changing p-value cutoff.
#' @param fc_column Column in the DE table defining the changing +/- log fold change.
#' @param bins Number of incremental changes in p-value/FC to examine.
#' @param constant_p When plotting changing FC, where should the p-value be held?
#' @param constant_fc When plotting changing p, where should the FC be held?
#' @return Plots and dataframes describing the changing definition of 'significant.'
#' @export
plot_num_siggenes <- function(table, p_column="limma_adjp", fc_column="limma_logfc",
                              bins=100, constant_p=0.05, constant_fc=0) {
    fc_column="limma_logfc"
    p_column="limma_adjp"
    bins = 100
    num_genes <- nrow(table)
    min_fc <- min(table[[fc_column]])
    neutral_fc <- 0.0
    max_fc <- max(table[[fc_column]])
    min_p <- 0.0
    max_p <- 1.0
    up_increments <- max_fc / bins
    down_increments <- min_fc / bins
    p_increments <- (max_p - min_p) / bins

    constant_up_fc <- constant_fc
    constant_down_fc <- constant_fc * -1.0
    start_up <- max_fc
    start_down <- min_fc
    start_p <- 0.00
    current_up_fc <- start_up
    current_down_fc <- start_down
    current_p <- start_p
    up_nums <- data.frame()
    down_nums <- data.frame()
    p_nums <- data.frame()
    for (inc in 1:bins) {
        current_up_fc <- current_up_fc - up_increments
        current_down_fc <- current_down_fc - down_increments
        current_p <- current_p + p_increments
        num_up <- sum(table[[fc_column]] >= current_up_fc & table[[p_column]] <= constant_p)
        num_down <- sum(table[[fc_column]] <= current_down_fc & table[[p_column]] <= constant_p)
        num_pup <- sum(table[[fc_column]] >= constant_up_fc & table[[p_column]] <= current_p)
        num_pdown <- sum(table[[fc_column]] <= constant_down_fc & table[[p_column]] <= current_p)
        up_nums <- rbind(up_nums, c(current_up_fc, num_up))
        down_nums <- rbind(down_nums, c(current_down_fc, num_down))
        p_nums <- rbind(p_nums, c(current_p, num_pup, num_pdown))
    }
    colnames(p_nums) <- c("p","up","down")
    colnames(up_nums) <- c("fc","num")
    colnames(down_nums) <- c("fc", "num")

    putative_up_inflection <- inflection::findiplist(x=as.matrix(up_nums[[1]]), y=as.matrix(up_nums[[2]]), 0)
    up_point_num <- putative_up_inflection[2,1]
    up_label <- paste0("At fc=", signif(up_nums[up_point_num, ][["fc"]], 4), " and p=", constant_p, ", ", up_nums[up_point_num, ][["num"]], " genes are de.")
    up_plot <- ggplot(data=up_nums, aes_string(x="fc", y="num")) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept=up_nums[[2]][[up_point_num]]) +
        ggplot2::geom_vline(xintercept=up_nums[[1]][[up_point_num]]) +
        ggplot2::geom_vline(xintercept=1.0, colour="red")

    putative_down_inflection <- inflection::findiplist(x=as.matrix(down_nums[[1]]), y=as.matrix(down_nums[[2]]), 0)
    down_point_num <- putative_down_inflection[1,2]
    down_plot <- ggplot(data=down_nums, aes_string(x="fc", y="num")) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept=down_nums[[2]][[down_point_num]]) +
        ggplot2::geom_vline(xintercept=down_nums[[1]][[down_point_num]]) +
        ggplot2::geom_vline(xintercept=-1.0, colour="red")

    putative_pup_inflection <- inflection::findiplist(x=as.matrix(p_nums[[1]]), y=as.matrix(p_nums[[2]]), 1)
    pup_point_num <- putative_pup_inflection[2,1]
    putative_pdown_inflection <- inflection::findiplist(x=as.matrix(p_nums[[1]]), y=as.matrix(p_nums[[2]]), 1)
    pdown_point_num <- putative_pdown_inflection[2,1]
    p_plot <- ggplot(data=p_nums) +
        ggplot2::geom_point(aes_string(x="p", y="up"), colour="darkred") +
        ggplot2::geom_point(aes_string(x="p", y="down"), colour="darkblue") +
        ggplot2::geom_vline(xintercept=0.05, colour="red") +
        ggplot2::geom_hline(yintercept=p_nums[[2]][[pup_point_num]], colour="darkred") +
        ggplot2::geom_vline(xintercept=p_nums[[1]][[pup_point_num]], colour="black") +
        ggplot2::geom_hline(yintercept=p_nums[[3]][[pdown_point_num]], colour="darkblue")

    retlist <- list(
        "up" = up_plot,
        "down" = down_plot,
        "p" = p_plot,
        "up_data" = up_nums,
        "down_data" = down_nums,
        "p_data" = p_nums)
    return(retlist)
}

#' Given the set of significant genes from combine_de_tables(), provide a view of how many are
#' significant up/down.
#'
#' These plots are pretty annoying, and I am certain that this function is not well written, but it
#' provides a series of bar plots which show the number of genes/contrast which are up and down
#' given a set of fold changes and p-value.
#'
#' @param combined  Result from combine_de_tables and/or extract_significant_genes().
#' @param fc_cutoffs  Choose 3 fold changes to define the queries.  0, 1, 2 mean greater/less than 0
#'     followed by 2 fold and 4 fold cutoffs.
#' @param fc_column  The column in the master-table to use for FC cutoffs.
#' @param p_type  Adjusted or not?
#' @param p  Chosen p-value cutoff.
#' @param z  Choose instead a z-score cutoff.
#' @param order  Choose a specific order for the plots.
#' @param maximum  Set a specific limit on the number of genes on the x-axis.
#' @param ...  More arguments are passed to arglist.
#' @export
significant_barplots <- function(combined, fc_cutoffs=c(0, 1, 2),
                                 fc_column="limma_logfc", p_type="adj",
                                 p=0.05, z=NULL, order=NULL, maximum=NULL, ...) {
    arglist <- list(...)
    sig_lists_up <- list(
        "limma" = list(),
        "edger" = list(),
        "deseq" = list())
    sig_lists_down <- list(
        "limma" = list(),
        "edger" = list(),
        "deseq" = list())
    plots <- list(
        "limma" = NULL,
        "edger" = NULL,
        "deseq" = NULL)
    tables_up <- list(
        "limma" = NULL,
        "edger" = NULL,
        "deseq" = NULL)
    tables_down <- list(
        "limma" = NULL,
        "edger" = NULL,
        "deseq" = NULL)
    table_length <- 0
    fc_names <- c()

    uplist <- list()
    downlist <- list()
    for (type in c("limma", "edger", "deseq")) {

        for (fc in fc_cutoffs) {
            ## This is a bit weird and circuituous
            ## The most common caller of this function is in fact extract_significant_genes
            fc_sig <- extract_significant_genes(combined, fc=fc,
                                                p=p, z=z, n=NULL, excel=FALSE,
                                                p_type=p_type, sig_bar=FALSE, ma=FALSE)
            table_length <- length(fc_sig[[type]][["ups"]])
            fc_name <- paste0("fc_", fc)
            fc_names <- append(fc_names, fc_name)

            for (tab in 1:table_length) { ## The table names are shared across methods and ups/downs
                table_name <- names(fc_sig[[type]][["ups"]])[tab]
                t_up <- nrow(fc_sig[[type]][["ups"]][[table_name]])
                t_down <- nrow(fc_sig[[type]][["downs"]][[table_name]])

                sig_lists_up[[type]][[fc_name]][[table_name]] <- t_up
                sig_lists_down[[type]][[fc_name]][[table_name]] <- t_down

            } ## End iterating through every table
        } ## End querying all fc cutoffs
        ## Now we need to collate the data and make the bars

        up_all <- list("limma" = numeric(),
                       "deseq" = numeric(),
                       "edger" = numeric())## The number of all genes FC > 0
        down_all <- list("limma" = numeric(),
                         "deseq" = numeric(),
                         "edger" = numeric())## The number of all genes FC < 0
        up_mid <- list("limma" = numeric(),
                       "deseq" = numeric(),
                       "edger" = numeric())## The number of genes 2<FC<4 (by default)
        down_mid <- list("limma" = numeric(),
                         "deseq" = numeric(),
                         "edger" = numeric())## The number of genes -2>FC>-4
        up_max <- list("limma" = numeric(),
                       "deseq" = numeric(),
                       "edger" = numeric())## The number of genes FC > 4
        down_max <- list("limma" = numeric(),
                         "deseq" = numeric(),
                         "edger" = numeric())  ## The number of genes FC < -4
        ##  The bar graph looks like
        ## ######### #### #  <-- Total width is the number of all >1FC genes
        ##         ^    ^------- Total >0FC - the set between 4FC and 2FC
        ##         |------------ Total >0FC - the smallest set >4FC

        papa_bear <- fc_names[[1]]  ## Because it is the largest grouping
        mama_bear <- fc_names[[2]]  ## The middle grouping
        baby_bear <- fc_names[[3]]  ## And the smallest grouping
        for (t in 1:table_length) {
            table_name <- names(sig_lists_up[[type]][[1]])[t]
            everything_up <- sig_lists_up[[type]][[papa_bear]][[table_name]] ## > 0 lfc
            mid_up <- sig_lists_up[[type]][[mama_bear]][[table_name]] ## > 1 lfc
            exclusive_up <- sig_lists_up[[type]][[baby_bear]][[table_name]] ## > 2 lfc
            ## Ah, I think the problem is that by calculating the numbers a,b,c
            ## It is stacking them and so I am getting a final bar of the sum of a,b,c
            up_all[[type]][[table_name]] <- everything_up
            up_mid[[type]][[table_name]] <- mid_up - exclusive_up
            up_max[[type]][[table_name]] <- exclusive_up
            up_all[[type]][[table_name]] <- up_all[[type]][[table_name]] - up_mid[[type]][[table_name]] - up_max[[type]][[table_name]]
            up_terminal <- up_all[[type]][[table_name]] + up_mid[[type]][[table_name]] + up_max[[type]][[table_name]]
            up_middle <- up_terminal - up_max[[type]][[table_name]]
            up_min <- up_terminal - up_mid[[type]][[table_name]]
            ## Now repeat for the set of down genes.
            everything_down <- sig_lists_down[[type]][[papa_bear]][[table_name]] ## > 0 lfc
            mid_down <- sig_lists_down[[type]][[mama_bear]][[table_name]] ## > 1 lfc
            exclusive_down <- sig_lists_down[[type]][[baby_bear]][[table_name]] ## > 2 lfc
            ## Ah, I think the problem is that by calculating the numbers a,b,c
            ## It is stacking them and so I am getting a final bar of the sum of a,b,c
            down_all[[type]][[table_name]] <- everything_down
            down_mid[[type]][[table_name]] <- mid_down - exclusive_down
            down_max[[type]][[table_name]] <- exclusive_down
            down_all[[type]][[table_name]] <- down_all[[type]][[table_name]] - down_mid[[type]][[table_name]] - down_max[[type]][[table_name]]
            down_terminal <- down_all[[type]][[table_name]] + down_mid[[type]][[table_name]] + down_max[[type]][[table_name]]
            down_middle <- down_terminal - down_max[[type]][[table_name]]
            down_min <- down_terminal - down_mid[[type]][[table_name]]
        } ## End for 1:table_length

        ## Prepare the tables for plotting.
        comparisons <- names(sig_lists_up[[type]][[1]])
        ## Once again, starting with only the up-stuff
        up <- cbind(comparisons, up_all[[type]], up_mid[[type]], up_max[[type]])
        up <- as.data.frame(up)
        colnames(up) <- c("comparisons","a_up_inner","b_up_middle","c_up_outer")
        uplist[[type]] <- up
        up <- reshape2::melt(up, id.var="comparisons")
        up[["comparisons"]] <- factor(up[["comparisons"]], levels=comparisons)
        up[["variable"]] <- factor(up[["variable"]],  levels=c("a_up_inner","b_up_middle","c_up_outer"))
        up[["value"]] <- as.numeric(up[["value"]])
        ## Repeat with the set of down materials
        down <- cbind(comparisons, down_all[[type]], down_mid[[type]], down_max[[type]])
        down <- as.data.frame(down)
        colnames(down) <- c("comparisons","a_down_inner","b_down_middle","c_down_outer")
        downlist[[type]] <- down
        colnames(down) <- c("comparisons","a_down_inner","b_down_middle","c_down_outer")
        down <- reshape2::melt(down, id.var="comparisons")
        down[["comparisons"]] <- factor(down[["comparisons"]], levels=comparisons)
        ##        down[["variable"]] <- factor(down[["variable"]],
        ##        levels=c("a_down_inner","b_down_middle","c_down_outer"))
        down[["variable"]] <- factor(down[["variable"]],
                                     levels=c("c_down_outer","b_down_middle","a_down_inner"))
        up[["variable"]] <- factor(up[["variable"]],
                                   levels=c("c_up_outer","b_up_middle","a_up_inner"))
        down[["value"]] <- as.numeric(down[["value"]]) * -1
        tables_up[[type]] <- up
        tables_down[[type]] <- down
        plots[[type]] <- plot_significant_bar(up, down, maximum=maximum)
        ##plots[[type]] <- plot_significant_bar(up, down, maximum=maximum, ...)
    } ## End iterating over the 3 types, limma/deseq/edger

    retlist <- list(
        "ups" = uplist,
        "downs" = downlist,
        "limma_up_table" = tables_up[["limma"]],
        "limma_down_table"= tables_down[["limma"]],
        "limma" = plots[["limma"]],
        "deseq_up_table" = tables_up[["deseq"]],
        "deseq_down_table"= tables_down[["deseq"]],
        "deseq" = plots[["deseq"]],
        "edger_up_table" = tables_up[["edger"]],
        "edger_down_table"= tables_down[["edger"]],
        "edger" = plots[["edger"]])
    return(retlist)
}

## EOF
