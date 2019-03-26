#' Make a MA plot of some limma output with pretty colors and shapes
#'
#' Yay pretty colors and shapes!
#'
#' @param pairwise The result from all_pairwise(), which should be changed to
#'   handle other invocations too.
#' @param type Type of table to use: deseq, edger, limma, basic.
#' @param table Result from edger to use, left alone it chooses the first.
#' @param logfc What logFC to use for the MA plot horizontal lines.
#' @param p_type Adjusted or raw pvalues?
#' @param p Cutoff to define 'significant' by p-value.
#' @param invert Invert the plot?
#' @param ... Extra arguments are passed to arglist.
#' @return a plot!
#' @seealso \code{\link{plot_ma_de}}
#' @examples
#' \dontrun{
#'  prettyplot <- edger_ma(all_aprwise) ## [sic, I'm witty! and can speel]
#' }
#' @export
extract_de_plots <- function(pairwise, type="edger", table=NULL, logfc=1,
                             p_type="adj", p=0.05, invert=FALSE, ...) {
  arglist <- list(...)
  table_source <- "all_pairwise"
  ## Possibilities include: all_pairwise, deseq_pairwise, limma_pairwise,
  ## edger_pairwise, basic_pairwise, combine_de_tables.
  if (!is.null(pairwise[["method"]])) {
    if (pairwise[["method"]] != type) {
      stop("The requested pairwise type and the provided input type do not match.")
    }
  }

  ## If the user did not ask for a specific table, assume the first one
  wanted_table <- NULL
  if (is.null(table)) {
    wanted_table <- 1
  } else {
    wanted_table <- table
  }

  ## if it is in fact all_pairwise, then there should be a set of
  ## slots 'limma', 'deseq', 'edger', 'basic' from which we can
  ## essentially convert the input by extracting the relevant type.
  if (class(pairwise)[1] == "all_pairwise") {
    table_source <- glue::glue("{type}_pairwise")
    pairwise <- pairwise[[type]]
  } else if (class(pairwise)[1] == "combined_de" |
             class(pairwise)[1] == "combined_table") {
    ## Then this came from combine_de...
    table_source <- "combined"
  } else if (!is.null(pairwise[["method"]])) {
    table_source <- glue::glue("{pairwise[['method']]}_pairwise")
  } else {
    stop("Unable to determine the source of this data.")
  }

  ## Depending on the source, choose appropriate column names.
  ## The expression column is the same across combined and _pairwise tables.
  ## The fc and p columns change depending on context.
  ## So does the set of possible tables.
  expr_col <- NULL
  fc_col <- NULL
  p_col <- NULL
  all_tables <- NULL
  if (table_source == "deseq_pairwise") {
    ## This column will need to be changed from base 10 to log scale.
    expr_col <- "baseMean"
    fc_col <- "logFC"  ## The most common
    if (p_type == "adj") {
      p_col <- "adj.P.Val"
    } else {
      p_col <- "P.Value"
    }
    all_tables <- pairwise[["all_tables"]]
  } else if (table_source == "edger_pairwise") {
    expr_col <- "logCPM"
    fc_col <- "logFC"  ## The most common
    if (p_type == "adj") {
      p_col <- "FDR"
    } else {
      p_col <- "PValue"
    }
    all_tables <- pairwise[["all_tables"]]
  } else if (table_source == "limma_pairwise") {
    expr_col <- "AveExpr"
    fc_col <- "logFC"  ## The most common
    if (p_type == "adj") {
      p_col <- "adj.P.Val"
    } else {
      p_col <- "P.Value"
    }
    all_tables <- pairwise[["all_tables"]]
  } else if (table_source == "basic_pairwise") {
    expr_col <- "numerator_median"
    fc_col <- "logFC"  ## The most common
    if (p_type == "adj") {
      p_col <- "adjp"
    } else {
      p_col <- "p"
    }
    all_tables <- pairwise[["all_tables"]]
  } else if (table_source == "ebseq_pairwise") {
    expr_col <- "ebseq_mean"
    fc_col <- "logFC"
    if (p_type == "adj") {
      p_col <- "ebseq_adjp"
    } else {
      p_col <- "ebseq_p"
    }
    all_tables <- pairwise[["all_tables"]]
  } else if (table_source == "combined") {
    if (type == "deseq") {
      expr_col <- "deseq_basemean"
    } else if (type == "edger") {
      expr_col <- "edger_logcpm"
    } else if (type == "limma") {
      expr_col <- "limma_ave"
    } else if (type == "basic") {
      expr_col <- "basic_nummed"
    } else if (type =="ebseq") {
      expr_col <- "ebseq_mean"
    }
    fc_col <- glue::glue("{type}_logfc")
    if (p_type == "adj") {
      p_col <- glue::glue("{type}_adjp")
    } else {
      p_col <- glue::glue("{type}_p")
    }
    all_tables <- pairwise[["data"]]
  } else {
    stop("Something went wrong, we should have only _pairwise and combined here.")
  }

  the_table <- NULL
  ## Now that we have the columns, figure out which table.
  if (class(all_tables) == "data.frame") {
    ## This came from the creation of combine_de_tables()
    the_table <- all_tables
  } else if (is.numeric(wanted_table)) {
    ## It is possible to just request the 1st, second, etc table
    ##the_table <- pairwise[["data"]][[table]]
    the_table <- all_tables[[wanted_table]]
  } else if (grepl(pattern="_vs_", x=wanted_table)) {
    ## The requested table might be a_vs_b, but sometimes a and b get flipped.
    ## Figure that out here and return the appropriate table.
    the_table <- wanted_table
    revname <- strsplit(x=the_table, split="_vs_")
    revname <- glue::glue("{revname[[1]][2]}_vs_{revname[[1]][1]}")
    possible_tables <- names(all_tables)
    if (!(the_table %in% possible_tables) & revname %in% possible_tables) {
      message("Trey you doofus, you reversed the name of the table.")
      the_table <- all_tables[[revname]]
    } else if (!(the_table %in% possible_tables) & !(revname %in% possible_tables)) {
      message("Unable to find the table in the set of possible tables.")
      message("The possible tables are: ", toString(possible_tables))
      stop()
    } else {
      the_table <- all_tables[[the_table]]
    }
  } else if (length(wanted_table) == 1) {
    ## One might request a name from the keepers list
    ## If so, figure that out here.
    table_parts <- pairwise[["keepers"]][[table]]
    if (is.null(table_parts)) {
      message("Unable to find the table in the set of possible tables.")
      message("The possible tables are: ", toString(possible_tables))
      stop()
    }
    the_table <- all_tables[[table]]
  } else if (length(wanted_table) == 2) {
    ## Perhaps one will ask for c(numerator, denominator)
    the_table <- glue::glue("{wanted_table[[1]]}_vs_{wanted_table[[2]]}")
    revname <- strsplit(x=the_table, split="_vs_")
    revname <- glue::glue("{revname[[1]][2]}_vs_{revname[[1]][1]}")
    possible_tables <- names(pairwise[["data"]])
    if (!(the_table %in% possible_tables) & revname %in% possible_tables) {
      message("Trey you doofus, you reversed the name of the table.")
      the_table <- all_tables[[revname]]
    } else if (!(the_table %in% possible_tables) & !(revname %in% possible_tables)) {
      stop("Unable to find the table in the set of possible tables.")
    } else {
      ##the_table <- pairwise[["data"]][[the_table]]
      the_tale <- all_tables[[the_table]]
    }
  } else {
    stop("Unable to discern the table requested.")
  }

  ## DESeq2 returns the median values as base 10, but we are using log2 (or log10?)
  if (type == "deseq") {
    the_table[["log_basemean"]] <- log2(x=the_table[[expr_col]] + 1.0)
    expr_col <- "log_basemean"
  }

  ma_material <- NULL
  vol_material <- NULL
  if (!is.null(the_table[[expr_col]]) &
      !is.null(the_table[[fc_col]]) &
      !is.null(the_table[[p_col]])) {
    ma_material <- plot_ma_de(
      table=the_table, expr_col=expr_col, fc_col=fc_col, p_col=p_col,
      logfc=logfc, p=p, invert=invert)
    ##...)
    vol_material <- plot_volcano_de(
      table=the_table, fc_col=fc_col, p_col=p_col,
      logfc=logfc, p=p)
    ##...)
  }

  retlist <- list(
    "ma" = ma_material,
    "volcano" = vol_material)
  return(retlist)
}

#' Perform a coefficient scatter plot of a limma/deseq/edger/basic table.
#'
#' Plot the gene abundances for two coefficients in a differential expression
#' comparison. By default, genes past 1.5 z scores from the mean are colored
#' red/green.
#'
#' @param output Result from the de_ family of functions, all_pairwise, or
#'   combine_de_tables().
#' @param toptable Chosen table to query for abundances.
#' @param type Query limma, deseq, edger, or basic outputs.
#' @param x The x-axis column to use, either a number of name.
#' @param y The y-axis column to use.
#' @param z Define the range of genes to color (FIXME: extend this to p-value
#'   and fold-change).
#' @param p Set a p-value cutoff for coloring the scatter plot (currently not
#'   supported).
#' @param lfc Set a fold-change cutoff for coloring points in the scatter plot
#'   (currently not supported.)
#' @param n Set a top-n fold-change for coloring the points in the scatter plot
#'   (this should work, actually).
#' @param loess Add a loess estimation (This is slow.)
#' @param alpha How see-through to make the dots.
#' @param color_low Color for the genes less than the mean.
#' @param color_high Color for the genes greater than the mean.
#' @param z_lines Add lines to show the z-score demarcations.
#' @param ... More arguments are passed to arglist.
#' @seealso \pkg{ggplot2}
#'  \code{\link{plot_linear_scatter}}
#' @examples
#' \dontrun{
#'  scatter_plot <- extract_coefficient_scatter(pairwise_output,
#'                                              type="deseq", x="uninfected", y="infected")
#' }
#' @export
extract_coefficient_scatter <- function(output, toptable=NULL, type="limma", x=1, y=2, z=1.5,
                                        p=NULL, lfc=NULL, n=NULL, loess=FALSE,
                                        alpha=0.4, color_low="#DD0000", z_lines=FALSE,
                                        color_high="#7B9F35", ...) {
  arglist <- list(...)
  ## This is an explicit test against all_pairwise() and reduces it to result from type.
  if (!is.null(output[[type]])) {
    output <- output[[type]]
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

  ## Extract the set of available names -- FIXME this should be standardized!!
  coefficients <- data.frame()
  thenames <- NULL
  if (type == "edger") {
    thenames <- names(output[["contrasts"]][["identities"]])
  } else if (type == "limma") {
    coefficients <- as.data.frame(output[["identity_comparisons"]][["coefficients"]])
    thenames <- colnames(coefficients)
  } else if (type == "deseq") {
    coefficients <- as.data.frame(output[["coefficients"]])
    thenames <- colnames(output[["coefficients"]])
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
  message("Actually comparing ", xname, " and ", yname, ".")

  ## Now extract the coefficent df
  if (type == "edger") {
    coefficient_df <- as.data.frame(output[["lrt"]][[1]][["coefficients"]])
    if (is.null(coefficient_df[[xname]]) || is.null(coefficient_df[[yname]])) {
      message("Did not find ", xname, " or ", yname, ".")
      return(NULL)
    }
    coefficient_df <- coefficient_df[, c(xname, yname)]
    coef_offset <- min(coefficient_df)
    ## This is dumb.
    coefficient_df <- coefficient_df + (coef_offset * -1.0)
  } else if (type == "limma") {
    coefficient_df <- as.data.frame(output[["pairwise_comparisons"]][["coefficients"]])
    if (is.null(coefficients[[x]]) || is.null(coefficients[[y]])) {
      message("Did not find ", x, " or ", y, ".")
      return(NULL)
    }
    coefficient_df <- coefficients[, c(x, y)]
  } else if (type == "deseq") {
    if (is.null(coefficients[[xname]]) || is.null(coefficients[[yname]])) {
      message("Did not find ", xname, " or ", yname, ".")
      return(NULL)
    }
    coefficient_df <- coefficients[, c(xname, yname)]
  } else if (type == "basic") {
    coefficient_df <- output[["medians"]]
    if (is.null(coefficients[[xname]]) || is.null(coefficients[[yname]])) {
      message("Did not find ", xname, " or ", yname, ".")
      return(NULL)
    }
    coefficient_df <- coefficient_df[, c(xname, yname)]
  }

  maxvalue <- max(coefficient_df) + 1.0
  minvalue <- min(coefficient_df) - 1.0
  plot <- sm(plot_linear_scatter(df=coefficient_df, loess=loess, gvis_filename=gvis_filename,
                                 gvis_trendline=gvis_trendline, first=xname, second=yname,
                                 tooltip_data=tooltip_data, base_url=base_url, alpha=alpha,
                                 pretty_colors=FALSE, color_low=color_low, color_high=color_high,
                                 p=p, lfc=lfc, n=n, z=z, z_lines=z_lines))
  plot[["scatter"]] <- plot[["scatter"]] +
    ggplot2::scale_x_continuous(limits=c(minvalue, maxvalue)) +
    ggplot2::scale_y_continuous(limits=c(minvalue, maxvalue))
  plot[["df"]] <- coefficient_df
  return(plot)
}

#' Create venn diagrams describing how well deseq/limma/edger agree.
#'
#' The sets of genes provided by limma and friends would ideally always agree,
#' but they do not. Use this to see out how much the (dis)agree.
#'
#' @param table Which table to query?
#' @param adjp Use adjusted p-values
#' @param p p-value cutoff, I forget what for right now.
#' @param lfc What fold-change cutoff to include?
#' @param ... More arguments are passed to arglist.
#' @return A list of venn plots
#' @seealso \pkg{venneuler} \pkg{Vennerable}
#' @examples
#' \dontrun{
#'  bunchovenns <- de_venn(pairwise_result)
#' }
#' @export
de_venn <- function(table, adjp=FALSE, p=0.05, lfc=0, ...) {
  arglist <- list(...)
  if (!is.null(table[["data"]])) {
    ## Then this is the result of combine_de
    retlist <- list()
    for (i in 1:length(names(table[["data"]]))) {
      a_table <- table[["data"]][[i]]
      retlist[[i]] <- de_venn(a_table, adjp=adjp, p=p, lfc=lfc, arglist)
    }
    return(retlist)
  }
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

  limma_sig <- sm(get_sig_genes(table, lfc=lfc,
                                column="limma_logfc", p_column=limma_p, p=p))
  edger_sig <- sm(get_sig_genes(table, lfc=lfc,
                                column="edger_logfc", p_column=edger_p, p=p))
  deseq_sig <- sm(get_sig_genes(table, lfc=lfc,
                                column="deseq_logfc", p_column=deseq_p, p=p))
  comp_up <- combine_tables(deseq_sig[["up_genes"]],
                            edger_sig[["up_genes"]],
                            limma_sig[["up_genes"]])
  comp_down <- combine_tables(deseq_sig[["down_genes"]],
                              edger_sig[["down_genes"]],
                              limma_sig[["down_genes"]])

  up_venn_lst <- list(
    "deseq" = comp_up[["deseq"]],
    "edger" = comp_up[["edger"]],
    "limma" = comp_up[["limma"]])
  down_venn_lst <- list(
    "deseq" = comp_down[["deseq"]],
    "edger" = comp_down[["edger"]],
    "limma" = comp_down[["limma"]])

  up_venn <- Vennerable::Venn(Sets=up_venn_lst)
  down_venn <- Vennerable::Venn(Sets=down_venn_lst)
  up_res <- Vennerable::plot(up_venn, doWeights=FALSE)
  up_venn_noweight <- grDevices::recordPlot()
  down_res <- Vennerable::plot(down_venn, doWeights=FALSE)
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

#' Given a DE table with p-values, plot them.
#'
#' Plot a multi-histogram containing (adjusted)p-values.
#'
#' @param combined Table to extract the values from.
#' @param type If provided, extract the {type}_p and {type}_adjp columns.
#' @param p_type Which type of pvalue to show (adjusted, raw, or all)?
#' @param columns Otherwise, extract whatever columns are provided.
#' @param ... Arguments passed through to the histogram plotter
#' @return Multihistogram of the result.
plot_de_pvals <- function(combined, type="limma", p_type="both", columns=NULL, ...) {
  if (is.null(type) & is.null(columns)) {
    stop("Some columns are required to extract p-values.")
  }
  if (is.null(columns)) {
    columns <- c(paste0(type, "_p"), paste0(type, "_adjp"))
  }
  plot_df <- combined[, columns]
  for (c in 1:ncol(plot_df)) {
    plot_df[[c]] <- as.numeric(plot_df[[c]])
  }

  if (p_type == "both") {
    p_stuff <- plot_multihistogram(plot_df, colors=c("darkred", "darkblue"), ...)
  } else if (p_type == "raw") {
    p_stuff <- plot_histogram(plot_df[[1]])
  } else {
    p_stuff <- plot_histogram(plot_df[[2]])
  }
  return(p_stuff)
}

#' Given a DE table with fold changes and p-values, show how 'significant'
#' changes with changing cutoffs.
#'
#' Sometimes one might want to know how many genes are deemed significant while
#' shifting the bars which define significant.  This provides that metrics as a
#' set of tables of numbers of significant up/down genes when p-value is held
#' constant, as well as number when fold-change is held constant.
#'
#' @param table DE table to examine.
#' @param methods List of methods to use when plotting.
#' @param bins Number of incremental changes in p-value/FC to examine.
#' @param constant_p When plotting changing FC, where should the p-value be held?
#' @param constant_fc When plotting changing p, where should the FC be held?
#' @return Plots and dataframes describing the changing definition of 'significant.'
#' @seealso \pkg{ggplot2}
#' @examples
#' \dontrun{
#'  crazy_sigplots <- plot_num_siggenes(pairwise_result)
#' }
#' @export
plot_num_siggenes <- function(table, methods=c("limma", "edger", "deseq", "ebseq"),
                              bins=100, constant_p=0.05, constant_fc=0) {
  min_fc <- -1
  max_fc <- 1
  lfc_columns <- c()
  p_columns <- c()
  kept_methods <- c()
  for (m in methods) {
    colname <- glue::glue("{m}_logfc")
    if (!is.null(table[[colname]])) {
      lfc_columns <- c(lfc_columns, colname)
      pcol <- glue::glue("{m}_adjp")
      kept_methods <- c(kept_methods, m)
      p_columns <- c(p_columns, pcol)
      test_fc <- min(table[[colname]])
      if (test_fc < min_fc) {
        min_fc <- test_fc
      }
      test_fc <- max(table[[colname]])
      if (test_fc > max_fc) {
        max_fc <- test_fc
      }
    }
  }
  num_genes <- nrow(table)
  neutral_fc <- 0.0
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
  pup_nums <- data.frame()
  pdown_nums <- data.frame()
  for (inc in 1:bins) {
    current_up_fc <- current_up_fc - up_increments
    current_down_fc <- current_down_fc - down_increments
    current_p <- current_p + p_increments
    up_nums_row <- c()
    down_nums_row <- c()
    pup_nums_row <- c()
    pdown_nums_row <- c()
    for (c in 1:length(lfc_columns)) {
      lfc_col <- lfc_columns[c]
      p_col <- p_columns[c]
      num_up <- sum(table[[lfc_col]] >= current_up_fc & table[[p_col]] <= constant_p)
      num_down <- sum(table[[lfc_col]] <= current_down_fc & table[[p_col]] <= constant_p)
      num_pup <- sum(table[[lfc_col]] >= constant_up_fc & table[[p_col]] <= current_p)
      num_pdown <- sum(table[[lfc_col]] <= constant_down_fc & table[[p_col]] <= current_p)
      up_nums_row <- c(up_nums_row, num_up)
      down_nums_row <- c(down_nums_row, num_down)
      pup_nums_row <- c(pup_nums_row, num_pup)
      pdown_nums_row <- c(pdown_nums_row, num_pdown)
    }

    up_nums <- rbind(up_nums, c(current_up_fc, up_nums_row))
    down_nums <- rbind(down_nums, c(current_down_fc, down_nums_row))
    pup_nums <- rbind(pup_nums, c(current_p, pup_nums_row))
    pdown_nums <- rbind(pdown_nums, c(current_p, pdown_nums_row))
  }
  colnames(pup_nums) <- c("p", kept_methods)
  pup_nums <- reshape2::melt(pup_nums, id.vars="p")
  colnames(pup_nums) <- c("p", "method", "value")
  colnames(pdown_nums) <- c("p", kept_methods)
  pdown_nums <- reshape2::melt(pdown_nums, id.vars="p")
  colnames(pdown_nums) <- c("p", "method", "value")
  colnames(up_nums) <- c("fc", kept_methods)
  up_nums <- reshape2::melt(up_nums, id.vars="fc")
  colnames(up_nums) <- c("fc", "method", "value")
  colnames(down_nums) <- c("fc", kept_methods)
  down_nums <- reshape2::melt(down_nums, id.vars="fc")
  colnames(down_nums) <- c("fc", "method", "value")

  up_plot <- ggplot(data=up_nums, aes_string(x="fc", y="value", color="method")) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(palette="Set1") +
    ggplot2::scale_color_brewer(palette="Set1") +
    ggplot2::geom_vline(xintercept=1.0, colour="red") +
    ggplot2::theme_bw(base_size=base_size)

  down_plot <- ggplot(data=down_nums, aes_string(x="fc", y="value", color="method")) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(palette="Set1") +
    ggplot2::scale_color_brewer(palette="Set1") +
    ggplot2::geom_vline(xintercept=-1.0, colour="red") +
    ggplot2::theme_bw(base_size=base_size)

  pup_plot <- ggplot(data=pup_nums, aes_string(x="p", y="value", color="method")) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(palette="Set1") +
    ggplot2::scale_color_brewer(palette="Set1") +
    ggplot2::geom_vline(xintercept=0.05, colour="red") +
    ggplot2::theme_bw(base_size=base_size)

  pdown_plot <- ggplot(data=pdown_nums, aes_string(x="p", y="value", color="method")) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(palette="Set1") +
    ggplot2::scale_color_brewer(palette="Set1") +
    ggplot2::geom_vline(xintercept=0.05, colour="red") +
    ggplot2::theme_bw(base_size=base_size)

  retlist <- list(
    "up" = up_plot,
    "down" = down_plot,
    "pup" = pup_plot,
    "pdown" = pdown_plot,
    "up_data" = up_nums,
    "down_data" = down_nums,
    "pup_data" = pup_nums,
    "pdown_data" = pdown_nums)
  return(retlist)
}

#' Plot the rank order of the data in two tables against each other.
#'
#' Steve Christensen has some neat plots showing the relationship between two
#' tables.  I though they were super-cool, so I co-opted the idea in this
#' function.
#'
#' @param first First table of values.
#' @param second Second table of values, if null it will use the first.
#' @param first_type Assuming this is from all_pairwise(), use this method.
#' @param second_type Ibid.
#' @param first_table Again, assuming all_pairwise(), use this to choose the
#'   table to extract.
#' @param alpha How see-through to make the dots?
#' @param second_table Ibid.
#' @param first_column What column to use to rank-order from the first table?
#' @param second_column What column to use to rank-order from the second table?
#' @param first_p_col Use this column for pretty colors from the first table.
#' @param second_p_col Use this column for pretty colors from the second table.
#' @param p_limit A p-value limit for coloring dots.
#' @param both_color If both columns are 'significant', use this color.
#' @param first_color If only the first column is 'significant', this color.
#' @param second_color If the second column is 'significant', this color.
#' @param no_color If neither column is 'significant', then this color.
#' @return a list with a plot and a couple summary statistics.
#' @export
rank_order_scatter <- function(first, second=NULL, first_type="limma",
                               second_type="limma", first_table=1, alpha=0.5,
                               second_table=2, first_column="logFC",
                               second_column="logFC", first_p_col="adj.P.Val",
                               second_p_col="adj.P.Val", p_limit=0.05,
                               both_color="red", first_color="green",
                               second_color="blue", no_color="black") {
  if (is.null(second)) {
    second <- first
  }
  if (!is.null(first[[first_type]])) {
    first <- first[[first_type]][["all_tables"]]
  }
  if (!is.null(second[[second_type]])) {
    second <- second[[second_type]][["all_tables"]]
  }
  table1 <- first[[first_table]]
  table2 <- second[[second_table]]
  merged <- merge(table1, table2, by="row.names")
  rownames(merged) <- merged[["Row.names"]]
  merged <- merged[, -1]

  if (first_column == second_column) {
    c1 <- glue::glue("{first_column}.x")
    c2 <- glue::glue("{first_column}.y")
  } else {
    c1 <- first_column
    c2 <- second_column
  }
  c1_idx <- order(merged[[c1]])
  merged <- merged[c1_idx, ]
  merged[["label"]] <- rownames(merged)
  c1_idx <- order(merged[[c1]])
  c2_idx <- order(merged[[c2]])
  merged[["x"]] <- as.numeric(c1_idx)
  merged[["y"]] <- as.numeric(c2_idx)

  merged[["state"]] <- "neither"
  if (first_p_col == second_p_col) {
    p1 <- glue::glue("{first_p_col}.x")
    p2 <- glue::glue("{first_p_col}.y")
  } else {
    p1 <- first_p_col
    p2 <- second_p_col
  }
  both_idx <- merged[[p1]] < p_limit & merged[[p2]] < p_limit
  merged[both_idx, "state"] <- "both"
  p1_idx <- merged[[p1]] < p_limit & merged[[p2]] >= p_limit
  merged[p1_idx, "state"] <- "first"
  p2_idx <- merged[[p2]] < p_limit & merged[[p1]] >= p_limit
  merged[p2_idx, "state"] <- "second"
  merged[["state"]] <- as.factor(merged[["state"]])

  first_table_colname <- glue::glue(
    "Table: {first_table}, Type: {first_type}, column: {first_column}")
  second_table_colname <- glue::glue(
    "Table: {second_table}, Type: {second_type}, column: {second_column}")

  plt <- ggplot(data=merged,
                aes_string(color="state", fill="state",
                           x="x", y="y", label="label")) +
    ggplot2::geom_point(size=1, alpha=alpha) +
    ggplot2::scale_color_manual(name="state",
                                values=c("both"=both_color,
                                         "first"=first_color,
                                         "second"=second_color,
                                         "neither"=no_color)) +
    ggplot2::geom_smooth(method="loess", color="lightblue") +
    ggplot2::ylab(glue::glue("Rank order of {second_table_colname}")) +
    ggplot2::xlab(glue::glue("Rank order of {first_table_colname}")) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=base_size, colour="black"))

  model_test <- try(lm(formula=y ~ x, data=merged), silent=TRUE)
  model_summary <- summary(model_test)
  cor <- cor.test(merged[[c1]], merged[[c2]], method="pearson")
  retlist <- list(
    "plot" = plt,
    "model" = model_test,
    "summary" = model_summary,
    "correlation" = cor)

  return(retlist)
}

#' Given the set of significant genes from combine_de_tables(), provide a view
#' of how many are significant up/down.
#'
#' These plots are pretty annoying, and I am certain that this function is not
#' well written, but it provides a series of bar plots which show the number of
#' genes/contrast which are up and down given a set of fold changes and
#' p-value.
#'
#' @param combined Result from combine_de_tables and/or extract_significant_genes().
#' @param lfc_cutoffs Choose 3 fold changes to define the queries.  0, 1, 2
#'   mean greater/less than 0 followed by 2 fold and 4 fold cutoffs.
#' @param invert Reverse the order of contrasts for readability?
#' @param p Chosen p-value cutoff.
#' @param z Choose instead a z-score cutoff.
#' @param p_type Adjusted or not?
#' @param according_to limma, deseq, edger, basic, or all of the above.
#' @param order Choose a specific order for the plots.
#' @param maximum Set a specific limit on the number of genes on the x-axis.
#' @param ... More arguments are passed to arglist.
#' @return list containing the significance bar plots and some information to
#'   hopefully help interpret them.
#' @seealso \pkg{ggplot2}
#' @examples
#' \dontrun{
#'  ## Damn I wish I were smrt enough to make this elegant, but I cannot.
#'  barplots <- significant_barplots(combined_result)
#' }
#' @export
significant_barplots <- function(combined, lfc_cutoffs=c(0, 1, 2), invert=FALSE,
                                 p=0.05, z=NULL, p_type="adj",
                                 according_to="all", order=NULL, maximum=NULL, ...) {
  arglist <- list(...)
  sig_lists_up <- list(
    "limma" = list(),
    "edger" = list(),
    "deseq" = list(),
    "ebseq" = list(),
    "basic" = list())
  sig_lists_down <- list(
    "limma" = list(),
    "edger" = list(),
    "deseq" = list(),
    "ebseq" = list(),
    "basic" = list())
  plots <- list(
    "limma" = NULL,
    "edger" = NULL,
    "deseq" = NULL,
    "ebseq" = NULL,
    "basic" = NULL)
  tables_up <- list(
    "limma" = NULL,
    "edger" = NULL,
    "deseq" = NULL,
    "ebseq" = NULL,
    "basic" = NULL)
  tables_down <- list(
    "limma" = NULL,
    "edger" = NULL,
    "deseq" = NULL,
    "ebseq" = NULL,
    "basic" = NULL)
  table_length <- 0
  fc_names <- c()

  uplist <- list()
  downlist <- list()

  types <- according_to
  if (according_to[[1]] == "all") {
    types <- c("limma", "edger", "deseq", "ebseq", "basic")
  }
  ##else if (according_to == c("limma", "edger", "deseq", "basic")) {
  ##  types <- c("limma", "edger", "deseq")
  ##}

  for (type in types) {
    test_column <- glue::glue("{type}_logfc")
    if (! test_column %in% colnames(combined[["data"]][[1]])) {
      message("We do not have the ", test_column, " in the data, skipping ", type, ".")
      next
    }
    for (fc in lfc_cutoffs) {
      ## This is a bit weird and circuituous
      ## The most common caller of this function is in fact extract_significant_genes
      fc_sig <- sm(extract_significant_genes(combined, lfc=fc, according_to=according_to,
                                             p=p, z=z, n=NULL, excel=FALSE,
                                             p_type=p_type, sig_bar=FALSE, ma=FALSE))
      table_length <- length(fc_sig[[type]][["ups"]])
      fc_name <- glue::glue("fc_{fc}")
      fc_names <- append(fc_names, fc_name)

      for (tab in 1:table_length) {
        ## The table names are shared across methods and ups/downs
        table_names <- names(fc_sig[[type]][["ups"]])
        if (isTRUE(invert)) {
          table_names <- rev(table_names)
        }
        table_name <- table_names[tab]
        t_up <- nrow(fc_sig[[type]][["ups"]][[table_name]])
        t_down <- nrow(fc_sig[[type]][["downs"]][[table_name]])

        sig_lists_up[[type]][[fc_name]][[table_name]] <- t_up
        sig_lists_down[[type]][[fc_name]][[table_name]] <- t_down
      } ## End iterating through every table
    } ## End querying all fc cutoffs
    ## Now we need to collate the data and make the bars

    up_all <- list("limma" = numeric(),
                   "deseq" = numeric(),
                   "edger" = numeric(),
                   "ebseq" = numeric(),
                   "basic" = numeric()
                   ) ## The number of all genes FC > 0
    down_all <- list("limma" = numeric(),
                     "deseq" = numeric(),
                     "edger" = numeric(),
                     "ebseq" = numeric(),
                     "basic" = numeric()
                     ) ## The number of all genes FC < 0
    up_mid <- list("limma" = numeric(),
                   "deseq" = numeric(),
                   "edger" = numeric(),
                   "ebseq" = numeric(),
                   "basic" = numeric()
                   ) ## The number of genes 2<FC<4 (by default)
    down_mid <- list("limma" = numeric(),
                     "deseq" = numeric(),
                     "edger" = numeric(),
                     "ebseq" = numeric(),
                     "basic" = numeric()
                     ) ## The number of genes -2>FC>-4
    up_max <- list("limma" = numeric(),
                   "deseq" = numeric(),
                   "edger" = numeric(),
                   "ebseq" = numeric(),
                   "basic" = numeric()
                   ) ## The number of genes FC > 4
    down_max <- list("limma" = numeric(),
                     "deseq" = numeric(),
                     "edger" = numeric(),
                     "ebseq" = numeric(),
                     "basic" = numeric()
                     ) ## The number of genes FC < -4
    ##  The bar graph looks like
    ## ######### #### #  <-- Total width is the number of all >1FC genes
    ##         ^    ^------- Total >0FC - the set between 4FC and 2FC
    ##         |------------ Total >0FC - the smallest set >4FC

    papa_bear <- fc_names[[1]]  ## Because it is the largest grouping
    mama_bear <- fc_names[[2]]  ## The middle grouping
    baby_bear <- fc_names[[3]]  ## And the smallest grouping
    for (t in 1:table_length) {
      table_names <- names(sig_lists_up[[type]][[1]])
      table_name <- table_names[t]
      ##table_names <- names(sig_lists_up[[type]][[1]])[t]
      everything_up <- sig_lists_up[[type]][[papa_bear]][[table_name]] ## > 0 lfc
      mid_up <- sig_lists_up[[type]][[mama_bear]][[table_name]] ## > 1 lfc
      exclusive_up <- sig_lists_up[[type]][[baby_bear]][[table_name]] ## > 2 lfc
      ## Ah, I think the problem is that by calculating the numbers a,b,c
      ## It is stacking them and so I am getting a final bar of the sum of a,b,c
      up_all[[type]][[table_name]] <- everything_up
      up_mid[[type]][[table_name]] <- mid_up - exclusive_up
      up_max[[type]][[table_name]] <- exclusive_up
      up_all[[type]][[table_name]] <- up_all[[type]][[table_name]] -
        up_mid[[type]][[table_name]] -
        up_max[[type]][[table_name]]
      up_terminal <- up_all[[type]][[table_name]] +
        up_mid[[type]][[table_name]] +
        up_max[[type]][[table_name]]
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
      down_all[[type]][[table_name]] <- down_all[[type]][[table_name]] -
        down_mid[[type]][[table_name]] -
        down_max[[type]][[table_name]]
      down_terminal <- down_all[[type]][[table_name]] +
        down_mid[[type]][[table_name]] +
        down_max[[type]][[table_name]]
      down_middle <- down_terminal - down_max[[type]][[table_name]]
      down_min <- down_terminal - down_mid[[type]][[table_name]]
    } ## End for 1:table_length

    ## Prepare the tables for plotting.
    comparisons <- names(sig_lists_up[[type]][[1]])
    ## Once again, starting with only the up-stuff
    up <- cbind(comparisons, up_all[[type]], up_mid[[type]], up_max[[type]])
    up <- as.data.frame(up)
    colnames(up) <- c("comparisons", "a_up_inner", "b_up_middle", "c_up_outer")
    uplist[[type]] <- up
    up <- suppressWarnings(reshape2::melt(up, id.var="comparisons"))
    up[["comparisons"]] <- factor(up[["comparisons"]], levels=comparisons)
    up[["variable"]] <- factor(up[["variable"]],
                               levels=c("a_up_inner", "b_up_middle", "c_up_outer"))
    up[["value"]] <- as.numeric(up[["value"]])
    ## Repeat with the set of down materials
    down <- cbind(comparisons, down_all[[type]], down_mid[[type]], down_max[[type]])
    down <- as.data.frame(down)
    colnames(down) <- c("comparisons", "a_down_inner", "b_down_middle", "c_down_outer")
    downlist[[type]] <- down
    colnames(down) <- c("comparisons", "a_down_inner", "b_down_middle", "c_down_outer")
    down <- suppressWarnings(reshape2::melt(down, id.var="comparisons"))
    down[["comparisons"]] <- factor(down[["comparisons"]], levels=comparisons)
    ##        down[["variable"]] <- factor(down[["variable"]],
    ##        levels=c("a_down_inner","b_down_middle","c_down_outer"))
    down[["variable"]] <- factor(down[["variable"]],
                                 levels=c("c_down_outer", "b_down_middle", "a_down_inner"))
    up[["variable"]] <- factor(up[["variable"]],
                               levels=c("c_up_outer", "b_up_middle", "a_up_inner"))
    down[["value"]] <- as.numeric(down[["value"]]) * -1
    tables_up[[type]] <- up
    tables_down[[type]] <- down
    plots[[type]] <- plot_significant_bar(up, down, maximum=maximum,
                                          ...)
    ## plots[[type]] <- plot_significant_bar(up, down, maximum=maximum) #, ...)
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
    "edger" = plots[["edger"]],
    "ebseq_up_table" = tables_up[["ebseq"]],
    "ebseq_down_table"= tables_down[["ebseq"]],
    "ebseq" = plots[["ebseq"]],
    "basic_up_table" = tables_up[["basic"]],
    "basic_down_table"= tables_down[["basic"]],
    "basic" = plots[["basic"]]
  )
  return(retlist)
}

## EOF
