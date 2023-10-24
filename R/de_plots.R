## de_plots.r: A series of plots which are in theory DE method agnostic.

#' Make a MA plot of some limma output with pretty colors and shapes.
#'
#' Yay pretty colors and shapes!
#' This function should be reworked following my rewrite of combine_de_tables().
#' It is certainly possible to make the logic here much simpler now.
#'
#' @param pairwise The result from all_pairwise(), which should be changed to
#'  handle other invocations too.
#' @param combined Result from one of the combine_de_table functions.
#' @param type Type of table to use: deseq, edger, limma, basic.
#' @param invert Invert the plot?
#' @param invert_colors vector of new colors.
#' @param numerator Use this factor as the numerator.
#' @param denominator Use this factor as the denominator.
#' @param alpha Use this transparency.
#' @param z z-score cutoff for coefficient significance.
#' @param logfc What logFC to use for the MA plot horizontal lines.
#' @param pval Cutoff to define 'significant' by p-value.
#' @param found_table Result from edger to use, left alone it chooses the first.
#' @param p_type Adjusted or raw pvalues?
#' @param color_high Color to use for the 'high' genes.
#' @param color_low Color to use for the 'low' genes.
#' @param loess Add a loess estimator to the coefficient plot?
#' @param z_lines Add the z-score lines?
#' @param label Label this number of top-diff genes.
#' @param label_columns Use this column for labelling genes.
#' @return a plot!
#' @seealso [plot_ma_de()] [plot_volcano_de()]
#' @examples
#' \dontrun{
#'  prettyplot <- edger_ma(all_aprwise) ## [sic, I'm witty! and can speel]
#' }
#' @export
extract_de_plots <- function(pairwise, combined = NULL, type = NULL,
                             invert = FALSE, invert_colors = c(),
                             numerator = NULL, denominator = NULL, alpha = 0.4, z = 1.5, n = NULL,
                             logfc = 1.0, pval = 0.05, found_table = NULL, p_type = "adj",
                             color_high = NULL, color_low = NULL, loess = FALSE,
                             z_lines = FALSE, label = 10, label_column = "hgncsymbol") {

  if (is.null(type)) {
    if (grepl(pattern = "pairwise", x = class(pairwise)[1])) {
      type <- gsub(x = class(pairwise)[1], pattern = "_pairwise$", replacement = "")
    }
    if (is.null(found_table)) {
      message("No table was provided, choosing the first.")
      found_table <- names(pairwise[["all_tables"]])[1]
    }
  }
  if (is.null(combined)) {
    combined = pairwise
  }
  if (is.null(numerator) && is.null(denominator)) {
    message("No numerator nor denominator was provided, arbitrarily choosing the first.")
    num_den_names <- get_num_den(names(pairwise[["all_tables"]])[1])
    numerator <- num_den_names[["numerator"]]
    denominator <- num_den_names[["denominator"]]
  }
  source_info <- get_plot_columns(combined, type, p_type = p_type)
  input <- source_info[["the_table"]]
  expr_col <- source_info[["expr_col"]]
  fc_col <- source_info[["fc_col"]]
  p_col <- source_info[["p_col"]]
  invert <- source_info[["invert"]]

  coef_result <- try(extract_coefficient_scatter(
    pairwise, type = type, x = denominator, y = numerator,
    z = z, n = n, loess = loess, alpha = alpha / 2.0,
    color_low = color_low, color_high = color_high,
    z_lines = z_lines))

  ma_material <- NULL
  vol_material <- NULL
  if (!is.null(input[[expr_col]]) &&
        !is.null(input[[fc_col]]) &&
         !is.null(input[[p_col]])) {
    ma_material <- plot_ma_condition_de(
      input = input, table_name = found_table,
      expr_col = expr_col, fc_col = fc_col, p_col = p_col,
      logfc = logfc, pval = pval, invert = invert,
      color_high = color_high, color_low = color_low,
      label = label, label_column = label_column)
    vol_material <- plot_volcano_condition_de(
      input = input, table_name = found_table,
      fc_col = fc_col, p_col = p_col,
      color_high = color_high, color_low = color_low,
      invert = invert, logfc = logfc, pval = pval,
      label = label, label_column = label_column)
  }

  retlist <- list(
    "coef" = coef_result,
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
#'  combine_de_tables().
#' @param toptable Chosen table to query for abundances.
#' @param type Query limma, deseq, edger, or basic outputs.
#' @param x The x-axis column to use, either a number of name.
#' @param y The y-axis column to use.
#' @param z Define the range of genes to color (FIXME: extend this to p-value
#'  and fold-change).
#' @param logfc Set a fold-change cutoff for coloring points in the scatter plot
#'  (currently not supported.)
#' @param n Set a top-n fold-change for coloring the points in the scatter plot
#'  (this should work, actually).
#' @param z_lines Add lines to show the z-score demarcations.
#' @param loess Add a loess estimation (This is slow.)
#' @param alpha How see-through to make the dots.
#' @param color_low Color for the genes less than the mean.
#' @param color_high Color for the genes greater than the mean.
#' @seealso [plot_linear_scatter()]
#' @examples
#' \dontrun{
#'  expt <- create_expt(metadata = "some_metadata.xlsx", gene_info = annotations)
#'  pairwise_output <- all_pairwise(expt)
#'  scatter_plot <- extract_coefficient_scatter(pairwise_output,
#'                                              type = "deseq", x = "uninfected", y = "infected")
#' }
#' @export
extract_coefficient_scatter <- function(output, toptable = NULL, type = "limma",
                                        x = 1, y = 2, z = 1.5, logfc = NULL, n = NULL,
                                        z_lines = FALSE, loess = FALSE, alpha = 0.4,
                                        color_low = "#DD0000", color_high = "#7B9F35") {
  ## This is an explicit test against all_pairwise() and reduces it to result from type.
  if (!is.null(output[[type]])) {
    output <- output[[type]]
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
  } else if (type == "ebseq") {
    thenames <- names(output[["conditions_table"]])
  } else {
    stop("I do not know what type you wish to query.")
  }

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
  } else if (type == "ebseq") {
    tables <- names(output[["all_tables"]])
    verted <- glue("{xname}_vs_{yname}")
    inverted <- glue("{yname}_vs_{xname}")
    coefficient_df <- data.frame()
    if (verted %in% tables) {
      table_idx <- verted == tables
      table <- output[["all_tables"]][[tables[table_idx]]]
      coefficient_df <- table[, c("ebseq_c1mean", "ebseq_c2mean")]
    } else if (inverted %in% tables) {
      table_idx <- inverted == tables
      table <- output[["all_tables"]][[tables[table_idx]]]
      coefficient_df <- table[, c("ebseq_c2mean", "ebseq_c1mean")]
    } else {
      stop("Did not find the table for ebseq.")
    }
    colnames(coefficient_df) <- c(xname, yname)
    coefficient_df[[1]] <- log2(coefficient_df[[1]])
    coefficient_df[[2]] <- log2(coefficient_df[[2]])
  } else if (type == "basic") {
    coefficient_df <- output[["medians"]]
    if (is.null(coefficient_df[[xname]]) || is.null(coefficient_df[[yname]])) {
      message("Did not find ", xname, " or ", yname, ".")
      return(NULL)
    }
    coefficient_df <- coefficient_df[, c(xname, yname)]
  }
  maxvalue <- max(coefficient_df) + 1.0
  minvalue <- min(coefficient_df) - 1.0
  plot <- plot_linear_scatter(df = coefficient_df, loess = loess, first = xname, second = yname,
                              alpha = alpha, pretty_colors = FALSE,
                              color_low = color_low, color_high = color_high,
                              n = n, z = z, logfc = logfc, z_lines = z_lines)
  plot[["scatter"]] <- plot[["scatter"]] +
    ggplot2::scale_x_continuous(limits = c(minvalue, maxvalue)) +
    ggplot2::scale_y_continuous(limits = c(minvalue, maxvalue))
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
#' @seealso [Vennerable] [get_sig_genes()]
#' @examples
#' \dontrun{
#'  bunchovenns <- de_venn(pairwise_result)
#' }
#' @export
de_venn <- function(table, adjp = FALSE, p = 0.05, lfc = 0, ...) {
  arglist <- list(...)
  if (!is.null(table[["data"]])) {
    ## Then this is the result of combine_de
    retlist <- list()
    for (i in seq_along(names(table[["data"]]))) {
      a_table <- table[["data"]][[i]]
      retlist[[i]] <- de_venn(a_table, adjp = adjp, p = p, lfc = lfc, arglist)
    }
    return(retlist)
  }
  limma_p <- "limma_p"
  deseq_p <- "deseq_p"
  edger_p <- "edger_p"
  if (isTRUE(adjp)) {
    limma_p <- "limma_adjp"
    deseq_p <- "deseq_adjp"
    edger_p <- "edger_adjp"
  }

  sig_data <- list()
  up_venn_lst <- list()
  down_venn_lst <- list()
  if (!is.null(table[["limma_logfc"]])) {
    sig_data[["limma"]] <- sm(get_sig_genes(
      table, lfc = lfc,
      column = "limma_logfc", p_column = limma_p, p = p))
    up_venn_lst[["limma"]] <- rownames(sig_data[["limma"]][["up_genes"]])
    down_venn_lst[["limma"]] <- rownames(sig_data[["limma"]][["down_genes"]])

  }
  if (!is.null(table[["deseq_logfc"]])) {
    sig_data[["deseq"]] <- sm(get_sig_genes(
      table, lfc = lfc,
      column = "deseq_logfc", p_column = deseq_p, p = p))
    up_venn_lst[["deseq"]] <- rownames(sig_data[["deseq"]][["up_genes"]])
    down_venn_lst[["deseq"]] <- rownames(sig_data[["deseq"]][["down_genes"]])
  }
  if (!is.null(table[["edger_logfc"]])) {
    sig_data[["edger"]] <- sm(get_sig_genes(
      table, lfc = lfc,
      column = "edger_logfc", p_column = edger_p, p = p))
    up_venn_lst[["edger"]] <- rownames(sig_data[["edger"]][["up_genes"]])
    down_venn_lst[["edger"]] <- rownames(sig_data[["edger"]][["down_genes"]])
  }

  up_venn <- Vennerable::Venn(Sets = up_venn_lst)
  down_venn <- Vennerable::Venn(Sets = down_venn_lst)
  tmp_file <- tmpmd5file(pattern = "venn", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  up_res <- Vennerable::plot(up_venn, doWeights = FALSE)
  up_venn_noweight <- grDevices::recordPlot()
  dev.off()
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  down_res <- Vennerable::plot(down_venn, doWeights = FALSE)
  down_venn_noweight <- grDevices::recordPlot()
  dev.off()
  removed <- file.remove(tmp_file)
  removed <- unlink(dirname(tmp_file))

  retlist <- list(
    "up_venn" = up_venn,
    "up_noweight" = up_venn_noweight,
    "down_venn" = down_venn,
    "down_noweight" = down_venn_noweight)
  return(retlist)
}

#' A small rat's nest of if statements intended to figure out what columns
#' are wanted to plot a MA/Volcano from any one of a diverse set of possible
#' input types.
#'
#' I split this function away from the main body of extract_de_plots()
#' so that I can come back to it and strip it down to something a bit
#' more legible.  Eventually I want to dispatch this logic off to
#' separate functions depending on the class of the input.
#'
#' This function should die in a fire.
#'
#' @param data Data structure in which to hunt columns/data.
#' @param type Type of method used to make the data.
#' @param p_type Use adjusted p-values?
get_plot_columns <- function(data, type, p_type = "adj") {
  ret <- list(
    "p_col" = "P.Val",
    "fc_col" = "logFC",
    "expr_col" = "baseMean",
    "wanted_table" = NULL,
    "invert" = FALSE,
    "input" = NULL,
    "table_source" = "data")

  ## Possibilities include: all_pairwise, deseq_pairwise, limma_pairwise,
  ## edger_pairwise, basic_pairwise, combine_de_tables.
  if (!is.null(data[["method"]])) {
    if (data[["method"]] != type) {
      stop("The requested pairwise type and the provided input type do not match.")
    }
  }

  ## If the user did not ask for a specific table, assume the first one
  wanted_table <- NULL
  found_table <- NULL
  if (is.null(found_table)) {
    wanted_table <- 1
  } else {
    wanted_table <- found_table
  }

  if ("combined_de" %in% class(data)) {
    wanted_tablename <- data[["kept"]][[wanted_table]]
    actual_tablenames <- data[["keepers"]][[wanted_table]]
    actual_numerator <- actual_tablenames[[1]]
    actual_denominator <- actual_tablenames[[2]]
    actual_tablename <- paste0(actual_numerator, "_vs_", actual_denominator)
    if (actual_tablename != wanted_tablename) {
      ret[["invert"]] <- TRUE
    }
    wanted_table <- wanted_tablename
    input <- data[["input"]]
  } else if ("combined_table" %in% class(data)) {
    table_source <- "combined_table"
  } else if (class(data)[1] == "all_pairwise") {
    ## if it is in fact all_pairwise, then there should be a set of
    ## slots 'limma', 'deseq', 'edger', 'basic' from which we can
    ## essentially convert the input by extracting the relevant type.
    table_source <- glue("{type}_pairwise")
    data <- data[[type]]
  } else if (!is.null(data[["method"]])) {
    table_source <- glue("{data[['method']]}_pairwise")
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
  if (table_source == "combined_table") {
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
    fc_col <- glue("{type}_logfc")
    if (p_type == "adj") {
      p_col <- glue("{type}_adjp")
    } else {
      p_col <- glue("{type}_p")
    }
    all_tables <- NULL
  } else if (table_source == "deseq_pairwise") {
    ## This column will need to be changed from base 10 to log scale.
    expr_col <- "baseMean"
    fc_col <- "logFC"  ## The most common
    if (p_type == "adj") {
      p_col <- "adj.P.Val"
    } else {
      p_col <- "P.Value"
    }
    all_tables <- data[["all_tables"]]
  } else if (table_source == "edger_pairwise") {
    expr_col <- "logCPM"
    fc_col <- "logFC"  ## The most common
    if (p_type == "adj") {
      p_col <- "FDR"
    } else {
      p_col <- "PValue"
    }
    all_tables <- data[["all_tables"]]
  } else if (table_source == "limma_pairwise") {
    expr_col <- "AveExpr"
    fc_col <- "logFC"  ## The most common
    if (p_type == "adj") {
      p_col <- "adj.P.Val"
    } else {
      p_col <- "P.Value"
    }
    all_tables <- data[["all_tables"]]
  } else if (table_source == "basic_pairwise") {
    ## basic_pairwise() may have columns which are 'numerator_median' or 'numerator_mean'.
    expr_col <- "numerator"
    if (!is.null(data[["all_tables"]][[1]][["numerator_mean"]])) {
      expr_col <- "numerator_mean"
    }
    fc_col <- "logFC"  ## The most common
    if (p_type == "adj") {
      p_col <- "adjp"
    } else {
      p_col <- "p"
    }
    all_tables <- data[["all_tables"]]
  } else if (table_source == "ebseq_pairwise") {
    expr_col <- "ebseq_mean"
    fc_col <- "logFC"
    if (p_type == "adj") {
      p_col <- "ebseq_adjp"
    } else {
      p_col <- "ebseq_ppde"
    }
    all_tables <- data[["all_tables"]]
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
    fc_col <- glue("{type}_logfc")
    if (p_type == "adj") {
      p_col <- glue("{type}_adjp")
    } else {
      p_col <- glue("{type}_p")
    }
    all_tables <- data[["data"]]
  } else {
    stop("Something went wrong, we should have only _pairwise and combined here.")
  }

  possible_tables <- names(all_tables)
  the_table <- NULL
  ## Now that we have the columns, figure out which table.
  if (is.null(all_tables)) {
    the_table <- data[["data"]]
  } else if ("data.frame" %in% class(all_tables)) {
    ## This came from the creation of combine_de_tables()
    the_table <- all_tables
  } else if (is.numeric(wanted_table)) {
    ## It is possible to just request the 1st, second, etc table
    ##the_table <- pairwise[["data"]][[table]]
    the_table <- all_tables[[wanted_table]]
  } else if (grepl(pattern = "_vs_", x = wanted_table)) {
    ## The requested table might be a_vs_b, but sometimes a and b get flipped.
    ## Figure that out here and return the appropriate table.
    the_table <- wanted_table
    revname <- strsplit(x = the_table, split = "_vs_")
    revname <- glue("{revname[[1]][2]}_vs_{revname[[1]][1]}")
    if (!(the_table %in% possible_tables) && revname %in% possible_tables) {
      mesg("Trey you doofus, you reversed the name of the table.")
      the_table <- all_tables[[revname]]
    } else if (!(the_table %in% possible_tables) && !(revname %in% possible_tables)) {
      message("Unable to find the table in the set of possible tables.")
      message("The possible tables are: ", toString(possible_tables))
      stop()
    } else {
      the_table <- all_tables[[the_table]]
    }
  } else if (length(wanted_table) == 1) {
    ## One might request a name from the keepers list
    ## If so, figure that out here.
    table_parts <- data[["keepers"]][[table]]
    if (is.null(table_parts)) {
      message("Unable to find the table in the set of possible tables.")
      message("The possible tables are: ", toString(possible_tables))
      stop()
    }
    the_table <- all_tables[[table]]
  } else if (length(wanted_table) == 2) {
    ## Perhaps one will ask for c(numerator, denominator)
    the_table <- glue("{wanted_table[[1]]}_vs_{wanted_table[[2]]}")
    revname <- strsplit(x = the_table, split = "_vs_")
    revname <- glue("{revname[[1]][2]}_vs_{revname[[1]][1]}")
    possible_tables <- names(data[["data"]])
    if (!(the_table %in% possible_tables) && revname %in% possible_tables) {
      mesg("Trey you doofus, you reversed the name of the table.")
      the_table <- all_tables[[revname]]
    } else if (!(the_table %in% possible_tables) && !(revname %in% possible_tables)) {
      stop("Unable to find the table in the set of possible tables.")
    } else {
      ##the_table <- pairwise[["data"]][[the_table]]
      the_table <- all_tables[[the_table]]
    }
  } else {
    stop("Unable to discern the table requested.")
  }

  ## DESeq2 returns the median values as base 10, but we are using log2 (or log10?)
  if (type == "deseq") {
    the_table[["log_basemean"]] <- log2(x = the_table[[expr_col]] + 1.0)
    expr_col <- "log_basemean"
  }

  ## Check that the wanted table is numeric
  if (is.numeric(wanted_table)) {
    wanted_table <- names(data[["all_tables"]])[wanted_table]
  }

  ret[["the_table"]] <- the_table
  ret[["expr_col"]] <- expr_col
  ret[["fc_col"]] <- fc_col
  ret[["p_col"]] <- p_col
  ret[["wanted_table"]] <- wanted_table
  return(ret)
}

#' Given a DE table with p-values, plot them.
#'
#' Plot a multi-histogram containing (adjusted)p-values.
#'
#' The assumption of this plot is that the adjustment will significantly
#' decrease the representation of genes in the 'highly significant' range of
#' p-values.  However, it is hoped that it will not utterly remove them.
#'
#' @param combined_data Table to extract the values from.
#' @param type If provided, extract the {type}_p and {type}_adjp columns.
#' @param p_type Which type of pvalue to show (adjusted, raw, or all)?
#' @param columns Otherwise, extract whatever columns are provided.
#' @param ... Arguments passed through to the histogram plotter
#' @return Multihistogram of the result.
#' @seealso [plot_histogram()]
plot_de_pvals <- function(combined_data, type = "limma", p_type = "both", columns = NULL, ...) {
  if (is.null(type) & is.null(columns)) {
    stop("Some columns are required to extract p-values.")
  }
  if (p_type == "all") {
    columns <- c(paste0(type, "_p"), paste0(type, "_adjp"), paste0(type, "_adjp_ihw"))
  } else if (p_type == "both") {
    columns <- c(paste0(type, "_p"), paste0(type, "_adjp"))
  } else if (is.null(columns)) {
    columns <- c(paste0(type, "_p"), paste0(type, "_adjp"))
  } else {
    columns <- c(paste0(type, "_p"), paste0(type, "_adjp"), paste0(type, "_adjp_", tolower(p_type)))
    p_type <- "all"
  }
  plot_df <- combined_data[, columns]
  for (c in seq_len(ncol(plot_df))) {
    plot_df[[c]] <- as.numeric(plot_df[[c]])
  }

  if (p_type == "all") {
    p_stuff <- plot_multihistogram(plot_df, colors = c("darkred", "darkblue", "forestgreen"),
                                   ...)
  } else if (p_type == "both") {
    p_stuff <- plot_multihistogram(plot_df, colors = c("darkred", "darkblue"), ...)
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
#' @seealso [ggplot2]
#' @examples
#' \dontrun{
#'  pairwise_result <- all_pairwise(expt)
#'  crazy_sigplots <- plot_num_siggenes(pairwise_result)
#' }
#' @export
plot_num_siggenes <- function(table, methods = c("limma", "edger", "deseq", "ebseq"),
                              bins = 100, constant_p = 0.05, constant_fc = 0) {
  min_fc <- -1
  max_fc <- 1
  lfc_columns <- c()
  p_columns <- c()
  kept_methods <- c()
  for (m in methods) {
    colname <- glue("{m}_logfc")
    if (!is.null(table[[colname]])) {
      lfc_columns <- c(lfc_columns, colname)
      pcol <- glue("{m}_adjp")
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
  for (inc in seq_len(bins)) {
    current_up_fc <- current_up_fc - up_increments
    current_down_fc <- current_down_fc - down_increments
    current_p <- current_p + p_increments
    up_nums_row <- c()
    down_nums_row <- c()
    pup_nums_row <- c()
    pdown_nums_row <- c()
    for (c in seq_along(lfc_columns)) {
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
  pup_nums <- reshape2::melt(pup_nums, id.vars = "p")
  colnames(pup_nums) <- c("p", "method", "value")
  colnames(pdown_nums) <- c("p", kept_methods)
  pdown_nums <- reshape2::melt(pdown_nums, id.vars = "p")
  colnames(pdown_nums) <- c("p", "method", "value")
  colnames(up_nums) <- c("fc", kept_methods)
  up_nums <- reshape2::melt(up_nums, id.vars = "fc")
  colnames(up_nums) <- c("fc", "method", "value")
  colnames(down_nums) <- c("fc", kept_methods)
  down_nums <- reshape2::melt(down_nums, id.vars = "fc")
  colnames(down_nums) <- c("fc", "method", "value")

  up_plot <- ggplot(data = up_nums,
                    aes(x = .data[["fc"]], y = .data[["value"]],
                        color = .data[["method"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::geom_vline(xintercept = 1.0, colour = "red") +
    ggplot2::theme_bw(base_size = base_size)

  down_plot <- ggplot(data = down_nums,
                      aes(x = .data[["fc"]], y = .data[["value"]],
                          color = .data[["method"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::geom_vline(xintercept=-1.0, colour = "red") +
    ggplot2::theme_bw(base_size = base_size)

  pup_plot <- ggplot(data = pup_nums,
                     aes(x = .data[["p"]], y = .data[["value"]],
                         color = .data[["method"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::geom_vline(xintercept = 0.05, colour = "red") +
    ggplot2::theme_bw(base_size = base_size)

  pdown_plot <- ggplot(data = pdown_nums,
                       aes(x = .data[["p"]], y = .data[["value"]],
                           color = .data[["method"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::geom_vline(xintercept = 0.05, colour = "red") +
    ggplot2::theme_bw(base_size = base_size)

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

#' Make a pretty MA plot from one of limma, deseq, edger, or basic.
#'
#' Because I can never remember, the following from wikipedia: "An MA plot is an
#' application of a Bland-Altman plot for visual representation of two channel
#' DNA microarray gene expression data which has been transformed onto the M
#' (log ratios) and A (mean average) scale."
#'
#' @param table Df of linear-modelling, normalized counts by sample-type,
#' @param expr_col Column showing the average expression across genes.
#' @param fc_col Column showing the logFC for each gene.
#' @param p_col Column containing the relevant p values.
#' @param pval Name of the pvalue column to use for cutoffs.
#' @param alpha How transparent to make the dots.
#' @param logfc Fold change cutoff.
#' @param label_numbers Show how many genes were 'significant', 'up', and 'down'?
#' @param size How big are the dots?
#' @param shapes Provide different shapes for up/down/etc?
#' @param invert Invert the ma plot?
#' @param label Label the top/bottom n logFC values?
#' @param ... More options for you
#' @return ggplot2 MA scatter plot.  This is defined as the rowmeans of the
#'  normalized counts by type across all sample types on the x axis, and the
#'  log fold change between conditions on the y-axis. Dots are colored
#'  depending on if they are 'significant.'  This will make a fun clicky
#'  googleVis graph if requested.
#' @seealso [limma_pairwise()] [deseq_pairwise()] [edger_pairwise()] [basic_pairwise()]
#' @examples
#'  \dontrun{
#'   plot_ma(voomed_data, table)
#'   ## Currently this assumes that a variant of toptable was used which
#'   ## gives adjusted p-values.  This is not always the case and I should
#'   ## check for that, but I have not yet.
#'  }
#' @export
plot_ma_de <- function(table, expr_col = "logCPM", fc_col = "logFC", p_col = "qvalue",
                       pval = 0.05, alpha = 0.4, logfc = 1.0, label_numbers = TRUE,
                       size = 2, shapes = TRUE, invert = FALSE, label = NULL,
                       label_column = "hgncsymbol", ...) {
  ## Set up the data frame which will describe the plot
  arglist <- list(...)
  ## I like dark blue and dark red for significant and insignificant genes respectively.
  ## Others may disagree and change that with sig_color, insig_color.
  sig_color <- "darkred"
  if (!is.null(arglist[["sig_color"]])) {
    sig_color <- arglist[["sig_color"]]
  }
  insig_color <- "darkblue"
  if (!is.null(arglist[["insig_color"]])) {
    insig_color <- arglist[["insig_color"]]
  }
  ## A recent request was to color gene families within these plots.
  ## Below there is a short function,  recolor_points() which handles this.
  ## The following 2 arguments are used for that.
  ## That function should work for other things like volcano or scatter plots.
  family <- NULL
  if (!is.null(arglist[["family"]])) {
    family <- arglist[["family"]]
  }
  family_color <- "red"
  if (!is.null(arglist[["family_color"]])) {
    family_color <- arglist[["family_color"]]
  }

  ## The data frame used for these MA plots needs to include a few aspects of the state
  ## Including: average expression (the y axis), log-fold change, p-value, a
  ## boolean of the p-value state, and a factor of the state which will be
  ## counted and provide some information on the side of the plot. One might
  ## note that I am pre-filling in this data frame with 4 blank entries. This is
  ## to make absolutely certain that ggplot will not re-order my damn
  ## categories.
  df <- data.frame(
    "avg" = c(0, 0, 0),
    "logfc" = c(0, 0, 0),
    "pval" = c(0, 0, 0),
    "pcut" = c(FALSE, FALSE, FALSE),
    "state" = c("a_upsig", "b_downsig", "c_insig"), stringsAsFactors = TRUE)

  ## Get rid of rows which will be annoying.
  ## If somehow a list got into the data table, this will fail, lets fix that now.
  tmp_table <- table
  for (c in seq_len(ncol(tmp_table))) {
    tmp_table[[c]] <- as.character(table[[c]])
  }
  rows_without_na <- complete.cases(tmp_table)
  rm(tmp_table)
  table <- table[rows_without_na, ]

  ## Extract the information of interest from my original table
  newdf <- data.frame("avg" = table[[expr_col]],
                      "logfc" = table[[fc_col]],
                      "pval" = table[[p_col]])
  if (!is.null(table[[label_column]])) {
    newdf[["label"]] <- table[[label_column]]
    rownames(newdf) <- make.names(rownames(table), unique = TRUE)
  } else {
    newdf[["label"]] <- rownames(table)
    rownames(newdf) <- rownames(table)
  }
  if (isTRUE(invert)) {
    newdf[["logfc"]] <- newdf[["logfc"]] * -1.0
  }
  ## Check if the data is on a log or base10 scale, if the latter, then convert it.
  if (max(newdf[["avg"]]) > 1000) {
    newdf[["avg"]] <- log(newdf[["avg"]])
  }

  ## Set up the state of significant/insiginificant vs. p-value and/or fold-change.
  newdf[["pval"]] <- as.numeric(format(newdf[["pval"]], scientific = FALSE))
  newdf[["pcut"]] <- newdf[["pval"]] <= pval
  newdf[["state"]] <- ifelse(newdf[["pval"]] > pval, "c_insig",
                             ifelse(newdf[["pval"]] <= pval &
                                      newdf[["logfc"]] >= logfc, "a_upsig",
                                    ifelse(newdf[["pval"]] <= pval &
                                             newdf[["logfc"]] <= (-1.0 * logfc),
                                           "b_downsig", "c_insig")))
  newdf[["state"]] <- as.factor(newdf[["state"]])
  df <- rbind(df, newdf)
  rm(newdf)

  ## Subtract one from each value because I filled in a fake value of each category to start.
  num_downsig <- sum(df[["state"]] == "b_downsig") - 1
  num_insig <- sum(df[["state"]] == "c_insig") - 1
  num_upsig <- sum(df[["state"]] == "a_upsig") - 1

  ## Make double-certain that my states are only factors or numbers where necessary.
  df[["avg"]] <- as.numeric(df[[1]])
  df[["logfc"]] <- as.numeric(df[[2]])
  df[["pval"]] <- as.numeric(df[[3]])
  df[["pcut"]] <- as.factor(df[[4]])
  df[["state"]] <- as.factor(df[[5]])
  df[["label"]] <- rownames(df)

  ## Set up the labels for the legend by significance.
  ## 4 states, 4 shapes -- these happen to be the 4 best shapes in R because they may be filled.
  ## shape 24 is the up arrow, 25 the down arrow, 21 the circle.
  state_shapes <- 21
  if (isTRUE(state_shapes)) {
    state_shapes <- c(24, 25, 21)
    names(state_shapes) <- c("a_upsig", "b_downsig", "c_insig")
  } else {
    state_shapes <- c(21, 21, 21)
    names(state_shapes) <- c("a_upsig", "b_downsig", "c_insig")
  }

  ## make the plot!
  plt <- ggplot(data = df,
                ## I am setting x, y, fill color, outline color, and the shape.
                aes(x = .data[["avg"]],
                    y = .data[["logfc"]],
                    label = .data[["label"]],
                    fill = as.factor(.data[["pcut"]]),
                    colour = as.factor(.data[["pcut"]]),
                    shape = as.factor(.data[["state"]]))) +
    ggplot2::geom_hline(yintercept = c((logfc * -1.0), logfc),
                        color = "red", size=(size / 3)) +
    ggplot2::geom_point(stat = "identity", size = size, alpha = alpha)
  if (isTRUE(label_numbers)) {
    plt <- plt +
      ## The following scale_shape_manual() sets the labels of the legend on the right side.
      ggplot2::scale_shape_manual(name = "State", values = state_shapes,
                                  labels = c(
                                    glue("Up Sig.: {num_upsig}"),
                                    glue("Down Sig.: {num_downsig}"),
                                    glue("Insig.: {num_insig}")),
                                  guide = ggplot2::guide_legend(override.aes = aes(size = 3,
                                                                                   fill = "grey")))
  } else {
    plt <- plt +
      ggplot2::scale_shape_manual(name = "State", values = state_shapes,
                                  guide = "none")
  }

  plt <- plt +
    ## Set the colors of the significant/insignificant points.
    ggplot2::scale_fill_manual(name = "as.factor(pcut)",
                               values = c("FALSE"=insig_color, "TRUE"=sig_color),
                               guide = "none") +
    ggplot2::scale_color_manual(name = "as.factor(pcut)",
                                values = c("FALSE"=insig_color, "TRUE"=sig_color),
                                guide = "none") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black")) +
    ggplot2::xlab("Average log2(Counts)") +
    ggplot2::ylab("log2(fold change)")

  ## Recolor a family of genes if requested.
  if (!is.null(family)) {
    plt <- recolor_points(plt, df, family, color = family_color)
  }

  if (!is.null(label)) {
    reordered_idx <- order(df[["logfc"]])
    reordered <- df[reordered_idx, ]
    top <- head(reordered, n = label)
    bottom <- tail(reordered, n = label)
    df_subset <- rbind(top, bottom)
    plt <- plt +
      ggrepel::geom_text_repel(
        data = df_subset,
        aes(label = .data[["label"]], x = .data[["avg"]], y = .data[["logfc"]]),
        colour = "black", box.padding = ggplot2::unit(0.5, "lines"),
        point.padding = ggplot2::unit(1.6, "lines"),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")))
  }

  ## Return the plot, some numbers, and the data frame used to make the plot so
  ## that I may check my work.
  retlist <- list(
    "num_upsig" = num_upsig,
    "num_downsig" = num_downsig,
    "num_insig" = num_insig,
    "plot" = plt,
    "df" = df)
  return(retlist)
}

plot_ma_condition_de <- function(input, table_name, expr_col = "logCPM",
                                 fc_col = "logFC", p_col = "qvalue",
                                 color_high = "red", color_low = "blue",
                                 pval = 0.05, alpha = 0.4, logfc = 1.0, label_numbers = TRUE,
                                 size = 2, shapes = TRUE, invert = FALSE,
                                 label = 10, label_column = "hgncsymbol", ...) {
  ## Set up the data frame which will describe the plot

  ## Example caller:
  ## ma_material <- plot_ma_condition_de(
  ##     pairwise, table = the_table, expr_col = expr_col, fc_col = fc_col, p_col = p_col,
  ##     logfc = logfc, pval = pval, invert = invert,

  arglist <- list(...)

  ## A recent request was to color gene families within these plots.
  ## Below there is a short function,  recolor_points() which handles this.
  ## The following 2 arguments are used for that.
  ## That function should work for other things like volcano or scatter plots.
  family <- NULL
  if (!is.null(arglist[["family"]])) {
    family <- arglist[["family"]]
  }
  family_color <- "red"
  if (!is.null(arglist[["family_color"]])) {
    family_color <- arglist[["family_color"]]
  }

  ## The data frame used for these MA plots needs to include a few aspects of the state
  ## Including: average expression (the y axis), log-fold change, p-value, a
  ## boolean of the p-value state, and a factor of the state which will be
  ## counted and provide some information on the side of the plot. One might
  ## note that I am pre-filling in this data frame with 4 blank entries. This is
  ## to make absolutely certain that ggplot will not re-order my damn
  ## categories.
  df <- data.frame(
    "avg" = c(0, 0, 0),
    "logfc" = c(0, 0, 0),
    "pval" = c(0, 0, 0),
    "pcut" = c(FALSE, FALSE, FALSE),
    "state" = c("a_upsig", "b_downsig", "c_insig"),
    "label" = c("", "", ""),
    stringsAsFactors = TRUE)

  ## Extract the information of interest from my original table
  newdf <- data.frame("avg" = input[[expr_col]],
                      "logfc" = input[[fc_col]],
                      "pval" = input[[p_col]])
  if (is.null(label_column)) {
    newdf[["label"]] <- rownames(input)
  } else if (label_column %in% colnames(input)) {
    newdf[["label"]] <- input[[label_column]]
  } else {
    message("The column: ", label_column, " is not in the data, using rownames.")
    newdf[["label"]] <- rownames(input)
  }
  rownames(newdf) <- rownames(input)
  rows_without_na <- complete.cases(newdf)
  newdf <- newdf[rows_without_na, ]

  if (isTRUE(invert)) {
    newdf[["logfc"]] <- newdf[["logfc"]] * -1.0
  }
  ## Check if the data is on a log or base10 scale, if the latter, then convert it.
  if (max(newdf[["avg"]]) > 1000) {
    newdf[["avg"]] <- log(newdf[["avg"]])
  }

  ## Set up the state of significant/insiginificant vs. p-value and/or fold-change.
  newdf[["pval"]] <- as.numeric(format(newdf[["pval"]], scientific = FALSE))
  newdf[["pcut"]] <- newdf[["pval"]] <= pval
  newdf[["state"]] <- "c_insig"
  newdf[["state"]] <- ifelse(newdf[["pval"]] > pval, "c_insig",
                             ifelse(newdf[["pval"]] <= pval &
                                      newdf[["logfc"]] >= logfc, "a_upsig",
                                    ifelse(newdf[["pval"]] <= pval &
                                             newdf[["logfc"]] <= (-1.0 * logfc),
                                           "b_downsig", "c_insig")))
  newdf[["state"]] <- as.factor(newdf[["state"]])
  df <- rbind(df, newdf)
  rm(newdf)

  ## Subtract one from each value because I filled in a fake value of each category to start.
  num_downsig <- sum(df[["state"]] == "b_downsig") - 1
  num_insig <- sum(df[["state"]] == "c_insig") - 1
  num_upsig <- sum(df[["state"]] == "a_upsig") - 1

  ## Make double-certain that my states are only factors or numbers where necessary.
  df[["avg"]] <- as.numeric(df[["avg"]])
  df[["logfc"]] <- as.numeric(df[["logfc"]])
  df[["pval"]] <- as.numeric(df[["pval"]])
  df[["pcut"]] <- as.factor(df[["pcut"]])
  df[["state"]] <- as.factor(df[["state"]])

  ## Set up the labels for the legend by significance.
  ## 4 states, 4 shapes -- these happen to be the 4 best shapes in R because they may be filled.
  ## shape 24 is the up arrow, 25 the down arrow, 21 the circle.
  state_shapes <- 21
  if (isTRUE(state_shapes)) {
    state_shapes <- c(24, 25, 21)
    names(state_shapes) <- c("a_upsig", "b_downsig", "c_insig")
  } else {
    state_shapes <- c(21, 21, 21)
    names(state_shapes) <- c("a_upsig", "b_downsig", "c_insig")
  }

  ## I am not sure why something is setting color high/low to NULL.
  if (is.null(color_high)) {
    color_high <- "red"
  }
  if (is.null(color_low)) {
    color_low <- "blue"
  }
  plot_colors <- c("#555555", color_high, color_low)
  names(plot_colors) <- c("c_insig", "a_upsig", "b_downsig")

  ## make the plot!
  plt <- ggplot(data = df,
                ## I am setting x, y, fill color, outline color, and the shape.
                aes(x = .data[["avg"]],
                    y = .data[["logfc"]],
                    label = .data[["label"]],
                    fill = as.factor(.data[["state"]]),
                    colour = as.factor(.data[["state"]]),
                    shape = as.factor(.data[["state"]]))) +
    ggplot2::geom_hline(yintercept = c((logfc * -1.0), logfc),
                        color = "red", size=(size / 3)) +
    ggplot2::geom_point(stat = "identity", size = size, alpha = alpha)
  if (isTRUE(label_numbers)) {
    plt <- plt +
      ## The following scale_shape_manual() sets the labels of the legend on the right side.
      ggplot2::scale_shape_manual(name = "State", values = state_shapes,
                                  labels = c(
                                    glue("Up Sig.: {num_upsig}"),
                                    glue("Down Sig.: {num_downsig}"),
                                    glue("Insig.: {num_insig}")),
                                  guide = ggplot2::guide_legend(override.aes = aes(size = 3,
                                                                                   fill = "grey")))
  } else {
    plt <- plt +
      ggplot2::scale_shape_manual(name = "State", values = state_shapes,
                                  guide = "none")
  }

  plt <- plt +
    ## Set the colors of the significant/insignificant points.
    ggplot2::scale_fill_manual(name = "state",
                               values = plot_colors,
                               guide = "none") +
    ggplot2::scale_color_manual(name = "state",
                                values = plot_colors,
                                guide = "none") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black")) +
    ggplot2::xlab("Average log2(Counts)") +
    ggplot2::ylab("log2(fold change)")

  ## Recolor a family of genes if requested.
  if (!is.null(family)) {
    plt <- recolor_points(plt, df, family, color = family_color)
  }

  if (!is.null(label)) {
    reordered_idx <- order(df[["logfc"]])
    reordered <- df[reordered_idx, ]
    ## I am now taking 1/2 of the number of desired labeled genes on each side.
    ## I think this more accurately fits the spirit of this idea.
    top <- head(reordered, n = (ceiling(label / 2)))
    bottom <- tail(reordered, n = (ceiling(label / 2)))
    df_subset <- rbind(top, bottom)
    plt <- plt +
      ggrepel::geom_text_repel(
        data = df_subset,
        aes(label = .data[["label"]], x = .data[["avg"]], y = .data[["logfc"]]),
        colour = "black", box.padding = ggplot2::unit(0.5, "lines"),
        point.padding = ggplot2::unit(1.6, "lines"),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")))
  }

  ## Return the plot, some numbers, and the data frame used to make the plot so
  ## that I may check my work.
  retlist <- list(
    "num_upsig" = num_upsig,
    "num_downsig" = num_downsig,
    "num_insig" = num_insig,
    "plot" = plt,
    "df" = df)
  return(retlist)
}

plot_sankey_de <- function(de_table, lfc = 1.0, p = 0.05,
                           lfc_column = "deseq_logfc", p_column = "deseq_adjp") {
  de_table <- t_cf_clinical_table_sva$data$outcome
  de_table[["start"]] <- "all"
  de_table[["lfc"]] <- "up"
  down_idx <- de_table[[lfc_column]] < 0
  de_table[down_idx, "lfc"] <- "down"
  up_idx <- de_table[[lfc_column]] > lfc
  de_table[up_idx, "lfc"] <- "sigup"
  down_idx <- de_table[[lfc_column]] < -1.0 * lfc
  de_table[down_idx, "lfc"] <- "sigdown"

  de_table[["adjusted_p"]] <- "insignificant"
  sig_idx <- de_table[[p_column]] <= p
  de_table[sig_idx, "adjusted_p"] <- "significant"

  meta <- de_table[, c("start", "lfc", "adjusted_p")]
  meta[["lfc"]] <- factor(meta[["lfc"]], levels = c("sigup", "up", "down", "sigdown"))
  meta[["adjusted_p"]] <- factor(meta[["adjusted_p"]], levels = c("significant", "insignificant"))

  #  color_choices <-
  test <- plot_meta_sankey(meta, drill_down = FALSE, factors = c("lfc", "adjusted_p"))
  return(test[["ggplot"]])
}

#' Make a pretty Volcano plot!
#'
#' Volcano plots and MA plots provide quick an easy methods to view the set of
#' (in)significantly differentially expressed genes.  In the case of a volcano
#' plot, it places the -log10 of the p-value estimate on the y-axis and the
#' fold-change between conditions on the x-axis.  Here is a neat snippet from
#' wikipedia: "The concept of volcano plot can be generalized to other
#' applications, where the x-axis is related to a measure of the strength of a
#' statistical signal, and y-axis is related to a measure of the statistical
#' significance of the signal."
#'
#' @param table Dataframe from limma's toptable which includes log(fold change) and an
#'  adjusted p-value.
#' @param alpha How transparent to make the dots.
#' @param color_by By p-value something else?
#' @param color_list List of colors for significance.
#' @param fc_col Which column contains the fc data?
#' @param fc_name Name of the fold-change to put on the plot.
#' @param line_color What color for the significance lines?
#' @param line_position Put the significance lines above or below the dots?
#' @param logfc Cutoff defining the minimum/maximum fold change for
#'  interesting.
#' @param p_col Which column contains the p-value data?
#' @param p_name Name of the p-value to put on the plot.
#' @param p Cutoff defining significant from not.
#' @param shapes_by_state Add fun shapes for the various significance states?
#' @param minimum_p If a pvalue is lower than this, then set it to
#'  this, thus artificially limiting the y-scale of a volcano plot.
#'  This is only valid if one thinks that the pvalues are artificially
#'  low and that is messing with the interpretation of the data.
#' @param size How big are the dots?
#' @param invert Flip the x-axis?
#' @param label Label the top/bottom n logFC values?
#' @param label_column Use this column of annotations for labels instead of rownames?
#' @param ... I love parameters!
#' @return Ggplot2 volcano scatter plot.  This is defined as the -log10(p-value)
#'   with respect to log(fold change).  The cutoff values are delineated with
#'   lines and mark the boundaries between 'significant' and not.  This will
#'   make a fun clicky googleVis graph if requested.
#' @seealso [all_pairwise()]
#' @examples
#' \dontrun{
#'  plot_volcano_de(table)
#'  ## Currently this assumes that a variant of toptable was used which
#'  ## gives adjusted p-values.  This is not always the case and I should
#'  ## check for that, but I have not yet.
#' }
#' @export
plot_volcano_de <- function(table, alpha = 0.5, color_by = "p",
                            color_list = c("FALSE" = "darkblue", "TRUE" = "darkred"),
                            fc_col = "logFC", fc_name = "log2 fold change",
                            line_color = "black", line_position = "bottom", logfc = 1.0,
                            p_col = "adj.P.Val", p_name = "-log10 p-value", p = 0.05,
                            shapes_by_state = FALSE, minimum_p = NULL,
                            size = 2, invert = FALSE, label = NULL,
                            label_column = "hgncsymbol", ...) {
  low_vert_line <- 0.0 - logfc
  horiz_line <- -1 * log10(p)

  if (! fc_col %in% colnames(table)) {
    stop("Column: ", fc_col, " is not in the table.")
  }
  if (! p_col %in% colnames(table)) {
    stop("Column: ", p_col, " is not in the table.")
  }

  df <- data.frame("xaxis" = as.numeric(table[[fc_col]]),
                   "yaxis" = as.numeric(table[[p_col]]))
  rownames(df) <- rownames(table)
  if (!is.null(minimum_p)) {
    low_idx <- df[["yaxis"]] < minimum_p
    df[low_idx, "yaxis"] <- minimum_p
  }
  ## Add the label column if it exists.
  if (!is.null(label_column) && !is.null(table[[label_column]])) {
    df[["label"]] <- table[[label_column]]
  } else {
    df[["label"]] <- rownames(table)
  }

  if (isTRUE(invert)) {
    df[["xaxis"]] <- df[["xaxis"]] * -1.0
  }


  ## This might have been converted to a string
  df[["logyaxis"]] <- -1.0 * log10(as.numeric(df[["yaxis"]]))
  df[["pcut"]] <- df[["yaxis"]] <= p
  df[["fccut"]] <- abs(df[["xaxis"]]) >= logfc

  df[["state"]] <- ifelse(table[[p_col]] > p, "pinsig",
                          ifelse(table[[p_col]] <= p &
                                   table[[fc_col]] >= logfc, "upsig",
                                 ifelse(table[[p_col]] <= p &
                                          table[[fc_col]] <= (-1 * logfc),
                                        "downsig", "fcinsig")))
  df[["pcut"]] <- as.factor(df[["pcut"]])
  df[["state"]] <- as.factor(df[["state"]])

  ## shape 25 is the down arrow, 22 is the square, 23 the diamond, 24 the up arrow
  state_shapes <- c(25, 22, 23, 24)
  names(state_shapes) <- c("downsig", "fcinsig", "pinsig", "upsig")

  color_column <- "pcut"
  color_column_number <- 2
  if (color_by != "p") {
    color_column <- "state"
    color_column_number <- 4
    color_list <- c("downsig" = "blue", "fcinsig" = "darkgrey",
                    "pinsig" = "darkgrey", "upsig" = "red")
  }
  ## Now make sure that the color column has the correct number of elements.
  if (length(color_list) != color_column_number) {
    mesg("The color list must have ", color_column_number,
         ", setting it to the default.")
  }

  ## Count the numbers in the categories
  num_downsig <- sum(df[["state"]] == "downsig")
  num_fcinsig <- sum(df[["state"]] == "fcinsig")
  num_pinsig <- sum(df[["state"]] == "pinsig")
  num_upsig <- sum(df[["state"]] == "upsig")

  plt <- NULL
  if (isTRUE(shapes_by_state)) {
    plt <- ggplot(data = df,
                  aes(x = .data[["xaxis"]], y = .data[["logyaxis"]],
                      label = .data[["label"]],
                      fill = color_column, colour = color_column, shape = "state"))
  } else {
    plt <- ggplot(data = df,
                  aes(x = .data[["xaxis"]], y = .data[["logyaxis"]],
                      label = .data[["label"]], fill = .data[[color_column]],
                      colour = .data[[color_column]]))
  }

  ## Now define when to put lines vs. points
  if (line_position == "bottom") {
    ## lines, then points.
    plt <- plt +
      ggplot2::geom_hline(yintercept = horiz_line, color = line_color, size=(size / 2)) +
      ggplot2::geom_vline(xintercept = logfc, color = line_color, size=(size / 2)) +
      ggplot2::geom_vline(xintercept = low_vert_line, color = line_color, size=(size / 2)) +
      ggplot2::geom_point(stat = "identity", size = size, alpha = alpha)
  } else {
    ## points, then lines
    plt <- plt +
      ggplot2::geom_point(stat = "identity", size = size, alpha = alpha) +
      ggplot2::geom_hline(yintercept = horiz_line, color = line_color, size=(size / 2)) +
      ggplot2::geom_vline(xintercept = logfc, color = line_color, size=(size / 2)) +
      ggplot2::geom_vline(xintercept = low_vert_line, color = line_color, size=(size / 2))
  }

  ## If shapes are being set by state,  add that to the legend now.
  if (isTRUE(shapes_by_state)) {
    plt <- plt +
      ggplot2::scale_shape_manual(
        name = "state", values = state_shapes,
        labels = c(
          glue("Down Sig.: {num_downsig}"),
          glue("FC Insig.: {num_fcinsig}"),
          glue("P Insig.: {num_pinsig}"),
          glue("Up Sig.: {num_upsig}")),
        guide = ggplot2::guide_legend(override.aes = aes(size = 3, fill = "grey")))
  }

  ## Now set the colors and axis labels
  plt <- plt +
    ggplot2::scale_fill_manual(name = color_column, values = color_list,
                               guide = "none") +
    ggplot2::scale_color_manual(name = color_column, values = color_list,
                                guide = "none") +
    ggplot2::xlab(label = fc_name) +
    ggplot2::ylab(label = p_name) +
    ## ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"))
  ##  axis.text.x = ggplot2::element_text(angle=-90))

  if (!is.null(label)) {
    if (is.numeric(label)) {
      reordered_idx <- order(df[["xaxis"]])
      reordered <- df[reordered_idx, ]
      sig_idx <- reordered[["logyaxis"]] > horiz_line
      reordered <- reordered[sig_idx, ]
      top <- head(reordered, n = label)
      bottom <- tail(reordered, n = label)
      df_subset <- rbind(top, bottom)
    } else if (is.character(label)) {
      sig_idx <- rownames(df) %in% label
      df_subset <- df[sig_idx, ]
    } else {
      stop("I do not understand this set of IDs to label.")
    }
    plt <- plt +
      ggrepel::geom_text_repel(data = df_subset,
                               aes(label = .data[["label"]], y = .data[["logyaxis"]],
                                   x = .data[["xaxis"]]),
                               colour = "black", box.padding = ggplot2::unit(0.5, "lines"),
                               point.padding = ggplot2::unit(1.6, "lines"),
                               arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")))
  }

  retlist <- list("plot" = plt,
                  "df" = df)
  return(retlist)
}

#' Theresa's volcano plots are objectively nicer because they are colored by condition.
#'
#' I therefore took a modified copy of her implementation and added it here.
#'
#' @param input Table of DE values, likely from combine_de_tables().
#' @param table_name Name the table!
#' @param alpha Make see-through.
#' @param fc_col Column containing the fold-change values.
#' @param fc_name Axis label.
#' @param line_color Color for the demarcation lines.
#' @param line_position Put the lines above or below the dots.
#' @param logfc Demarcation line for fold-change significance.
#' @param p_col Column containing the significance information.
#' @param p_name Axis label for the significance.
#' @param pval Demarcation for (in)significance.
#' @param shapes_by_state Change point shapes according to their states?
#' @param color_high Color for the ups.
#' @param color_low and the downs.
#' @param size Point size
#' @param invert Flip the plot?
#' @param label Label some points?
#' @param label_column Using this column in the data.
#' @export
plot_volcano_condition_de <- function(input, table_name, alpha = 0.5,
                                      fc_col = "logFC", fc_name = "log2 fold change",
                                      line_color = "black", line_position = "bottom", logfc = 1.0,
                                      p_col = "adj.P.Val", p_name = "-log10 p-value", pval = 0.05,
                                      shapes_by_state = FALSE,
                                      color_high = "darkred", color_low = "darkblue",
                                      size = 2, invert = FALSE, label = NULL,
                                      label_column = "hgncsymbol", label_size = 6) {

  low_vert_line <- 0.0 - logfc
  horiz_line <- -1 * log10(pval)

  if (! fc_col %in% colnames(input)) {
    stop("Column: ", fc_col, " is not in the table.")
  }
  if (! p_col %in% colnames(input)) {
    stop("Column: ", p_col, " is not in the table.")
  }
  df <- data.frame("xaxis" = as.numeric(input[[fc_col]]),
                   "yaxis" = as.numeric(input[[p_col]]))
  rownames(df) <- rownames(input)

  if (isTRUE(invert)) {
    df[["xaxis"]] <- df[["xaxis"]] * -1.0
  }
  ## Add the label column if it exists.
  if (!is.null(label_column) && !is.null(input[[label_column]])) {
    df[["label"]] <- input[[label_column]]
  } else {
    df[["label"]] <- rownames(input)
  }

  ## This might have been converted to a string
  df[["logyaxis"]] <- -1.0 * log10(as.numeric(df[["yaxis"]]))
  df[["pcut"]] <- df[["yaxis"]] <= pval
  df[["fccut"]] <- abs(df[["xaxis"]]) >= logfc

  df[["state"]] <- "insignificant"
  numerator_sig <- df[["xaxis"]] >= logfc & df[["yaxis"]] <= pval
  df[numerator_sig, "state"] <- "up"
  denominator_sig <- df[["xaxis"]] <= -1.0 * logfc & df[["yaxis"]] <= pval
  df[denominator_sig, "state"] <- "down"
  df[["state"]] <- as.factor(df[["state"]])

  ## I am not sure why something is setting color high/low to NULL.
  if (is.null(color_high)) {
    color_high <- "red"
  }
  if (is.null(color_low)) {
    color_low <- "blue"
  }
  plot_colors <- c("#555555", color_high, color_low)
  names(plot_colors) <- c("insignificant", "up", "down")

  plt <- ggplot(data = df,
                aes(x = .data[["xaxis"]], y = .data[["logyaxis"]],
                    label = .data[["label"]], fill = .data[["state"]],
                    colour = .data[["state"]]))

  ## Now define when to put lines vs. points
  if (is.null(line_position)) {
    plt <- plt +
      ggplot2::geom_point(stat = "identity", size = size, alpha = alpha)
  } else if (line_position == "bottom") {
    ## lines, then points.
    plt <- plt +
      ggplot2::geom_hline(yintercept = horiz_line, color = line_color, size = (size / 2)) +
      ggplot2::geom_vline(xintercept = logfc, color = line_color, size = (size / 2)) +
      ggplot2::geom_vline(xintercept = low_vert_line, color = line_color, size = (size / 2)) +
      ggplot2::geom_point(stat = "identity", size = size, alpha = alpha)

  } else if (line_position == "top") {
    ## points, then lines
    plt <- plt +
      ggplot2::geom_point(stat = "identity", size = size, alpha = alpha) +
      ggplot2::geom_hline(yintercept = horiz_line, color = line_color, size = (size / 2)) +
      ggplot2::geom_vline(xintercept = logfc, color = line_color, size = (size / 2)) +
      ggplot2::geom_vline(xintercept = low_vert_line, color = line_color, size = (size / 2))
  } else {
    mesg("Not printing volcano demarcation lines.")
    plt <- plt +
      ggplot2::geom_point(stat = "identity", size = size, alpha = alpha)
  }

  ## Now set the colors and axis labels
  plt <- plt +
    ggplot2::scale_fill_manual(name = "state", values = plot_colors,
                               guide = "none") +
    ggplot2::scale_color_manual(name = "state", values = plot_colors,
                                guide = "none") +
    ggplot2::xlab(label = fc_name) +
    ggplot2::ylab(label = p_name) +
    ## ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"))
  ##  axis.text.x = ggplot2::element_text(angle=-90))

  if (!is.null(label)) {
    if (is.numeric(label)) {
      reordered_idx <- order(df[["xaxis"]])
      reordered <- df[reordered_idx, ]
      sig_idx <- reordered[["logyaxis"]] > horiz_line
      reordered <- reordered[sig_idx, ]
      top <- head(reordered, n = label)
      bottom <- tail(reordered, n = label)
      df_subset <- rbind(top, bottom)
    } else if (is.character(label)) {
      sig_idx <- df[["label"]] %in% label
      mesg("Found ", sum(sig_idx), " of the labeled genes.")
      df_subset <- df[sig_idx, ]
    } else {
      stop("I do not understand this set of IDs to label.")
    }
    plt <- plt +
      ggrepel::geom_text_repel(data = df_subset,
                               aes(label = .data[["label"]], y = .data[["logyaxis"]],
                                   x = .data[["xaxis"]]),
                               colour = "black", box.padding = ggplot2::unit(0.5, "lines"),
                               point.padding = ggplot2::unit(1.6, "lines"),
                               size = label_size,
                               arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")))
  }

  retlist <- list("plot" = plt,
                  "df" = df)
  return(retlist)
}

#' Plot the rank order of the data in two tables against each other.
#'
#' Steve Christensen has some neat plots showing the relationship between two
#' tables.  I thought they were cool, so I co-opted the idea in this
#' function.
#'
#' @param first First table of values.
#' @param second Second table of values, if null it will use the first.
#' @param first_type Assuming this is from all_pairwise(), use this method.
#' @param second_type Ibid.
#' @param first_table Again, assuming all_pairwise(), use this to choose the
#'  table to extract.
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
rank_order_scatter <- function(first, second = NULL, first_type = "limma",
                               second_type = "limma", first_table = NULL, alpha = 0.5,
                               second_table = NULL, first_column = "logFC",
                               second_column = "logFC", first_p_col = "adj.P.Val",
                               second_p_col = "adj.P.Val", p_limit = 0.05,
                               both_color = "red", first_color = "green",
                               second_color = "blue", no_color = "black") {
  if (is.null(second)) {
    second <- first
  }

  ## If the two elements are equal, then default to two different tables
  ## If the two elements are not equal, default to the same table.
  test <- try(testthat::expect_equal(first, second), silent = TRUE)
  if (class(test)[1] == "try-error") {
    ## They are not equal.
    if (is.null(first_table)) {
      message("No first table was provided, setting it to the first table.")
      first_table <- 1
    }
    if (is.null(second_table)) {
      message("No second table was provided, setting it to the first table.")
      second_table <- 1
    }
  } else {
    ## Two different de results.
    if (is.null(first_table)) {
      first_table <- 1
    }
    if (is.null(second_table)) {
      second_table <- 2
    }
  }

  if (!is.null(first[[first_type]])) {
    first <- first[[first_type]][["all_tables"]]
  }
  if (!is.null(second[[second_type]])) {
    second <- second[[second_type]][["all_tables"]]
  }
  table1 <- first[[first_table]]
  table2 <- second[[second_table]]
  merged <- merge(table1, table2, by = "row.names")
  rownames(merged) <- merged[["Row.names"]]
  merged <- merged[, -1]

  if (first_column == second_column) {
    c1 <- glue("{first_column}.x")
    c2 <- glue("{first_column}.y")
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
    p1 <- glue("{first_p_col}.x")
    p2 <- glue("{first_p_col}.y")
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

  first_table_colname <- glue(
    "Table: {first_table}, Type: {first_type}, column: {first_column}")
  second_table_colname <- glue(
    "Table: {second_table}, Type: {second_type}, column: {second_column}")

  plt <- ggplot(data = merged,
                aes(color = .data[["state"]], fill = .data[["state"]],
                    x = .data[["x"]], y = .data[["y"]], label = .data[["label"]])) +
    ggplot2::geom_point(size = 1, alpha = alpha) +
    ggplot2::scale_color_manual(name = "state",
                                values = c("both"=both_color,
                                           "first"=first_color,
                                           "second"=second_color,
                                           "neither"=no_color)) +
    ggplot2::geom_smooth(method = "loess", color = "lightblue") +
    ggplot2::ylab(glue("Rank order of {second_table_colname}")) +
    ggplot2::xlab(glue("Rank order of {first_table_colname}")) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = "none",
                   axis.text = ggplot2::element_text(size = base_size, colour = "black"))

  model_test <- try(lm(formula = y ~ x, data = merged), silent = TRUE)
  model_summary <- summary(model_test)
  cor <- cor.test(merged[[c1]], merged[[c2]], method = "pearson")
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
#'  mean greater/less than 0 followed by 2 fold and 4 fold cutoffs.
#' @param invert Reverse the order of contrasts for readability?
#' @param p Chosen p-value cutoff.
#' @param z Choose instead a z-score cutoff.
#' @param p_type Adjusted or not?
#' @param according_to limma, deseq, edger, basic, or all of the above.
#' @param order Choose a specific order for the plots.
#' @param maximum Set a specific limit on the number of genes on the x-axis.
#' @param ... More arguments are passed to arglist.
#' @return list containing the significance bar plots and some information to
#'  hopefully help interpret them.
#' @examples
#' \dontrun{
#'  expt <- create_expt(metadata = "some_metadata.xlsx", gene_info = annotations)
#'  pairwise_result <- all_pairwise(expt)
#'  combined_result <- combine_de_tables(pairwise_result)
#'  ## Damn I wish I were smrt enough to make this elegant, but I cannot.
#'  barplots <- significant_barplots(combined_result)
#' }
#' @export
significant_barplots <- function(combined, lfc_cutoffs = c(0, 1, 2), invert = FALSE,
                                 p = 0.05, z = NULL, p_type = "adj",
                                 according_to = "all", order = NULL, maximum = NULL, ...) {
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
    test_column <- glue("{type}_logfc")
    if (! test_column %in% colnames(combined[["data"]][[1]])) {
      message("We do not have the ", test_column, " in the data, skipping ", type, ".")
      next
    }
    for (fc in lfc_cutoffs) {
      ## This is a bit weird and circuituous
      ## The most common caller of this function is in fact extract_significant_genes
      fc_sig <- sm(extract_significant_genes(combined, lfc = fc, according_to = according_to,
                                             p = p, z = z, n = NULL, excel = FALSE,
                                             p_type = p_type, sig_bar = FALSE, ma = FALSE))
      table_length <- length(fc_sig[[type]][["ups"]])
      fc_name <- glue("fc_{fc}")
      fc_names <- append(fc_names, fc_name)

      for (tab in seq_len(table_length)) {
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
    for (t in seq_len(table_length)) {
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
    up <- suppressWarnings(reshape2::melt(up, id.var = "comparisons"))
    up[["comparisons"]] <- factor(up[["comparisons"]], levels = comparisons)
    up[["variable"]] <- factor(up[["variable"]],
                               levels = c("a_up_inner", "b_up_middle", "c_up_outer"))
    up[["value"]] <- as.numeric(up[["value"]])
    ## Repeat with the set of down materials
    down <- cbind(comparisons, down_all[[type]], down_mid[[type]], down_max[[type]])
    down <- as.data.frame(down)
    colnames(down) <- c("comparisons", "a_down_inner", "b_down_middle", "c_down_outer")
    downlist[[type]] <- down
    colnames(down) <- c("comparisons", "a_down_inner", "b_down_middle", "c_down_outer")
    down <- suppressWarnings(reshape2::melt(down, id.var = "comparisons"))
    down[["comparisons"]] <- factor(down[["comparisons"]], levels = comparisons)
    ##        down[["variable"]] <- factor(down[["variable"]],
    ##        levels = c("a_down_inner","b_down_middle","c_down_outer"))
    down[["variable"]] <- factor(down[["variable"]],
                                 levels = c("c_down_outer", "b_down_middle", "a_down_inner"))
    up[["variable"]] <- factor(up[["variable"]],
                               levels = c("c_up_outer", "b_up_middle", "a_up_inner"))
    down[["value"]] <- as.numeric(down[["value"]]) * -1
    tables_up[[type]] <- up
    tables_down[[type]] <- down
    plots[[type]] <- plot_significant_bar(up, down, maximum = maximum,
                                          ...)
    ## plots[[type]] <- plot_significant_bar(up, down, maximum = maximum) #, ...)
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

#' Extract overlapping groups from an upset
#'
#' Taken from: https://github.com/hms-dbmi/UpSetR/issues/85
#' and lightly modified to match my style and so I could more
#' easily understand what it is doing.
#'
#' @param lst upset data structure.
#' @param sort Sort the result?
#' @export
overlap_groups <- function(input, sort = TRUE) {
  ## FIXME: Make use of S4 here
  input_mtrx <- NULL
  element_names <- NULL
  if ("list" %in% class(input)) {
    input_mtrx <- UpSetR::fromList(input) == 1
    element_names <- unlist(input)
  } else if ("upset" %in% class(input)) {
    stop("The upsetR fromList seems to strip out the gene names, don't use it until I figure out what is up.")
  }

  ## lst could look like this:
  ## $one
  ## [1] "a" "b" "c" "e" "g" "h" "k" "l" "m"
  ## $two
  ## [1] "a" "b" "d" "e" "j"
  ## $three
  ## [1] "a" "e" "f" "g" "h" "i" "j" "l" "m"

  ##     one   two three
  ## a  TRUE  TRUE  TRUE
  ## b  TRUE  TRUE FALSE
  ##...
  ## condensing matrix to unique combinations elements
  combination_mtrx <- unique(input_mtrx)
  groups <- list()
  num_combinations <- nrow(combination_mtrx)
  ## going through all unique combinations and collect elements for each in a list
  for (i in seq_len(num_combinations)) {
    combination <- combination_mtrx[i, ]
    my_elements <- which(apply(input_mtrx, 1, function(x) all(x == combination)))
    attr(my_elements, "groups") <- combination
    groups[[paste(colnames(combination_mtrx)[combination], collapse = ":")]] <- my_elements
    #my_elements
    ## attr(,"groups")
    ##   one   two three
    ## FALSE FALSE  TRUE
    ##  f  i
    ## 12 13
  }
  if (sort) {
    groups <- groups[order(sapply(groups, function(x) length(x)), decreasing = TRUE)]
  }
  attr(groups, "elements") <- element_names
  return(groups)
  ## save element list to facilitate access using an index in case rownames are not named
}

#' Mostly as a reminder of how to get the gene IDs from a specific group in an upset plot.
#'
#' Given a set of groups from upsetr, extract the elements from one of them.
#' @param overlapping_groups Result from overlap_groups, which just makes an indexed
#'  version of the genes by venn/upset group.
#' @param group Name of the subset of interest, something like 'a:b' for the union of a:b.
#' @export
overlap_geneids <- function(overlapping_groups, group) {
  gene_ids <- attr(overlapping_groups, "elements")[overlapping_groups[[group]]]
  return(gene_ids)
}

#' Make an upset plot of all up/down genes in a set of contrasts.
#'
#' This is intended to give a quick and dirty view of the genes
#' observed in a series of de comparisons.
#'
#' @param combined Result from combine_de_tables.
#' @param according_to Choose the lfc column to use.
#' @param lfc Choose the logFC
#' @param adjp and the p-value.
upsetr_all <- function(combined, according_to = "deseq",
                       lfc = 1.0, adjp = 1.0) {
  ud_list <- list()
  for (t in names(combined[["data"]])) {
    t_data <- combined[["data"]][[t]]
    fc_col <- paste0(according_to, "_logfc")
    p_col <- paste0(according_to, "_adjp")
    up_name <- glue("{t}_up")
    up_idx <- t_data[[fc_col]] >= lfc &
      t_data[[p_col]] <= adjp
    up_ids <- rownames(t_data)[up_idx]
    ud_list[[up_name]] <- up_ids
    down_name <- glue("{t}_down")
    down_idx <- t_data[[fc_col]] <= (-1.0 * lfc) &
      t_data[[p_col]] <= adjp
    down_ids <- rownames(t_data)[down_idx]
    ud_list[[down_name]] <- down_ids
  }

  upset_combined <- UpSetR::upset(data = UpSetR::fromList(ud_list),
                                  nsets = length(ud_list))
  return(upset_combined)
}

#' Use UpSetR to compare significant gene lists.
#'
#' @param sig datastructure of significantly DE genes.
#' @param according_to Choose your favorite method.
#' @param contrasts Choose a specific contrast(s)
#' @param up Make a plot of the up genes?
#' @param down Make a plot of the down genes?
#' @param both Make a plot of the up+down genes?
#' @param scale Make the numbers larger and easier to read?
#' @param ... Other parameters to pass to upset().
#' @export
upsetr_sig <- function(sig, according_to = "deseq", contrasts = NULL, up = TRUE,
                       down = TRUE, both = FALSE, scale = 2, ...) {

  ## Start by pulling the gene lists from the significant gene sets.
  start <- sig[[according_to]]
  ups <- NULL
  downs <- NULL
  if (isTRUE(up)) {
    ups <- start[["ups"]]
  }
  if (isTRUE(down)) {
    downs <- start[["downs"]]
  }

  ## Figure out which contrasts (if not all of them) are desired to plot.
  wanted_contrasts <- names(start[["ups"]])
  if (!is.null(contrasts)) {
    wanted_idx <- wanted_contrasts %in% contrasts
    wanted_contrasts <- wanted_contrasts[wanted_idx]
  }

  ## Setup the lists of significant genes to plot.
  upsetr_up_list <- list()
  upsetr_down_list <- list()
  upsetr_both_list <- list()
  for (entry in wanted_contrasts) {
    if (!is.null(ups)) {
      upsetr_up_list[[entry]] <- rownames(ups[[entry]])
    }
    if (!is.null(downs)) {
      upsetr_down_list[[entry]] <- rownames(downs[[entry]])
    }
    if (isTRUE(both)) {
      upsetr_both_list[[entry]] <- c(rownames(ups[[entry]]),
                                     rownames(downs[[entry]]))
    }
  } ## End looking for things to list

  ## Do the plots.
  retlist <- list()
  if (isTRUE(up)) {
    retlist[["up"]] <- UpSetR::upset(UpSetR::fromList(upsetr_up_list),
                                     text.scale = scale, ...)
    retlist[["up_groups"]] <- overlap_groups(upsetr_up_list)
  }
  if (isTRUE(down)) {
    retlist[["down"]] <- UpSetR::upset(UpSetR::fromList(upsetr_down_list),
                                       text.scale = scale, ...)
    retlist[["down_groups"]] <- overlap_groups(upsetr_down_list)
  }
  if (isTRUE(both)) {
    retlist[["both"]] <- UpSetR::upset(UpSetR::fromList(upsetr_both_list),
                                       text.scale = scale, ...)
    retlist[["both_groups"]] <- overlap_groups(upsetr_both_list)
  }

  return(retlist)
}

## EOF
