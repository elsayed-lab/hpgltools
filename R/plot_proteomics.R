#' Plot mzXML peak intensities with respect to m/z.
#'
#' I want to have a pretty plot of peak intensities and m/z.  The plot provided
#' by this function is interesting, but suffers from some oddities; notably that
#' it does not currently separate the MS1 and MS2 data.
#'
#' @param mzxml_data The data structure from extract_mzxml or whatever it is.
#' @param loess Do a loess smoothing from which to extract a function
#'  describing the data?  This is terribly slow, and in the data I have
#'  examined so far, not very helpful, so it is FALSE by default.
#' @param alpha Make the plotted dots opaque to this degree.
#' @param ms1 Include MS1 data in the plot?
#' @param ms2 Include MS2 data in the plot?
#' @param x_scale Plot the x-axis on a non linear scale?
#' @param y_scale Plot the y-axis on a non linear scale?
#' @param ... Extra arguments for the downstream functions.
#' @return ggplot2 goodness.
#' @export
plot_intensity_mz <- function(mzxml_data, loess = FALSE, alpha = 0.5, ms1=TRUE, ms2=TRUE,
                              x_scale = NULL, y_scale = NULL, ...) {
  arglist <- list(...)
  metadata <- mzxml_data[["metadata"]]
  colors <- mzxml_data[["colors"]]
  sample_data <- mzxml_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)
  keepers <- c()
  for (i in 1:samples) {
    name <- metadata[i, "sampleid"]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    keepers <- c(keepers, i)
    ## message("Adding ", name)
    plotted_table <- sample_data[[i]][["scans"]]
    if (!isTRUE(ms1)) {
      kept_idx <- plotted_table[["level"]] != "MS1"
      plotted_table <- plotted_table[kept_idx, ]
    }
    if (!isTRUE(ms2)) {
      kept_idx <- plotted_table[["level"]] != "MS2"
      plotted_table <- plotted_table[kept_idx, ]
    }
    plotted_data <- plotted_table[, c("basepeakmz", "basepeakintensity")]
    plotted_data[["sample"]] <- name
    plotted_data <- plotted_data[, c("sample", "basepeakmz", "basepeakintensity")]
    colnames(plotted_data) <- c("sample", "mz", "intensity")
    ## Re-order the columns because I like sample first.
    plot_df <- rbind(plot_df, plotted_data)
  }

  ## Drop rows from the metadata and colors which had errors.
  metadata <- metadata[keepers, ]
  colors <- colors[keepers]

  chosen_palette <- "Dark2"
  sample_colors <- sm(
      grDevices::colorRampPalette(
                     RColorBrewer::brewer.pal(samples, chosen_palette))(samples))

  ## Randomize the rows of the df so we can see if any sample is actually overrepresented
  plot_df <- plot_df[sample(nrow(plot_df)), ]

  if (!is.null(x_scale)) {
    plot_df[["mz"]] <- check_plot_scale(plot_df[["mz"]], scale)[["data"]]
  }
  if (!is.null(y_scale)) {
    plot_df[["intensity"]] <- check_plot_scale(plot_df[["intensity"]], scale)[["data"]]
  }

  int_vs_mz <- ggplot(data = plot_df, aes_string(x = "mz", y = "intensity",
                                                 fill = "sample", colour = "sample")) +
    ggplot2::geom_point(alpha = alpha, size = 0.5) +
    ggplot2::scale_fill_manual(
                 name = "Sample", values = sample_colors,
                 guide = ggplot2::guide_legend(override.aes = aes(size = 3))) +
    ggplot2::scale_color_manual(
                 name = "Sample", values = sample_colors,
                 guide = ggplot2::guide_legend(override.aes = aes(size = 3))) +
    ggplot2::theme_bw(base_size = base_size)

  if (!is.null(x_scale)) {
    int_vs_mz <- int_vs_mz + ggplot2::scale_x_continuous(trans = scales::log2_trans())
  }
  if (!is.null(y_scale)) {
    int_vs_mz <- int_vs_mz + ggplot2::scale_y_continuous(trans = scales::log2_trans())
  }

  if (isTRUE(lowess)) {
    int_vs_mz <- int_vs_mz +
      ggplot2::geom_smooth(method = "loess", size = 1.0)
  }
  retlist <- list(
      "data" = plotted_data,
      "plot" = int_vs_mz)
  return(retlist)
}

#' Make a boxplot out of some of the various data available in the mzxml data.
#'
#' There are a few data within the mzXML raw data files which are likely
#' candidates for simple summary via a boxplot/densityplot/whatever.  For the
#' moment I am just doing boxplots of a few of them.  Since my metadata
#' extractor dumps a couple of tables, one must choose a desired table and
#' column from it to plot.
#'
#' @param mzxml_data Provide a list of mzxml data, one element for each sample.
#' @param table One of precursors or scans
#' @param column One of the columns from the table; if 'scans' is chosen, then
#'  likely choices include: 'peakscount', 'basepeakmz', 'basepeakintensity'; if
#'  'precursors' is chosen, then the only likely choice for the moment is
#'  'precursorintensity'.
#' @param violin Print the samples as violins rather than only box/whiskers?
#' @param names Names for the x-axis of the plot.
#' @param plot_title Title the plot?
#' @param scale Put the data on a specific scale?
#' @param ... Further arguments, presumably for colors or some such.
#' @return Boxplot describing the requested column of data in the set of mzXML files.
#' @export
plot_mzxml_boxplot <- function(mzxml_data, table = "precursors", column = "precursorintensity",
                               violin = FALSE, names = NULL, plot_title = NULL, scale = NULL, ...) {
  arglist <- list(...)
  metadata <- mzxml_data[["metadata"]]
  colors <- mzxml_data[["colors"]]
  sample_data <- mzxml_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)
  keepers <- c()
  for (i in 1:samples) {
    name <- metadata[i, "sampleid"]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    keepers <- c(keepers, i)
    ## message("Adding ", name)
    names(colors)[i] <- name
    plotted_table <- sample_data[[i]][[table]]
    plotted_data <- as.data.frame(plotted_table[[column]])
    plotted_data[["sample"]] <- name
    plotted_data[["color"]] <- colors[[name]]
    colnames(plotted_data) <- c(column, "sample", "color")
    ## Re-order the columns because I like sample first.
    plotted_data <- plotted_data[, c("sample", column, "color")]
    plot_df <- rbind(plot_df, plotted_data)
  }

  ## Drop rows from the metadata and colors which had errors.
  if (length(keepers) > 0) {
    metadata <- metadata[keepers, ]
    colors <- colors[keepers]
  } else {
    stop("Something bad happened to the set of kept samples.")
  }

  scale_data <- check_plot_scale(plot_df[[column]], scale)
  if (is.null(scale)) {
    scale <- scale_data[["scale"]]
  }
  plot_df[[column]] <- scale_data[["data"]]
  plot_df[["color"]] <- as.factor(plot_df[["color"]])

  boxplot <- ggplot2::ggplot(data = plot_df, ggplot2::aes_string(x = "sample", y = column))
  if (isTRUE(violin)) {
    boxplot <- boxplot +
      ggplot2::geom_violin(aes_string(fill = "sample"),
                           width = 1, scale = "area", show.legend = FALSE) +
      ggplot2::geom_boxplot(na.rm = TRUE, alpha = 0.3, color = "black", size = 0.5,
                            outlier.alpha = 0.01, width = 0.2)
  } else {
    boxplot <- boxplot +
      sm(ggplot2::geom_boxplot(aes_string(fill = "sample"),
                               na.rm = TRUE, fill = colors, size = 0.5,
                               outlier.size = 1.5,
                               outlier.colour = ggplot2::alpha("black", 0.2)))
  }
  boxplot <- boxplot +
    ggplot2::scale_fill_manual(values = as.character(colors)) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::xlab("Sample") + ggplot2::ylab(column)
  if (!is.null(plot_title)) {
    boxplot <- boxplot + ggplot2::ggtitle(plot_title)
  }
  if (!is.null(names)) {
    boxplot <- boxplot + ggplot2::scale_x_discrete(labels = names)
  }
  scale <- "log"
  if (scale == "log") {
    boxplot <- boxplot + ggplot2::scale_y_continuous(trans = scales::log2_trans())
  } else if (scale == "logdim") {
    boxplot <- boxplot + ggplot2::coord_trans(y = "log2")
  } else if (isTRUE(scale)) {
    boxplot <- boxplot + ggplot2::scale_y_log10()
  }

  return(boxplot)
}

#' Count some aspect(s) of the pyprophet data and plot them.
#'
#' This function is mostly redundant with the plot_mzxml_boxplot above.
#' Unfortunately, the two data types are subtly different enough that I felt it
#' not worth while to generalize the functions.
#'
#' @param pyprophet_data List containing the pyprophet results.
#' @param type What to count/plot?
#' @param keep_real Do we keep the real data when plotting the data? (perhaps
#'  we only want the decoys)
#' @param keep_decoys Do we keep the decoys when plotting the data?
#' @param expt_names Names for the x-axis of the plot.
#' @param label_chars Maximum number of characters before abbreviating sample names.
#' @param plot_title Title the plot?
#' @param scale Put the data on a specific scale?
#' @param ... Further arguments, presumably for colors or some such.
#' @return Boxplot describing the desired column from the data.
#' @export
plot_pyprophet_counts <- function(pyprophet_data, type = "count", keep_real = TRUE,
                                  keep_decoys = TRUE, expt_names = NULL, label_chars = 10,
                                  plot_title = NULL, scale = NULL, ...) {
  arglist <- list(...)
  metadata <- pyprophet_data[["metadata"]]
  colors <- pyprophet_data[["colors"]]
  sample_data <- pyprophet_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)

  ## Reset the sample names if one wants a specific column from the metadata.
  if (!is.null(expt_names) & class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      names(sample_data) <- make.names(metadata[[expt_names]], unique = TRUE)
    } else {
      names(sample_data) <- expt_names
    }
  }

  keepers <- c()
  plotted_data <- data.frame()
  for (i in 1:samples) {
    name <- names(sample_data)[i]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    keepers <- c(keepers, i)
    ## message("Adding ", name)
    plotted_table <- sample_data[[i]]
    if (!isTRUE(keep_decoys)) {
      good_idx <- plotted_table[["decoy"]] != 1
      plotted_table <- plotted_table[good_idx, ]
    }
    if (!isTRUE(keep_real)) {
      good_idx <- plotted_table[["decoy"]] != 0
      plotted_table <- plotted_table[good_idx, ]
    }

    row_condition <- as.character(metadata[i, "condition"])
    row_color <- colors[i]
    if (type == "protein_count") {
      proteins <- length(unique(plotted_table[["proteinname"]]))
      row <- c(name, proteins, row_condition, row_color)
    } else if (type == "count") {
      ## Taking nrow is a simplistic way to count the number of identifications.
      row <- c(name, nrow(plotted_table), row_condition, row_color)
    } else if (type == "intensity") {
      row <- c(name, sum(as.numeric(plotted_table[["intensity"]])), row_condition, row_color)
    } else {
      row <- c(name, sum(as.numeric(plotted_table[[type]])), row_condition, row_color)
    }

    plotted_data <- rbind(plotted_data, row)
  } ## End the for loop.

  y_label <- "Identifications per sample"
  if (type == "count") {
    colnames(plotted_data) <- c("id", "sum", "condition", "colors")
    plotted_data[["sum"]] <- as.numeric(plotted_data[["sum"]])
  } else if (type == "intensity") {
    y_label <- "Sum of intensities per sample."
    colnames(plotted_data) <- c("id", "sum", "condition", "colors")
    plotted_data[["sum"]] <- as.numeric(plotted_data[["sum"]])
  } else {
    y_label <- glue::glue("Sum of {type} per sample.")
    colnames(plotted_data) <- c("id", "sum", "condition", "colors")
    plotted_data[["sum"]] <- as.numeric(plotted_data[["sum"]])
  }

  if (!is.null(label_chars) & is.numeric(label_chars)) {
    plot_df[["id"]] <- abbreviate(plot_df[["id"]], minlength = label_chars)
  }
  our_plot <- plot_sample_bars(plotted_data, integerp = TRUE,
                               text = TRUE, yscale = "log2",
                               ylabel = y_label)

  retlist <- list(
      "df" = plotted_data,
      "plot" = our_plot)
  return(retlist)
}

#' Invoked plot_pyprophet_counts() twice, once for the x-axis, and once for the y.
#'
#' Then plot the result, hopefully adding some new insights into the state of
#' the post-pyprophet results.  By default, this puts the number of
#' identifications (number of rows) on the x-axis for each sample, and the sum
#' of intensities on the y.  Currently missing is the ability to change this
#' from sum to mean/median/etc.  That should trivially be possible via the
#' addition of arguments for the various functions of interest.
#'
#' @param pyprophet_data List of pyprophet matrices by sample.
#' @param keep_real Use the real identifications (as opposed to the decoys)?
#' @param size Size of the glyphs used in the plot.
#' @param label_size Set the label sizes.
#' @param keep_decoys Use the decoy identifications (vs. the real)?
#' @param expt_names Manually change the labels to some other column than sample.
#' @param label_chars Maximum number of characters in the label before
#'  shortening.
#' @param x_type Column in the data to put on the x-axis.
#' @param y_type Column in the data to put on the y-axis.
#' @param plot_title Plot title.
#' @param scale Put the data onto the log scale?
#' @param ... Extra arguments passed along.
#' @export
plot_pyprophet_xy <- function(pyprophet_data, keep_real = TRUE, size = 6, label_size = 4,
                              keep_decoys = TRUE, expt_names = NULL, label_chars = 10,
                              x_type = "count", y_type = "intensity",
                              plot_title = NULL, scale = NULL, ...) {
  arglist <- list(...)

  x_data <- plot_pyprophet_counts(pyprophet_data,
                                  type = x_type,
                                  keep_real = keep_real,
                                  keep_decoys = keep_decoys,
                                  expt_names = expt_names,
                                  label_chars = label_chars,
                                  plot_title = plot_title,
                                  scale = scale,
                                  ...)
  y_data <- plot_pyprophet_counts(pyprophet_data,
                                  type = y_type,
                                  keep_real = keep_real,
                                  keep_decoys = keep_decoys,
                                  expt_names = expt_names,
                                  label_chars = label_chars,
                                  plot_title = plot_title,
                                  scale = scale,
                                  ...)
  the_df <- x_data[["df"]]
  y_df <- y_data[["df"]]
  colnames(the_df)[2] <- x_type
  the_df[[y_type]] <- y_df[[2]]


  color_listing <- the_df[, c("condition", "colors")]
  color_listing <- unique(color_listing)
  color_list <- as.character(color_listing[["colors"]])
  names(color_list) <- as.character(color_listing[["condition"]])

  sc_plot <- ggplot(data = the_df,
                    aes_string(x = x_type, y = y_type, label = "id")) +
    ggplot2::geom_point(size = size, shape = 21,
                        aes_string(colour = "as.factor(condition)",
                                   fill = "as.factor(condition)")) +
    ggplot2::geom_point(size = size, shape = 21, colour = "black", show.legend = FALSE,
                        aes_string(fill = "as.factor(condition)")) +
    ggplot2::scale_color_manual(name = "Condition",
                                guide = "legend",
                                values = color_list) +
    ggplot2::scale_fill_manual(name = "Condition",
                               guide = "legend",
                               values = color_list) +
    ggrepel::geom_text_repel(aes_string(label = "id"),
                             size = label_size, box.padding = ggplot2::unit(0.5, "lines"),
                             point.padding = ggplot2::unit(1.6, "lines"),
                             arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")))
  return(sc_plot)
}

#' Make a boxplot out of some of the various data available in the pyprophet
#' data.
#'
#' This function is mostly redundant with the plot_mzxml_boxplot above.
#' Unfortunately, the two data types are subtly different enough that I felt it
#' not worth while to generalize the functions.
#'
#' @param pyprophet_data List containing the pyprophet results.
#' @param column What column of the pyprophet scored data to plot?
#' @param keep_real Do we keep the real data when plotting the data? (perhaps
#'  we only want the decoys)
#' @param keep_decoys Do we keep the decoys when plotting the data?
#' @param expt_names Names for the x-axis of the plot.
#' @param label_chars Maximum number of characters before abbreviating sample names.
#' @param plot_title Title the plot?
#' @param scale Put the data on a specific scale?
#' @param ... Further arguments, presumably for colors or some such.
#' @return Boxplot describing the desired column from the data.
#' @export
plot_pyprophet_distribution <- function(pyprophet_data, column = "delta_rt", keep_real = TRUE,
                                        keep_decoys = TRUE, expt_names = NULL, label_chars = 10,
                                        plot_title = NULL, scale = NULL, ...) {
  arglist <- list(...)
  metadata <- pyprophet_data[["metadata"]]
  colors <- pyprophet_data[["colors"]]
  sample_data <- pyprophet_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)

  ## Reset the sample names if one wants a specific column from the metadata.
  if (!is.null(expt_names) & class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      names(sample_data) <- make.names(metadata[[expt_names]], unique = TRUE)
    } else {
      names(sample_data) <- expt_names
    }
  }

  keepers <- c()
  for (i in 1:samples) {
    name <- names(sample_data)[i]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    keepers <- c(keepers, i)
    ## message("Adding ", name)
    plotted_table <- sample_data[[i]]
    if (!isTRUE(keep_decoys)) {
      good_idx <- plotted_table[["decoy"]] != 1
      plotted_table <- plotted_table[good_idx, ]
    }
    if (!isTRUE(keep_real)) {
      good_idx <- plotted_table[["decoy"]] != 0
      plotted_table <- plotted_table[good_idx, ]
    }
    plotted_data <- as.data.frame(plotted_table[c("sequence", "proteinname",
                                                  "aggr_fragment_annotation", column)])
    plotted_data[["sample"]] <- name
    colnames(plotted_data) <- c("sequence", "proteinname", "fragment", column, "sample")
    ## Re-order the columns because I like sample first.
    plotted_data <- plotted_data[, c("sample", column, "sequence", "proteinname", "fragment")]
    plot_df <- rbind(plot_df, plotted_data)
  }
  ## I am not certain this is valid.
  plot_df[[column]] <- abs(plot_df[[column]])

  ##  testing <- data.table::as.data.table(plot_df)
  ##  recast_dt <- data.table::dcast.data.table(data = testing,
  ##                                            formula = sequence+proteinname~sample,
  ##                                            fun.aggregate = mean,
  ##                                            value.var = "intensity")
  ##  names <- recast_dt[["proteinname"]]
  ##  sequences <- recast_dt[["sequence"]]
  ##  recast_dt[, c("proteinname", "sequence") := NULL]
  ##  nan_idx <- is.na(recast_dt)
  ##  recast_dt[nan_idx] <- 0
  ##  recast_norm <- log2(
  ##    1 + preprocessCore::normalize.quantiles.robust(as.matrix(recast_dt)))
  ##  remelt <- as.data.table(recast_norm)
  ##  remelt[["proteinname"]] <- names
  ##  remelt[["sequence"]] <- sequences
  ##  remelted <- data.table::melt(data = remelt, value.name = "intensity")
  ##  colnames(remelted) <- c("proteinname", "sequence", "sample", "intensity")

  ## Drop rows from the metadata and colors which had errors.
  if (length(keepers) > 0) {
    metadata <- metadata[keepers, ]
    colors <- colors[keepers]
  } else {
    stop("Something bad happened to the set of kept samples.")
  }

  scale_data <- check_plot_scale(plot_df[[column]], scale)
  if (is.null(scale)) {
    scale <- scale_data[["scale"]]
  }
  plot_df[[column]] <- scale_data[["data"]]

  if (!is.null(label_chars) & is.numeric(label_chars)) {
    plot_df[["sample"]] <- abbreviate(plot_df[["sample"]], minlength = label_chars)
  }
  boxplot <- ggplot2::ggplot(data = plot_df, ggplot2::aes_string(x = "sample", y = column)) +
    sm(ggplot2::geom_boxplot(na.rm = TRUE,
                             ggplot2::aes_string(fill = "sample"),
                             fill = colors,
                             size = 0.5,
                             outlier.size = 1.5,
                             outlier.colour = ggplot2::alpha("black", 0.2))) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::xlab("Sample") + ggplot2::ylab(column)
  if (!is.null(plot_title)) {
    boxplot <- boxplot + ggplot2::ggtitle(plot_title)
  }

  density <- ggplot(data = plot_df, ggplot2::aes_string(x = column, colour = "sample")) +
    ggplot2::geom_density(aes_string(x = column, y = "..count..", fill = "sample"),
                          position = "identity", na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = as.character(colors)) +
    ggplot2::scale_fill_manual(values = ggplot2::alpha(as.character(colors), 0.1)) +
    ggplot2::ylab("Number of genes.") +
    ggplot2::xlab("Number of hits/gene.") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   legend.key.size = ggplot2::unit(0.3, "cm"))
  density <- directlabels::direct.label(density)

  violin <- ggplot(data = plot_df, aes_string(x = "sample", y = column)) +
    ggplot2::geom_violin(aes_string(fill = "sample"), width = 1, scale = "area") +
    ggplot2::geom_boxplot(aes_string(fill = "sample"), outlier.alpha = 0.01, width = 0.2) +
    ggplot2::scale_fill_manual(values = as.character(colors)) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   legend.position = "none")

  dotboxplot <- boxplot +
    ggplot2::geom_jitter(shape = 16, position = ggplot2::position_jitter(0.1),
                         size = 2, alpha = 0.2)

  if (scale == "log") {
    boxplot <- boxplot + ggplot2::scale_y_continuous(trans = scales::log2_trans())
    dotboxplot <- dotboxplot + ggplot2::scale_y_continuous(trans = scales::log2_trans())
    violin <- violin + ggplot2::scale_y_continuous(trans = scales::log2_trans())
  } else if (scale == "logdim") {
    boxplot <- boxplot + ggplot2::coord_trans(y = "log2")
    dotboxplot <- dotboxplot + ggplot2::coord_trans(y = "log2")
    violin <- violin + ggplot2::coord_trans(y = "log2")
  } else if (isTRUE(scale)) {
    boxplot <- boxplot + ggplot2::scale_y_log10()
    dotboxplot <- dotboxplot + ggplot2::scale_y_log10()
    violin <- violin + ggplot2::scale_y_log10()
  }

  retlist <- list(
      "violin" = violin,
      "boxplot" = boxplot,
      "dotboxplot" = dotboxplot,
      "density" = density)
  return(retlist)
}

#' Read data from pyprophet and plot columns from it.
#'
#' More proteomics diagnostics!  Now that I am looking more closely, I think
#' this should be folded into plot_pyprophet_distribution().
#'
#' @param pyprophet_data Data from extract_pyprophet_data()
#' @param column Chosen column to plot.
#' @param keep_real FIXME: This should be changed to something like 'data_type'
#'   here and in plot_pyprophet_distribution.
#' @param keep_decoys Do we keep the decoys when plotting the data?
#' @param expt_names Names for the x-axis of the plot.
#' @param label_chars Maximum number of characters before abbreviating sample
#'   names.
#' @param protein chosen protein(s) to plot.
#' @param plot_title Title the plot?
#' @param scale Put the data on a specific scale?
#' @param legend Include the legend?
#' @param order_by Reorder the samples by some factor, presumably condition.
#' @param show_all Skip samples for which no observations were made.
#' @param ... Further arguments, presumably for colors or some such.
#' @return Boxplot describing the desired column from the data.
#' @export
plot_pyprophet_protein <- function(pyprophet_data, column = "intensity", keep_real = TRUE,
                                   keep_decoys = FALSE, expt_names = NULL, label_chars = 10,
                                   protein = NULL, plot_title = NULL, scale = NULL, legend = NULL,
                                   order_by = "condition", show_all = TRUE,
                                   ...) {
  arglist <- list(...)

  metadata <- pyprophet_data[["metadata"]]
  colors <- pyprophet_data[["colors"]]
  sample_data <- pyprophet_data[["sample_data"]]
  plot_df <- data.frame()

  samples <- metadata[["sampleid"]]
  sample_order <- metadata[["sampleid"]]
  if (!is.null(order_by)) {
    new_order <- order(metadata[[order_by]])
    metadata <- metadata[new_order, ]
    sample_order <- metadata[["sampleid"]]
    samples <- metadata[["sampleid"]]
    colors <- colors[new_order]
  }

  ## Reset the sample names if one wants a specific column from the metadata.
  if (!is.null(expt_names) & class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      names(sample_data) <- make.names(metadata[[expt_names]], unique = TRUE)
    } else {
      names(sample_data) <- expt_names
    }
    names(colors) <- names(sample_data)
  }

  slength <- length(samples)
  keepers <- c()
  for (i in 1:slength) {
    name <- samples[i]
    if (class(sample_data[[name]])[1] == "try-error") {
      next
    }
    if (is.null(sample_data[[name]])) {
      next
    }
    keepers <- c(keepers, i)
    ## message("Adding ", name)
    plotted_table <- sample_data[[i]]
    if (!isTRUE(keep_decoys)) {
      good_idx <- plotted_table[["decoy"]] != 1
      plotted_table <- plotted_table[good_idx, ]
    }
    if (!isTRUE(keep_real)) {
      good_idx <- plotted_table[["decoy"]] != 0
      plotted_table <- plotted_table[good_idx, ]
    }
    plotted_data <- as.data.frame(plotted_table[c("sequence", "proteinname",
                                                  "aggr_fragment_annotation", column)])
    plotted_data[["sample"]] <- name
    colnames(plotted_data) <- c("sequence", "proteinname", "fragment", column, "sample")
    ## Re-order the columns because I like sample first.
    plotted_data <- plotted_data[, c("sample", column, "sequence", "proteinname", "fragment")]
    plot_df <- rbind(plot_df, plotted_data)
  }
  ## Drop rows from the metadata and colors which had errors.
  if (length(keepers) > 0) {
    metadata <- metadata[keepers, ]
    colors <- colors[keepers]
  } else {
    stop("Something bad happened to the set of kept samples.")
  }

  ## Now we have a data frame of intensities by sample.
  ## We want them to be ordered as specified in order_by.
  ## Therefore we need to add in rows for those samples which are missing identifications.
  null_samples <- rep(FALSE, length(sample_order))
  names(null_samples) <- sample_order
  data_minimum <- 0
  final_df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(final_df) <- c("sample", column, "sequence", "proteinname", "fragment")
  if (is.null(protein)) {
    stop("This requires a protein ID to search.")
  } else {
    kept_prot_idx <- grepl(pattern = protein, x = plot_df[["proteinname"]])
    plot_df <- plot_df[kept_prot_idx, ]
    data_minimum <- min(as.numeric(plot_df[[column]])) / 2
    for (s in 1:length(sample_order)) {
      sample <- sample_order[s]
      found <- sample %in% plot_df[["sample"]]
      if (isTRUE(found)) {
        rows_idx <- plot_df[["sample"]] == sample
        added_rows <- plot_df[rows_idx, ]
        colnames(added_rows) <- colnames(plot_df)
        final_df <- rbind(final_df, added_rows)
      } else {
        if (isTRUE(show_all)) {
          ## Cut the minimum of the data by half and use it for the blank.
          null_samples[sample] <- TRUE
          blank_row <- c(sample, data_minimum, "", protein, NA)
          final_df <- rbind(final_df, blank_row)
          colnames(final_df) <- colnames(plot_df)
        } ## End checking if we are going to show all samples.
      } ## End checking if a given sample is not in the df at this time.
    } ## End checking if we are to show all samples even if they have no observations.
  } ## End checking if we want to look at a single protein.
  final_df[["sequence"]] <- as.factor(final_df[["sequence"]])
  final_df[[column]] <- as.numeric(final_df[[column]])

  ## Fix the darn colors!
  ## The problem occurs when a sample does not have sufficient identifications
  ## to make a violin. When that happens, the colors get out of sync
  final_df[["colors"]] <- ""
  kept_colors <- c()
  for (sample in names(colors)) {
    sample_color <- colors[sample]
    sample_idx <- final_df[["sample"]] == sample
    sample_sum <- sum(sample_idx)
    if (sample_sum >= 3) {
      kept_colors <- c(kept_colors, sample_color)
    }
    final_df[sample_idx, "colors"] <- sample_color
  }

  scale_data <- check_plot_scale(final_df[[column]], scale = scale)
  ## ...)
  if (is.null(scale)) {
    scale <- scale_data[["scale"]]
  }
  final_df[[column]] <- scale_data[["data"]]

  color_names <- names(colors)
  if (!is.null(label_chars) & is.numeric(label_chars)) {
    final_df[["sample"]] <- abbreviate(final_df[["sample"]], minlength = label_chars)
    color_names <- abbreviate(color_names, minlength = label_chars)
    names(null_samples) <- abbreviate(names(null_samples), minlength = label_chars)
  }

  final_df[["sample"]] <- factor(final_df[["sample"]], levels = color_names)
  observations_by_sample <- table(final_df[["sample"]])
  obs <- as.numeric(observations_by_sample)
  names(obs) <- names(observations_by_sample)
  for (c in 1:length(null_samples)) {
    sample_id <- names(null_samples)[c]
    if (isTRUE(null_samples[sample_id])) {
      obs[sample_id] <- 0
    }
  }
  if (!isTRUE(show_all)) {
    obs_idx <- obs != 0
    obs <- obs[obs_idx]
  }

  sum_obs <- sum(obs)
  violin <- ggplot2::ggplot(data = final_df, ggplot2::aes_string(x = "sample", y = column)) +
    ggplot2::geom_violin(aes_string(fill = "sample"), width = 1, scale = "area") +
    ggplot2::scale_fill_manual(values = as.character(kept_colors)) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::xlab("Sample") +
    ggplot2::ylab(column) +
    ggplot2::geom_jitter(shape = 16, position = ggplot2::position_jitter(0.1),
                         size = 2, alpha = 0.5) +
    ggplot2::annotate("text", x = 1:length(obs),
                      y = max(final_df[[column]] + (0.2 * max(final_df[[column]]))),
                      label = as.character(obs)) +
    ggplot2::labs(caption = glue::glue("Number observed peptides in all samples: {sum_obs}"))
  if (!is.null(plot_title)) {
    violin <- violin + ggplot2::ggtitle(plot_title)
  }

  if (scale == "log") {
    violin <- violin + ggplot2::scale_y_continuous(
                                    labels = scales::scientific,
                                    trans = scales::log2_trans())
  } else if (scale == "logdim") {
    violin <- violin + ggplot2::coord_trans(y = "log2")
  } else if (isTRUE(scale)) {
    violin <- violin + ggplot2::scale_y_log10()
  }
  if (is.null(legend)) {
    violin <- violin + ggplot2::theme(legend.position = "none")
  }

  return(violin)
}

#' Plot some data from the result of extract_pyprophet_data()
#'
#' extract_pyprophet_data() provides a ridiculously large data table of a scored
#' openswath data after processing by pyprophet.
#'
#' @param pyprophet_data List of pyprophet data, one element for each sample,
#'  taken from extract_peprophet_data()
#' @param xaxis Column to plot on the x-axis
#' @param xscale Change the scale of the x-axis?
#' @param yaxis guess!
#' @param yscale Change the scale of the y-axis?
#' @param alpha How see-through to make the dots?
#' @param color_by Change the colors of the points either by sample or condition?
#' @param legend Include a legend of samples?
#' @param sample Which sample(s) to include?
#' @param size_column Use a column for scaling the sizes of dots in the plot?
#' @param rug Add a distribution rug to the axes?
#' @param ... extra options which may be used for plotting.
#' @return a plot!
#' @export
plot_pyprophet_points <- function(pyprophet_data, xaxis = "mass", xscale = NULL, sample = NULL,
                                  yaxis = "leftwidth", yscale = NULL, alpha = 0.4, color_by = "sample",
                                  legend = TRUE, size_column = "mscore", rug = TRUE, ...) {
  arglist <- list(...)

  metadata <- pyprophet_data[["metadata"]]
  colors <- pyprophet_data[["colors"]]
  sample_data <- pyprophet_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)
  keepers <- c()
  for (i in 1:samples) {
    name <- metadata[i, "sampleid"]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    ## message("Adding ", name)
    plotted_table <- sample_data[[i]]
    plotted_table[["sample"]] <- name
    if (is.null(sample)) {
      keepers <- c(keepers, i)
      plot_df <- rbind(plot_df, plotted_table)
    } else {
      if (name %in% sample) {
        keepers <- c(keepers, i)
        plot_df <- rbind(plot_df, plotted_table)
      }
    }
  }
  samples <- length(keepers)

  ## Drop rows from the metadata and colors which had errors.
  metadata <- metadata[keepers, ]
  ## These colors are by condition.
  colors <- colors[keepers]
  if (color_by == "sample") {
    chosen_palette <- "Dark2"
    sample_colors <- sm(
        grDevices::colorRampPalette(
                       RColorBrewer::brewer.pal(samples, chosen_palette))(samples))
  } else {
    sample_colors <- colors
  }

  ## Randomize the rows of the df so we can see if any sample is actually overrepresented
  plot_df <- plot_df[sample(nrow(plot_df)), ]

  if (is.null(plot_df[[xaxis]])) {
    stop("The x axis data seems to be missing.")
  }
  if (is.null(plot_df[[yaxis]])) {
    stop("The y axis data seems to be missing.")
  }

  if (!is.null(xscale)) {
    plot_df[[xaxis]] <- check_plot_scale(plot_df[[xaxis]], scale)[["data"]]
  }
  if (!is.null(yscale)) {
    plot_df[[yaxis]] <- check_plot_scale(plot_df[[yaxis]], scale)[["data"]]
  }

  x_vs_y <- ggplot(data = plot_df, aes_string(x = xaxis, y = yaxis,
                                              fill = "sample", colour = "sample")) +
    ggplot2::geom_point(alpha = alpha, size = 0.5) +
    ggplot2::scale_fill_manual(
                 name = "Sample", values = sample_colors,
                 guide = ggplot2::guide_legend(override.aes = aes(size = 3))) +
    ggplot2::scale_color_manual(
                 name = "Sample", values = sample_colors,
                 guide = ggplot2::guide_legend(override.aes = aes(size = 3))) +
    ggplot2::theme_bw(base_size = base_size)

  if (!is.null(xscale)) {
    x_vs_y <- x_vs_y + ggplot2::scale_x_continuous(trans = scales::log2_trans())
  }
  if (!is.null(yscale)) {
    x_vs_y <- x_vs_y + ggplot2::scale_y_continuous(trans = scales::log2_trans())
  }
  if (isTRUE(lowess)) {
    x_vs_y <- x_vs_y +
      ggplot2::geom_smooth(method = "loess", size = 1.0)
  }
  if (!isTRUE(legend)) {
    x_vs_y <- x_vs_y +
      ggplot2::theme(legend.position = "none")
  }
  if (isTRUE(rug)) {
    x_vs_y <- x_vs_y + ggplot2::geom_rug(colour = "gray50", alpha = alpha)
  }
  retlist <- list(
      "data" = plot_df,
      "plot" = x_vs_y)
  return(retlist)
}

#' Plot some data from the result of extract_peprophet_data()
#'
#' extract_peprophet_data() provides a ridiculously large data table of a comet
#' result after processing by RefreshParser and xinteract/peptideProphet.
#' This table has some 37-ish columns and I am not entirely certain which ones
#' are useful as diagnostics of the data.  I chose a few and made options to
#' pull some/most of the rest.  Lets play!
#'
#' @param table Big honking data table from extract_peprophet_data()
#' @param xaxis Column to plot on the x-axis
#' @param xscale Change the scale of the x-axis?
#' @param yaxis guess!
#' @param yscale Change the scale of the y-axis?
#' @param size_column Use a column for scaling the sizes of dots in the plot?
#' @param ... extra options which may be used for plotting.
#' @return a plot!
#' @export
plot_peprophet_data <- function(table, xaxis = "precursor_neutral_mass", xscale = NULL,
                                yaxis = "num_matched_ions", yscale = NULL,
                                size_column = "prophet_probability", ...) {
  arglist <- list(...)
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["chosen_palette"]])) {
    chosen_palette <- arglist[["chosen_palette"]]
  }
  color_column <- "decoy"
  if (!is.null(arglist[["color_column"]])) {
    color_column <- arglist[["color_column"]]
  }
  if (is.null(table[[color_column]])) {
    table[["color"]] <- "black"
  } else {
    table[["color"]] <- as.factor(table[[color_column]])
  }
  color_list <- NULL
  num_colors <- nlevels(as.factor(table[["color"]]))
  if (num_colors == 2) {
    color_list <- c("darkred", "darkblue")
  } else {
    color_list <- sm(
        grDevices::colorRampPalette(
                       RColorBrewer::brewer.pal(num_colors, chosen_palette))(num_colors))
  }

  if (is.null(table[[xaxis]])) {
    stop(glue("The x-axis column: {xaxis} does not appear in the data."))
  }
  if (is.null(table[[yaxis]])) {
    stop(glue("The y-axis column: {yaxis} does not appear in the data."))
  }

  table <- as.data.frame(table)
  if (is.null(table[[size_column]])) {
    table[["size"]] <- 1
  } else {
    if (class(table[[size_column]]) == "numeric") {
      ## quants <- as.numeric(quantile(unique(table[[size_column]])))
      ## size_values <- c(4, 8, 12, 16, 20)
      ## names(size_values) <- quants
      table[["size"]] <- table[[size_column]]
    } else {
      table[["size"]] <- 1
    }
  }
  ##min_val <- min(table[[size_column]])
  ##max_val <- max(table[[size_column]])
  range <- as.numeric(quantile(unique(table[["size"]])))
  table[table[["size"]] >= range[[5]], "size"] <- "06biggest"
  table[table[["size"]] >= range[[4]] &
        table[["size"]] < range[[5]], "size"] <- "05big"
  table[table[["size"]] >= range[[3]] &
        table[["size"]] < range[[4]], "size"] <- "04medium_big"
  table[table[["size"]] >= range[[2]] &
        table[["size"]] < range[[3]], "size"] <- "03medium_small"
  table[table[["size"]] >= range[[1]] &
        table[["size"]] < range[[2]], "size"] <- "02small"
  table[table[["size"]] < range[[1]], "size"] <- "01smallest"

  ## Setting the factor/vector of sizes is a bit confusing to me.
  table[["size"]] <- as.factor(table[["size"]])
  levels(table[["size"]]) <- c("01smallest", "02small", "03medium_small",
                               "04medium_big", "05big", "06biggest")
  my_sizes <- c("01smallest"=0.4, "02small"=8, "03medium_small"=1.2,
                "04medium_big"=1.6, "05big"=2.0, "06biggest"=2.4)

  scale_x_cont <- "raw"
  if (!is.null(xscale)) {
    if (is.numeric(xscale)) {
      table[[xaxis]] <- log(table[[xaxis]] + 1) / log(xscale)
    } else if (xscale == "log2") {
      scale_x_cont <- "log2"
    } else if (xscale == "log10") {
      scale_x_cont <- "log10"
    } else {
      message("I do not understand your scale.")
    }
  }
  scale_y_cont <- "raw"
  if (!is.null(yscale)) {
    if (is.numeric(yscale)) {
      table[[xaxis]] <- log(table[[yaxis]] + 1) / log(yscale)
    } else if (yscale == "log2") {
      scale_y_cont <- "log2"
    } else if (yscale == "log10") {
      scale_y_cont <- "log10"
    } else {
      message("I do not understand your scale.")
    }
  }

  table[["text"]] <- glue("{table[['protein']]}:{table[['peptide']]}")

  a_plot <- ggplot(data = table, aes_string(x = xaxis, y = yaxis, text = "text",
                                            color = "color", size = "size")) +
    ggplot2::geom_point(alpha = 0.4, aes_string(fill = "color", color = "color")) +
    ggplot2::scale_color_manual(name = "color", values = color_list) +
    ggplot2::geom_rug() +
    ggplot2::scale_size_manual(values = c(0.2, 0.6, 1.0, 1.4, 1.8, 2.2))
  if (scale_x_cont == "log2") {
    a_plot <- ggplot2::scale_x_continuous(trans = scales::log2_trans())
  } else if (scale_x_cont == "log10") {
    a_plot <- ggplot2::scale_x_continuous(trans = scales::log10_trans())
  }
  if (scale_y_cont == "log2") {
    a_plot <- ggplot2::scale_y_continuous(trans = scales::log2_trans())
  } else if (scale_y_cont == "log10") {
    a_plot <- ggplot2::scale_y_continuous(trans = scales::log10_trans())
  }

  return(a_plot)
}

#' Plot the average mass and expected intensity of a set of sequences given an
#' enzyme.
#'
#' This uses the cleaver package to generate a plot of expected intensities
#' vs. weight for a list of protein sequences.
#'
#' @param pep_sequences Set of protein sequences.
#' @param enzyme One of the allowed enzymes for cleaver.
#' @param start Limit the set of fragments from this point
#' @param end to this point.
#' @return List containing the distribution of weights and the associated plot.
#' @export
plot_cleaved <- function(pep_sequences, enzyme = "trypsin", start = 600, end = 1500) {
  products <- cleaver::cleave(pep_sequences, enzym = enzyme)
  ## old_par <- par(pin = c(6,3))
  pep_sizes <- data.frame()
  plot(NA, xlim = c(start, end), ylim = c(0, 1),
       xlab = "mass in Daltons", ylab = "relative intensity",
       main = glue("Digested sequences with: {enzyme}"))
  for (pep in 1:length(products)) {
    seq <- names(pep_sequences)[[pep]]
    prod <- products[[pep]]
    for (i in 1:length(prod)) {
      atoms <- try(BRAIN::getAtomsFromSeq(prod[[i]]), silent = TRUE)
      if (class(atoms) != "try-error") {
        d <- BRAIN::useBRAIN(atoms)
        avg_mass <- d[["avgMass"]]
        most_likely <- max(d[["isoDistr"]])
        lines(avg_mass, most_likely, type = "h", col = 2)
        row <- c(seq, avg_mass, most_likely)
        pep_sizes <- rbind(pep_sizes, row)
      }
    }
  }
  colnames(pep_sizes) <- c("cds", "average_mass", "highest_likelihood")
  plot <- grDevices::recordPlot()
  retlist <- list("plot" = plot, "sizes" = pep_sizes)
  return(retlist)
}

#' Make a histogram of how many peptides are expected at every integer dalton
#' from a given start to end size for a given enzyme digestion.
#'
#' This is very similar to plot_cleaved() above, but tries to be a little bit smarter.
#'
#' @param pep_sequences Protein sequences as per plot_cleaved().
#' @param enzyme Compatible enzyme name from cleaver.
#' @param start Print histogram from here
#' @param end to here.
#' @param color Make the bars this color.
#' @return List containing the plot and size distribution.
#' @export
cleavage_histogram <- function(pep_sequences, enzyme = "trypsin",
                               start = 600, end = 1500, color = "black") {
  products <- cleaver::cleave(pep_sequences, enzym = enzyme)
  prod_df <- as.data.frame(products)
  prod_df <- dplyr::as.tbl(prod_df[, c("group_name", "value")])
  colnames(prod_df) <- c("group_name", "sequence")

  new_df <- prod_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(mass = gather_masses(sequence))

  plot <- ggplot2::ggplot(data = new_df, ggplot2::aes_string(x = "mass")) +
    ggplot2::geom_histogram(binwidth = 1, colour = color) +
    ggplot2::scale_x_continuous(limits = c(start, end))

  retlist <- list(
      "plot" = plot,
      "masses" = new_df)
  return(retlist)
}

## EOF
