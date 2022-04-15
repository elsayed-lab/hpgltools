## plot_distribution.r: A few plots to describe data distributions
## Currently this includes boxplots, density plots, and qq plots.

#' Make a ggplot boxplot of a set of samples.
#'
#' Boxplots and density plots provide complementary views of data distributions.
#' The general idea is that if the box for one sample is significantly shifted
#' from the others, then it is likely an outlier in the same way a density plot
#' shifted is an outlier.
#'
#' @param data Expt or data frame set of samples.
#' @param colors Color scheme, if not provided will make its own.
#' @param plot_title A title!
#' @param order Set the order of boxen.
#' @param violin  Print this as a violin rather than a just box/whiskers?
#' @param scale Whether to log scale the y-axis.
#' @param expt_names Another version of the sample names for printing.
#' @param label_chars  Maximum number of characters for abbreviating sample names.
#' @param ... More parameters are more fun!
#' @return Ggplot2 boxplot of the samples.  Each boxplot
#' contains the following information: a centered line describing the
#' median value of counts of all genes in the sample, a box around the
#' line describing the inner-quartiles around the median (quartiles 2
#' and 3 for those who are counting), a vertical line above/below the
#' box which shows 1.5x the inner quartile range (a common metric of
#' the non-outliers), and single dots for each gene which is outside
#' that range.  A single dot is transparent.
#' @seealso [ggplot2]
#' @examples
#' \dontrun{
#'  a_boxplot <- plot_boxplot(expt)
#'  a_boxplot  ## ooo pretty boxplot look at the lines
#' }
#' @export
plot_boxplot <- function(data, colors = NULL, plot_title = NULL, order = NULL,
                         violin = FALSE, scale = NULL, expt_names = NULL, label_chars = 10,
                         ...) {
  arglist <- list(...)
  data_class <- class(data)[1]
  if (data_class == "expt") {
    design <- pData(data)
    colors <- data[["colors"]]
    data <- as.data.frame(exprs(data))
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
    design <- pData(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    ## some functions prefer matrix, so I am keeping this explicit for the moment
    data <- as.data.frame(data)
  } else {
    stop("This function understands types: expt, ExpressionSet, data.frame, and matrix.")
  }

  ## I am now using this check of the data in a few places, so I function'd it.
  scale_data <- check_plot_scale(data, scale)
  scale <- scale_data[["scale"]]
  data <- scale_data[["data"]]

  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(dim(data)[2])
  }
  data_matrix <- as.matrix(data)
  ## Likely only needed when using quantile norm/batch correction and it sets a value to < 0
  data[data < 0] <- 0

  if (!is.null(expt_names) & class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      colnames(data) <- make.names(design[[expt_names]], unique = TRUE)
    } else {
      colnames(data) <- expt_names
    }
  }
  if (!is.null(label_chars) & is.numeric(label_chars)) {
    colnames(data) <- abbreviate(colnames(data), minlength = label_chars)
  }

  data[["id"]] <- rownames(data)
  dataframe <- reshape2::melt(data, id = c("id"))
  colnames(dataframe) <- c("gene", "sample", "reads")

  dataframe[["sample"]] <- factor(dataframe[["sample"]])
  if (!is.null(order)) {
    if (order == "lexical") {
      ## Order it by sample names lexically
      lexical <- order(levels(dataframe[["sample"]]))
      new_levels <- levels(dataframe[["sample"]])[lexical]
      levels(dataframe[["sample"]]) <- new_levels
    } else {
      new_df <- data.frame()
      for (o in order) {
        matches <- grep(pattern = o, x = dataframe[["sample"]])
        adders <- dataframe[matches, ]
        new_df <- rbind(new_df, adders)
      }
      dataframe <- new_df
      dataframe[["sample"]] <- factor(dataframe[["sample"]], order)
    }
  }

  ## The use of data= and aes() leads to no visible binding for global variable warnings
  ## I am not sure what to do about them in this context.
  boxplot <- ggplot2::ggplot(data = dataframe, aes_string(x = "sample", y = "reads"))
  if (isTRUE(violin)) {
    boxplot <- boxplot +
      ggplot2::geom_violin(aes_string(fill = "sample"), width = 1, scale = "area",
                           show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = as.character(colors), guide = "none") +
      ggplot2::geom_boxplot(aes_string(fill = "sample"), outlier.alpha = 0.01,
                            width = 0.1)
  } else {
    boxplot <- boxplot +
      sm(ggplot2::geom_boxplot(aes_string(fill = "sample"),
                               na.rm = TRUE, fill = colors, size = 0.5,
                               outlier.size = 1.5,
                               guide = "none",
                               outlier.colour = ggplot2::alpha("black", 0.2)))
  }
  boxplot <- boxplot +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::xlab("Sample") +
    ggplot2::ylab("Per-gene (pseudo)count distribution")
  if (!is.null(plot_title)) {
    boxplot <- boxplot + ggplot2::ggtitle(plot_title)
  }

  if (scale == "log") {
    boxplot <- boxplot +
      ggplot2::scale_y_continuous(labels = scales::scientific,
                                  trans = "log2")
  } else if (scale == "logdim") {
    boxplot <- boxplot +
      ggplot2::coord_trans(y = "log2", labels = scales::scientific)
  } else if (isTRUE(scale)) {
    boxplot <- boxplot +
      ggplot2::scale_y_continuous(trans = "log10",
                                  labels = scales::scientific)
  }
  return(boxplot)
}

#' Create a density plot, showing the distribution of each column of data.
#'
#' Density plots and boxplots are cousins and provide very similar views of data
#' distributions. Some people like one, some the other.  I think they are both
#' colorful and fun!
#'
#' @param data Expt, expressionset, or data frame.
#' @param colors Color scheme to use.
#' @param expt_names Names of the samples.
#' @param position How to place the lines, either let them overlap (identity), or stack them.
#' @param direct Use direct.labels for labeling the plot?
#' @param fill Fill the distributions?  This might make the plot unreasonably colorful.
#' @param plot_title Title for the plot.
#' @param scale Plot on the log scale?
#' @param colors_by Factor for coloring the lines
#' @param label_chars Maximum number of characters in sample names before abbreviation.
#' @param ... sometimes extra arguments might come from graph_metrics()
#' @return ggplot2 density plot!
#' @seealso [ggplot2]
#' @examples
#' \dontrun{
#'  funkytown <- plot_density(data)
#' }
#' @export
plot_density <- function(data, colors = NULL, expt_names = NULL, position = "identity", direct = TRUE,
                         fill = NULL, plot_title = NULL, scale = NULL, colors_by = "condition",
                         label_chars = 10, ...) {
  ## also position='stack'
  data_class <- class(data)[1]
  design <- NULL
  if (data_class == "expt") {
    design <- pData(data)
    colors <- data[["colors"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
    design <- pData(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    ## some functions prefer matrix, so I am keeping this explicit for the moment
    data <- as.matrix(data)
  } else {
    stop("This function understands types: expt, ExpressionSet, data.frame, and matrix.")
  }

  if (is.null(scale)) {
    if (max(data) > 10000) {
      mesg("This data will benefit from being displayed on the log scale.")
      mesg("If this is not desired, set scale='raw'")
      scale <- "log"
      negative_idx <- data < 0
      if (sum(negative_idx) > 0) {
        mesg("Some data are negative.  We are on a log scale, setting them to 0.5.")
        data[negative_idx] <- 0.5
        message("Changed ", sum(negative_idx), " negative features.")
      }
      zero_idx <- data == 0
      if (sum(zero_idx) > 0) {
        mesg("Some entries are 0.  We are on a log scale, setting them to 0.5.")
        data[zero_idx] <- 0.5
        message("Changed ", sum(zero_idx), " zero count features.")
      }
    } else {
      scale <- "raw"
    }
  }

  if (!is.null(expt_names)) {
    if (class(expt_names) == "character" & length(expt_names) == 1) {
      ## Then this refers to an experimental metadata column.
      colnames(data) <- design[[expt_names]]
    } else {
      colnames(data) <- expt_names
    }
  }
  if (!is.null(label_chars) & is.numeric(label_chars)) {
    colnames(data) <- abbreviate(colnames(data), minlength = label_chars)
  }

  ## If the columns lose the connectivity between the sample and values, then
  ## the ggplot below will fail with env missing.
  melted <- data.table::as.data.table(reshape2::melt(data))
  if (dim(melted)[2] == 3) {
    colnames(melted) <- c("id", "sample", "counts")
  } else if (dim(melted)[2] == 2) {
    colnames(melted) <- c("sample", "counts")
  } else {
    stop("Could not properly melt the data.")
  }
  densityplot <- NULL
  if (is.null(fill)) {
    densityplot <- ggplot2::ggplot(data = melted,
                                   ggplot2::aes_string(x = "counts", colour = "sample"))
  } else {
    fill <- "sample"
    densityplot <- ggplot2::ggplot(data = melted,
                                   ggplot2::aes_string(x = "counts", colour = "sample", fill = "fill"))
  }
  densityplot <- densityplot +
    ggplot2::geom_density(ggplot2::aes_string(x = "counts", y = "..count..", fill = "sample"),
                          position = position, na.rm = TRUE) +
    ggplot2::ylab("Number of genes.") + ggplot2::xlab("Number of hits/gene.") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   legend.key.size = ggplot2::unit(0.3, "cm"))
  if (!is.null(plot_title)) {
    densityplot <- densityplot + ggplot2::ggtitle(plot_title)
  }

  if (scale == "log") {
    densityplot <- densityplot + ggplot2::scale_x_continuous(trans = scales::log2_trans(),
                                                             labels = scales::scientific)
  } else if (scale == "logdim") {
    densityplot <- densityplot +
      ggplot2::coord_trans(x = "log2") +
      ggplot2::scale_x_continuous(labels = scales::scientific)
  } else if (isTRUE(scale)) {
    densityplot <- densityplot +
      ggplot2::scale_x_log10(labels = scales::scientific)
  }

  if (!is.null(colors_by)) {
    densityplot <- densityplot + ggplot2::scale_colour_manual(values = as.character(colors)) +
      ggplot2::scale_fill_manual(values = ggplot2::alpha(as.character(colors), 0.1))
  }

  if (isTRUE(direct)) {
    densityplot <- directlabels::direct.label(densityplot)
  }

  condition_summary <- data.table::data.table()
  batch_summary <- data.table::data.table()
  counts <- NULL
  if (!is.null(design)) {
    if (!is.null(design[["condition"]])) {
      melted[, "condition" := design[sample, "condition"]]
      condition_summary <- data.table::setDT(melted)[, list("min"=min(counts),
                                                            "1st"=quantile(x = counts, probs = 0.25),
                                                            "median"=median(x = counts),
                                                            "mean"=mean(counts),
                                                            "3rd"=quantile(x = counts, probs = 0.75),
                                                            "max"=max(counts)),
                                                     by = "condition"]
    }
    if (!is.null(design[["batch"]])) {
      melted[, "batch" := design[sample, "batch"]]
      batch_summary <- data.table::setDT(melted)[, list("min"=min(counts),
                                                        "1st"=quantile(x = counts, probs = 0.25),
                                                        "median"=median(x = counts),
                                                        "mean"=mean(counts),
                                                        "3rd"=quantile(x = counts, probs = 0.75),
                                                        "max"=max(counts)),
                                                 by = "batch"]
    }
  }

  sample_summary <- data.table::setDT(melted)[, list("min"=min(counts),
                                                     "1st"=quantile(x = counts, probs = 0.25),
                                                     "median"=median(x = counts),
                                                     "mean"=mean(counts),
                                                     "3rd"=quantile(x = counts, probs = 0.75),
                                                     "max"=max(counts)),
                                              by = "sample"]
  retlist <- list(
      "plot" = densityplot,
      "condition_summary" = condition_summary,
      "batch_summary" = batch_summary,
      "sample_summary" = sample_summary,
      "table" = melted)
  return(retlist)
}

#' Quantile/quantile comparison of the mean of all samples vs. each sample.
#'
#' This allows one to visualize all individual data columns against the mean of
#' all columns of data in order to see if any one is significantly different
#' than the cloud.
#'
#' @param data Expressionset, expt, or dataframe of samples.
#' @param labels What kind of labels to print?
#' @param ... Arguments passed presumably from graph_metrics().
#' @return List containing:
#'  logs = a recordPlot() of the pairwise log qq plots.
#'  ratios = a recordPlot() of the pairwise ratio qq plots.
#'  means = a table of the median values of all the summaries of the qq plots.
#' @seealso [Biobase]
#' @export
plot_qq_all <- function(data, labels = "short", ...) {
  arglist <- list(...)
  data_class <- class(data)[1]
  if (data_class == "expt") {
    design <- pData(data)
    colors <- data[["colors"]]
    data <- as.data.frame(exprs(data))
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
    design <- pData(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    data <- as.data.frame(data)
  } else {
    stop("This understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  sample_data <- data[, c(1, 2)]
  means <- rowMeans(data)
  sample_data[["mean"]] <- means
  logs <- list()
  ratios <- list()
  means <- list()
  comparisons <- length(colnames(data))
  row_columns <- ceiling(sqrt(comparisons))
  rows <- nrow(data)
  ## I want to make a square containing the graphs.
  count <- 1
  for (i in 1:comparisons) {
    ith <- colnames(data)[i]
    message("Making plot of ", ith, "(", i, ") vs. a sample distribution.")
    tmpdf <- data.frame("ith"=data[, i],
                        "mean"=sample_data[["mean"]])
    colnames(tmpdf) <- c(ith, "mean")
    tmpqq <- plot_single_qq(tmpdf, x = 1, y = 2, labels = labels)
    logs[[count]] <- tmpqq[["log"]]
    ratios[[count]] <- tmpqq[["ratio"]]
    means[[count]] <- tmpqq[["summary"]][["Median"]]
    count <- count + 1
  }

  tmp_file <- tempfile(pattern = "multi", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  result <- plot_multiplot(logs)
  log_plots <- grDevices::recordPlot()
  dev.off()
  file.remove(tmp_file)

  tmp_file <- tempfile(pattern = "multi", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  plot_multiplot(ratios)
  ratio_plots <- grDevices::recordPlot()
  dev.off()
  file.remove(tmp_file)

  plots <- list(logs = log_plots, ratios = ratio_plots, medians = means)
  return(plots)
}

#' Perform a qqplot between two columns of a matrix.
#'
#' Given two columns of data, how well do the distributions match one another?
#' The answer to that question may be visualized through a qq plot!
#'
#' @param data Data frame/expt/expressionset.
#' @param x First column to compare.
#' @param y Second column to compare.
#' @param labels Include the lables?
#' @return a list of the logs, ratios, and mean between the plots as ggplots.
#' @seealso [Biobase]
#' @export
plot_single_qq <- function(data, x = 1, y = 2, labels = TRUE) {
  data_class <- class(data)[1]
  if (data_class == "expt") {
    design <- pData(data)
    colors <- data[["colors"]]
    names <- data[["names"]]
    data <- as.data.frame(exprs(data))
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
    design <- pData(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    data <- as.data.frame(data)
  } else {
    stop("This understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  xlabel <- colnames(data)[x]
  ylabel <- colnames(data)[y]
  xvector <- as.vector(data[, x])
  yvector <- as.vector(data[, y])
  ratio_df <- data.frame("sorted_x" = sort(xvector),
                         "sorted_y" = sort(yvector))
  ratio_df[["ratio"]] <- ratio_df[["sorted_x"]] / ratio_df[["sorted_y"]]
  ## Drop elements with zero
  ## First by checking for NAN, then Inf.
  na_idx <- is.na(ratio_df[["ratio"]])
  ratio_df <- ratio_df[!na_idx, ]
  nan_idx <- is.infinite(ratio_df[["ratio"]])
  ratio_df <- ratio_df[!nan_idx, ]
  ratio_df[["increment"]] <- as.vector(1:nrow(ratio_df))

  if (labels == "short") {
    y_string <- glue("{xlabel} : {ylabel}")
  } else {
    y_string <- glue("Ratio of sorted {xlabel}  and {ylabel}.")
  }
  ratio_plot <- ggplot2::ggplot(ratio_df,
                                ggplot2::aes_string(x = "increment", y = "ratio")) +
    ggplot2::geom_point(colour = sm(grDevices::densCols(ratio_df[["ratio"]])), stat = "identity",
                        size = 1, alpha = 0.2, na.rm = TRUE) +
    ggplot2::scale_y_continuous(limits = c(0, 2))
  if (isTRUE(labels)) {
    ratio_plot <- ratio_plot +
      ggplot2::xlab("Sorted gene") +
      ggplot2::ylab(y_string) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(legend.position = "none")
  } else if (labels == "short") {
    ratio_plot <- ratio_plot +
      ggplot2::ylab(y_string) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     legend.position = "none",
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank())
  } else {
    ratio_plot <- ratio_plot +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     legend.position = "none",
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank())
  }

  log_df <- data.frame(cbind(log(ratio_df[["sorted_x"]] + 1.5),
                             log(ratio_df[["sorted_y"]] + 1.5)))
  gg_max <- max(log_df)
  colnames(log_df) <- c(xlabel, ylabel)
  log_df[["sub"]] <- log_df[, 1] - log_df[, 2]
  sorted_x <- as.vector(log_df[, 1])
  log_ratio_plot <- ggplot2::ggplot(log_df,
                                    ggplot2::aes_string(x = "get(xlabel)", y = "get(ylabel)")) +
    ggplot2::geom_point(colour = grDevices::densCols(x = sorted_x), stat = "identity") +
    ggplot2::scale_y_continuous(limits = c(0, gg_max)) +
    ggplot2::scale_x_continuous(limits = c(0, gg_max))
  if (isTRUE(labels)) {
    log_ratio_plot <- log_ratio_plot +
      ggplot2::xlab(glue("log sorted {xlabel}")) +
      ggplot2::ylab(glue("log sorted {ylabel}")) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(legend.position = "none")
  } else if (labels == "short") {
    log_ratio_plot <- log_ratio_plot +
      ggplot2::xlab("gene") +
      ggplot2::ylab(y_string) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     legend.position = "none",
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank())
  } else {
    log_ratio_plot <- log_ratio_plot +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     legend.position = "none",
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank())
  }
  log_summary <- summary(log_df[["sub"]])
  qq_plots <- list(
      "ratio" = ratio_plot,
      "log" = log_ratio_plot,
      "summary" = log_summary)
  return(qq_plots)
}

#' Plot the representation of the top-n genes in the total counts / sample.
#'
#' One question we might ask is: how much do the most abundant genes in a
#' samples comprise the entire sample?  This plot attempts to provide a visual
#' hint toward answering this question.  It does so by rank-ordering all the
#' genes in every sample and dividing their counts by the total number of reads
#' in that sample.  It then smooths the points to provide the resulting trend.
#' The steeper the resulting line, the more over-represented these top-n genes
#' are.  I suspect, but haven't tried yet, that the inflection point of the
#' resulting curve is also a useful diagnostic in this question.
#'
#' @param data Dataframe/matrix/whatever for performing topn-plot.
#' @param plot_title A title for the plot.
#' @param num The N in top-n genes, if null, do them all.
#' @param expt_names Column or character list of sample names.
#' @param plot_labels Method for labelling the lines.
#' @param label_chars Maximum number of characters before abbreviating samples.
#' @param plot_legend Add a legend to the plot?
#' @param ... Extra arguments, currently unused.
#' @return List containing the ggplot2
#' @export
plot_topn <- function(data, plot_title = NULL, num = 100, expt_names = NULL,
                      plot_labels = "direct", label_chars = 10, plot_legend = FALSE, ...) {
  arglist <- list(...)
  data_class <- class(data)
  if (data_class == "expt") {
    design <- pData(data)
    colors <- data[["colors"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
    design <- pData(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    data <- as.matrix(data)
  } else {
    stop("This understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  columns <- colSums(data)
  testing <- data / columns
  newdf <- data.frame(row.names = 1:nrow(testing))
  for (col in colnames(testing)) {
    ranked <- order(testing[, col], decreasing = TRUE)
    newdf[, col] <- testing[ranked, col]
  }

  if (num > 0) {
    newdf <- head(newdf, n = num)
  }
  if (num < 0) {
    newdf <- tail(newdf, n=-1 * num)
  }
  newdf[["rank"]] <- rownames(newdf)

  ## Choose the smoothing algorithm
  smoother <- "loess"
  if (is.null(arglist[["smoother"]])) {
    if (nrow(newdf) > 5000) {
      smoother <- "gam"
    }
  } else {
    smoother <- arglist[["smoother"]]
  }

  if (!is.null(expt_names) & class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      colnames(newdf) <- make.names(design[[expt_names]], unique = TRUE)
    } else {
      colnames(newdf) <- expt_names
    }
    colnames(newdf)[ncol(newdf)] <- "rank"
  }
  if (!is.null(label_chars) & is.numeric(label_chars)) {
    colnames(newdf) <- abbreviate(colnames(newdf), minlength = label_chars)
  }

  tmpdf <- reshape2::melt(newdf, id.vars = "rank")
  colnames(tmpdf) <- c("rank", "sample", "pct")
  tmpdf[["rank"]] <- as.numeric(tmpdf[["rank"]])
  tmpdf[["value"]] <- as.numeric(tmpdf[["pct"]])
  tmpdf[["sample"]] <- as.factor(tmpdf[["sample"]])

  topn_plot <- ggplot(tmpdf, aes_string(x = "rank", y = "pct", color = "sample")) +
    ggplot2::geom_smooth(method = smoother, level = 0.5, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = base_size)

  if (!is.null(plot_title)) {
    topn_plot <- topn_plot +
      ggplot2::ggtitle(plot_title)
  }

  if (plot_labels == "direct") {
    topn_plot <- topn_plot +
      directlabels::geom_dl(aes_string(label = "sample"), method = "smart.grid")
  }
  if (isFALSE(plot_legend)) {
    topn_plot <- topn_plot +
      ggplot2::theme(legend.position = "none")
  }

  retlist <- list(
      "plot" = topn_plot,
      "table" = tmpdf)
  return(retlist)
}

#' Look at the (biological)coefficient of variation/quartile coefficient of dispersion
#' with respect to an experimental factor.
#'
#' I want to look at the (B)CV of some data with respect to condition/batch/whatever.
#' This function should make that possible, with some important caveats.  The
#' most appropriate metric  is actually the biological coefficient of variation
#' as calculated by DESeq2/EdgeR; but the metrics I am currently taking are the
#' simpler and less appropriate CV(sd/mean) and QCD(q3-q1/q3+q1).
#'
#' @param data Expressionset/epxt to poke at.
#' @param x_axis Factor in the experimental design we may use to group the data
#'  and calculate the dispersion metrics.
#' @param colors Set of colors to use when making the violins
#' @param plot_title Optional title to include with the plot.
#' @param ... Extra arguments to pass along.
#' @return List of plots showing the coefficients vs. genes along with the data.
#' @export
plot_variance_coefficients <- function(data, x_axis = "condition", colors = NULL,
                                       plot_title = NULL, ...) {
  arglist <- list(...)
  plot_legend <- FALSE
  if (!is.null(arglist[["plot_legend"]])) {
    plot_legend <- arglist[["plot_legend"]]
  }

  data_class <- class(data)
  if (data_class == "expt") {
    design <- pData(data)

    colors <- data[["colors"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
    design <- pData(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    data <- as.matrix(data)
  } else {
    stop("This understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  melted <- data.table::as.data.table(reshape2::melt(data))
  if (ncol(melted) == 3) {
    colnames(melted) <- c("gene", "sample", "exprs")
  } else if (dim(melted)[2] == 2) {
    colnames(melted) <- c("sample", "exprs")
  } else {
    stop("Could not properly melt the data.")
  }

  for (add in x_axis) {
    if (add %in% colnames(design)) {
      tmp_df <- as.data.frame(design[[add]])
      colnames(tmp_df) <- add
      tmp_df[["sample"]] <- rownames(design)
      rownames(tmp_df) <- rownames(design)
      melted <- merge(melted, tmp_df, by = "sample")
    } else if (add %in% colnames(melted)) {
      message("Skipping variable: ", add, " because it is in the data.")
    } else {
      stop("Could not find the metadata variable: ", add)
    }
  }

  ## The various forms of evaluation in the hadleyverse is getting ridiculous.
  .data <- NULL
  message("Naively calculating coefficient of variation/dispersion with respect to ",
          x_axis, ".")
  cv_data <- melted %>%
    dplyr::group_by(.data[["gene"]], .data[[x_axis]]) %>%
    dplyr::summarize(
               "mean_exprs" = mean(.data[["exprs"]], na.rm = TRUE),
               "sd_exprs" = sd(.data[["exprs"]], na.rm = TRUE),
               "q1" = quantile(.data[["exprs"]], probs = 0.25),
               "q3" = quantile(.data[["exprs"]], probs = 0.75))
  cv_data[["cv"]] <- cv_data[["sd_exprs"]] / cv_data[["mean_exprs"]]
  cv_data[["disp"]] <- (cv_data[["q3"]] - cv_data[["q1"]]) / (cv_data[["q3"]] + cv_data[["q1"]])
  na_idx <- is.na(cv_data[["cv"]])
  cv_data[na_idx, "cv"] <- 0
  na_idx <- is.na(cv_data[["disp"]])
  cv_data[na_idx, "disp"] <- 0
  ## The metrics of dispersion taken so far are not really appropriate for
  ## RNASeq distributed data. Ideally, I would like to subset the expressionset
  ## according to the x_axis factor and perform a DESeq2/edgeR dispersion
  ## estimate for the remaining pile of data, then add the results to cv_data.
  ## The following piece of code is a simple way to get the normal, pooled
  ## dispersion information. In theory I should be able to refactor this to do
  ## what I want.
  ## message("Using edgeR to calculate dispersions with respect to: ", x_axis)
  ## test <- import_edger(data, design[[x_axis]])
  ## disp_model <- model.matrix(object = as.formula(paste0("~", x_axis)), data = design)
  ## disp_data <- edgeR::estimateDisp(test, design = disp_model, group = design[[x_axis]])
  ## ## The problem here is that the tagwise.dispersions do not have names, I think we
  ## ## we can assume that they are in the same order as the original rownames.
  ## cv_data[["bcv"]] <- disp_data[["tagwise.dispersion"]]
  message("Finished calculating dispersion estimates.")

  color_list <- NULL
  if (x_axis == "condition") {
    names(colors) <- design[, "condition"]
    color_list <- colors
  } else if (!is.null(arglist[["colors"]])) {
    color_list <- arglist[["colors"]]
  } else {
    num_colors <- length(levels(as.factor(cv_data[[x_axis]])))
    color_list <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(num_colors)
    names(color_list) <- levels(cv_data[[x_axis]])
  }

  get_mean_cv <- function(x) {
    data.frame("y" = mean(x),
               "label" = signif(mean(x, na.rm = TRUE), digits = 2))
  }

  ## Add the number of samples of each type to the top of the plot with this.
  sample_numbers <- list()
  for (l in levels(as.factor(cv_data[[x_axis]]))) {
    sample_numbers[[l]] <- sum(design[[x_axis]] == l)
  }
  cv_data[["x_axis"]] <- cv_data[[x_axis]]
  y_labels <- list(
      "bcv" = "Biological coefficient of variation",
      "cv" = "Coefficient of variation",
      "disp" = "Quartile coefficient of dispersion")
  retlst <- list()
  for (type in c("cv", "disp")) {
    retlst[[type]] <- ggplot(cv_data, aes_string(x = "x_axis", y = type)) +
      ggplot2::geom_violin(aes_string(fill = "x_axis"), width = 1, scale = "area") +
      ggplot2::scale_fill_manual(values = color_list, name = x_axis) +
      ggplot2::stat_summary(fun.data = get_mean_cv, geom = "text",
                            position = ggplot2::position_fill(vjust=-0.1)) +
      ggplot2::annotate("text", x = 1:length(sample_numbers),
                        y = max(cv_data[[type]] + (0.05 * max(cv_data[[type]]))),
                        label = as.character(sample_numbers)) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::ylab(as.character(y_labels[type])) +
      ggplot2::xlab("")
    if (!is.null(plot_title)) {
      retlst[[type]] <- retlst[[type]] +
        ggplot2::ggtitle(plot_title)
    }
    if (isFALSE(plot_legend)) {
      retlst[[type]] <- retlst[[type]] +
        ggplot2::theme(legend.position = "none")
    }
  }
  retlst[["data"]] <- cv_data
  retlst[["plot"]] <- retlst[["cv"]]
  return(retlst)
}

## EOF
