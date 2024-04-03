## plot_dotplot.r: Dotplots in various contexts, currently just smc/smd

#' Make a dotplot of some categorised factors and a set of SVs (for other factors).
#'
#' This should make a quick df of the factors and surrogates and plot them.
#'
#' @param expt Experiment from which to acquire the design, counts, etc.
#' @param svest Set of surrogate variable estimations from sva/svg
#'  or batch estimates.
#' @param sv Which surrogate to plot?
#' @param chosen_factor Factor to compare against.
#' @param factor_type This may be a factor or range, it is intended to plot
#'  a scatterplot if it is a range, a dotplot if a factor.
#' @return surrogate variable plot as per Leek's work
#' @seealso [ggplot2]
#' @examples
#' \dontrun{
#'  estimate_vs_snps <- plot_svfactor(start, surrogate_estimate, "snpcategory")
#' }
#' @export
plot_svfactor <- function(expt, svest, sv = 1, chosen_factor = "batch", factor_type = "factor") {
  meta <- pData(expt)
  chosen <- meta[[chosen_factor]]
  names <- sampleNames(expt)
  my_colors <- colors(expt)
  sv_df <- data.frame(
    "adjust" = svest[, sv],  ## Take a single estimate from compare_estimates()
    "factors" = chosen,
    "samplenames" = names)

  sv_melted <- reshape2::melt(sv_df, idvars = "factors")
  minval <- min(sv_df[["adjust"]])
  maxval <- max(sv_df[["adjust"]])
  my_binwidth <- (maxval - minval) / 40
  sv_plot <- ggplot(sv_melted, aes(x = .data[["factors"]], y = .data[["value"]])) +
    ggplot2::geom_dotplot(binwidth = my_binwidth, binaxis = "y",
                          stackdir = "center", binpositions = "all",
                          colour = "black", fill = my_colors) +
    ggplot2::xlab(glue("Experimental factor: {chosen_factor}")) +
    ggplot2::ylab("1st surrogate variable estimation") +
    ggplot2::geom_text(aes(x = .data[["factors"]], y = .data[["value"]],
                           label = .data[["samplenames"]]),
                       angle = 45, size = 3, vjust = 2) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black")) +
    ggplot2::theme_bw(base_size = base_size)
  return(sv_plot)
}

#' Make a dotplot of known batches vs. SVs.
#'
#' This should make a quick df of the factors and surrogates and plot them.
#' Maybe it should be folded into plot_svfactor?  Hmm, I think first I will
#' write this and see if it is better.
#'
#' @param expt Experiment from which to acquire the design, counts, etc.
#' @param svs Set of surrogate variable estimations from sva/svg
#'  or batch estimates.
#' @param sv Which surrogate variable to show?
#' @param batch_column Which experimental design column to use?
#' @param factor_type This may be a factor or range, it is intended to plot
#'  a scatterplot if it is a range, a dotplot if a factor.
#' @param id_column Use this column for the sample IDs.
#' @return Plot of batch vs surrogate variables as per Leek's work.
#' @seealso [sva] [ggplot2]
#' @examples
#' \dontrun{
#'  estimate_vs_snps <- plot_batchsv(start, surrogate_estimate, "snpcategory")
#' }
#' @export
plot_batchsv <- function(expt, svs, sv = 1, batch_column = "batch", factor_type = "factor",
                         id_column = "sampleid") {
  meta <- pData(expt)
  chosen <- meta[[batch_column]]
  names(chosen) <- sampleNames(expt)
  num_batches <- length(unique(chosen))
  samples <- meta[[id_column]]
  if (is.null(samples)) {
    samples <- rownames(meta)
  }

  factor_df <- data.frame(
    "sample" = samples,
    "factor" = as.integer(as.factor(pData(expt)[[batch_column]])),
    "fill" = expt[["colors"]],
    "condition" = expt[["conditions"]],
    "batch" = expt[["batches"]],
    "color" = "black",
    "svs" = svs[, sv])
  if (num_batches <= 5) {
    factor_df[["shape"]] <- 20 + as.numeric(as.factor(factor_df[["batch"]]))
  } else {
    factor_df[["shape"]] <- 21
  }
  factor_df[["shape"]] <- as.factor(factor_df[["shape"]])

  color_list <- as.character(factor_df[["fill"]])
  names(color_list) <- as.character(factor_df[["condition"]])

  sample_factor <- ggplot(factor_df,
                          aes(x = .data[["sample"]],
                              y = .data[["factor"]],
                              shape = .data[["batch"]],
                              fill = .data[["condition"]])) +
    ggplot2::geom_point(size = 5,
                        aes(shape = .data[["batch"]],
                            colour = .data[["condition"]],
                            fill = .data[["condition"]])) +
    ggplot2::geom_point(size = 5,
                        colour = "black",
                        show.legend = FALSE,
                        aes(shape = .data[["batch"]],
                            fill = .data[["condition"]])) +
    ggplot2::scale_shape_manual(name = "Batch",
                                labels = levels(as.factor(factor_df[["batch"]])),
                                guide = ggplot2::guide_legend(
                                  override.aes = list(size = 5, fill = "grey")),
                                values = 21:25) +
    ggplot2::scale_color_manual(name = "Condition",
                                guide = "legend",
                                values = color_list) +
    ggplot2::scale_fill_manual(name = "Condition",
                               guide = "legend",
                               values = color_list) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  ##, hjust = 1.5, vjust = 0.5))
  factor_svs <- ggplot(factor_df,
                       aes(x = .data[["factor"]],
                           y = .data[["svs"]],
                           fill = .data[["condition"]],
                           colour = .data[["condition"]],
                           shape = .data[["shape"]])) +
    ggplot2::geom_point(size = 5,
                        aes(shape = .data[["batch"]],
                            colour = .data[["condition"]],
                            fill = .data[["condition"]])) +
    ggplot2::geom_point(size = 5,
                        colour = "black",
                        show.legend = FALSE,
                        aes(shape = .data[["batch"]],
                            fill = .data[["condition"]])) +
    ggplot2::scale_shape_manual(name = "Batch",
                                labels = levels(as.factor(factor_df[["batch"]])),
                                guide = ggplot2::guide_legend(
                                  override.aes = list(size = 5, fill = "grey")),
                                values = 21:25) +
    ggplot2::scale_color_manual(name = "Condition",
                                guide = "legend",
                                values = color_list) +
    ggplot2::scale_fill_manual(name = "Condition",
                               guide = "legend",
                               values = color_list) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  ##, hjust = 1.5, vjust = 0.5))

  svs_sample <- ggplot(factor_df,
                       aes(x = .data[["sample"]],
                           y = .data[["svs"]],
                           fill = .data[["condition"]],
                           colour = .data[["condition"]],
                           shape = .data[["shape"]])) +
    ggplot2::geom_point(size = 5,
                        aes(shape = .data[["batch"]],
                            colour = .data[["condition"]],
                            fill = .data[["condition"]])) +
    ggplot2::geom_point(size = 5,
                        colour = "black",
                        show.legend = FALSE,
                        aes(shape = .data[["batch"]],
                            fill = .data[["condition"]])) +
    ggplot2::scale_shape_manual(name = "Batch",
                                labels = levels(as.factor(factor_df[["batch"]])),
                                guide = ggplot2::guide_legend(
                                  override.aes = list(size = 5, fill = "grey")),
                                values = 21:25) +
    ggplot2::scale_color_manual(name = "Condition",
                                guide = "legend",
                                values = color_list) +
    ggplot2::scale_fill_manual(name = "Condition",
                               guide = "legend",
                               values = color_list) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  ## An alternate possibility:  hjust = 1.5, vjust = 0.5))

  plots <- list(
    "sample_factor" = sample_factor,
    "factor_svs" = factor_svs,
    "svs_sample" = svs_sample)
  class(plots) <- "sv_plots"
  return(plots)
}

#' make a dotplot of some categorised factors and a set of principle components.
#'
#' This should make a quick df of the factors and PCs and plot them.
#'
#' @param pc_df Df of principle components.
#' @param expt Expt containing counts, metadata, etc.
#' @param exp_factor Experimental factor to compare against.
#' @param component Which principal component to compare against?
#' @return Plot of principle component vs factors in the data
#' @seealso [ggplot2]
#' @examples
#' \dontrun{
#'  estimate_vs_pcs <- plot_pcfactor(pcs, times)
#' }
#' @export
plot_pcfactor <- function(pc_df, expt, exp_factor = "condition", component = "PC1") {
  meta <- pData(expt)
  samplenames <- sampleNames(expt)
  my_colors <- colors(expt)
  minval <- min(pc_df[[component]])
  maxval <- max(pc_df[[component]])
  my_binwidth <- (maxval - minval) / 40
  sv_plot <- ggplot(pc_df, aes(x = .data[[exp_factor]], y = .data[[component]])) +
    ggplot2::geom_dotplot(binwidth = my_binwidth, binaxis = "y",
                          stackdir = "center", binpositions = "all",
                          colour = "black", fill = my_colors) +
    ggplot2::xlab(glue("Experimental factor: {exp_factor}")) +
    ggplot2::ylab("1st surrogate variable estimation") +
    ##ggplot2::geom_text(
    ##           ggplot2::aes(x = .data[[exp_factor]], y = .data[[component]], label = "strains"),
    ##           angle = 45, size = 3, vjust = 2) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  return(sv_plot)
}

#' Make an R plot of the standard median correlation or distance among samples.
#'
#' This was written by a mix of Kwame Okrah <kokrah at gmail dot com>, Laura
#' Dillon <dillonl at umd dot edu>, and Hector Corrada Bravo <hcorrada at umd dot edu>
#' I reimplemented it using ggplot2 and tried to make it a little more flexible.
#' The general idea is to take the pairwise correlations/distances of the
#' samples, then take the medians, and plot them.  This version of the plot is
#' no longer actually a dotplot, but a point plot, but who is counting?
#'
#' @param data Expt, expressionset, or data frame.
#' @param design Specify metadata if desired.
#' @param colors Color scheme if data is not an expt.
#' @param method Correlation or distance method to use.
#' @param plot_legend Include a legend on the side?
#' @param expt_names Use pretty names for the samples?
#' @param label_chars Maximum number of characters before abbreviating sample names.
#' @param plot_title Title for the graph.
#' @param dot_size How large should the glyphs be?
#' @param ... More parameters to make you happy!
#' @return ggplot of the standard median something
#'  among the samples.  This will also write to an
#'  open device.  The resulting plot measures the median correlation of
#'  each sample among its peers.  It notes 1.5* the interquartile range
#'  among the samples and makes a horizontal line at that correlation
#'  coefficient.  Any sample which falls below this line is considered
#'  for removal because it is much less similar to all of its peers.
#' @seealso [matrixStats] [ggplot2]
#' @examples
#' \dontrun{
#'  smc_plot = hpgl_smc(expt = expt)
#' }
#' @export
plot_sm <- function(data, design = NULL, colors = NULL, method = "pearson", plot_legend = FALSE,
                    expt_names = NULL, label_chars = 10, plot_title = NULL, dot_size = 5,
                    ...) {
  arglist <- list(...)
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }

  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(ncol(data), chosen_palette))(ncol(data))
  }
  colors <- as.character(colors)

  properties <- NULL
  if (method == "pearson" || method == "spearman" || method == "robust") {
    mesg("Performing ", method, " correlation.")
    properties <- hpgl_cor(data, method = method)
  } else {
    mesg("Performing ", method, " distance.")
    properties <- as.matrix(dist(t(data), method = method))
  }

  prop_median <- matrixStats::rowMedians(properties)
  prop_spread <- stats::quantile(prop_median, p = c(1, 3) / 4)
  prop_iqr <- diff(prop_spread)
  log_scale <- FALSE
  outer_limit <- NULL
  ylimit <- NULL
  type <- "unknown"
  if (method == "pearson" || method == "spearman" || method == "robust") {
    outer_limit <- prop_spread[1] - (1.5 * prop_iqr)
    ylimit <- c(pmin(min(prop_median), outer_limit), max(prop_median))
    ylimit <- ylimit[[1]]
    type <- "correlation"
  } else {
    log_scale <- TRUE
    outer_limit <- prop_spread[2] + (1.5 * prop_iqr)
    ylimit <- c(pmin(min(prop_median), outer_limit), max(prop_median))
    ylimit <- ylimit[[2]]
    type <- "distance"
  }

  sm_df <- data.frame(
    "sample" = rownames(properties),
    "sm" = prop_median,
    "color" = colors)

  if (!is.null(design)) {
    if (!is.null(design[["condition"]])) {
      sm_df[["condition"]] <- as.factor(design[["condition"]])
    } else {
      sm_df[["condition"]] <- "undefined"
    }
    if (is.null(design[["batch"]])) {
      sm_df[["batch"]] <- "undefined"
    } else {
      sm_df[["batch"]] <- as.factor(design[["batch"]])
    }
    if (class(expt_names) == "character" && length(expt_names) == 1) {
      ## Then this refers to an experimental metadata column.
      sm_df[["sample"]] <- design[[expt_names]]
    }
  } else {
    sm_df[["condition"]] <- "undefined"
    sm_df[["batch"]] <- "undefined"
  }

  color_listing <- sm_df[, c("condition", "color")]
  color_listing <- unique(color_listing)
  color_list <- as.character(color_listing[["color"]])
  names(color_list) <- as.character(color_listing[["condition"]])
  sm_df[["num"]] <- seq_len(nrow(sm_df))
  num_batches <- nlevels(sm_df[["batch"]])

  if (!is.null(label_chars) && is.numeric(label_chars)) {
    sm_df[["sample"]] <- abbreviate(sm_df[["sample"]], minlength = label_chars)
  }

  legend_position <- "right"
  if (isFALSE(plot_legend)) {
    legend_position <- "none"
  }

  minval <- min(sm_df[["sm"]])
  maxval <- max(sm_df[["sm"]])
  my_binwidth <- (maxval - minval) / 40
  if (num_batches <= 5) {
    sm_plot <- ggplot(sm_df,
                      aes(x = .data[["num"]], y = .data[["sm"]],
                          shape = .data[["batch"]], fill = .data[["condition"]])) +
      ggplot2::geom_hline(colour = "red", yintercept = ylimit, linewidth = 1) +
      ggplot2::geom_point(size = dot_size,
                          aes(shape = .data[["batch"]],
                              colour = .data[["condition"]],
                              fill = .data[["condition"]])) +
      ggplot2::geom_point(size = dot_size, colour = "black", show.legend = FALSE,
                          aes(shape = .data[["batch"]], fill = .data[["condition"]])) +
      ggplot2::scale_color_manual(name = "Condition",
                                  guide = "legend",
                                  values = color_list) +
      ggplot2::scale_fill_manual(name = "Condition",
                                 guide = "legend",
                                 values = color_list) +
      ggplot2::scale_shape_manual(name = "Batch",
                                  labels = levels(as.factor(sm_df[["batch"]])),
                                  guide = ggplot2::guide_legend(
                                    override.aes = list(size = 5, fill = "grey")),
                                  values = 21:25) +
      ggplot2::scale_x_continuous(labels = sm_df[["sample"]],
                                  breaks = 1:nrow(sm_df),
                                  limits = c(1, nrow(sm_df))) +
      ggplot2::ylab(glue("Standard Median {method}")) +
      ggplot2::xlab("Sample") +
      ggplot2::ggtitle(plot_title) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                     legend.position = legend_position,
                     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
    ## Perhaps instead: hjust = 1.5, vjust = 0.5))

  } else {
    sm_plot <- ggplot(
      sm_df,
      aes(x = .data[["sample"]], y = .data[["sm"]],
          shape = .data[["batch"]], fill = .data[["condition"]])) +
      ggplot2::geom_hline(color = "red", yintercept = ylimit, linewidth = 1) +
      ggplot2::geom_dotplot(binwidth = my_binwidth,
                            binaxis = "y",
                            stackdir = "center",
                            binpositions = "all",
                            colour = "black",
                            dotsize = 1,
                            aes(fill = as.factor(.data[["condition"]]))) +
      ggplot2::ylab(glue("Standard Median {method}")) +
      ggplot2::xlab("Sample") +
      ggplot2::ggtitle(plot_title) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        legend.position = legend_position,
        axis.text.x = ggplot2::element_text(size = base_size, colour = "black",
                                            angle = 90, hjust = 1))

  }
  if (type == "distance") {
    sm_plot <- sm_plot +
      ggplot2::scale_y_continuous(labels = scales::scientific)
  }

  retlist <- list(
    "measurement" = properties,
    "medians" = prop_median,
    "quantile" = prop_spread,
    "plot" = sm_plot)
  class(retlist) <- "standardmedian_plot"
  return(retlist)
}
setGeneric("plot_sm")

## EOF
