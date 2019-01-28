## plot_dotplot.r: Dotplots in various contexts, currently just smc/smd

#' Make a dotplot of some categorised factors and a set of SVs (for other factors).
#'
#' This should make a quick df of the factors and surrogates and plot them.
#'
#' @param expt Experiment from which to acquire the design, counts, etc.
#' @param svest Set of surrogate variable estimations from sva/svg
#'  or batch estimates.
#' @param chosen_factor Factor to compare against.
#' @param factor_type This may be a factor or range, it is intended to plot
#'  a scatterplot if it is a range, a dotplot if a factor.
#' @return surrogate variable plot as per Leek's work
#' @seealso \pkg{ggplot2}
#' @examples
#' \dontrun{
#'  estimate_vs_snps <- plot_svfactor(start, surrogate_estimate, "snpcategory")
#' }
#' @export
plot_svfactor <- function(expt, svest, sv=1, chosen_factor="batch", factor_type="factor") {
  chosen <- expt[["design"]][[chosen_factor]]
  sv_df <- data.frame(
    "adjust" = svest[, sv],  ## Take a single estimate from compare_estimates()
    "factors" = chosen,
    "samplenames" = rownames(expt[["design"]])
  )
  samplenames <- rownames(expt[["design"]])
  my_colors <- expt[["colors"]]

  sv_melted <- reshape2::melt(sv_df, idvars="factors")
  minval <- min(sv_df[["adjust"]])
  maxval <- max(sv_df[["adjust"]])
  my_binwidth <- (maxval - minval) / 40
  sv_plot <- ggplot2::ggplot(sv_melted, ggplot2::aes_string(x="factors", y="value")) +
    ggplot2::geom_dotplot(binwidth=my_binwidth, binaxis="y",
                          stackdir="center", binpositions="all",
                          colour="black", fill=my_colors) +
    ggplot2::xlab(glue("Experimental factor: {chosen_factor}")) +
    ggplot2::ylab("1st surrogate variable estimation") +
    ggplot2::geom_text(ggplot2::aes_string(x="factors", y="value", label="samplenames"),
                       angle=45, size=3, vjust=2) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black")) +
    ggplot2::theme_bw(base_size=base_size)
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
#' @param batch_column Which experimental design column to use?
#' @param factor_type This may be a factor or range, it is intended to plot
#'  a scatterplot if it is a range, a dotplot if a factor.
#' @return Plot of batch vs surrogate variables as per Leek's work.
#' @seealso \pkg{sva} \pkg{ggplot2}
#' @examples
#' \dontrun{
#'  estimate_vs_snps <- plot_batchsv(start, surrogate_estimate, "snpcategory")
#' }
#' @export
plot_batchsv <- function(expt, svs, sv=1, batch_column="batch", factor_type="factor") {
  chosen <- pData(expt)[, batch_column]
  names(chosen) <- sampleNames(expt)
  num_batches <- length(unique(chosen))

  factor_df <- data.frame(
    "sample" = expt[["design"]][["sampleid"]],
    "factor" = as.integer(as.factor(expt[["design"]][[batch_column]])),
    "fill" = expt[["colors"]],
    "condition" = expt[["conditions"]],
    "batch" = expt[["batches"]],
    "color" = "black",
    "svs" = svs[, sv])
  if (num_batches <= 5) {
    factor_df[["shape"]] <- 20 + as.numeric(factor_df[["batch"]])
  } else {
    factor_df[["shape"]] <- 21
  }
  factor_df[["shape"]] <- as.factor(factor_df[["shape"]])

  color_list <- as.character(factor_df[["fill"]])
  names(color_list) <- as.character(factor_df[["condition"]])

  sample_factor <- ggplot(factor_df,
                          aes_string(x="sample",
                                     y="factor",
                                     shape="batch",
                                     fill="condition")) +
    ggplot2::geom_point(size=5,
                        aes_string(shape="batch",
                                   colour="condition",
                                   fill="condition")) +
    ggplot2::geom_point(size=5,
                        colour="black",
                        show.legend=FALSE,
                        aes_string(shape="batch",
                                   fill="condition")) +
    ggplot2::scale_shape_manual(name="Batch",
                                labels=levels(as.factor(factor_df[["batch"]])),
                                guide=ggplot2::guide_legend(
                                                 override.aes=list(size=5, fill="grey")),
                                values=21:25) +
    ggplot2::scale_color_manual(name="Condition",
                                guide="legend",
                                values=color_list) +
    ggplot2::scale_fill_manual(name="Condition",
                               guide="legend",
                               values=color_list) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   axis.text.x=ggplot2::element_text(angle=90, vjust=0.5))
  ##, hjust=1.5, vjust=0.5))

  factor_svs <- ggplot2::ggplot(factor_df,
                                aes_string(x="factor",
                                           y="svs",
                                           fill="condition",
                                           colour="condition",
                                           shape="shape")) +
    ggplot2::geom_point(size=5,
                        aes_string(shape="batch",
                                   colour="condition",
                                   fill="condition")) +
    ggplot2::geom_point(size=5,
                        colour="black",
                        show.legend=FALSE,
                        aes_string(shape="batch",
                                   fill="condition")) +
    ggplot2::scale_shape_manual(name="Batch",
                                labels=levels(as.factor(factor_df[["batch"]])),
                                guide=ggplot2::guide_legend(
                                                 override.aes=list(size=5, fill="grey")),
                                values=21:25) +
    ggplot2::scale_color_manual(name="Condition",
                                guide="legend",
                                values=color_list) +
    ggplot2::scale_fill_manual(name="Condition",
                               guide="legend",
                               values=color_list) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   axis.text.x=ggplot2::element_text(angle=90, vjust=0.5))
  ##, hjust=1.5, vjust=0.5))

  svs_sample <- ggplot2::ggplot(factor_df,
                                aes_string(x="sample",
                                           y="svs",
                                           fill="condition",
                                           colour="condition",
                                           shape="shape")) +
    ggplot2::geom_point(size=5,
                        aes_string(shape="batch",
                                   colour="condition",
                                   fill="condition")) +
    ggplot2::geom_point(size=5,
                        colour="black",
                        show.legend=FALSE,
                        aes_string(shape="batch",
                                   fill="condition")) +
    ggplot2::scale_shape_manual(name="Batch",
                                labels=levels(as.factor(factor_df[["batch"]])),
                                guide=ggplot2::guide_legend(
                                                 override.aes=list(size=5, fill="grey")),
                                values=21:25) +
    ggplot2::scale_color_manual(name="Condition",
                                guide="legend",
                                values=color_list) +
    ggplot2::scale_fill_manual(name="Condition",
                               guide="legend",
                               values=color_list) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   axis.text.x=ggplot2::element_text(angle=90, vjust=0.5))
  ## An alternate possibility:  hjust=1.5, vjust=0.5))

  plots <- list(
    "sample_factor" = sample_factor,
    "factor_svs" = factor_svs,
    "svs_sample" = svs_sample)
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
#' @seealso \pkg{ggplot2}
#' @examples
#' \dontrun{
#'  estimate_vs_pcs <- plot_pcfactor(pcs, times)
#' }
#' @export
plot_pcfactor <- function(pc_df, expt, exp_factor="condition", component="PC1") {
  samplenames <- rownames(expt[["design"]])

  my_colors <- expt[["colors"]]
  minval <- min(pc_df[[component]])
  maxval <- max(pc_df[[component]])
  my_binwidth <- (maxval - minval) / 40
  sv_plot <- ggplot2::ggplot(pc_df, aes_string(x=exp_factor, y=component)) +
    ggplot2::geom_dotplot(binwidth=my_binwidth, binaxis="y",
                          stackdir="center", binpositions="all",
                          colour="black", fill=my_colors) +
    ggplot2::xlab(glue("Experimental factor: {exp_factor}")) +
    ggplot2::ylab("1st surrogate variable estimation") +
    ##ggplot2::geom_text(
    ##           ggplot2::aes_string(x=exp_factor, y=component, label="strains"),
    ##           angle=45, size=3, vjust=2) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   axis.text.x=ggplot2::element_text(angle=90, vjust=0.5))
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
#' @param colors Color scheme if data is not an expt.
#' @param method Correlation or distance method to use.
#' @param plot_legend  Include a legend on the side?
#' @param expt_names Use pretty names for the samples?
#' @param label_chars Maximum number of characters before abbreviating sample names.
#' @param title Title for the graph.
#' @param dot_size  How large should the glyphs be?
#' @param ... More parameters to make you happy!
#' @return ggplot of the standard median something
#'  among the samples.  This will also write to an
#'  open device.  The resulting plot measures the median correlation of
#'  each sample among its peers.  It notes 1.5* the interquartile range
#'  among the samples and makes a horizontal line at that correlation
#'  coefficient.  Any sample which falls below this line is considered
#'  for removal because it is much less similar to all of its peers.
#' @seealso \pkg{matrixStats} \pkg{grDevices}
#'  \code{\link{hpgl_cor}} \code{\link[matrixStats]{rowMedians}}
#'  \code{\link[stats]{quantile}} \code{\link{diff}} \code{\link[grDevices]{recordPlot}}
#' @examples
#' \dontrun{
#'  smc_plot = hpgl_smc(expt=expt)
#' }
#' @export
plot_sm <- function(data, colors=NULL, method="pearson", plot_legend=FALSE,
                    expt_names=NULL, label_chars=10, title=NULL, dot_size=5, ...) {
  arglist <- list(...)
  data_class <- class(data)[1]
  conditions <- NULL
  if (data_class == "expt") {
    design <- data[["design"]]
    colors <- data[["colors"]]
    conditions <- data[["conditions"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    design <- pData(data)
    conditions <- pData(data)[["conditions"]]
    data <- exprs(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    data <- as.data.frame(data)
    ## some functions prefer matrix, so I am keeping this explicit for the moment
  } else {
    stop("This function currently only understands classes of type:
 expt, ExpressionSet, data.frame, and matrix.")
  }

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
  if (method == "pearson" | method == "spearman" | method == "robust") {
    message("Performing correlation.")
    properties <- hpgl_cor(data, method=method)
  } else {
    message("Performing distance.")
    properties <- as.matrix(dist(t(data)))
  }

  prop_median <- matrixStats::rowMedians(properties)
  prop_spread <- stats::quantile(prop_median, p=c(1, 3) / 4)
  prop_iqr <- diff(prop_spread)
  log_scale <- FALSE

  outer_limit <- NULL
  ylimit <- NULL
  type <- "unknown"
  if (method == "pearson" | method == "spearman" | method == "robust") {
    outer_limit <- prop_spread[1] - 1.5 * prop_iqr
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

  if (is.null(conditions)) {
    conditions <- colors
  }

  sm_df <- data.frame(
    "sample" = rownames(properties),
    "sm" = prop_median,
    "condition" = conditions,
    "color" = colors)
  color_listing <- sm_df[, c("condition", "color")]
  color_listing <- unique(color_listing)
  color_list <- as.character(color_listing[["color"]])
  names(color_list) <- as.character(color_listing[["condition"]])
  sm_df[["num"]] <- 1:nrow(sm_df)
  if (is.null(design[["batch"]])) {
    sm_df[["batch"]] <- "undefined"
  } else {
    sm_df[["batch"]] <- as.factor(design[["batch"]])
  }
  num_batches <- nlevels(sm_df[["batch"]])

  if (class(expt_names) == "character" && length(expt_names) == 1) {
    ## Then this refers to an experimental metadata column.
    sm_df[["sample"]] <- design[[expt_names]]
  }
  if (!is.null(label_chars) && is.numeric(label_chars)) {
    sm_df[["sample"]] <- abbreviate(sm_df[["sample"]], minlength=label_chars)
  }

  legend_position <- "right"
  if (isFALSE(plot_legend)) {
    legend_position <- "none"
  }

  minval <- min(sm_df[["sm"]])
  maxval <- max(sm_df[["sm"]])
  my_binwidth <- (maxval - minval) / 40
  if (num_batches <= 5) {
    sm_plot <- ggplot(sm_df, aes_string(
                               x="num", y="sm", shape="batch", fill="condition")) +
      ggplot2::geom_hline(colour="red", yintercept=ylimit, size=1) +
      ggplot2::geom_point(size=dot_size,
                          aes_string(shape="batch",
                                     colour="condition",
                                     fill="condition")) +
      ggplot2::geom_point(size=dot_size, colour="black", show.legend=FALSE,
                          aes_string(shape="batch", fill="condition")) +
      ggplot2::scale_color_manual(name="Condition",
                                  guide="legend",
                                  values=color_list) +
      ggplot2::scale_fill_manual(name="Condition",
                                 guide="legend",
                                 values=color_list) +
      ggplot2::scale_shape_manual(name="Batch",
                                  labels=levels(as.factor(sm_df[["batch"]])),
                                  guide=ggplot2::guide_legend(
                                                   override.aes=list(size=5, fill="grey")),
                                  values=21:25) +
      ggplot2::scale_x_continuous(labels=sm_df[["sample"]],
                                  breaks=1:nrow(sm_df),
                                  limits=c(1, nrow(sm_df))) +
    ggplot2::ylab(glue("Standard Median {method}")) +
    ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   legend.position=legend_position,
                   axis.text.x=ggplot2::element_text(angle=90, vjust=0.5))
    ## Perhaps instead: hjust=1.5, vjust=0.5))

  } else {
    sm_plot <- ggplot2::ggplot(
                          sm_df,
                          aes_string(x="sample", y="sm", shape="batch", fill="condition")) +
      ggplot2::geom_hline(color="red", yintercept=ylimit, size=1) +
      ggplot2::geom_dotplot(binwidth=my_binwidth,
                            binaxis="y",
                            stackdir="center",
                            binpositions="all",
                            colour="black",
                            dotsize=1,
                            aes_string(fill="as.factor(condition)")) +
      ggplot2::ylab(glue("Standard Median {method}")) +
      ggplot2::xlab("Sample") +
      ggplot2::ggtitle(title) +
      ggplot2::theme_bw(base_size=base_size) +
      ggplot2::theme(
                 legend.position=legend_position,
                 axis.text.x=ggplot2::element_text(size=base_size, colour="black",
                                                   angle=90, hjust=1))

  }
  if (type == "distance") {
    sm_plot <- sm_plot +
      ggplot2::scale_y_continuous(labels=scales::scientific)
  }
  return(sm_plot)
}

## EOF
