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
plot_svfactor <- function(expt, svest, chosen_factor="batch", factor_type="factor") {
  chosen <- expt[["design"]][[chosen_factor]]
  sv_df <- data.frame(
    "adjust" = svest[, 1],  ## Take a single estimate from compare_estimates()
    "factors" = chosen,
    "samplenames" = rownames(expt[["design"]])
  )
  samplenames <- rownames(expt[["design"]])
  my_colors <- expt[["colors"]]

  sv_melted <- reshape2::melt(sv_df, idvars="factors")
  minval <- min(sv_df$adjust)
  maxval <- max(sv_df$adjust)
  my_binwidth <- (maxval - minval) / 40
  sv_plot <- ggplot2::ggplot(sv_melted, ggplot2::aes_string(x="factors", y="value")) +
    ggplot2::geom_dotplot(binwidth=my_binwidth, binaxis="y",
                          stackdir="center", binpositions="all",
                          colour="black", fill=my_colors) +
    ggplot2::xlab(paste0("Experimental factor: ", chosen_factor)) +
    ggplot2::ylab(paste0("1st surrogate variable estimation")) +
    ##ggplot2::geom_text(ggplot2::aes_string(x="factors", y="value", label="strains"), angle=45, size=3, vjust=2) +
    ggplot2::geom_text(ggplot2::aes_string(x="factors", y="value", label="samplenames"),
                       angle=45, size=3, vjust=2) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=10, colour="black")) +
    ggplot2::theme_bw()
  return(sv_plot)
}

## I want to make some notes on getting ggplot plots to look like I want:
## 1. Order of geom_thingie()/scale_thingie() matters
## 2. The initial df needs to have a column for every color, fill, shape, etc.
## 3. The initial aes needs to have a mapping for the same color, fill, shape, etc.
## 4. Lay down a geom/scale for every element to change.

###        ggplot2::geom_point(size=5,
###                            aes_string(shape="as.factor(shape)",
###                                       colour="condition",
###                                       fill="condition")) +
###        ggplot2::geom_point(size=5,
###                            colour="black",
###                            show.legend=FALSE,
###                            aes_string(shape="as.factor(shape)",fill="condition")) +
###        ggplot2::scale_color_manual(values=color_list) +
###        ggplot2::scale_fill_manual(values=color_list) +
###        ggplot2::scale_shape_manual(values=21, guide=FALSE)

#' Make a dotplot of known batches vs. SVs.
#'
#' This should make a quick df of the factors and surrogates and plot them.  Maybe it should be
#' folded into plot_svfactor?  Hmm, I think first I will write this and see if it is better.
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
plot_batchsv <- function(expt, svs, batch_column="batch", factor_type="factor") {
  chosen <- expt[["design"]][[batch_column]]
  num_batches <- length(unique(chosen))

  factor_df <- data.frame(
    "sample" = expt[["design"]][["sampleid"]],
    "factor" = as.integer(as.factor(expt[["design"]][[batch_column]])),
    "fill" = expt[["colors"]],
    "condition" = expt[["conditions"]],
    "batch" = expt[["batches"]],
    "shape" = 21,
    "color" = "black",
    "svs" = svs[, 1])
  color_list <- as.character(factor_df[["fill"]])
  names(color_list) <- as.character(factor_df[["condition"]])

  sample_factor <- ggplot(factor_df, aes_string(x="sample", y="factor")) +
    ggplot2::geom_dotplot(binaxis="y", stackdir="center",
                          binpositions="all", colour="black",
                          ##fill=factor_df[["fill"]])
                          fill=color_list)

  factor_svs <- ggplot2::ggplot(data=as.data.frame(factor_df),
                                aes_string(x="factor",
                                           y="svs",
                                           fill="condition",
                                           colour="condition",
                                           shape="shape")) +
    ggplot2::geom_point(size=5,
                        aes_string(shape="as.factor(shape)",
                                   colour="condition",
                                   fill="condition")) +
    ggplot2::geom_point(size=5,
                        colour="black",
                        show.legend=FALSE,
                        aes_string(shape="as.factor(shape)",
                                   fill="condition")) +
    ggplot2::scale_color_manual(values=color_list) +
    ggplot2::scale_fill_manual(values=color_list) +
    ggplot2::scale_shape_manual(values=21, guide=FALSE)

  svs_sample <- ggplot(factor_df, aes_string(x="sample", y="svs")) +
    ggplot2::geom_dotplot(binaxis="y", stackdir="center",
                          binpositions="all", colour="black",
                          fill=factor_df[["fill"]]) +
    ggplot2::geom_text(aes_string(x="sample", y="svs", label="batch"),
                       size=4, vjust=2) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(size=10, colour="black",
                                                     angle=90, hjust=1)) +
    ggplot2::theme_bw()

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
  sv_plot <- ggplot2::ggplot(pc_df, aes_string(x="factors", y="value")) +
    ggplot2::geom_dotplot(binwidth=my_binwidth, binaxis="y",
                          stackdir="center", binpositions="all",
                          colour="black", fill=my_colors) +
    ggplot2::xlab(paste0("Experimental factor: ", exp_factor)) +
    ggplot2::ylab(paste0("1st surrogate variable estimation")) +
    ggplot2::geom_text(ggplot2::aes_string(x="factors", y="value", label="strains"), angle=45, size=3, vjust=2) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(size=10, colour="black")) +
    ggplot2::theme_bw()
  return(sv_plot)
}

#' Make an R plot of the standard median correlation or distance among samples.
#'
#' This was written by a mix of Kwame Okrah <kokrah at gmail dot com>, Laura
#' Dillon <dillonl at umd dot edu>, and Hector Corrada Bravo <hcorrada at umd dot edu>
#' I reimplemented it using ggplot2 and tried to make it a little more flexible.
#' The general idea is to take the pairwise correlations/distances of the samples, then take the
#' medians, and plot them.
#'
#' @param data Expt, expressionset, or data frame.
#' @param colors Color scheme if data is not an expt.
#' @param method Correlation or distance method to use.
#' @param names Use pretty names for the samples?
#' @param title Title for the graph.
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
plot_sm <- function(data, colors=NULL, method="pearson", names=NULL, title=NULL, ...) {
  data_class <- class(data)[1]
  arglist <- list(...)
  conditions <- NULL
  if (data_class == "expt") {
    design <- data[["design"]]
    colors <- data[["colors"]]
    names <- data[["names"]]
    conditions <- data[["conditions"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    design <- pData(data)
    conditions <- pData(data)[["conditions"]]
    data <- exprs(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
  } else {
    stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }

  if (is.null(names)) {
    names <- colnames(data)
  }

  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(ncol(data), chosen_palette))(ncol(data))
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

  outer_limit <- NULL
  ylimit <- NULL
  if (method == "pearson" | method == "spearman" | method == "robust") {
    outer_limit <- prop_spread[1] - 1.5 * prop_iqr
    ylimit <- c(pmin(min(prop_median), outer_limit), max(prop_median))
    ylimit <- ylimit[[1]]
  } else {
    outer_limit <- prop_spread[2] + (1.5 * prop_iqr)
    ylimit <- c(pmin(min(prop_median), outer_limit), max(prop_median))
    ylimit <- ylimit[[2]]
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

  minval <- min(sm_df[["sm"]])
  maxval <- max(sm_df[["sm"]])
  my_binwidth <- (maxval - minval) / 40

  sm_plot <- ggplot2::ggplot(sm_df, aes_string(x="sample", y="sm", fill="condition")) +
    ggplot2::geom_hline(color="red", yintercept=ylimit, size=1) +
    ggplot2::geom_dotplot(binwidth=my_binwidth,
                          binaxis="y",
                          stackdir="center",
                          binpositions="all",
                          colour="black",
                          dotsize=1,
                          aes_string(fill="as.factor(condition)")) +
    ggplot2::scale_fill_manual(name="Condition",
                               guide="legend",
                               values=color_list) +
    ggplot2::ylab(paste0("Standard Median ", method)) +
    ggplot2::xlab(paste0("Sample")) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(size=10, colour="black",
                                                     angle=90, hjust=1)) +
    ggplot2::theme_bw()

  return(sm_plot)
}

## EOF
