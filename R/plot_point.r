## plot_scatter.r: Various scatter plots

#' Steal edgeR's plotBCV() and make it a ggplot2.
#'
#' This was written primarily to understand what that function is doing in edgeR.
#'
#' @param data  A dataframe/expt/exprs with count data
#' @return a plot! of the BCV a la ggplot2.
#' @seealso \pkg{edgeR}
#'  \code{\link[edgeR]{plotBCV}}
#' @examples
#' \dontrun{
#'  bcv <- plot_bcv(expt)
#'  summary(bcv$data)
#'  bcv$plot
#' }
#' @export
plot_bcv <- function(data) {
  data_class <- class(data)[1]
  if (data_class == "expt") {
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    data <- as.data.frame(data)
    ## some functions prefer matrix, so I am keeping this explicit for the moment
  } else {
    stop("This function only understands types: expt, ExpressionSet, data.frame, and matrix.")
  }
  data <- edgeR::DGEList(counts=data)
  edisp <- edgeR::estimateDisp(data)
  avg_log_cpm <- edisp[["AveLogCPM"]]
  if (is.null(avg_log_cpm)) {
    avg_log_cpm <- edgeR::aveLogCPM(edisp[["counts"]], offset=edgeR::getOffset(edisp))
  }
  disper <- edgeR::getDispersion(edisp)
  if (is.null(disper)) {
    stop("No dispersions to plot")
  }
  if (attr(disper, "type") == "common") {
    disper <- rep(disper, length=length(avg_log_cpm))
  }
  disp_df <- data.frame("A" = avg_log_cpm,
                        "disp" = sqrt(disper))
  fitted_disp <- gplots::lowess(disp_df[["A"]], disp_df[["disp"]], f=0.5)
  f <- stats::approxfun(fitted_disp, rule=2)
  disp_df[["label"]] <- rownames(disp_df)
  disp_plot <- ggplot(disp_df, aes_string(x="A", y="disp", label="label")) +
    ggplot2::geom_point() +
    ggplot2::xlab("Average log(CPM)") +
    ggplot2::ylab("Dispersion of Biological Variance") +
    ggplot2::stat_density2d(geom="tile", aes_string(fill="..density..^0.25"),
                            contour=FALSE, show.legend=FALSE) +
    ggplot2::scale_fill_gradientn(
               colours=grDevices::colorRampPalette(c("white", "black"))(256)) +
    ggplot2::geom_smooth(method="loess") +
    ggplot2::stat_function(fun=f, colour="red") +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=base_size, colour="black"))
  ret <- list("data"=disp_df, "plot"=disp_plot)
  return(ret)
}

#' Make a scatter plot between two sets of numbers with a cheesy distance metric
#' and some statistics of the two sets.
#'
#' The distance metric should be codified and made more intelligent.
#' Currently it creates a dataframe of distances which are absolute
#' distances from each axis, multiplied by each other, summed by axis,
#' then normalized against the maximum.
#'
#' @param df Dataframe likely containing two columns.
#' @param tooltip_data Df of tooltip information for gvis graphs.
#' @param gvis_filename Filename to write a fancy html graph.
#' @param size Size of the dots.
#' @return Ggplot2 scatter plot.  This plot provides a "bird's eye"
#' view of two data sets.  This plot assumes the two data structures
#' are not correlated, and so it calculates the median/mad of each
#' axis and uses these to calculate a stupid, home-grown distance
#' metric away from both medians.  This distance metric is used to
#' color dots which are presumed the therefore be interesting because
#' they are far from 'normal.'  This will make a fun clicky googleVis
#' graph if requested.
#' @seealso \pkg{ggplot2}
#'  \code{\link{plot_gvis_scatter}} \code{\link[ggplot2]{geom_point}}
#'  \code{\link{plot_linear_scatter}}
#' @examples
#' \dontrun{
#'  dist_scatter(lotsofnumbers_intwo_columns, tooltip_data=tooltip_dataframe,
#'                    gvis_filename="html/fun_scatterplot.html")
#' }
#' @export
plot_dist_scatter <- function(df, tooltip_data=NULL, gvis_filename=NULL, size=2) {
  df <- data.frame(df[, c(1, 2)])
  df <- df[complete.cases(df), ]
  df_columns <- colnames(df)
  df_x_axis <- df_columns[1]
  df_y_axis <- df_columns[2]
  colnames(df) <- c("first", "second")
  first_median <- summary(df[, 1])["Median"]
  second_median <- summary(df[, 2])["Median"]
  first_mad <- stats::mad(df[, 1])
  second_mad <- stats::mad(df[, 2])
  mydist <- sillydist(df[, 1], df[, 2], first_median, second_median)
  mydist[["x"]] <- abs((mydist[, 1] - first_median) / abs(first_median))
  mydist[["y"]] <- abs((mydist[, 2] - second_median) / abs(second_median))
  mydist[["x"]] <- mydist[["x"]] / max(mydist[["x"]])
  mydist[["y"]] <- mydist[["y"]] / max(mydist[["y"]])
  mydist[["dist"]] <- mydist[["x"]] * mydist[["y"]]
  mydist[["dist"]] <- mydist[["dist"]] / max(mydist[["dist"]])
  line_size <- size / 2
  df[["label"]] <- rownames(df)
  first_vs_second <- ggplot(df, aes_string(x="first", y="second", label="label")) +
    ggplot2::xlab(glue("Expression of {df_x_axis}")) +
    ggplot2::ylab(glue("Expression of {df_y_axis}")) +
    ggplot2::geom_vline(
               color="grey", xintercept=(first_median - first_mad), size=line_size) +
    ggplot2::geom_vline(
               color="grey", xintercept=(first_median + first_mad), size=line_size) +
    ggplot2::geom_vline(
               color="darkgrey", xintercept=first_median, size=line_size) +
    ggplot2::geom_hline(
               color="grey", yintercept=(second_median - second_mad), size=line_size) +
    ggplot2::geom_hline(
               color="grey", yintercept=(second_median + second_mad), size=line_size) +
    ggplot2::geom_hline(color="darkgrey", yintercept=second_median, size=line_size) +
    ggplot2::geom_point(
               colour=grDevices::hsv(mydist[["dist"]], 1, mydist[["dist"]]),
               alpha=0.6, size=size) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=base_size, colour="black"))
  if (!is.null(gvis_filename)) {
    plot_gvis_scatter(df, tooltip_data=tooltip_data, filename=gvis_filename)
  }
  return(first_vs_second)
}

#' Make a scatter plot between two groups with a linear model superimposed and
#' some supporting statistics.
#'
#' @param df Dataframe likely containing two columns.
#' @param tooltip_data Df of tooltip information for gvis graphs.
#' @param gvis_filename  Filename to write a fancy html graph.
#' @param cormethod What type of correlation to check?
#' @param size Size of the dots on the plot.
#' @param identity Add the identity line?
#' @param loess Add a loess estimation?
#' @param gvis_trendline Add a trendline to the gvis plot?  There are a couple
#'   possible types, I think linear is the most common.
#' @param z_lines  Include lines defining the z-score boundaries.
#' @param first First column to plot.
#' @param second Second column to plot.
#' @param base_url Base url to add to the plot.
#' @param pretty_colors Colors!
#' @param color_high Chosen color for points significantly above the mean.
#' @param color_low Chosen color for points significantly below the mean.
#' @param alpha  Choose an alpha channel to define how see-through the dots are.
#' @param ... Extra args likely used for choosing significant genes.
#' @return List including a ggplot2 scatter plot and some histograms.  This plot
#'   provides a "bird's eye" view of two data sets.  This plot assumes a
#'   (potential) linear correlation between the data, so it calculates the
#'   correlation between them.  It then calculates and plots a robust linear
#'   model of the data using an 'SMDM' estimator (which I don't remember how to
#'   describe, just that the document I was reading said it is good).  The
#'   median/mad of each axis is calculated and plotted as well.  The distance
#'   from the linear model is finally used to color the dots on the plot.
#'   Histograms of each axis are plotted separately and then together under a
#'   single cdf to allow tests of distribution similarity.  This will make a fun
#'   clicky googleVis graph if requested.
#' @seealso \pkg{robust} \pkg{stats} \pkg{ggplot2}
#'  \code{\link[robust]{lmRob}} \code{\link[stats]{weights}} \code{\link{plot_histogram}}
#' @examples
#' \dontrun{
#'  plot_linear_scatter(lotsofnumbers_intwo_columns, tooltip_data=tooltip_dataframe,
#'                      gvis_filename="html/fun_scatterplot.html")
#' }
#' @export
plot_linear_scatter <- function(df, tooltip_data=NULL, gvis_filename=NULL,
                                cormethod="pearson", size=2, loess=FALSE,
                                identity=FALSE, gvis_trendline=NULL,
                                z_lines=FALSE, first=NULL, second=NULL,
                                base_url=NULL, pretty_colors=TRUE,
                                color_high=NULL, color_low=NULL, alpha=0.4, ...) {
  ## At this time, one might expect arglist to contain
  ## z, p, fc, n and these will therefore be passed to get_sig_genes()
  arglist <- list(...)
  if (isTRUE(color_high)) {
    color_high <- "#FF0000"
  }
  if (isTRUE(color_low)) {
    color_low <- "#7B9F35"
  }

  df <- data.frame(df[, c(1, 2)])
  df <- df[complete.cases(df), ]
  correlation <- try(cor.test(df[, 1], df[, 2], method=cormethod, exact=FALSE))
  if (class(correlation) == "try-error") {
    correlation <- NULL
  }
  df_columns <- colnames(df)
  df_x_axis <- df_columns[1]
  df_y_axis <- df_columns[2]
  colnames(df) <- c("first", "second")
  model_test <- try(robustbase::lmrob(formula=second ~ first,
                                      data=df, method="SMDM"), silent=TRUE)
  linear_model <- NULL
  linear_model_summary <- NULL
  linear_model_rsq <- NULL
  linear_model_weights <- NULL
  linear_model_intercept <- NULL
  linear_model_slope <- NULL
  if (class(model_test) == "try-error") {
    model_test <- try(lm(formula=second ~ first, data=df), silent=TRUE)
  } else {
    linear_model <- model_test
    linear_model_summary <- summary(linear_model)
    linear_model_rsq <- linear_model_summary[["r.squared"]]
    linear_model_weights <- stats::weights(linear_model, type="robustness", na.action=NULL)
    linear_model_intercept <- stats::coef(linear_model_summary)[1]
    linear_model_slope <- stats::coef(linear_model_summary)[2]
  }
  if (class(model_test) == "try-error") {
    model_test <- try(glm(formula=second ~ first, data=df), silent=TRUE)
  } else {
    linear_model <- model_test
    linear_model_summary <- summary(linear_model)
    linear_model_rsq <- linear_model_summary[["r.squared"]]
    linear_model_weights <- stats::weights(linear_model, type="robustness", na.action=NULL)
    linear_model_intercept <- stats::coef(linear_model_summary)[1]
    linear_model_slope <- stats::coef(linear_model_summary)[2]
  }

  if (class(model_test) == "try-error") {
    message("Could not create a linear model of the data.")
    message("Going to perform a scatter plot without linear model.")
    plot <- plot_scatter(df)
    ret <- list(data=df, scatter=plot)
    return(ret)
  }

  first_median <- summary(df[["first"]])[["Median"]]
  second_median <- summary(df[["second"]])[["Median"]]
  first_mad <- stats::mad(df[["first"]], na.rm=TRUE)
  second_mad <- stats::mad(df[["second"]], na.rm=TRUE)
  line_size <- size / 2
  df[["label"]] <- rownames(df)
  first_vs_second <- ggplot(df, aes_string(x="first", y="second", label="label")) +
    ggplot2::xlab(glue("Expression of {df_x_axis}")) +
    ggplot2::ylab(glue("Expression of {df_y_axis}")) +
    ggplot2::geom_vline(
               color="grey", xintercept=(first_median - first_mad), size=line_size) +
    ggplot2::geom_vline(
               color="grey", xintercept=(first_median + first_mad), size=line_size) +
    ggplot2::geom_hline(
               color="grey", yintercept=(second_median - second_mad), size=line_size) +
    ggplot2::geom_hline(
               color="grey", yintercept=(second_median + second_mad), size=line_size) +
    ggplot2::geom_hline(
               color="darkgrey", yintercept=second_median, size=line_size) +
    ggplot2::geom_vline(
               color="darkgrey", xintercept=first_median, size=line_size) +
    ggplot2::geom_abline(
               colour="grey", slope=linear_model_slope,
               intercept=linear_model_intercept, size=line_size)
  ## The axes and guide-lines are set up, now add the points

  low_df <- high_df <- NULL
  if (!is.null(color_low) | !is.null(color_high)) {
    ## If you want to color the above or below identity line points, then you
    ## will need subsets to define them
    tmpdf <- df
    tmpdf[["ratio"]] <- tmpdf[, 2] - tmpdf[, 1]
    subset_points <- sm(get_sig_genes(tmpdf, column="ratio", ...))
    high_subset <- subset_points[["up_genes"]]
    low_subset <- subset_points[["down_genes"]]
    original_df <- tmpdf
    high_index <- rownames(original_df) %in% rownames(high_subset)
    high_df <- original_df[high_index, ]
    low_index <- rownames(original_df) %in% rownames(low_subset)
    low_df <- original_df[low_index, ]
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(colour="black", size=size, alpha=alpha)
  }

  if (isTRUE(z_lines)) {
    z <- 1.5
    if (!is.null(arglist[["z"]])) {
      z <- arglist[["z"]]
    }
    first_vs_second <- first_vs_second +
      ggplot2::geom_abline(colour="grey", slope=linear_model_slope,
                           intercept=linear_model_intercept + z, size=line_size / 3) +
      ggplot2::geom_abline(colour="grey", slope=linear_model_slope,
                           intercept=linear_model_intercept - z, size=line_size / 3)
  }

  ## Add a color to the dots which are lower than the identity line by some amount
  if (!is.null(color_low)) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(data=low_df, colour=color_low)
  }
  if (!is.null(color_high)) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(data=high_df, colour=color_high)
  }

  if (isTRUE(pretty_colors)) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(size=size, alpha=alpha,
                          colour=grDevices::hsv(linear_model_weights * 9/20,
                                                linear_model_weights/20 + 19/20,
                                                (1.0 - linear_model_weights)))
  } else {
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(colour="black", size=size, alpha=alpha)
  }

  if (loess == TRUE) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_smooth(method="loess")
  }

  if (identity == TRUE) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_abline(colour="darkgreen", slope=1, intercept=0, size=1)
  }

  first_vs_second <- first_vs_second +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=base_size, colour="black"))

  if (!is.null(gvis_filename)) {
    plot_gvis_scatter(df, tooltip_data=tooltip_data, filename=gvis_filename,
                      trendline=gvis_trendline, base_url=base_url)
  }
  if (!is.null(first) & !is.null(second)) {
    colnames(df) <- c(first, second)
  } else if (!is.null(first)) {
    colnames(df) <- c(first, "second")
  } else if (!is.null(second)) {
    colnames(df) <- c("first", second)
  }
  x_histogram <- plot_histogram(data.frame(df[, 1]), fillcolor="lightblue", color="blue")
  y_histogram <- plot_histogram(data.frame(df[, 2]), fillcolor="pink", color="red")
  both_histogram <- plot_multihistogram(df)
  plots <- list(
    "data" = df,
    "scatter" = first_vs_second,
    "x_histogram" = x_histogram,
    "y_histogram" = y_histogram,
    "both_histogram" = both_histogram,
    "correlation" = correlation,
    "lm_model" = linear_model,
    "lm_summary" = linear_model_summary,
    "lm_weights" = linear_model_weights,
    "lm_rsq" = linear_model_rsq,
    "first_median" = first_median,
    "first_mad" = first_mad,
    "second_median" = second_median,
    "second_mad" = second_mad)
  return(plots)
}

#' Make a pretty MA plot from one of limma, deseq, edger, or basic.
#'
#' Because I can never remember, the following from wikipedia: "An MA plot is an
#' application of a Bland-Altman plot for visual representation of two channel
#' DNA microarray gene expression data which has been transformed onto the M
#' (log ratios) and A (mean average) scale."
#'
#' @param table  Df of linear-modelling, normalized counts by sample-type,
#' @param expr_col  Column showing the average expression across genes.
#' @param fc_col  Column showing the logFC for each gene.
#' @param p_col  Column containing the relevant p values.
#' @param pval_cutoff  Name of the pvalue column to use for cutoffs.
#' @param alpha  How transparent to make the dots.
#' @param logfc_cutoff  Fold change cutoff.
#' @param label_numbers  Show how many genes were 'significant', 'up', and 'down'?
#' @param size  How big are the dots?
#' @param tooltip_data  Df of tooltip information for gvis.
#' @param gvis_filename  Filename to write a fancy html graph.
#' @param invert  Invert the ma plot?
#' @param ...  More options for you
#' @return  ggplot2 MA scatter plot.  This is defined as the rowmeans of the
#'   normalized counts by type across all sample types on the x axis, and the
#'   log fold change between conditions on the y-axis. Dots are colored
#'   depending on if they are 'significant.'  This will make a fun clicky
#'   googleVis graph if requested.
#' @seealso \pkg{limma} \pkg{googleVis} \pkg{DESeq2} \pkg{edgeR}
#'  \code{\link{plot_gvis_ma}} \code{\link[limma]{toptable}}
#'  \code{\link[limma]{voom}} \code{\link{hpgl_voom}}
#'  \code{\link[limma]{lmFit}} \code{\link[limma]{makeContrasts}}
#'  \code{\link[limma]{contrasts.fit}}
#' @examples
#'  \dontrun{
#'   plot_ma(voomed_data, table, gvis_filename="html/fun_ma_plot.html")
#'   ## Currently this assumes that a variant of toptable was used which
#'   ## gives adjusted p-values.  This is not always the case and I should
#'   ## check for that, but I have not yet.
#'  }
#' @export
plot_ma_de <- function(table, expr_col="logCPM", fc_col="logFC", p_col="qvalue",
                       pval_cutoff=0.05, alpha=0.4, logfc_cutoff=1, label_numbers=TRUE,
                       size=2, tooltip_data=NULL, gvis_filename=NULL, invert=FALSE, ...) {
  ## Set up the data frame which will describe the plot
  arglist <- list(...)
  ## I like dark blue and dark red for significant and insignificant genes respectively.
  ## Others may disagree and change that with sig_color, insig_color.
  sig_color <- "darkblue"
  if (!is.null(arglist[["sig_color"]])) {
    sig_color <- arglist[["sig_color"]]
  }
  insig_color <- "darkred"
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
    "state" = c("a_upsig", "b_downsig", "c_insig"), stringsAsFactors=TRUE)

  ## Get rid of rows which will be annoying.
  ## If somehow a list got into the data table, this will fail, lets fix that now.
  tmp_table <- table
  for (c in 1:ncol(tmp_table)) {
    tmp_table[[c]] <- as.character(table[[c]])
  }
  rows_without_na <- complete.cases(tmp_table)
  rm(tmp_table)
  table <- table[rows_without_na, ]

  ## Extract the information of interest from my original table
  newdf <- data.frame("avg" = table[[expr_col]],
                      "logfc" = table[[fc_col]],
                      "pval" = table[[p_col]])
  rownames(newdf) <- rownames(table)
  if (isTRUE(invert)) {
    newdf[["logfc"]] <- newdf[["logfc"]] * -1.0
  }
  ## Check if the data is on a log or base10 scale, if the latter, then convert it.
  if (max(newdf[["avg"]]) > 1000) {
    newdf[["avg"]] <- log(newdf[["avg"]])
  }

  ## Set up the state of significant/insiginificant vs. p-value and/or fold-change.
  newdf[["pval"]] <- as.numeric(format(newdf[["pval"]], scientific=FALSE))
  newdf[["pcut"]] <- newdf[["pval"]] <= pval_cutoff
  newdf[["state"]] <- ifelse(newdf[["pval"]] > pval_cutoff, "c_insig",
                      ifelse(newdf[["pval"]] <= pval_cutoff &
                             newdf[["logfc"]] >= logfc_cutoff, "a_upsig",
                      ifelse(newdf[["pval"]] <= pval_cutoff &
                             newdf[["logfc"]] <= (-1.0 * logfc_cutoff),
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
  state_shapes <- c(24, 25, 21)
  names(state_shapes) <- c("a_upsig", "b_downsig", "c_insig")

  ## make the plot!
  plt <- ggplot(data=df,
                ## I am setting x, y, fill color, outline color, and the shape.
                aes_string(x="avg",
                           y="logfc",
                           label="label",
                           fill="as.factor(pcut)",
                           colour="as.factor(pcut)",
                           shape="as.factor(state)")) +
    ggplot2::geom_hline(yintercept=c((logfc_cutoff * -1.0), logfc_cutoff),
                        color="red", size=(size / 3)) +
    ggplot2::geom_point(stat="identity", size=size, alpha=alpha)
  if (isTRUE(label_numbers)) {
    plt <- plt +
      ## The following scale_shape_manual() sets the labels of the legend on the right side.
      ggplot2::scale_shape_manual(name="State", values=state_shapes,
                                  labels=c(
                                    glue("Up Sig.: {num_upsig}"),
                                    glue("Down Sig.: {num_downsig}"),
                                    glue("Insig.: {num_insig}")),
                                  guide=ggplot2::guide_legend(override.aes=aes(size=3,
                                                                               fill="grey")))
  } else {
    plt <- plt +
      ggplot2::scale_shape_manual(name="State", values=state_shapes,
                                  guide="none")
  }

  plt <- plt +
    ## Set the colors of the significant/insignificant points.
    ggplot2::scale_fill_manual(name="as.factor(pcut)",
                               values=c("FALSE"=insig_color, "TRUE"=sig_color),
                               guide=FALSE) +
    ggplot2::scale_color_manual(name="as.factor(pcut)",
                                values=c("FALSE"=insig_color, "TRUE"=sig_color),
                                guide=FALSE) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black")) +
    ggplot2::xlab("Average log2(Counts)") +
    ggplot2::ylab("log2(fold change)")

  ## Make a gvis plot if requested.
  if (!is.null(gvis_filename)) {
    plot_gvis_ma(df, tooltip_data=tooltip_data, filename=gvis_filename, ...)
  }

  ## Recolor a family of genes if requested.
  if (!is.null(family)) {
    plt <- recolor_points(plt, df, family, color=family_color)
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

#' Quick point-recolorizer given an existing plot, df, list of rownames to
#' recolor, and a color.
#'
#' This function should make it easy to color a family of genes in any of the
#' point plots.
#'
#' @param plot  Geom_point based plot
#' @param df  Data frame used to create the plot
#' @param ids  Set of ids which must be in the rownames of df to recolor
#' @param color  Chosen color for the new points.
#' @param ...  Extra arguments are passed to arglist.
#' @return prettier plot.
recolor_points <- function(plot, df, ids, color="red", ...) {
  arglist <- list(...)
  alpha <- 0.3
  if (!is.null(arglist[["alpha"]])) {
    alpha <- arglist[["alpha"]]
  }

  point_index <- rownames(df) %in% ids
  newdf <- df[point_index, ]
  newplot <- plot +
    ggplot2::geom_point(data=newdf,  colour=color, fill=color, alpha=alpha)
  return(newplot)
}

#' Make a ggplot graph of the number of non-zero genes by sample.
#'
#' This puts the number of genes with > 0 hits on the y-axis and CPM on the
#' x-axis. Made by Ramzi Temanni <temanni at umd dot edu>.
#'
#' @param data Expt, expressionset, or dataframe.
#' @param design Eesign matrix.
#' @param colors Color scheme.
#' @param plot_labels How do you want to label the graph? 'fancy' will use
#'   directlabels() to try to match the labels with the positions without
#'   overlapping anything else will just stick them on a 45' offset next to the
#'   graphed point.
#' @param expt_names  Column or character list of preferred sample names.
#' @param label_chars  How many characters for sample names before abbreviation.
#' @param plot_legend  Print a legend for this plot?
#' @param title Add a title?
#' @param ... rawr!
#' @return a ggplot2 plot of the number of non-zero genes with respect to each
#'   library's CPM.
#' @seealso \pkg{ggplot2}
#'  \code{\link[ggplot2]{geom_point}} \code{\link[directlabels]{geom_dl}}
#' @examples
#' \dontrun{
#'  nonzero_plot = plot_nonzero(expt=expt)
#'  nonzero_plot  ## ooo pretty
#' }
#' @export
plot_nonzero <- function(data, design=NULL, colors=NULL, plot_labels=NULL,
                         expt_names=NULL, label_chars=10, plot_legend=FALSE,
                         title=NULL, ...) {
  arglist <- list(...)
  hpgl_env <- environment()
  names <- NULL
  data_class <- class(data)[1]
  if (data_class == "expt") {
    design <- data[["design"]]
    colors <- data[["colors"]]
    names <- data[["samplenames"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    ## some functions prefer matrix, so I am keeping this explicit for the moment
    data <- as.data.frame(data)
  } else {
    stop("This function understands types: expt, ExpressionSet, data.frame, and matrix.")
  }

  shapes <- as.integer(as.factor(design[["batch"]]))
  condition <- design[["condition"]]
  batch <- design[["batch"]]

  if (!is.null(expt_names) && class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      colnames(data) <- make.names(design[[expt_names]], unique=TRUE)
    } else {
      colnames(data) <- expt_names
    }
  }
  if (!is.null(label_chars) && is.numeric(label_chars)) {
    colnames(data) <- abbreviate(colnames(data), minlength=label_chars)
  }
  nz_df <- data.frame(
    "id" = colnames(data),
    "nonzero_genes" = colSums(data >= 1),
    "cpm" = colSums(data) * 1e-6,
    "condition" = condition,
    "batch" = batch,
    "color" = as.character(colors))

  color_listing <- nz_df[, c("condition", "color")]
  color_listing <- unique(color_listing)
  color_list <- as.character(color_listing[["color"]])
  names(color_list) <- as.character(color_listing[["condition"]])
  nz_df[["label"]] <- rownames(nz_df)

  non_zero_plot <- ggplot(data=nz_df,
                          aes_string(x="cpm", y="nonzero_genes", label="label"),
                          environment=hpgl_env) +
    ggplot2::geom_point(size=3, shape=21,
                        aes_string(colour="as.factor(condition)",
                                   fill="as.factor(condition)")) +
    ggplot2::geom_point(size=3, shape=21, colour="black", show.legend=FALSE,
                        aes_string(fill="as.factor(condition)")) +
    ggplot2::scale_color_manual(name="Condition",
                                guide="legend",
                                values=color_list) +
    ggplot2::scale_fill_manual(name="Condition",
                               guide="legend",
                               values=color_list) +
    ggplot2::ylab("Number of non-zero genes observed.") +
    ggplot2::xlab("Observed CPM") +
    ggplot2::theme_bw(base_size=base_size)

  if (is.null(plot_labels)) {
    plot_labels <- "repel"
  }
  if (plot_labels == FALSE) {
    message("Not putting labels on the plot.")
  } else if (plot_labels == "normal") {
    non_zero_plot <- non_zero_plot +
      ggplot2::geom_text(ggplot2::aes_string(x="cpm", y="nonzero_genes", label="id",
                                             angle=45, size=4, vjust=2))
  } else if (plot_labels == "repel") {
    non_zero_plot <- non_zero_plot +
      ggrepel::geom_text_repel(ggplot2::aes_string(label="id"),
                               size=5, box.padding=ggplot2::unit(0.5, "lines"),
                               point.padding=ggplot2::unit(1.6, "lines"),
                               arrow=ggplot2::arrow(length=ggplot2::unit(0.01, "npc")))
  } else if (plot_labels == "dlsmart") {
    non_zero_plot <- non_zero_plot +
      directlabels::geom_dl(ggplot2::aes_string(label="id"), method="smart.grid")
  } else {
    non_zero_plot <- non_zero_plot +
      directlabels::geom_dl(ggplot2::aes_string(label="id"), method="first.qp")
  }

  if (!is.null(title)) {
    non_zero_plot <- non_zero_plot + ggplot2::ggtitle(title)
  }
  non_zero_plot <- non_zero_plot +
    ggplot2::theme(axis.ticks=ggplot2::element_blank(),
                   axis.text=ggplot2::element_text(size=base_size, colour="black"))
  if (isFALSE(plot_legend)) {
    non_zero_plot <- non_zero_plot +
      ggplot2::theme(legend.position="none")
  }

  retlist <- list(
    "plot" = non_zero_plot,
    "table" = nz_df)
  return(retlist)
}

#' Plot all pairwise MA plots in an experiment.
#'
#' Use affy's ma.plot() on every pair of columns in a data set to help diagnose
#' problematic samples.
#'
#' @param data Expt expressionset or data frame.
#' @param log Is the data in log format?
#' @param ... Options are good and passed to arglist().
#' @return List of affy::maplots
#' @seealso \pkg{affy}
#'  \code{\link[affy]{ma.plot}}
#' @examples
#' \dontrun{
#'  ma_plots = plot_pairwise_ma(expt=some_expt)
#' }
#' @export
plot_pairwise_ma <- function(data, log=NULL, ...) {
  data_class <- class(data)[1]
  if (data_class == "expt") {
    design <- data[["design"]]
    colors <- data[["colors"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    ## some functions prefer matrix, so I am keeping this explicit for the moment
    data <- as.data.frame(data)
  } else {
    stop("This function understands types: expt, ExpressionSet, data.frame, and matrix.")
  }
  plot_list <- list()
  for (c in 1:(length(colnames(data)) - 1)) {
    nextc <- c + 1
    for (d in nextc:length(colnames(data))) {
      first <- as.numeric(data[, c])
      second <- as.numeric(data[, d])
      if (max(first) > 1000) {
        if (is.null(log)) {
          message("I suspect you want to set log=TRUE for this.")
          message("In fact, I am so sure, I am doing it now.")
          message("If I am wrong, set log=FALSE, but I'm not.")
          log <- TRUE
        }
      } else if (max(first) < 80) {
        if (!is.null(log)) {
          message("I suspect you want to set log=FALSE for this.")
          message("In fact, I am so  sure, I am doing it now.")
          message("If I am wrong, set log=TRUE.")
          log <- FALSE
        }
      }
      firstname <- colnames(data)[c]
      secondname <- colnames(data)[d]
      name <- glue("{firstname}_{secondname}")
      if (isTRUE(log)) {
        first <- log2(first + 1.0)
        second <- log2(second + 1.0)
      }
      m <- first - second
      a <- (first + second) / 2
      affy::ma.plot(A=a, M=m, plot.method="smoothScatter",
                    show.statistics=TRUE, add.loess=TRUE)
      title(glue("MA of {firstname} vs {secondname}."))
      plot_list[[name]] <- grDevices::recordPlot()
    }
  }
  return(plot_list)
}

#' Make a pretty scatter plot between two sets of numbers.
#'
#' This function tries to supplement a normal scatterplot with some information
#' describing the relationship between the columns of data plotted.
#'
#' @param df Dataframe likely containing two columns.
#' @param gvis_filename Filename to write a fancy html graph.
#' @param tooltip_data Df of tooltip information for gvis.
#' @param size Size of the dots on the graph.
#' @param color Color of the dots on the graph.
#' @return Ggplot2 scatter plot.
#' @seealso \pkg{ggplot2} \pkg{googleVis}
#'  \code{\link{plot_gvis_scatter}} \code{\link[ggplot2]{geom_point}}
#'  \code{\link{plot_linear_scatter}}
#' @examples
#' \dontrun{
#'  plot_scatter(lotsofnumbers_intwo_columns, tooltip_data=tooltip_dataframe,
#'               gvis_filename="html/fun_scatterplot.html")
#' }
#' @export
plot_scatter <- function(df, tooltip_data=NULL, color="black",
                         gvis_filename=NULL, size=2) {
  df <- data.frame(df[, c(1, 2)])
  df <- df[complete.cases(df), ]
  df_columns <- colnames(df)
  df_x_axis <- df_columns[1]
  df_y_axis <- df_columns[2]
  colnames(df) <- c("first", "second")
  df[["label"]] <- rownames(df)
  first_vs_second <- ggplot(df, aes_string(x="first", y="second",
                                           label="label")) +
    ggplot2::xlab(glue("Expression of {df_x_axis}")) +
    ggplot2::ylab(glue("Expression of {df_y_axis}")) +
    ggplot2::geom_point(colour=color, alpha=0.6, size=size) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=10, colour="black"))
  if (!is.null(gvis_filename)) {
    plot_gvis_scatter(df, tooltip_data=tooltip_data, filename=gvis_filename)
  }
  return(first_vs_second)
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
#' @param table  Dataframe from limma's toptable which includes log(fold change) and an
#'     adjusted p-value.
#' @param alpha  How transparent to make the dots.
#' @param color_by  By p-value something else?
#' @param color_list  A list of colors for significance.
#' @param fc_col  Which column contains the fc data?
#' @param fc_name  Name of the fold-change to put on the plot.
#' @param gvis_filename  Filename to write a fancy html graph.
#' @param line_color  What color for the significance lines?
#' @param line_position  Put the significance lines above or below the dots?
#' @param logfc_cutoff  Cutoff defining the minimum/maximum fold change for
#'   interesting.  This is log, so I went with +/- 0.8 mostly arbitrarily as the
#'   default.
#' @param p_col  Which column contains the p-value data?
#' @param p_name  Name of the p-value to put on the plot.
#' @param pval_cutoff  Cutoff defining significant from not.
#' @param shapes_by_state  Add fun shapes for the various significance states?
#' @param size  How big are the dots?
#' @param tooltip_data  Df of tooltip information for gvis.
#' @param ...  I love parameters!
#' @return Ggplot2 volcano scatter plot.  This is defined as the -log10(p-value)
#'   with respect to log(fold change).  The cutoff values are delineated with
#'   lines and mark the boundaries between 'significant' and not.  This will
#'   make a fun clicky googleVis graph if requested.
#' @seealso \pkg{limma}
#'  \code{\link{plot_gvis_ma}} \code{\link[limma]{toptable}}
#'  \code{\link[limma]{voom}} \code{\link{hpgl_voom}} \code{\link[limma]{lmFit}}
#'  \code{\link[limma]{makeContrasts}} \code{\link[limma]{contrasts.fit}}
#' @examples
#' \dontrun{
#'  plot_volcano_de(table, gvis_filename="html/fun_ma_plot.html")
#'  ## Currently this assumes that a variant of toptable was used which
#'  ## gives adjusted p-values.  This is not always the case and I should
#'  ## check for that, but I have not yet.
#' }
#' @export
plot_volcano_de <- function(table, alpha=0.6, color_by="p",
                            color_list=c("FALSE"="darkred", "TRUE"="darkblue"),
                            fc_col="logFC", fc_name="log2 fold change", gvis_filename=NULL,
                            line_color="black", line_position="bottom", logfc_cutoff=1.0,
                            p_col="adj.P.Val", p_name="-log10 p-value", pval_cutoff=0.05,
                            shapes_by_state=TRUE, size=2, tooltip_data=NULL, ...) {
  low_vert_line <- 0.0 - logfc_cutoff
  horiz_line <- -1 * log10(pval_cutoff)

  df <- data.frame("xaxis" = as.numeric(table[[fc_col]]),
                   "yaxis" = as.numeric(table[[p_col]]), stringsAsFactors=TRUE)
  ## This might have been converted to a string
  df[["logyaxis"]] <- -1.0 * log10(as.numeric(df[["yaxis"]]))
  df[["pcut"]] <- df[["yaxis"]] <= pval_cutoff
  df[["state"]] <- ifelse(table[[p_col]] > pval_cutoff, "pinsig",
                   ifelse(table[[p_col]] <= pval_cutoff &
                          table[[fc_col]] >= logfc_cutoff, "upsig",
                   ifelse(table[[p_col]] <= pval_cutoff &
                          table[[fc_col]] <= (-1 * logfc_cutoff),
                          "downsig", "fcinsig")))
  df[["pcut"]] <- as.factor(df[["pcut"]])
  df[["state"]] <- as.factor(df[["state"]])
  df[["label"]] <- rownames(df)

  ## shape 25 is the down arrow, 22 is the square, 23 the diamond, 24 the up arrow
  state_shapes <- c(25, 22, 23, 24)
  names(state_shapes) <- c("downsig", "fcinsig", "pinsig", "upsig")

  color_column <- "pcut"
  color_column_number <- 2
  default_color_list <- c("FALSE"="darkred", "TRUE"="darkblue")
  if (color_by != "p") {
    color_column <- "state"
    color_column_number <- 4
    default_color_list <- c("downsig"="blue", "fcinsig"="darkgrey",
                            "pinsig"="darkgrey", "upsig"="red")
  }
  ## Now make sure that the color column has the correct number of elements.
  if (length(color_list) != color_column_number) {
    message("The color list must have ", color_column_number,
            ", setting it to the default.")
    color_list <- default_color_list
  }

  ## Count the numbers in the categories
  num_downsig <- sum(df[["state"]] == "downsig")
  num_fcinsig <- sum(df[["state"]] == "fcinsig")
  num_pinsig <- sum(df[["state"]] == "pinsig")
  num_upsig <- sum(df[["state"]] == "upsig")

  plt <- NULL
  if (isTRUE(shapes_by_state)) {
    plt <- ggplot(data=df,
                  aes_string(x="xaxis", y="logyaxis", label="label",
                             fill=color_column, colour=color_column, shape="state"))
  } else {
    plt <- ggplot(data=df,
                  aes_string(x="xaxis", y="logyaxis", label="label",
                             fill=color_column, colour=color_column))
  }

  ## Now define when to put lines vs. points
  if (line_position == "bottom") {
    ## lines, then points.
    plt <- plt +
      ggplot2::geom_hline(yintercept=horiz_line, color=line_color, size=(size / 2)) +
      ggplot2::geom_vline(xintercept=logfc_cutoff, color=line_color, size=(size / 2)) +
      ggplot2::geom_vline(xintercept=low_vert_line, color=line_color, size=(size / 2)) +
      ggplot2::geom_point(stat="identity", size=size, alpha=alpha)
  } else {
    ## points, then lines
    plt <- plt +
      ggplot2::geom_point(stat="identity", size=size, alpha=alpha) +
      ggplot2::geom_hline(yintercept=horiz_line, color=line_color, size=(size / 2)) +
      ggplot2::geom_vline(xintercept=logfc_cutoff, color=line_color, size=(size / 2)) +
      ggplot2::geom_vline(xintercept=low_vert_line, color=line_color, size=(size / 2))
  }

  ## If shapes are being set by state,  add that to the legend now.
  if (isTRUE(shapes_by_state)) {
    plt <- plt +
      ggplot2::scale_shape_manual(
                 name="state", values=state_shapes,
                 labels=c(
                   glue("Down Sig.: {num_downsig}"),
                   glue("FC Insig.: {num_fcinsig}"),
                   glue("P Insig.: {num_pinsig}"),
                   glue("Up Sig.: {num_upsig}")),
                 guide=ggplot2::guide_legend(override.aes=aes(size=3, fill="grey")))
      }

  ## Now set the colors and axis labels
  plt <- plt +
    ggplot2::scale_fill_manual(name=color_column, values=color_list,
                               guide=FALSE) +
    ggplot2::scale_color_manual(name=color_column, values=color_list,
                                guide=FALSE) +
    ggplot2::xlab(label=fc_name) +
    ggplot2::ylab(label=p_name) +
    ## ggplot2::guides(shape=ggplot2::guide_legend(override.aes=list(size=3))) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"))
  ##  axis.text.x=ggplot2::element_text(angle=-90))

  gvis_result <- NULL
  if (!is.null(gvis_filename)) {
    gvis_result <- plot_gvis_volcano(table, fc_col=fc_col, p_col=p_col,
                                     logfc_cutoff=logfc_cutoff, fc_cutoff=logfc_cutoff,
                                     pval_cutoff=pval_cutoff, tooltip_data=tooltip_data,
                                     filename=gvis_filename)
  }
  retlist <- list("plot" = plt,
                  "df" = df,
                  "gvis" = gvis_result)
  return(retlist)
}

## EOF
