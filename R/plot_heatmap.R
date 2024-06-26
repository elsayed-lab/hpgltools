## plot_heatmap.r: Heatmaps separated by usage

#' Make a heatmap.3 description of the correlation between samples.
#'
#' Given a set of count tables and design, this will calculate the pairwise
#' correlations and plot them as a heatmap.  It attempts to standardize the
#' inputs and eventual output.
#'
#' @param expt_data Dataframe, expt, or expressionset to work with.
#' @param expt_colors Color scheme for the samples, not needed if this is an expt.
#' @param expt_design Design matrix describing the experiment, not needed if this is an expt.
#' @param method Correlation statistic to use. (pearson, spearman, kendall, robust).
#' @param expt_names Alternate names to use for the samples.
#' @param batch_row Name of the design row used for 'batch' column colors.
#' @param plot_title Title for the plot.
#' @param label_chars  Limit on the number of label characters.
#' @param ... More options are wonderful!
#' @return Gplots heatmap describing describing how the samples are clustering
#'   vis a vis pairwise correlation.
#' @seealso [grDevice] [gplot2::heatmap.2()]
#' @examples
#' \dontrun{
#'  corheat_plot <- hpgl_corheat(expt = expt, method = "robust")
#' }
#' @export
plot_corheat <- function(expt_data, expt_colors = NULL, expt_design = NULL,
                         method = "pearson", expt_names = NULL,
                         batch_row = "batch", plot_title = NULL, label_chars = 10, ...) {
  map_list <- plot_heatmap(expt_data, expt_colors = expt_colors, expt_design = expt_design,
                           method = method, expt_names = expt_names, type = "correlation",
                           batch_row = batch_row, plot_title = plot_title, label_chars = label_chars, ...)
  class(map_list) <- "correlation_heatmap"
  return(map_list)
}

#' Make a heatmap.3 of the distances (euclidean by default) between samples.
#'
#' Given a set of count tables and design, this will calculate the pairwise
#' distances and plot them as a heatmap.  It attempts to standardize the inputs
#' and eventual output.
#'
#' @param expt_data Dataframe, expt, or expressionset to work with.
#' @param expt_colors Color scheme (not needed if an expt is provided).
#' @param expt_design Design matrix (not needed if an expt is provided).
#' @param method Distance metric to use.
#' @param expt_names Alternate names to use for the samples.
#' @param batch_row Name of the design row used for 'batch' column colors.
#' @param plot_title Title for the plot.
#' @param label_chars Limit on the number of label characters.
#' @param ... More parameters!
#' @return a recordPlot() heatmap describing the distance between samples.
#' @seealso [gplots::heatmap.2()]
#' @examples
#' \dontrun{
#'  disheat_plot = plot_disheat(expt = expt, method = "euclidean")
#' }
#' @export
plot_disheat <- function(expt_data, expt_colors = NULL, expt_design = NULL,
                         method = "euclidean", expt_names = NULL,
                         batch_row = "batch",  plot_title = NULL, label_chars = 10, ...) {
  map_list <- plot_heatmap(expt_data, expt_colors = expt_colors, expt_design = expt_design,
                           method = method, expt_names = expt_names, type = "distance",
                           batch_row = batch_row, plot_title = plot_title,
                           label_chars = label_chars, ...)
  class(map_list) <- "distance_heatmap"
  return(map_list)
}

#' Make a heatmap.3 plot, does the work for plot_disheat and plot_corheat.
#'
#' This does what is says on the tin.  Sets the colors for correlation or
#' distance heatmaps, handles the calculation of the relevant metrics, and plots
#' the heatmap.
#'
#' @param expt_data Dataframe, expt, or expressionset to work with.
#' @param expt_colors Color scheme for the samples.
#' @param expt_design Design matrix describing the experiment vis a vis
#'  conditions and batches.
#' @param method Distance or correlation metric to use.
#' @param expt_names Alternate names to use for the samples.
#' @param type Defines the use of correlation, distance, or sample heatmap.
#' @param batch_row Name of the design row used for 'batch' column colors.
#' @param plot_title Title for the plot.
#' @param label_chars Limit on the number of label characters.
#' @param ... I like elipses!
#' @return a recordPlot() heatmap describing the distance between samples.
#' @seealso [gplots::heatmap.2()]
#' @export
plot_heatmap <- function(expt_data, expt_colors = NULL, expt_design = NULL,
                         method = "pearson", expt_names = NULL,
                         type = "correlation", batch_row = "batch", plot_title = NULL,
                         label_chars = 10, ...) {
  arglist <- list(...)

  margin_list <- c(12, 9)
  if (!is.null(arglist[["margin_list"]])) {
    margin_list <- arglist[["margin_list"]]
  }
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }
  keysize <- 2
  if (!is.null(arglist[["keysize"]])) {
    keysize <- arglist[["keysize"]]
  }
  remove_equal <- FALSE
  if (!is.null(arglist[["remove_equal"]])) {
    remove_equal <- arglist[["remove_equal"]]
  }

  if (is.null(expt_colors)) {
    num_cols <- ncol(expt_data)
    expt_colors <- sm(grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(num_cols, chosen_palette))(num_cols))
  }
  expt_colors <- as.character(expt_colors)

  if (is.null(expt_names)) {
    expt_names <- colnames(expt_data)
  } else {
    if (class(expt_names) == "character" && length(expt_names) == 1) {
      ## Then this refers to an experimental metadata column.
      expt_names <- expt_design[[expt_names]]
    }
  }
  if (!is.null(label_chars) && is.numeric(label_chars)) {
    expt_names <- abbreviate(expt_names, minlength = label_chars)
  }

  if (isTRUE(remove_equal)) {
    cv_min <- 1
    if (!is.null(arglist[["cv_min"]])) {
      cv_min <- arglist[["cv_min"]]
    }
    cv_max <- Inf
    if (!is.null(arglist[["cv_max"]])) {
      cv_min <- arglist[["cv_max"]]
    }
    test <- genefilter::cv(cv_min, cv_max)
    filter_list <- genefilter::filterfun(test)
    answer <- genefilter::genefilter(expt_data, filter_list)
    expt_data <- expt_data[answer, ]
  }

  heatmap_data <- NULL
  heatmap_colors <- NULL
  if (type == "correlation") {
    heatmap_data <- hpgl_cor(expt_data, method = method)
    heatmap_colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(100)
    if (method == "cordist") {
      heatmap_colors <- grDevices::colorRampPalette(
        c("yellow2", "goldenrod", "darkred"),
        bias = 0.5)(100)
    }
  } else if (type == "distance") {
    heatmap_data <- as.matrix(dist(t(expt_data)), method = method)
    heatmap_colors <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(9, "GnBu"))(100)
  } else {
    heatmap_colors <- gplots::redgreen(75)
    heatmap_data <- as.matrix(expt_data)
  }
  ## Set the batch colors depending on # of observed batches
  if (is.null(expt_design)) {
    row_colors <- rep("white", length(expt_colors))
  } else if (is.null(expt_design[[batch_row]])) {
    row_colors <- rep("white", length(expt_colors))
  } else if (length(levels(as.factor(expt_design[[batch_row]]))) >= 2) {
    ## We have >= 2 batches, and so will fill in the column colors
    num_batch_colors <- length(levels(as.factor(expt_design[[batch_row]])))
    batch_color_assignments <- as.integer(as.factor(expt_design[[batch_row]]))
    row_colors <- RColorBrewer::brewer.pal(12, "Set3")[batch_color_assignments]
  } else {
    ## If we just have 1 batch, make it... white (to disappear).
    row_colors <- rep("white", length(expt_colors))
  }

  map <- NULL
  na_idx <- is.na(heatmap_data)
  heatmap_data[na_idx] <- 0
  tmp_file <- tmpmd5file(pattern = "heat", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  if (type == "correlation") {
    map <- heatmap.3(heatmap_data, keysize = keysize, labRow = expt_names,
                     labCol = expt_names, ColSideColors = expt_colors,
                     RowSideColors = row_colors, margins = margin_list,
                     scale = "none", trace = "none",
                     linewidth = 0.5, main = plot_title, ...)
  } else {
    map <- heatmap.3(heatmap_data, keysize = keysize, labRow = expt_names,
                     labCol = expt_names, ColSideColors = expt_colors,
                     RowSideColors = row_colors, margins = margin_list,
                     scale = "none", trace = "none",
                     linewidth = 0.5, main = plot_title,
                     col = rev(heatmap_colors), ...)
  }
  recorded_heatmap_plot <- grDevices::recordPlot()
  dev.off()
  removed <- suppressWarnings(file.remove(tmp_file))
  removed <- unlink(dirname(tmp_file))

  retlist <- list("map" = map,
                  "plot" = recorded_heatmap_plot,
                  "data" = heatmap_data)
  return(retlist)
}
setGeneric("plot_heatmap")

#' Potential replacement for heatmap.2 based plots.
#'
#' Heatplus is an interesting tool, I have a few examples of using it and intend
#' to include them here.
#'
#' @param expt Experiment to try plotting.
#' @param type What comparison method to use on the data (distance or correlation)?
#' @param method What distance/correlation method to perform?
#' @param annot_columns Set of columns to include as terminal columns next to the heatmap.
#' @param annot_rows Set of columns to include as terminal rows below the heatmap.
#' @param cutoff Cutoff used to define color changes in the annotated clustering.
#' @param cluster_colors Choose colors for the clustering?
#' @param scale Scale the heatmap colors?
#' @param cluster_width How much space to include between clustering?
#' @param cluster_function Choose an alternate clustering function than hclust()?
#' @param heatmap_colors Choose your own heatmap cluster palette?
#' @return List containing the returned heatmap along with some parameters used to create it.
#' @seealso [Heatplus] [fastcluster]
#' @export
plot_heatplus <- function(expt, type = "correlation", method = "pearson", annot_columns = "batch",
                          annot_rows = "condition", cutoff = 1.0,
                          cluster_colors = NULL, scale = "none",
                          cluster_width = 2.0, cluster_function = NULL, heatmap_colors = NULL) {
  data <- exprs(expt)
  if (type == "correlation") {
    data <- hpgl_cor(data, method = method)
  } else {
    data <- hpgl_dist(data, method = "euclidean")
  }

  if (is.null(cluster_function)) {
    cluster_function <- fastcluster::hclust
  }

  mydendro <- list(
      "clustfun" = cluster_function,
      "lwd" = cluster_width)
  des <- pData(expt)
  col_data <- as.data.frame(des[, annot_columns])
  colnames(col_data) <- annot_columns
  row_data <- as.data.frame(des[, annot_rows])
  colnames(row_data) <- annot_rows

  myannot <- list(
      "inclRef" = FALSE,
      "Col" = list("data" = col_data),
      "Row" = list("data" = row_data))

  if (is.null(cluster_colors)) {
    cluster_colors <- Heatplus::BrewerClusterCol
  }
  myclust <- list("cuth" = cutoff,
                  "col" = cluster_colors)
  mylabs <- list(
      "Row" = list("nrow" = 4),
      "Col" = list("nrow" = 4))

  first_map <- Heatplus::annHeatmap2(data, dendrogram = mydendro, annotation = myannot,
                                     cluster = myclust, labels = mylabs)

  number_colors <- length(levels(as.factor(first_map[["data"]][["x"]])))
  if (is.null(heatmap_colors)) {
    heatmap_colors <- colorRampPalette(c("darkblue", "beige"))(number_colors)
  }

  if (is.null(cluster_colors)) {
    num_clusters <- max(first_map[["cluster"]][["Row"]][["grp"]])
    chosen_palette <- "Dark2"
    new_colors <- sm(grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(
        num_clusters, chosen_palette))(num_clusters))
  }
  myclust <- list("cuth" = 1.0,
                  "col" = new_colors)

  final_map <- Heatplus::annHeatmap2(
    data, dendrogram = mydendro, annotation = myannot,
    cluster = myclust, labels = mylabs, scale = scale, col = heatmap_colors)

  tmp_file <- tmpmd5file(pattern = "heat", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  plot(final_map)
  rec_plot <- grDevices::recordPlot()
  dev.off()
  removed <- suppressWarnings(file.remove(tmp_file))
  removed <- unlink(dirname(tmp_file))

  retlist <- list(
    "annotations" = myannot,
    "clusters" = myclust,
    "labels" = mylabs,
    "colors" = heatmap_colors,
    "first_map" = first_map,
    "map" = final_map,
    "plot" = rec_plot)
  return(retlist)
}

#' Make a heatmap.3 description of the similarity of the genes among samples.
#'
#' Sometimes you just want to see how the genes of an experiment are related to
#' each other.  This can handle that.  These heatmap functions should probably
#' be replaced with neatmaps or heatplus or whatever it is, as the annotation
#' dataframes in them are pretty awesome.
#'
#' @param data Expt/expressionset/dataframe set of samples.
#' @param colors Color scheme of the samples (not needed if input is an expt).
#' @param design Design matrix describing the experiment (gotten for free if an expt).
#' @param heatmap_colors Specify a colormap.
#' @param expt_names Alternate samples names.
#' @param dendrogram Where to put dendrograms?
#' @param row_label Passed through to heatmap.2.
#' @param plot_title Title of the plot!
#' @param Rowv Reorder the rows by expression?
#' @param Colv Reorder the columns by expression?
#' @param label_chars Maximum number of characters before abbreviating sample names.
#' @param filter Filter the data before performing this plot?
#' @param ... More parameters for a good time!
#' @return a recordPlot() heatmap describing the samples.
#' @seealso [gplots::heatmap.2()]
#' @export
plot_sample_heatmap <- function(data, colors = NULL, design = NULL, heatmap_colors = NULL,
                                expt_names = NULL, dendrogram = "column",
                                row_label = NA, plot_title = NULL, Rowv = TRUE,
                                Colv = TRUE, label_chars = 10, filter = TRUE, ...) {
  data <- as.data.frame(data)
  if (is.null(heatmap_colors)) {
    heatmap_colors <- gplots::redgreen(75)
  }
  if (is.null(names)) {
    names <- colnames(data)
  }
  data <- as.matrix(data)

  if (is.null(expt_names)) {
    expt_names <- colnames(data)
  } else {
    if (class(expt_names) == "character" && length(expt_names) == 1) {
      ## Then this refers to an experimental metadata column.
      expt_names <- design[[expt_names]]
    }
  }
  if (!is.null(label_chars) && is.numeric(label_chars)) {
    expt_names <- abbreviate(expt_names, minlength = label_chars)
  }

  ## drop NAs to help hclust()
  na_idx <- is.na(data)
  data[na_idx] <- -20
  tmp_file <- tmpmd5file(pattern = "heat", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  heatmap.3(data, keysize = 0.8, labRow = row_label, col = heatmap_colors, dendrogram = dendrogram,
            labCol = expt_names, margins = c(12, 8), trace = "none", ColSideColors = colors,
            linewidth = 0.5, main = plot_title, Rowv = Rowv, Colv = Colv)
  hpgl_heatmap_plot <- grDevices::recordPlot()
  dev.off()
  removed <- suppressWarnings(file.remove(tmp_file))
  removed <- unlink(dirname(tmp_file))

  return(hpgl_heatmap_plot)
}
setGeneric("plot_sample_heatmap")

#' An experiment to see if I can visualize the genes with the highest variance.
#'
#' @param expt ExpressionSet
#' @param fun mean or median
#' @param fact Which factor to slice/dice the data?
#' @param row_label Label the rows?
#' @param plot_title Title for the plot
#' @param Rowv Row vs (yeah I forgot what this does.)
#' @param Colv Col vs
#' @param label_chars Maximum number of characters in the sample IDs.
#' @param dendrogram Make a tree of the samples?
#' @param min_delta Minimum delta value for filtering
#' @param x_factor When plotting two factors against each other, which is x?
#' @param y_factor When plotting two factors against each other, which is y?
#' @param min_cvsd Include only those with a minimal CV?
#' @param cv_min Minimum cv to examine (I think this should be slightly lower)
#' @param cv_max Maximum cV to examine (I think this should be limited to ~ 0.7?)
#' @param remove_equal Filter uninteresting genes.
#' @export
plot_sample_cvheatmap <- function(expt, fun = "mean", fact = "condition",
                                  row_label = NA, plot_title = NULL, Rowv = TRUE,
                                  Colv = TRUE, label_chars = 10, dendrogram = "column",
                                  min_delta = 0.5, x_factor = 1, y_factor = 2, min_cvsd = NULL,
                                  cv_min = 1, cv_max = Inf, remove_equal = TRUE) {
  ## I am certain there is a better way to do this, but I am tired and not thinking well.
  colors <- c()
  expt_colors <- expt[["colors"]]
  fact_info <- pData(expt)[[fact]]
  names(expt_colors) <- fact_info
  for (i in seq_along(expt_colors)) {
    name <- names(expt_colors)[i]
    this_color <- as.character(expt_colors[i])
    colors[name] <- this_color
  }

  if (isTRUE(remove_equal)) {
    expt <- normalize_expt(expt, filter = "cv", cv_min = 1, cv_max = Inf)
  }
  cvs <- as.matrix(median_by_factor(expt, fun = fun, fact = fact)[["cvs"]])
  if (!is.null(min_cvsd)) {
    cv_sds <- matrixStats::rowSds(cvs)
    if (is.numeric(min_cvsd)) {
      keepers <- cv_sds >= min_cvsd
    } else if (is.character(min_cvsd)) {
      ## Min. 1st Qu. Median ...
      min_cvsd <- summary(cv_sds)[[min_cvsd]]
      keepers <- cv_sds >= min_cvsd
    } else {
      keepers <- cv_sds >= summary(cv_sds)[5]
    }
    cvs <- cvs[keepers, ]
  }

  heatmap_colors <- gplots::redgreen(75)
  if (is.null(names)) {
    names <- colnames(data)
  }

  tmp_file <- tmpmd5file(pattern = "heat", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  heatmap.3(cvs, keysize = 0.8, labRow = rownames(cvs), col = heatmap_colors, dendrogram = dendrogram,
            margins = c(12, 8), trace = "none", ColSideColors = colors,
            linewidth = 0.5, main = plot_title, Rowv = Rowv, Colv = Colv)
  cv_heatmap_plot <- grDevices::recordPlot()
  dev.off()
  removed <- suppressWarnings(file.remove(tmp_file))
  removed <- unlink(dirname(tmp_file))

  point_df <- cvs[, c(x_factor, y_factor)]

  return(cv_heatmap_plot)
}

#' a minor change to heatmap.2 makes heatmap.3
#'
#' heatmap.2 is the devil.
#'
#' @param x data
#' @param Rowv add rows?
#' @param Colv add columns?
#' @param distfun distance function to use
#' @param hclustfun clustering function to use
#' @param dendrogram which axes to put trees on
#' @param reorderfun reorder the rows/columns?
#' @param symm symmetrical?
#' @param scale add the scale?
#' @param na.rm remove nas from the data?
#' @param revC reverse the columns?
#' @param add.expr no clue
#' @param breaks also no clue
#' @param symbreaks still no clue
#' @param col colors!
#' @param colsep column separator
#' @param rowsep row separator
#' @param sepcolor color to put between columns/rows
#' @param sepwidth how much to separate
#' @param cellnote mur?
#' @param notecex size of the notes
#' @param notecol color of the notes
#' @param na.color a parameter call to bg
#' @param trace do a trace for rows/columns?
#' @param tracecol color of the trace
#' @param hline the hline
#' @param vline the vline
#' @param linecol the line color
#' @param margins margins are good
#' @param ColSideColors colors for the columns as annotation
#' @param RowSideColors colors for the rows as annotation
#' @param cexRow row size
#' @param cexCol column size
#' @param labRow hmmmm
#' @param labCol still dont know
#' @param srtRow srt the row?
#' @param srtCol srt the column?
#' @param adjRow adj the row?
#' @param adjCol adj the column?
#' @param offsetRow how far to place the text from the row
#' @param offsetCol how far to place the text from the column
#' @param key add a key?
#' @param keysize if so, how big?
#' @param density.info for the key, what information to add
#' @param denscol tracecol hmm ok
#' @param symkey I like keys
#' @param densadj adj the dens?
#' @param key.title title for the key
#' @param key.xlab text for the x axis of the key
#' @param key.ylab text for the y axis of the key
#' @param key.xtickfun add text to the ticks of the key x axis
#' @param key.ytickfun add text to the ticks of the key y axis
#' @param key.par parameters for the key
#' @param main the main title of the plot
#' @param xlab main x label
#' @param ylab main y label
#' @param lmat the lmat
#' @param lhei the lhei
#' @param lwid the lwid
#' @param extrafun I do enjoy me some extra fun
#' @param linewidth the width of lines
#' @param ... because this function did not already have enough options
#' @return a heatmap!
#' @seealso \code{\link[gplots]{heatmap.2}}
#' @export
heatmap.3 <- function(x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist, hclustfun = fastcluster::hclust,
                      dendrogram = c("both", "row", "column", "none"),
                      reorderfun = function(d, w) reorder(d, w),
                      symm = FALSE, scale = c("none", "row", "column"),
                      na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks,
                      symbreaks = min(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors", colsep, rowsep, sepcolor = "white",
                      sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan",
                      na.color = par("bg"), trace = c("column", "row", "both", "none"),
                      tracecol = "cyan", hline = median(breaks), vline = median(breaks),
                      linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors,
                      cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
                      labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, NA),
                      adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5,
                      key = TRUE, keysize = 1.5, density.info = c("histogram", "density", "none"),
                      denscol = tracecol, symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25,
                      key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL,
                      key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL, ylab = NULL,
                      lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, linewidth = 1.0, ...) {
  if (!is.null(main)) {
    if (main == FALSE) {
      main <- NULL
    }
  }
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low) / (high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
             "none"
           else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.",
            "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dendrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
    if (length(rowInd) > nr || any(rowInd < 1 | rowInd >
                                   nr))
      stop("Rowv dendrogram doesn't match size of x")
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
    if (length(colInd) > nc || any(colInd < 1 | colInd >
                                   nc))
      stop("Colv dendrogram doesn't match size of x")
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
                               x
                             else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
                               x
                             else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
                (1:nr)[rowInd]
              else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
                (1:nc)[colInd]
              else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) <
      1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || length(ColSideColors) !=
          nc)
        stop("'ColSideColors' must be a character vector of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
                                             1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || length(RowSideColors) !=
          nr)
        stop("'RowSideColors' must be a character vector of length nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
                                             1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr),
        axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!gtools::invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  if (is.null(srtCol))
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + offsetCol,
         tick = 0, cex.axis = cexCol, hadj = adjCol[1], padj = adjCol[2])
  else {
    if (is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol))
        adjCol <- c(1, NA)
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2,
                   tick = 0)
      text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * strheight("M"),
           labels = labCol, adj = adjCol, cex = cexCol, srt = srtCol)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtCol ignored.")
  }
  if (is.null(srtRow)) {
    axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow,
         tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
  }
  else {
    if (is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, labels = rep("", nr), las = 2,
                   line = -0.5, tick = 0)
      text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"),
           y = ypos, labels = labRow, adj = adjRow, cex = cexRow,
           srt = srtRow)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtRow ignored.")
  }
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0,
                              xright = csep + 0.5 + sepwidth[1],
                              ytop = ncol(x) + 1, lty = 1,
                              lwd = linewidth, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5,
                              xright = nrow(x) + 1,
                              ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2],
                              lty = 1, lwd = linewidth,
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = linewidth, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i - 0.5 + hline.vals, col = linecol,
               lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = linewidth, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    mar <- c(5, 4, 2, 1)
    if (!is.null(key.xlab) && is.na(key.xlab))
      mar[1] <- 2
    if (!is.null(key.ylab) && is.na(key.ylab))
      mar[2] <- 2
    if (!is.null(key.title) && is.na(key.title))
      mar[3] <- 1
    par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
    if (length(key.par) > 0)
      do.call(par, key.par)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    if (is.null(key.xtickfun)) {
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      xargs <- list(at = xv, labels = lv)
    }
    else {
      xargs <- key.xtickfun()
    }
    xargs$side <- 1
    do.call(axis, xargs)
    if (is.null(key.xlab)) {
      if (scale == "row")
        key.xlab <- "Row Z-Score"
      else if (scale == "column")
        key.xlab <- "Column Z-Score"
      else key.xlab <- "Value"
    }
    if (!is.na(key.xlab)) {
      mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5)
    }
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = linewidth)
      if (is.null(key.ytickfun)) {
        yargs <- list(at = pretty(dens$y)/max(dens$y) *
                        0.95, labels = pretty(dens$y))
      }
      else {
        yargs <- key.ytickfun()
      }
      yargs$side <- 2
      do.call(axis, yargs)
      if (is.null(key.title))
        key.title <- "Color Key\nand Density Plot"
      if (!is.na(key.title))
        title(key.title)
      par(cex = 0.5)
      if (is.null(key.ylab))
        key.ylab <- "Density"
      if (!is.na(key.ylab))
        mtext(side = 2, key.ylab, line = par("mgp")[1],
              padj = 0.5)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = linewidth, type = "s",
            col = denscol)
      if (is.null(key.ytickfun)) {
        yargs <- list(at = pretty(hy)/max(hy) * 0.95,
                      labels = pretty(hy))
      }
      else {
        yargs <- key.ytickfun()
      }
      yargs$side <- 2
      do.call(axis, yargs)
      if (is.null(key.title))
        key.title <- "Color Key\nand Histogram"
      if (!is.na(key.title))
        title(key.title)
      par(cex = 0.5)
      if (is.null(key.ylab))
        key.ylab <- "Count"
      if (!is.na(key.ylab))
        mtext(side = 2, key.ylab, line = par("mgp")[1],
              padj = 0.5)
    }
    else title("Color Key")
    if (trace %in% c("both", "column")) {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if (!is.null(vline)) {
        abline(v = vline.vals, col = linecol, lty = 2)
      }
    }
    if (trace %in% c("both", "row")) {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if (!is.null(hline)) {
        abline(v = hline.vals, col = linecol, lty = 2)
      }
    }
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  if (!is.null(extrafun))
    extrafun()
  invisible(retval)
}

## EOF
