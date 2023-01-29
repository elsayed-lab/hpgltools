## Note to self, I think for future ggplot2 plots, I must start by creating the data frame
## Then cast every column in it explicitly, and only then invoke ggplot(data = df ...)

## If I see something like:
## 'In sample_data$mean = means : Coercing LHS to a list'
## That likely means that I was supposed to have data in the
## data.frame() format, but instead it is a matrix.  In functions
## where this is a danger, it is a likely good idea to cast it as a
## data frame.

#' Look at the range of the data for a plot and use it to suggest if a plot
#' should be on log scale.
#'
#' There are a bunch of plots which often-but-not-always benefit from being
#' displayed on a log scale rather than base 10.  This is a quick and dirty
#' heuristic which suggests the appropriate scale.  If the data 'should' be on
#' the log scale and it has 0s, then they are moved to 1 so that when logged
#' they will return to 0.  Similarly, if there are negative numbers and the
#' intended scale is log, then this will set values less than 0 to zero to avoid
#' imaginary numbers.
#'
#' @param data Data to plot.
#' @param scale If known, this will be used to define what (if any) values to
#'   change.
#' @param max_data Define the upper limit for the heuristic.
#' @param min_data Define the lower limit for the heuristic.
check_plot_scale <- function(data, scale = NULL, max_data = 10000, min_data = 10) {
  if (is.null(scale)) {
    if (max(data) > max_data & min(data) < min_data) {
      mesg("This data will benefit from being displayed on the log scale.")
      mesg("If this is not desired, set scale='raw'")
      scale <- "log"
      negative_idx <- data < 0
      if (sum(negative_idx) > 0) {
        data[negative_idx] <- 0
        message("Changed ", sum(negative_idx), " negative features to 0.")
      }
      zero_idx <- data == 0
      if (sum(zero_idx) > 0) {
        message(sum(zero_idx), " entries are 0.  We are on a log scale, adding 1 to the data.")
        data <- data + 1
      }
    } else {
      scale <- "raw"
    }
  } else {
    mesg("An explicit scale was requested: ", scale, ".")
  }
  retlist <- list(
      "data" = data,
      "scale" = scale)
  return(retlist)
}

#' Translate the hexadecimal color codes to three decimal numbers.
#'
#' @param rgb hexadecimal color input.
color_int <- function(rgb) {
  hex <- gsub(pattern = "^\\#", replacement = "", x = rgb)
  red <- as.integer(as.hexmode(gsub(pattern = "^(.{2}).{4}$", replacement = "\\1", x = hex)))
  green <- as.integer(as.hexmode(gsub(pattern = "^.{2}(.{2}).{2}$", replacement = "\\1", x = hex)))
  blue <- as.integer(as.hexmode(gsub(pattern = "^.{4}(.{2})$", replacement = "\\1", x = hex)))
  retlist <- list(
      "red" = red,
      "green" = green,
      "blue" = blue)
  return(retlist)
}

#' Simplify plotly ggplot conversion so that there are no shenanigans.
#'
#' I am a fan of ggplotly, but its conversion to an html file is not perfect.
#' This hopefully will get around the most likely/worst problems.
#'
#' @param gg Plot from ggplot2.
#' @param filename Output filename.
#' @param selfcontained htmlwidgets: Return the plot as a self-contained file
#'  with images re-encoded base64.
#' @param libdir htmlwidgets: Directory into which to put dependencies.
#' @param background htmlwidgets: String for the background of the image.
#' @param plot_title htmlwidgets: Title of the page!
#' @param knitrOptions htmlwidgets: I am not a fan of camelCase, but
#'  nonetheless, options from knitr for htmlwidgets.
#' @param ... Any remaining elipsis options are passed to ggplotly.
#' @return The final output filename
#' @seealso [htmlwidgets] [plotly] [ggplot2]
#' @export
ggplt <- function(gg, filename = "ggplot.html",
                  selfcontained = TRUE, libdir = NULL, background = "white",
                  plot_title = class(gg)[[1]], knitrOptions = list(), ...) {
  base <- basename(filename)
  dir <- dirname(filename)

  ## 202210: There is a deprecated function call in plotly, which is out
  ## of the scope of my interest.
  out <- suppressWarnings(plotly::ggplotly(gg,
                                           ...))
  widget <- htmlwidgets::saveWidget(
                             plotly::as_widget(out), base, selfcontained, libdir = libdir,
                             background = background, title = plot_title, knitrOptions = knitrOptions)
  final <- base
  if (dir != ".") {
    final <- file.path(dir, base)
    moved <- file.rename(base, final)
  }
  return(final)
}

#' Make lots of graphs!
#'
#' Plot out a set of metrics describing the state of an experiment
#' including library sizes, # non-zero genes, heatmaps, boxplots,
#' density plots, pca plots, standard median distance/correlation, and
#' qq plots.
#'
#' @param expt an expt to process
#' @param cormethod The correlation test for heatmaps.
#' @param distmethod define the distance metric for heatmaps.
#' @param title_suffix Text to add to the titles of the plots.
#' @param qq Include qq plots?
#' @param ma Include pairwise ma plots?
#' @param gene_heat Include a heatmap of the gene expression data?
#' @param ... Extra parameters optionally fed to the various plots
#' @return a loooong list of plots including the following:
#' \enumerate{
#'   \item nonzero = a ggplot2 plot of the non-zero genes vs library size
#'   \item libsize = a ggplot2 bar plot of the library sizes
#'   \item boxplot = a ggplot2 boxplot of the raw data
#'   \item corheat = a recordPlot()ed pairwise correlation heatmap of the raw data
#'   \item smc = a recordPlot()ed view of the standard median pairwise correlation of the raw data
#'   \item disheat = a recordPlot()ed pairwise euclidean distance heatmap of the raw data
#'   \item smd = a recordPlot()ed view of the standard median pairwise distance of the raw data
#'   \item pcaplot = a recordPlot()ed PCA plot of the raw samples
#'   \item pcatable = a table describing the relative contribution of condition/batch of the raw data
#'   \item pcares =  a table describing the relative contribution of condition/batch of the raw data
#'   \item pcavar = a table describing the variance of the raw data
#'   \item qq = a recordPlotted() view comparing the quantile/quantiles between
#'      the mean of all data and every raw sample
#'   \item density = a ggplot2 view of the density of each raw sample (this is
#'      complementary but more fun than a boxplot)
#' }
#' @seealso [plot_nonzero()] [plot_legend()] [plot_libsize()] [plot_disheat()]
#'  [plot_corheat()] [plot_topn()] [plot_pca()] [plot_sm()] [plot_boxplot()]
#' @examples
#' \dontrun{
#'  toomany_plots <- graph_metrics(expt)
#'  toomany_plots$pcaplot
#'  norm <- normalize_expt(expt, convert = "cpm", batch = TRUE, filter_low = TRUE,
#'                         transform = "log2", norm = "rle")
#'  holy_asscrackers <- graph_metrics(norm, qq = TRUE, ma = TRUE)
#' }
#' @export
graph_metrics <- function(expt, cormethod = "pearson", distmethod = "euclidean",
                          title_suffix = NULL, qq = FALSE, ma = FALSE, gene_heat = FALSE,
                          ...) {
  arglist <- list(...)
  if (!exists("expt", inherits = FALSE)) {
    stop("The input data does not exist.")
  }
  dev_length <- length(dev.list())
  if (dev_length > 1) {
    message("Hey! You have ", dev_length,
            " plotting devices open, this might be in error.")
  }
  ## First gather the necessary data for the various plots.
  old_options <- options(scipen = 10)
  nonzero_title <- "Non zero genes"
  libsize_title <- "Library sizes"
  boxplot_title <- "Boxplot"
  corheat_title <- "Correlation heatmap"
  smc_title <- "Standard Median Correlation"
  disheat_title <- "Distance heatmap"
  smd_title <- "Standard Median Distance"
  pca_title <- "Principle Component Analysis"
  tsne_title <- "T-SNE Analysis"
  dens_title <- "Density plot"
  cv_title <- "Coefficient of variance plot"
  topn_title <- "Top-n representation"
  pc_loading_title <- "Expression of top-n PC loading-genes"
  if (!is.null(title_suffix)) {
    nonzero_title <- glue("{nonzero_title}: {title_suffix}")
    libsize_title <- glue("{libsize_title}: {title_suffix}")
    boxplot_title <- glue("{boxplot_title}: {title_suffix}")
    corheat_title <- glue("{corheat_title}: {title_suffix}")
    smc_title <- glue("{smc_title}: {title_suffix}")
    disheat_title <- glue("{disheat_title}: {title_suffix}")
    smd_title <- glue("{smd_title}: {title_suffix}")
    pca_title <- glue("{pca_title}: {title_suffix}")
    tsne_title <- glue("{tsne_title}:  {title_suffix}")
    dens_title <- glue("{dens_title}: {title_suffix}")
    cv_title <- glue("{cv_title}: {title_suffix}")
    topn_title <- glue("{topn_title}: {title_suffix}")
    pc_loading_title <- glue("{pc_loading_title}: {title_suffix}")
  }

  ## I am putting the ... arguments on a separate line so that I can check that
  ## each of these functions is working properly in an interactive session.
  mesg("Graphing number of non-zero genes with respect to CPM by library.")
  nonzero <- try(plot_nonzero(expt, plot_title = nonzero_title,
                              ...))
  if ("try-error" %in% class(nonzero)) {
    nonzero <- list()
  }
  mesg("Graphing library sizes.")
  libsize <- try(plot_libsize(expt, plot_title = libsize_title,
                              ...))
  if ("try-error" %in% class(libsize)) {
    libsize <- list()
  }
  mesg("Graphing a boxplot.")
  boxplot <- try(plot_boxplot(expt, plot_title = boxplot_title,
                              ...))
  if ("try-error" %in% class(boxplot)) {
    boxplot <- NULL
  }
  mesg("Graphing a correlation heatmap.")
  corheat <- try(plot_corheat(expt, method = cormethod, plot_title = corheat_title,
                              ...))
  if ("try-error" %in% class(corheat)) {
    corheat <- list()
  }
  mesg("Graphing a standard median correlation.")
  smc <- try(plot_sm(expt, method = cormethod, plot_title = smc_title,
                     ...))
  if ("try-error" %in% class(smc)) {
    smc <- NULL
  }
  mesg("Graphing a distance heatmap.")
  disheat <- try(plot_disheat(expt, method = distmethod, plot_title = disheat_title,
                              ...))
  if ("try-error" %in% class(disheat)) {
    disheat <- list()
  }
  mesg("Graphing a standard median distance.")
  smd <- try(plot_sm(expt, method = distmethod, plot_title = smd_title,
                     ...))
  if ("try-error" %in% class(smd)) {
    smd <- NULL
  }
  mesg("Graphing a PCA plot.")
  pca <- try(plot_pca(expt, plot_title = pca_title,
                      ...))
  if ("try-error" %in% class(pca)) {
    pca <- list()
  }
  mesg("Graphing a T-SNE plot.")
  tsne <- try(plot_tsne(expt, plot_title = tsne_title,
                        ...))
  if ("try-error" %in% class(tsne)) {
    tsne <- list()
  }
  mesg("Plotting a density plot.")
  density <- try(plot_density(expt, plot_title = dens_title,
                              ...))
  if ("try-error" %in% class(density)) {
    density <- list()
  }
  mesg("Plotting a CV plot.")
  cv <- try(plot_variance_coefficients(expt, plot_title = dens_title,
                                       ...))
  if ("try-error" %in% class(cv)) {
    cv <- list()
  }
  mesg("Plotting the representation of the top-n genes.")
  topn <- try(plot_topn(expt, plot_title = topn_title,
                        ...))
  if ("try-error" %in% class(topn)) {
    topn <- list()
  }
  tmp_expt <- sm(normalize_expt(expt, filter = TRUE))
  pcload <- list()
  if (nrow(exprs(tmp_expt)) > ncol(exprs(tmp_expt))) {
    mesg("Plotting the expression of the top-n PC loaded genes.")
    pcload <- try(plot_pcload(tmp_expt, plot_title = pc_loading_title))
    if ("try-error" %in% class(pcload)) {
      pcload <- list()
    }
  }

  mesg("Printing a color to condition legend.")
  legend <- try(plot_legend(expt))
  if ("try-error" %in% class(legend)) {
    legend <- list()
  }
  qq_logs <- NULL
  qq_ratios <- NULL
  if (isTRUE(qq)) {
    mesg("QQ plotting!")
    qq_plots <- try(sm(suppressWarnings(plot_qq_all(tmp_expt,
                                                    ...))))
    if ("try-error" %in% class(qq_plots)) {
      qq_plots <- list()
    }
    qq_logs <- qq_plots[["logs"]]
    qq_ratios <- qq_plots[["ratios"]]
  }

  ma_plots <- NULL
  if (isTRUE(ma)) {
    mesg("Many MA plots!")
    ma_plots <- try(suppressWarnings(plot_pairwise_ma(expt,
                                                      ...)))
    if ("try-error" %in% class(ma_plots)) {
      ma_plots <- list()
    }
  }

  gene_heatmap <- NULL
  if (isTRUE(gene_heat)) {
    mesg("gene heatmap!")
    gene_heatmap <- try(suppressWarnings(plot_sample_heatmap(tmp_expt,
                                                             ...)))
    if ("try-error" %in% class(gene_heatmap)) {
      gene_heatmap <- NULL
    }
  }

  ret_data <- list(
      "boxplot" = boxplot,
      "corheat" = corheat[["plot"]],
      "cvplot" = cv[["plot"]],
      "density" = density[["plot"]],
      "density_table" = density[["table"]],
      "disheat" = disheat[["plot"]],
      "gene_heatmap" = gene_heatmap,
      "legend" = legend[["plot"]],
      "legend_colors" = legend[["colors"]],
      "libsize" = libsize[["plot"]],
      "libsizes" = libsize[["table"]],
      "libsize_summary" = libsize[["summary"]],
      "ma" = ma_plots,
      "nonzero" = nonzero[["plot"]],
      "nonzero_table" = nonzero[["table"]],
      "pc_loadplot" = pcload[["plot"]],
      "pc_summary" = pca[["residual_df"]],
      "pc_propvar" = pca[["prop_var"]],
      "pc_plot" = pca[["plot"]],
      "pc_table" = pca[["table"]],
      "qqlog" = qq_logs,
      "qqrat" = qq_ratios,
      "smc" = smc[["plot"]],
      "smd" = smd[["plot"]],
      "topnplot" = topn[["plot"]],
      "tsne_summary" = tsne[["residual_df"]],
      "tsne_propvar" = tsne[["prop_var"]],
      "tsne_plot" = tsne[["plot"]],
      "tsne_table" = tsne[["table"]]
  )
  new_options <- options(old_options)
  return(ret_data)
}

#' Scab the legend from a PCA plot and print it alone
#'
#' This way I can have a legend object to move about.
#'
#' @param stuff This can take either a ggplot2 pca plot or some data from which to make one.
#' @return A legend!
#' @export
plot_legend <- function(stuff) {
  plot <- NULL
  if (class(stuff)[[1]] == "gg") {
    ## Then assume it is a pca plot
    plot <- stuff
  } else {
    plot <- plot_pca(stuff)[["plot"]]
  }

  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot))
  leg <- which(sapply(tmp[["grobs"]], function(x) x[["name"]]) == "guide-box")
  legend <- tmp[["grobs"]][[leg]]
  tmp_file <- tempfile(pattern = "legend", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  grid::grid.newpage()
  grid::grid.draw(legend)
  legend_plot <- grDevices::recordPlot()
  dev.off()
  removed <- file.remove(tmp_file)
  ret <- list(
      colors = plot[["data"]][, c("condition", "batch", "colors")],
      plot = legend_plot)
  return(ret)
}

## I thought multiplot() was a part of ggplot(), but no, weird:
## http://stackoverflow.com/questions/24387376/r-wired-error-could-not-find-function-multiplot
## Also found at:
## http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
#' Make a grid of plots.
#'
#' @param plots List of plots
#' @param file File to write to
#' @param cols Number of columns in the grid
#' @param layout Set the layout specifically
#' @return a multiplot!
#' @export
plot_multiplot <- function(plots, file, cols = NULL, layout = NULL) {
  ## Make a list from the ... arguments and plotlist
  ##  plots <- c(list(...), plotlist)
  num_plots <- length(plots)
  if (is.null(cols)) {
    cols <- ceiling(sqrt(length(plots)))
  }
  ## If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    ## Make the panel
    ## ncol: Number of columns of plots
    ## nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(num_plots / cols)),
                     ncol = cols, nrow = ceiling(num_plots / cols))
  }

  if (num_plots==1) {
    print(plots[[1]])
  } else {
    ## Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
                                 layout = grid::grid.layout(nrow(layout), ncol(layout))))
    ## Make each plot, in the correct location
    for (i in seq_len(num_plots)) {
      ## Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx[["row"]],
                                            layout.pos.col = matchidx[["col"]]))
    }
  }
}

## EOF
