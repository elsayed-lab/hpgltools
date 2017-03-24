## Note to self, I think for future ggplot2 plots, I must start by creating the data frame
## Then cast every column in it explicitly, and only then invoke ggplot(data=df ...)

## If I see something like:
## 'In sample_data$mean = means : Coercing LHS to a list'
## That likely means that I was supposed to have data in the
## data.frame() format, but instead it is a matrix.  In functions
## where this is a danger, it is a likely good idea to cast it as a
## data frame.

#' Make lots of graphs!
#'
#' Plot out a set of metrics describing the state of an experiment
#' including library sizes, # non-zero genes, heatmaps, boxplots,
#' density plots, pca plots, standard median distance/correlation, and
#' qq plots.
#'
#' @param expt  an expt to process
#' @param cormethod   the correlation test for heatmaps.
#' @param distmethod define the distance metric for heatmaps.
#' @param title_suffix   text to add to the titles of the plots.
#' @param qq   include qq plots
#' @param ma   include pairwise ma plots
#' @param ... extra parameters optionally fed to the various plots
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
#'   \item qq = a recordPlotted() view comparing the quantile/quantiles between the mean of all data and every raw sample
#'   \item density = a ggplot2 view of the density of each raw sample (this is complementary but more fun than a boxplot)
#' }
#' @seealso \pkg{Biobase} \pkg{ggplot2} \pkg{grDevices} \pkg{gplots}
#'  \code{\link[Biobase]{exprs}} \code{\link{hpgl_norm}} \code{\link{plot_nonzero}} \code{\link{plot_libsize}}
#'  \code{\link{plot_boxplot}} \code{\link{plot_corheat}} \code{\link{plot_sm}} \code{\link{plot_disheat}}
#'  \code{\link{plot_pca}} \code{\link{plot_qq_all}} \code{\link{plot_pairwise_ma}}
#' @examples
#' \dontrun{
#'  toomany_plots <- graph_metrics(expt)
#'  toomany_plots$pcaplot
#'  norm <- normalize_expt(expt, convert="cpm", batch=TRUE, filter_low=TRUE,
#'                         transform="log2", norm="rle")
#'  holy_asscrackers <- graph_metrics(norm, qq=TRUE, ma=TRUE)
#'  ## good luck, you are going to be waiting a while for the ma plots to print!
#' }
#' @export
graph_metrics <- function(expt, cormethod="pearson", distmethod="euclidean", title_suffix=NULL,
                          qq=NULL, ma=NULL, ...) {
    ## First gather the necessary data for the various plots.
    old_options <- options(scipen=10)
    ##old_options <- options(device = function(...) {
    ##    .Call("R_GD_nullDevice", PACKAGE = "grDevices")
    ##})
    nonzero_title <- "Non zero genes"
    libsize_title <- "Library sizes"
    boxplot_title <- "Boxplot"
    corheat_title <- "Correlation heatmap"
    smc_title <- "Standard Median Correlation"
    disheat_title <- "Distance heatmap"
    smd_title <- "Standard Median Distance"
    pca_title <- "Principle Component Analysis"
    dens_title <- "Density plot"
    if (!is.null(title_suffix)) {
        nonzero_title <- paste0(nonzero_title, ": ", title_suffix)
        libsize_title <- paste0(libsize_title, ": ", title_suffix)
        boxplot_title <- paste0(boxplot_title, ": ", title_suffix)
        corheat_title <- paste0(corheat_title, ": ", title_suffix)
        smc_title <- paste0(smc_title, ": ", title_suffix)
        disheat_title <- paste0(disheat_title, ": ", title_suffix)
        smd_title <- paste0(smd_title, ": ", title_suffix)
        pca_title <- paste0(pca_title, ": ", title_suffix)
        dens_title <- paste0(dens_title, ": ", title_suffix)
    }
    message("Graphing number of non-zero genes with respect to CPM by library.")
    nonzero_plot <- try(plot_nonzero(expt, title=nonzero_title, ...))
    message("Graphing library sizes.")
    libsize_plot <- try(plot_libsize(expt, title=libsize_title, ...))
    message("Graphing a boxplot.")
    boxplot <- try(plot_boxplot(expt, title=boxplot_title, ...))
    message("Graphing a correlation heatmap.")
    corheat <- try(plot_corheat(expt, method=cormethod, title=corheat_title, ...))
    message("Graphing a standard median correlation.")
    smc <- try(plot_sm(expt, method=cormethod, title=smc_title, ...))
    message("Graphing a distance heatmap.")
    disheat <- try(plot_disheat(expt, method=distmethod, title=disheat_title, ...))
    message("Graphing a standard median distance.")
    smd <- try(plot_sm(expt, method=distmethod, title=smd_title, ...))
    message("Graphing a PCA plot.")
    pca <- try(plot_pca(expt, title=pca_title, ...))
    message("Plotting a density plot.")
    density <- try(plot_density(expt, title=dens_title))
    message("Printing a color to condition legend.")
    legend <- try(plot_legend(pca$plot))

    qq_logs <- NULL
    qq_ratios <- NULL
    if (isTRUE(qq)) {
        message("QQ plotting!")
        qq_plots <- try(suppressWarnings(plot_qq_all(expt)))
        qq_logs <- qq_plots$logs
        qq_ratios <- qq_plots$ratios
    }

    ma_plots <- NULL
    if (isTRUE(ma)) {
        message("Many MA plots!")
        ma_plots <- try(suppressWarnings(plot_pairwise_ma(expt)))
    }

    ret_data <- list(
        "nonzero" = nonzero_plot,
        "libsize" = libsize_plot,
        "boxplot" = boxplot,
        "corheat" = corheat$plot,
        "smc" = smc,
        "disheat" = disheat$plot,
        "smd" = smd,
        "pcaplot" = pca$plot,
        "pcatable" = pca$table,
        "pcares" = pca$res,
        "pcavar" = pca$variance,
        "density" = density,
        "legend" = legend,
        "qqlog" = qq_logs,
        "qqrat" = qq_ratios,
        "ma" = ma_plots)
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
    grid::grid.newpage()
    grid::grid.draw(legend)
    legend_plot <- grDevices::recordPlot()
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
#' @param plots  a list of plots
#' @param file  a file to write to
#' @param cols   the number of columns in the grid
#' @param layout  set the layout specifically
#' @return a multiplot!
#' @export
plot_multiplot <- function(plots, file, cols=NULL, layout=NULL) {
  ## Make a list from the ... arguments and plotlist
  ##  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  if (is.null(cols)) {
      cols <- ceiling(sqrt(length(plots)))
  }
  ## If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
      ## Make the panel
      ## ncol: Number of columns of plots
      ## nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                       ncol=cols, nrow=ceiling(numPlots / cols))
  }

  if (numPlots==1) {
      print(plots[[1]])
  } else {
      ## Set up the page
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow(layout), ncol(layout))))
      ## Make each plot, in the correct location
      for (i in 1:numPlots) {
          ## Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind=TRUE))
          print(plots[[i]], vp=grid::viewport(layout.pos.row=matchidx$row,
                                              layout.pos.col=matchidx$col))
      }
  }
}

## EOF
