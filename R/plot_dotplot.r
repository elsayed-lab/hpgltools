## Time-stamp: <Mon Apr 25 16:32:25 2016 Ashton Trey Belew (abelew@gmail.com)>

## plot_dotplot.r: Dotplots in various contexts, currently just smc/smd

#'   Make an R plot of the standard median correlation among samples.
#'
#' This was written by a mix of Kwame Okrah <kokrah at gmail dot com>, Laura
#' Dillon <dillonl at umd dot edu>, and Hector Corrada Bravo <hcorrada at umd dot edu>
#'
#' @param data  an expt, expressionset, or data frame.
#' @param colors   a color scheme
#' @param method   a correlation method to use.
#' @param names   use pretty names for the samples?
#' @param title   title for the graph.
#' @param ... more parameters to make you happy
#' @return a recordPlot() of the standard median pairwise correlation
#' among the samples.  This will also write to an
#' open device.  The resulting plot measures the median correlation of
#' each sample among its peers.  It notes 1.5* the interquartile range
#' among the samples and makes a horizontal line at that correlation
#' coefficient.  Any sample which falls below this line is considered
#' for removal because it is much less similar to all of its peers.
#' @seealso \link{hpgl_cor} \link[matrixStats]{rowMedians}
#' \link[stats]{quantile} \link{diff} \link[grDevices]{recordPlot}
#' @examples
#' \dontrun{
#'  smc_plot = hpgl_smc(expt=expt)
#' }
#' @export
hpgl_smc <- function(data, colors=NULL, method="pearson", names=NULL, title=NULL, ...) {
    data_class <- class(data)[1]
    arglist <- list(...)
    if (data_class == "expt") {
        design <- data$design
        colors <- data$colors
        names <- data$names
        data <- Biobase::exprs(data$expressionset)
    } else if (data_class == 'ExpressionSet') {
        data <- Biobase::exprs(data)
    } else if (data_class == 'matrix' | data_class == 'data.frame') {
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }

    chosen_palette <- "Dark2"
    if (!is.null(arglist$palette)) {
        chosen_palette <- arglist$palette
    }

    if (is.null(names)) {
        names <- colnames(data)
    }

    if (is.null(colors)) {
        colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(ncol(data), chosen_palette))(ncol(data))
    }
    colors <- as.character(colors)
    correlations <- hpgl_cor(data, method=method)
    cor_median <- matrixStats::rowMedians(correlations)
    cor_spread <- stats::quantile(cor_median, p=c(1,3)/4)
    cor_iqr <- diff(cor_spread)
    outer_limit <- cor_spread[1] - 1.5 * cor_iqr
    ylimit <- c(pmin(min(cor_median), outer_limit), max(cor_median))
    plot(cor_median, xaxt="n", ylim=ylimit,
         xlab="", main="", ylab="Median pairwise correlation",
         ## col=hpgl_colors, pch=16, cex=1.5)
         bg=colors, col="black", pch=21, cex=1.5)
    title(title)
    axis(side=1, at=seq(along=cor_median), labels=names, las=2)
    abline(h=outer_limit, lty=2)
    abline(v=1:length(names), lty=3, col="black")
    hpgl_smc_plot <- grDevices::recordPlot()
    return(hpgl_smc_plot)
}

#'   Make an R plot of the standard median distance among samples.
#'
#' @param data an expt/expressionset/data frame of samples.
#' @param colors   a color scheme
#' @param method   a distance metric to use.
#' @param names   use pretty names for the samples?
#' @param title   title for the graph.
#' @param ... parameters make me happy
#' @return smd_plot a recordPlot of plot.  This will also write to an
#' open device.  This plot takes the median distance of each sample
#' with all of its peers.  It then calculates 1.5* the interquartile
#' range of distances.  Any sample which has a median distance greater
#' than this is considered for removal.
#' @seealso \code{\link{dist}}, \code{\link{quantile}},
#' \code{\link{diff}}, \code{\link{recordPlot}}
#' @examples
#' \dontrun{
#'  smd_plot = hpgl_smd(expt=expt)
#' }
#' @export
hpgl_smd <- function(data, colors=NULL, names=NULL, method="euclidean", title=NULL, ...) {
    arglist <- list(...)
    data_class <- class(data)[1]
    if (data_class == "expt") {
        design <- data$design
        colors <- data$colors
        names <- data$names
        data <- Biobase::exprs(data$expressionset)
    } else if (data_class == "ExpressionSet") {
        data <- Biobase::exprs(data)
    } else if (data_class == "matrix" | data_class == "data.frame") {
        data <- as.matrix(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }

    chosen_palette <- "Dark2"
    if (!is.null(arglist$palette)) {
        chosen_palette <- arglist$palette
    }

    if (is.null(names)) {
        names <- colnames(data)
    }
    if (is.null(colors)) {
        colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(ncol(data), chosen_palette))(ncol(data))
    }
    colors <- as.character(colors)
    dists <- as.matrix(dist(t(data)), method=method)
    dist_median <- matrixStats::rowMedians(dists)
    dist_spread <- stats::quantile(dist_median, p=c(1,3)/4)
    dist_iqr <- diff(dist_spread)
    outer_limit <- dist_spread[2] + (1.5 * dist_iqr)
##    ylimit = c(min(dist_median), max(dist_median))
    ylimit <- c(min(dist_median), pmax(max(dist_median), outer_limit))
    plot(dist_median, xaxt="n", ylim=ylimit,
         xlab="", main="", ylab="Median pairwise distance",
         ## col=hpgl_colors, pch=16, cex=1.5)
         bg=colors, col="black", pch=21, cex=1.5)
    title(title)
    axis(side=1, at=seq(along=dist_median), labels=names, las=2)
    abline(h=outer_limit, lty=2)
    abline(v=1:length(names), lty=3, col="black")
    hpgl_smd_plot <- grDevices::recordPlot()
    return(hpgl_smd_plot)
}

#' make a dotplot of some categorised factors and a set of SVs (for other factors)
#'
#' This should make a quick df of the factors and surrogates and plot them.
#'
#' @param expt an experiment from which to acquire the design, counts, etc
#' @param svest a set of surrogate variable estimations from sva/svg
#'        or batch estimates
#' @param chosen_factor a factor to compare against
#' @param factor_type this may be a factor or range, it is intended to plot
#'        a scatterplot if it is a range, a dotplot if a factor
#' @examples
#' \dontrun{
#' estimate_vs_snps <- plot_svfactor(start, surrogate_estimate, "snpcategory")
#' }
plot_svfactor <- function(expt, svest, chosen_factor="snpcategory", factor_type="factor") {
    chosen <- expt[["design"]][[chosen_factor]]
    sv_df <- data.frame(
        "adjust" = svest$model_adjust[, 1],  ## Take a single estimate from compare_estimates()
        "factors" = chosen,
        "samplenames" = rownames(expt[["design"]])
    )
    samplenames <- rownames(expt[["design"]])
    my_colors <- expt[["colors"]]

    sv_melted <- reshape2::melt(sv_df, idvars="factors")
    sv_plot <- ggplot2::ggplot(sv_df, ggplot2::aes_string(x="factors", y="value")) +
        ggplot2::geom_dotplot(binaxis="y", stackdir="center", binpositions="all", colour="black", fill=my_colors) +
        ggplot2::xlab(paste0("Experimental factor: ", chosen_factor)) +
        ggplot2::ylab(paste0("1st surrogate variable estimation")) +
        ggplot2::theme_bw()

    return(sv_plot)
}
