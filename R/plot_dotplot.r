## Time-stamp: <Fri Apr 29 23:01:14 2016 Ashton Trey Belew (abelew@gmail.com)>

## plot_dotplot.r: Dotplots in various contexts, currently just smc/smd

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
#' @export
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
    minval <- min(sv_df$adjust)
    maxval <- max(sv_df$adjust)
    my_binwidth <- (maxval - minval) / 40
    sv_plot <- ggplot2::ggplot(sv_melted, ggplot2::aes_string(x="factors", y="value")) +
        ggplot2::geom_dotplot(binwidth=my_binwidth, binaxis="y", stackdir="center", binpositions="all", colour="black", fill=my_colors) +
        ggplot2::xlab(paste0("Experimental factor: ", chosen_factor)) +
        ggplot2::ylab(paste0("1st surrogate variable estimation")) +
        ggplot2::geom_text(ggplot2::aes_string(x="factors", y="value", label="strains"), angle=45, size=3, vjust=2) +
        ggplot2::theme_bw()
    return(sv_plot)
}

#' make a dotplot of some categorised factors and a set of principle components.
#'
#' This should make a quick df of the factors and PCs and plot them.
#'
#' @param pc_df Df of principle components.
#' @param exp_factor Experimental factor to compare against.
#' @examples
#' \dontrun{
#' estimate_vs_pcs <- plot_pcfactor(pcs, times)
#' }
#' @export
plot_pcfactor <- function(pc_df, exp_factpr) {
    chosen <- expt[["design"]][[chosen_factor]]
    sv_df <- data.frame(
        "adjust" = svest$model_adjust[, 1],  ## Take a single estimate from compare_estimates()
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
        ggplot2::geom_dotplot(binwidth=my_binwidth, binaxis="y", stackdir="center", binpositions="all", colour="black", fill=my_colors) +
        ggplot2::xlab(paste0("Experimental factor: ", chosen_factor)) +
        ggplot2::ylab(paste0("1st surrogate variable estimation")) +
        ggplot2::geom_text(ggplot2::aes_string(x="factors", y="value", label="strains"), angle=45, size=3, vjust=2) +
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
plot_sm <- function(data, colors=NULL, method="pearson", names=NULL, title=NULL, ...) {
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

    properties <- NULL
    if (method == "pearson" | method == "spearman" | method == "robust") {
        message("Performing correlation.")
        properties <- hpgl_cor(data, method=method)
    } else {
        message("Performing distance.")
        properties <- as.matrix(dist(t(data)))
    }

    prop_median <- matrixStats::rowMedians(properties)
    prop_spread <- stats::quantile(prop_median, p=c(1,3)/4)
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

    sm_df <- data.frame(
        "sample" = rownames(properties),
        "sm" = prop_median,
        "color" = colors)

    minval <- min(sm_df[["sm"]])
    maxval <- max(sm_df[["sm"]])
    my_binwidth <- (maxval - minval) / 40

    sm_plot <- ggplot2::ggplot(sm_df, ggplot2::aes_string(x="sample", y="sm")) +
        ggplot2::geom_hline(color="red", yintercept=ylimit, size=2) +
        ggplot2::geom_dotplot(binwidth=my_binwidth, binaxis="y", stackdir="center", binpositions="all", colour="black", fill=sm_df[["color"]]) +
        ggplot2::ylab(paste0("Standard Median ", method)) +
        ggplot2::xlab(paste0("Sample")) +
        ggplot2::ggtitle(title) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

    return(sm_plot)
}

## EOF
