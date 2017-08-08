## plot_distribution.r: A few plots to describe data distributions
## Currently this includes boxplots, density plots, and qq plots.

#' Make a ggplot boxplot of a set of samples.
#'
#' Boxplots and density plots provide complementary views of data distributions.  The general idea
#' is that if the box for one sample is significantly shifted from the others, then it is likely an
#' outlier in the same way a density plot shifted is an outlier.
#'
#' @param data Expt or data frame set of samples.
#' @param colors Color scheme, if not provided will make its own.
#' @param names Another version of the sample names for printing.
#' @param scale Whether to log scale the y-axis.
#' @param title A title!
#' @param ... More parameters are more fun!
#' @return Ggplot2 boxplot of the samples.  Each boxplot
#' contains the following information: a centered line describing the
#' median value of counts of all genes in the sample, a box around the
#' line describing the inner-quartiles around the median (quartiles 2
#' and 3 for those who are counting), a vertical line above/below the
#' box which shows 1.5x the inner quartile range (a common metric of
#' the non-outliers), and single dots for each gene which is outside
#' that range.  A single dot is transparent.
#' @seealso \pkg{ggplot2} \pkg{reshape2} \link[ggplot2]{geom_boxplot}
#'  \code{\link[reshape2]{melt}} \code{\link[ggplot2]{scale_x_discrete}}
#' @examples
#' \dontrun{
#'  a_boxplot <- plot_boxplot(expt)
#'  a_boxplot  ## ooo pretty boxplot look at the lines
#' }
#' @export
plot_boxplot <- function(data, colors=NULL, names=NULL, title=NULL, scale=NULL, ...) {
    plot_env <- environment()
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
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }

    if (is.null(scale)) {
        if (max(data) > 10000) {
            message("This data will benefit from being displayed on the log scale.")
            message("If this is not desired, set scale='raw'")
            scale <- "log"
            negative_idx <- data < 0
            if (sum(negative_idx) > 0) {
                message("Some data are negative.  We are on log scale, setting them to 0.5.")
                data[negative_idx] <- 0.5
                message(paste0("Changed ", sum(negative_idx), " negative features."))
            }
            zero_idx <- data == 0
            if (sum(zero_idx) > 0) {
                message("Some entries are 0.  We are on log scale, setting them to 0.5.")
                data[zero_idx] <- 0.5
                message(paste0("Changed ", sum(zero_idx), " zero count features."))
            }
        } else {
            scale <- "raw"
        }
    }

    if (is.null(colors)) {
        colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(dim(data)[2])
    }
    data_matrix <- as.matrix(data)
    data[data < 0] <- 0 ## Likely only needed when using quantile norm/batch correction and it sets a value to < 0

    data[["id"]] <- rownames(data)
    dataframe <- reshape2::melt(data, id=c("id"))
    colnames(dataframe) <- c("gene", "variable", "value")
    ## The use of data= and aes() leads to no visible binding for global variable warnings
    ## I am not sure what to do about them in this context.
    boxplot <- ggplot2::ggplot(data=dataframe, ggplot2::aes_string(x="variable", y="value")) +
        sm(ggplot2::geom_boxplot(na.rm=TRUE,
                                 ggplot2::aes_string(fill="variable"),
                                 fill=colors, size=0.5,
                                 outlier.size=1.5,
                                 outlier.colour=ggplot2::alpha("black", 0.2))) +
        ggplot2::theme_bw() + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1)) +
        ggplot2::xlab("Sample") + ggplot2::ylab("Per-gene (pseudo)count distribution")
    if (!is.null(title)) {
        boxplot <- boxplot + ggplot2::ggtitle(title)
    }
    if (!is.null(names)) {
        boxplot <- boxplot + ggplot2::scale_x_discrete(labels=names)
    }

    if (scale == "log") {
        boxplot <- boxplot + ggplot2::scale_y_continuous(trans=scales::log2_trans())
    } else if (scale == "logdim") {
        boxplot <- boxplot + ggplot2::coord_trans(y="log2")
    } else if (isTRUE(scale)) {
        boxplot <- boxplot + ggplot2::scale_y_log10()
    }
    return(boxplot)
}

#' Create a density plot, showing the distribution of each column of data.
#'
#' Density plots and boxplots are cousins and provide very similar views of data distributions.
#' Some people like one, some the other.  I think they are both colorful and fun!
#'
#' @param data Expt, expressionset, or data frame.
#' @param colors Color scheme to use.
#' @param sample_names Names of the samples.
#' @param position How to place the lines, either let them overlap (identity), or stack them.
#' @param fill Fill the distributions?  This might make the plot unreasonably colorful.
#' @param scale Plot on the log scale?
#' @param title Title for the plot.
#' @param colors_by Factor for coloring the lines
#' @return Ggplot2 density plot!
#' @seealso \pkg{ggplot2}
#'  \code{\link[ggplot2]{geom_density}}
#' @examples
#' \dontrun{
#'  funkytown <- plot_density(data)
#' }
#' @export
plot_density <- function(data, colors=NULL, sample_names=NULL, position="identity", direct=TRUE,
                         fill=NULL, title=NULL, scale=NULL, colors_by="condition") {
    ## also position='stack'
    plot_env <- environment()
    data_class <- class(data)[1]
    if (data_class == "expt") {
        design <- pData(data)
        colors <- data[["colors"]]
        names <- data[["names"]]
        data <- exprs(data)
    } else if (data_class == "ExpressionSet") {
        data <- exprs(data)
        design <- pData(data)
    } else if (data_class == "matrix" | data_class == "data.frame") {
        data <- as.matrix(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }

    if (is.null(scale)) {
        if (max(data) > 10000) {
            message("This data will benefit from being displayed on the log scale.")
            message("If this is not desired, set scale='raw'")
            scale <- "log"
            negative_idx <- data < 0
            if (sum(negative_idx) > 0) {
                message("Some data are negative.  We are on log scale, setting them to 0.5.")
                data[negative_idx] <- 0.5
                message(paste0("Changed ", sum(negative_idx), " negative features."))
            }
            zero_idx <- data == 0
            if (sum(zero_idx) > 0) {
                message("Some entries are 0.  We are on log scale, setting them to 0.5.")
                data[zero_idx] <- 0.5
                message(paste0("Changed ", sum(zero_idx), " zero count features."))
            }
        } else {
            scale <- "raw"
        }
    }

    if (!is.null(sample_names)) {
        colnames(data) <- make.names(sample_names, unique=TRUE)
    }
    ## If the columns lose the connectivity between the sample and values, then
    ## the ggplot below will fail with env missing.
    melted <- reshape2::melt(data)
    if (dim(melted)[2] == 3) {
        colnames(melted) <- c("id", "sample", "counts")
    } else if (dim(melted)[2] == 2) {
        colnames(melted) <- c("sample", "counts")
    } else {
        stop("Could not properly melt the data.")
    }
    densityplot <- NULL
    if (is.null(fill)) {
        densityplot <- ggplot2::ggplot(data=melted,
                                       ggplot2::aes_string(x="counts", colour="sample"),
                                       environment=plot_env)
    } else {
        fill <- "sample"
        densityplot <- ggplot2::ggplot(data=melted,
                                       ggplot2::aes_string(x="counts", colour="sample", fill="fill"),
                                       environment=plot_env)
    }

    densityplot <- densityplot +
        ggplot2::geom_density(ggplot2::aes_string(x="counts", y="..count..", fill="sample"),
                              position=position, na.rm=TRUE) +
        ggplot2::ylab("Number of genes.") + ggplot2::xlab("Number of hits/gene.") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.key.size=ggplot2::unit(0.3, "cm"))
    if (!is.null(title)) {
        densityplot <- densityplot + ggplot2::ggtitle(title)
    }

    if (scale == "log") {
        densityplot <- densityplot + ggplot2::scale_x_continuous(trans=scales::log2_trans())
    } else if (scale == "logdim") {
        densityplot <- densityplot + ggplot2::coord_trans(x="log2")
    } else if (isTRUE(scale)) {
        densityplot <- densityplot + ggplot2::scale_x_log10()
    }

    if (!is.null(colors_by)) {
        densityplot <- densityplot + ggplot2::scale_colour_manual(values=as.character(colors)) +
            ggplot2::scale_fill_manual(values=ggplot2::alpha(as.character(colors), 0.1))
    }

    if (isTRUE(direct)) {
        densityplot <- directlabels::direct.label(densityplot)
    }
    return(densityplot)
}

#' Quantile/quantile comparison of the mean of all samples vs. each sample.
#'
#' This allows one to visualize all individual data columns against the mean of all columns of data
#' in order to see if any one is significantly different than the cloud.
#'
#' @param data Expressionset, expt, or dataframe of samples.
#' @param labels What kind of labels to print?
#' @return List containing:
#'  logs = a recordPlot() of the pairwise log qq plots.
#'  ratios = a recordPlot() of the pairwise ratio qq plots.
#'  means = a table of the median values of all the summaries of the qq plots.
#' @seealso \pkg{Biobase}
#' @export
plot_qq_all <- function(data, labels="short") {
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
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }

    sample_data <- data[, c(1, 2)]
    ## This is bizarre, performing this operation with transform fails when called from a function
    ## but works fine when called interactively, wtf indeed?
    ##    sample_data = transform(sample_data, mean=rowMeans(plot_df))
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
        message(paste("Making plot of ", ith, "(", i, ") vs. a sample distribution.", sep=""))
        tmpdf <- data.frame("ith"=data[, i],
                            "mean"=sample_data[["mean"]])
        colnames(tmpdf) <- c(ith, "mean")
        tmpqq <- plot_single_qq(tmpdf, x=1, y=2, labels=labels)
        logs[[count]] <- tmpqq[["log"]]
        ratios[[count]] <- tmpqq[["ratio"]]
        means[[count]] <- tmpqq[["summary"]][["Median"]]
        count <- count + 1
    }
    result <- plot_multiplot(logs)
    log_plots <- grDevices::recordPlot()
    plot_multiplot(ratios)
    ratio_plots <- grDevices::recordPlot()
    plots <- list(logs=log_plots, ratios=ratio_plots, medians=means)
    return(plots)
}

#' Perform a qqplot between two columns of a matrix.
#'
#' Given two columns of data, how well do the distributions match one another?  The answer to that
#' question may be visualized through a qq plot!
#'
#' @param data Data frame/expt/expressionset.
#' @param x First column to compare.
#' @param y Second column to compare.
#' @param labels Include the lables?
#' @return a list of the logs, ratios, and mean between the plots as ggplots.
#' @seealso \pkg{Biobase}
#' @export
plot_single_qq <- function(data, x=1, y=2, labels=TRUE) {
    plot_env <- environment()
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
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
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
        y_string <- paste(xlabel, " : ", ylabel, sep="")
    } else {
        y_string <- paste("Ratio of sorted ", xlabel, " and ", ylabel, ".", sep="")
    }
    ratio_plot <- ggplot2::ggplot(ratio_df,
                                  ggplot2::aes_string(x="increment", y="ratio"),
                                  environment=plot_env) +
        ggplot2::geom_point(colour=sm(grDevices::densCols(ratio_df[["ratio"]])), stat="identity",
                            size=1, alpha=0.2, na.rm=TRUE) +
        ggplot2::scale_y_continuous(limits=c(0, 2))
    if (isTRUE(labels)) {
        ratio_plot <- ratio_plot +
            ggplot2::xlab("Sorted gene") +
            ggplot2::ylab(y_string) +
            ggplot2::theme(legend.position="none")
    } else if (labels == "short") {
        ratio_plot <- ratio_plot +
            ggplot2::ylab(y_string) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks=ggplot2::element_blank(),
                           axis.title.x=ggplot2::element_blank(),
                           legend.position="none",
                           panel.background=ggplot2::element_blank(),
                           panel.border=ggplot2::element_blank(),
                           panel.grid.major=ggplot2::element_blank(),
                           panel.grid.minor=ggplot2::element_blank(),
                           plot.background=ggplot2::element_blank())
    } else {
        ratio_plot <- ratio_plot + ggplot2::theme_bw() +
            ggplot2::theme(axis.line=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks=ggplot2::element_blank(),
                           axis.title.x=ggplot2::element_blank(),
                           axis.title.y=ggplot2::element_blank(),
                           legend.position="none",
                           panel.background=ggplot2::element_blank(),
                           panel.border=ggplot2::element_blank(),
                           panel.grid.major=ggplot2::element_blank(),
                           panel.grid.minor=ggplot2::element_blank(),
                           plot.background=ggplot2::element_blank())
    }

    log_df <- data.frame(cbind(log(ratio_df[["sorted_x"]] + 1.5),
                               log(ratio_df[["sorted_y"]] + 1.5)))
    gg_max <- max(log_df)
    colnames(log_df) <- c(xlabel, ylabel)
    log_df[["sub"]] <- log_df[, 1] - log_df[, 2]
    sorted_x <- as.vector(log_df[, 1])
    log_ratio_plot <- ggplot2::ggplot(log_df,
                                      ggplot2::aes_string(x="get(xlabel)", y="get(ylabel)"),
                                      environment=plot_env) +
    ggplot2::geom_point(colour=grDevices::densCols(x=sorted_x), stat="identity") +
    ggplot2::scale_y_continuous(limits=c(0, gg_max)) +
    ggplot2::scale_x_continuous(limits=c(0, gg_max))
    if (isTRUE(labels)) {
        log_ratio_plot <- log_ratio_plot +
            ggplot2::xlab(paste("log sorted ", xlabel)) +
            ggplot2::ylab(paste("log sorted ", ylabel)) +
            ggplot2::theme(legend.position="none")
    } else if (labels == "short") {
        log_ratio_plot <- log_ratio_plot +
            ggplot2::xlab("gene") +
            ggplot2::ylab(y_string) +
            ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks=ggplot2::element_blank(),
                           axis.title.x=ggplot2::element_blank(),
                           legend.position="none",
                           panel.background=ggplot2::element_blank(),
                           panel.border=ggplot2::element_blank(),
                           panel.grid.major=ggplot2::element_blank(),
                           panel.grid.minor=ggplot2::element_blank(),
                           plot.background=ggplot2::element_blank())
    } else {
        log_ratio_plot <- log_ratio_plot +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.line=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks=ggplot2::element_blank(),
                           axis.title.x=ggplot2::element_blank(),
                           axis.title.y=ggplot2::element_blank(),
                           legend.position="none",
                           panel.background=ggplot2::element_blank(),
                           panel.border=ggplot2::element_blank(),
                           panel.grid.major=ggplot2::element_blank(),
                           panel.grid.minor=ggplot2::element_blank(),
                           plot.background=ggplot2::element_blank())
    }
    ratio_plot <- ratio_plot + ggplot2::theme_bw()
    log_ratio_plot <- log_ratio_plot + ggplot2::theme_bw()
    log_summary <- summary(log_df$sub)
    qq_plots <- list(ratio=ratio_plot, log=log_ratio_plot, summary=log_summary)
    return(qq_plots)
}

#' Perform qq plots of every column against every other column of a dataset.
#'
#' This function is stupid, don't use it.  It makes more sense to just use plot_qq, however I am not
#' quite read to delete this function yet.
#'
#' @param data Dataframe to perform pairwise qqplots with.
#' @return List containing the recordPlot() output of the ratios, logs, and means among samples.
#' @seealso \pkg{Biobase}
#' @export
plot_qq_all_pairwise <- function(data) {
    data_class <- class(data)[1]
    names <- NULL
    if (data_class == "expt") {
        design <- pData(data)
        colors <- data[["colors"]]
        names <- data[["names"]]
        data <- exprs(data)
    } else if (data_class == "ExpressionSet") {
        data <- exprs(data)
        design <- pData(data)
    } else if (data_class == "matrix" | data_class == "data.frame") {
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }

    logs <- list()
    ratios <- list()
    rows <- length(colnames(data))
    means <- matrix(nrow=rows, ncol=rows)
    count <- 1
    for (i in 1:rows) {
        for (j in 1:rows) {
            ith <- colnames(data)[i]
            jth <- colnames(data)[j]
            message(paste("Making plot of ", ith, "(", i, ") vs. ", jth, "(", j, ") as element: ", count, ".", sep=""))
            tmp <- plot_qq_plot(data, x=i, y=j, labels=names)
            logs[[count]] <- tmp$log
            ratios[[count]] <- tmp$ratio
            means[i, j] <- tmp$summary[["Mean"]]
            count <- count + 1
        }
    }
    plot_multiplot(logs)
    log_plots <- grDevices::recordPlot()
    plot_multiplot(ratios)
    ratio_plots <- grDevices::recordPlot()
    heatmap.3(means, trace="none")
    means_heatmap <- grDevices::recordPlot()
    plots <- list(logs=log_plots, ratios=ratio_plots, means=means_heatmap)
    return(plots)
}

## EOF
