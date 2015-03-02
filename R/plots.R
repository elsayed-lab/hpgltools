## Time-stamp: <Sun Mar  1 17:05:42 2015 Ashton Trey Belew (abelew@gmail.com)>

#' Make a bunch of graphs describing the state of an experiment
#' before/after normalization.
#'
#' @param expt an expt experiment containing the data, design, and colors
#' @param norm_type normalization strategy for the data.  Defaults to
#' quantile.
#' @param filter whether to log2/10 filter the data.  Defaults to
#' log2.
#' @param out_type whether to cpm/rpkm the data.  Defaults to cpm.
#' @param filter_low whether to low-count filter the data.  Defaults
#' to TRUE.
#' @param cormethod define the correlation test for heatmaps.
#' Defaults to pearson (Available: pearson, spearman, kendal, robust)
#' @param distmethod define the distance metric for heatmaps.
#' Defaults to euclidean (Lots are available, I don't understand them.)
#' @param ... extra parameters optionally fed to the various plots
#' 
#' @return a loooong list of plots including the following:
#'   nonzero = a ggplot2 plot of the non-zero genes vs library size
#'   libsize = a ggplot2 bar plot of the library sizes
#'   raw_boxplot = a ggplot2 boxplot of the raw data
#'   norm_boxplot = a ggplot2 boxplot of the normalized data
#'   raw_corheat = a recordPlot()ed pairwise correlation heatmap of the raw data
#'   norm_corheat = a recordPlot()ed pairwise correlation heatmap of the normalized data
#'   raw_smc = a recordPlot()ed view of the standard median pairwise correlation of the raw data
#'   norm_smc = a recordPlot()ed view of the standard median pairwise correlation of the normalized data
#'   raw_disheat = a recordPlot()ed pairwise euclidean distance heatmap of the raw data
#'   norm_disheat = a recordPlot()ed pairwise euclidean distance heatmap of the normalized data
#'   raw_smd = a recordPlot()ed view of the standard median pairwise distance of the raw data
#'   norm_smd = a recordPlot()ed view of the standard median pairwise distance of the normalized data
#'   raw_pcaplot = a recordPlot()ed PCA plot of the raw samples
#'   norm_pcaplot = a recordPlot()ed PCA plot of the normalized samples
#'   raw_pcatable = a table describing the relative contribution of condition/batch of the raw data
#'   norm_pcatable = a table describing the relative contribution of condition/batch of the normalized data
#'   raw_pcares =  a table describing the relative contribution of condition/batch of the raw data
#'   norm_pcares = a table describing the relative contribution of condition/batch of the normalized data
#'   raw_pcavar = a table describing the variance of the raw data
#'   norm_pcavar = a table describing the variance of the normalized data
#'   raw_qq = a recordPlotted() view comparing the quantile/quantiles between the mean of all data and every raw sample
#'   norm_qq = a recordPlotted() view comparing the quantile/quantiles between the mean of all data and every normalized sample
#'   raw_density = a ggplot2 view of the density of each raw sample (this is complementary but more fun than a boxplot)
#'   norm_density = a ggplot2 view of the density of each normalized
#'   batch_boxplot = a ggplot2 boxplot of the data after calling limma's removeBatchEffect()
#'   batch_disheat = a recordPlot() of a distance heatmap after calling limma's removeBatchEffect()
#'   batch_corheat = a recordPlot() of a correlation heatmap after calling limma's removeBatchEffect()
#'   batch_pcaplot = a ggplot2 PCA plot after removeBatchEffect()
#'   batch_pcatable = a pca table of the removeBatchEffect()'d data
#'   batch_pcares = a pcares table of the removeBatchEffect()'d data
#'   batch_pcavar = a pcavar table of the removeBatchEffect()'d data
#' 
#' @seealso \code{\link{exprs}}, \code{\link{hpgl_norm}},
#' \code{\link{graph_nonzero}}, \code{\link{hpgl_libsize}},
#' \code{\link{hpgl_boxplot}}, \code{\link{hpgl_corheat}},
#' \code{\link{hpgl_smc}}, \code{\link{hpgl_disheat}},
#' \code{\link{hpgl_smd}}, \code{\link{hpgl_pca}},
#' \code{\link{replayPlot}}, \code{\link{recordPlot}}
#' 
#' @export
#' @examples
#' ## toomany_plots = graph_metrics(expt)
#' ## testnorm = graph_metrics(expt, norm_type="tmm", filter="log2", out_type="rpkm", cormethod="robust")
graph_metrics = function(expt, transform="log2", norm="quant", convert="cpm", filter_low=TRUE, cormethod="pearson", distmethod="euclidean", ...) {
    expt_design = expt$design
    expt_colors = expt$colors
    expt_names = expt$names
    expt_raw_data = Biobase::exprs(expt$expressionset)
    expt_norm_data = hpgltools::hpgl_norm(expt=expt, transform=transform, norm=norm, convert=convert, filter_low=filter_low, ...)$counts
    message("Graphing number of non-zero genes with respect to CPM by library.")
    nonzero_plot = hpgltools::graph_nonzero(expt=expt, title="Non zero genes.", ...)
    message("Graphing library sizes.")
    libsize_plot = hpgltools::hpgl_libsize(expt=expt, title="Library sizes.", ...)
    message("Graphing a raw data boxplot on log scale.")
    raw_boxplot = hpgltools::hpgl_boxplot(expt=expt, title="Boxplot of log(raw data).", scale="log", ...)
    message("Graphing a normalized boxplot.")
    norm_boxplot = hpgltools::hpgl_boxplot(df=expt_norm_data, names=expt$names, colors=expt_colors, title="Boxplot of normalized data.", ...)
    message("Graphing a raw-data correlation heatmap.")
    raw_corheat = hpgltools::hpgl_corheat(expt=expt, method=cormethod, title="Correlation heatmap of raw data.", ...)
    message("Graphing a raw-data standard median correlation.")
    raw_smc = hpgltools::hpgl_smc(expt=expt, method=cormethod, title="Standard Median Correlation, raw data.", ...)
    message("Graphing a normalized correlation heatmap.")
    norm_corheat = hpgltools::hpgl_corheat(df=expt_norm_data, names=expt$names, colors=expt_colors, design=expt_design, method=cormethod, title="Corrleation heatmap of normalized data.", ...)
    message("Graphing a normalized standard median correlation.")
    norm_smc = hpgltools::hpgl_smc(df=expt_norm_data, names=expt$names, colors=expt_colors, method=cormethod, title="Standard Median Correlation, norm. data.", ...)
    message("Graphing a raw-data distance heatmap.")
    raw_disheat = hpgltools::hpgl_disheat(expt=expt, method=distmethod, title="Distance heatmap, raw data.", ...)
    message("Graphing a raw-data standard median distance.")
    raw_smd = hpgltools::hpgl_smd(expt=expt, method=distmethod, title="Standard Median Distance, raw data.", ...)
    message("Graphing a normalized distance heatmap.")
    norm_disheat = hpgltools::hpgl_disheat(df=expt_norm_data, names=expt$names, colors=expt_colors, design=expt_design, method=distmethod, title="Distance heatmap, norm. data.", ...)
    message("Graphing a normalized standard median distance.")
    norm_smd = hpgltools::hpgl_smd(df=expt_norm_data, names=expt$names, colors=expt_colors, method=distmethod, title="Standard Median Distance, norm. data.", ...)
    message("Graphing a PCA plot of the raw data.")
    raw_pca = try(hpgltools::hpgl_pca(expt=expt, fancy_labels=FALSE, title="PCA plot of raw data.", ...))
    message("Printing a qqplot of the raw data.")
    raw_qq = try(suppressWarnings(hpgltools::hpgl_qq_all(df=data.frame(exprs(expt$expressionset)))))
    raw_density = try(hpgltools::hpgl_density_plot(expt=expt, title="Density plot of raw data."))
    norm_density = try(hpgltools::hpgl_density_plot(df=expt_norm_data, title="Density plot of normalized data."))    
    message("Graphing a PCA plot of the normalized data.")
    norm_pca = try(hpgltools::hpgl_pca(df=expt_norm_data, names=expt$names, fancy_labels=FALSE, colors=expt_colors, design=expt_design, title="PCA plot of norm. data.", ...))
    message("Printing a qqplot of the normalized data.")
    norm_qq = try(suppressWarnings(hpgltools::hpgl_qq_all(df=expt_norm_data)))
    batch_removed = limma::removeBatchEffect(Biobase::exprs(expt$expressionset), batch=expt$batches)
    batch_boxplot = hpgltools::hpgl_boxplot(df=batch_removed, names=expt$names, colors=expt_colors, title="Boxplot of batch removed data.", scale="log", ...)
    batch_disheat = hpgltools::hpgl_disheat(df=batch_removed, names=expt$names, colors=expt_colors, design=expt_design, method=distmethod, title="Distance heatmap of batch removed data.", ...)
    batch_corheat = hpgltools::hpgl_corheat(df=batch_removed, names=expt$names, colors=expt_colors, design=expt_design, method=cormethod, title="Correlation heatmap of batch removed data.", ...)
    batch_pca = try(hpgltools::hpgl_pca(df=batch_removed, names=expt$names, fancy_labels=FALSE, colors=expt_colors, design=expt_design, title="PCA plot of batch removed data.", ...))
    
    ret_data = list(
        nonzero=nonzero_plot, libsize=libsize_plot, raw_boxplot=raw_boxplot,
        norm_boxplot=norm_boxplot, raw_corheat=raw_corheat, raw_smc=raw_smc,
        norm_corheat=norm_corheat, norm_smc=norm_smc, raw_disheat=raw_disheat,
        raw_smd=raw_smd, norm_disheat=norm_disheat, norm_smd=norm_smd,
        raw_pcaplot=raw_pca$plot, norm_pcaplot=norm_pca$plot,
        raw_pcatable=raw_pca$table, norm_pcatable=norm_pca$table,
        raw_pcares=raw_pca$res, norm_pcares=norm_pca$res,
        raw_pcavar=raw_pca$variance, norm_pcavar=norm_pca$variance,
        raw_qq=raw_qq, ##norm_qq=norm_qq,
        raw_density=raw_density, norm_density=norm_density,
        batch_boxplot=batch_boxplot, batch_disheat=batch_disheat, batch_corheat=batch_corheat,
        batch_pcaplot=batch_pca$plot, batch_pcatable=batch_pca$table,
        batch_pcares=batch_pca$res, batch_pcavar=batch_pca$variance
    )
    return(ret_data)
}

#' Make a ggplot graph of the number of non-zero genes by sample.  Made by Ramzi Temanni.
#'
#' @param expt an expt set of samples
#' @param df alternately a data frame which must be accompanied by
#' @param design a design matrix and
#' @param colors a color scheme
#' 
#' @return a ggplot2 plot of the number of non-zero genes with respect to each library's CPM
#' @seealso \code{\link{geom_point}}, \code{\link{geom_dl}}
#' 
#' @export
#' @examples
#' ## nonzero_plot = graph_nonzero(expt=expt)
#' ## nonzero_plot  ## ooo pretty
graph_nonzero = function(df=NULL, design=NULL, colors=NULL, expt=NULL, title=NULL, ...) {
    hpgl_env = environment()
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    if (is.null(expt)) {
        hpgl_design = design
        hpgl_colors = colors
        hpgl_df = df
    } else if (is.null(df)) {
        hpgl_design = expt$design
        hpgl_colors = expt$colors
        hpgl_df = exprs(expt$expressionset)
    } else {
        stop("Both df and expt are defined, choose one.")
    }
    if (is.null(expt$names)) {
        hpgl_labels = colnames(hpgl_df)
    } else {
        hpgl_labels = expt$names
    }
    hpgl_shapes = as.integer(hpgl_design$batch)
    non_zero = data.frame(id=colnames(hpgl_df),
        nonzero_genes=colSums(hpgl_df > 1),
        cpm=colSums(hpgl_df) * 1e-6,
        condition=hpgl_design$condition,
        batch=hpgl_design$batch)
    non_zero_plot = ggplot2::ggplot(data=non_zero, ggplot2::aes(x=cpm, y=nonzero_genes), environment=hpgl_env,
        colour=hpgl_colors, shape=hpgl_shapes) +
        geom_point(stat="identity", size=3, colour=hpgl_colors, shape=hpgl_shapes) +
        directlabels::geom_dl(aes(label=hpgl_labels), method="smart.grid", colour=hpgl_colors) +
        ylab("Number of non-zero genes observed.") +
        xlab("Observed CPM")
    if (!is.null(title)) {
        non_zero_plot = non_zero_plot + ggplot2::ggtitle(title)
    }
    non_zero_plot = non_zero_plot + ggplot2::theme(axis.ticks=element_blank(), axis.text.x=element_text(angle=90))
    return(non_zero_plot)
}

#' Make a ggplot graph of library sizes.
#'
#' @param expt an expt set of samples
#' @param df alternately a data frame which must be accompanied by
#' @param design a design matrix and
#' @param colors a color scheme
#' @param scale whether or not to log10 the y-axis
#' 
#' @return a ggplot2 bar plot of every sample's size
#' @seealso \code{\link{geom_bar}}, \code{\link{geom_text}},
#' \code{\link{prettyNum}}, \code{\link{scale_y_log10}}
#' 
#' @export
#' @examples
#' ## libsize_plot = hpgl_libsize(expt=expt)
#' ## libsize_plot  ## ooo pretty bargraph
hpgl_libsize = function(df=NULL, colors=NULL, expt=NULL, scale=TRUE, names=NULL, title=NULL, text=TRUE, ...) {
    hpgl_env = environment()
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    if (is.null(expt)) {
        hpgl_colors = colors
        hpgl_df = df
        hpgl_names = names
    } else if (is.null(df)) {
        hpgl_colors = expt$colors
        hpgl_df = Biobase::exprs(expt$expressionset)
        hpgl_names = expt$names
    } else {
        stop("Both df and expt are defined, choose one.")
    }
    if (is.null(hpgl_colors)) {
        hpgl_colors = colorRampPalette(brewer.pal(ncol(hpgl_df),"Dark2"))(ncol(hpgl_df))
    }
    hpgl_colors = as.character(hpgl_colors)
    tmp = data.frame(id=colnames(hpgl_df),
        sum=colSums(hpgl_df),
        colors=factor(hpgl_colors))
    tmp$order = factor(tmp$id, as.character(tmp$id))
    libsize_plot = ggplot2::ggplot(data=tmp, ggplot2::aes(x=order, y=sum),
        environment=hpgl_env, colour=tmp$colors) +
            geom_bar(aes(x=order), stat="identity", colour="black", fill=tmp$colors) +
            xlab("Sample ID") +
            ylab("Library size in (pseudo)counts.") +
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1.5, vjust=0.5))
    if (text == TRUE) {
        libsize_plot = libsize_plot + geom_text(ggplot2::aes(reorder(order), label=prettyNum(tmp$sum, big.mark=",")), angle=90, size=3, color="white", hjust=1.2)
    }
    if (!is.null(title)) {
        libsize_plot = libsize_plot + ggtitle(title)
    }
    if (scale == TRUE) {
        message("Adding log10")
        libsize_plot = libsize_plot + scale_y_log10()
    }
    if (!is.null(hpgl_names)) {
        libsize_plot = libsize_plot + scale_x_discrete(labels=hpgl_names)
    }
    return(libsize_plot)
}

#' Make a ggplot boxplot of a set of samples.
#'
#' @param expt an expt set of samples
#' @param df alternately a data frame which must be accompanied by
#' @param design a design matrix and
#' @param colors a color scheme
#' @param names a nicer version of the sample names
#' @param scale whether to log scale the y-axis
#' 
#' @return a ggplot2 boxplot of the samples.  Each boxplot
#' contains the following information: a centered line describing the
#' median value of counts of all genes in the sample, a box around the
#' line describing the inner-quartiles around the median (quartiles 2
#' and 3 for those who are counting), a vertical line above/below the
#' box which shows 1.5x the inner quartile range (a common metric of
#' the non-outliers), and single dots for each gene which is outside
#' that range.  A single dot is transparent.
#' @seealso \code{\link{geom_boxplot}}, \code{\link{melt}},
#' \code{\link{scale_x_discrete}}
#' 
#' @export
#' @examples
#' ## a_boxplot = hpgl_boxplot(expt=expt)
#' ## a_boxplot  ## ooo pretty boxplot look at the lines
hpgl_boxplot = function(df=NULL, colors_fill=NULL, names=NULL, expt=NULL, title=NULL, scale="raw", ...) {
    hpgl_env = environment()
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    if (is.null(expt)) {
        hpgl_design = design
        hpgl_colors = colors_fill
        hpgl_names = names
        hpgl_df = df
    } else if (is.null(df)) {
        hpgl_design = expt$design
        hpgl_colors = expt$colors
        hpgl_df = exprs(expt$expressionset)
        hpgl_names = expt$names
    } else {
        stop("Both df and expt are defined, choose one.")
    }
    if (is.null(hpgl_colors)) {
        hpgl_colors = colorRampPalette(brewer.pal(9,"Blues"))(dim(df)[2])
    }
    hpgl_df = data.frame(hpgl_df)
    if (scale == "log") {
        hpgl_df = hpgl_df + 1
        hpgl_df = log2(hpgl_df)
    }    
    hpgl_df$id = rownames(hpgl_df)
    dataframe = melt(hpgl_df, id=c("id"))
    colnames(dataframe) = c("gene","variable","value")
    boxplot = ggplot2::ggplot(data=dataframe, aes(x=variable, y=value)) +
        suppressWarnings(geom_boxplot(aes(fill=variable),
                     fill=hpgl_colors,
                     size=0.5,
                     outlier.size=1.5,
                     outlier.colour=alpha("black", 0.2))) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90, hjust=1)) +
        xlab("Sample") +
        ylab("Per-gene log(counts)")
    if (!is.null(title)) {
        boxplot = boxplot + ggtitle(title)
    }
    if (!is.null(hpgl_names)) {
        boxplot = boxplot + scale_x_discrete(labels=hpgl_names)
    }
    return(boxplot)
}

#' quantile/quantile comparison of all samples (in this case the mean of all samples, and each sample)
#'
#' @param df the counts
#' @param expt or the expt class
#' @param verbose be chatty while running?
#' @param against either the mean of all samples, or a sample dataset generated on the fly
#' @param norm normalize the data? -- this is no longer needed I think
#'
#' @return a list containing:
#'   logs = a recordPlot() of the pairwise log qq plots
#'   ratios = a recordPlot() of the pairwise ratio qq plots
#'   means = a table of the median values of all the summaries of the qq plots
#'
#' @export
hpgl_qq_all = function(df=NULL, expt=NULL, verbose=FALSE, labels="short") {
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    hpgl_df = NULL
    if (is.null(df)) {
        hpgl_df = exprs(expt$expressionset)
    } else {
        hpgl_df = df
    }
    hpgl_df = as.data.frame(hpgl_df)
    sample_data = hpgl_df[,c(1,2)]
    ## This is bizarre, performing this operation with transform fails when called from a function
    ## but works fine when called interactively, wtf indeed?
    ##    sample_data = transform(sample_data, mean=rowMeans(hpgl_df))
    means = rowMeans(hpgl_df)
    sample_data$mean=means
    logs = list()
    ratios = list()
    means = list()
    comparisons = length(colnames(hpgl_df))
    row_columns = ceiling(sqrt(comparisons))
    rows = nrow(hpgl_df)
    ## I want to make a square containing the graphs.
    count = 1
    for (i in 1:comparisons) {
        ith = colnames(hpgl_df)[i]
        if (verbose) {
            message(paste("Making plot of ", ith, "(", i, ") vs. a sample distribution.", sep=""))
        }
        tmpdf = data.frame(ith=hpgl_df[,i], mean=sample_data$mean)
        colnames(tmpdf) = c(ith, "mean")
        tmpqq = hpgl_qq_plot(df=tmpdf, x=1, y=2, labels=labels)
        logs[[count]] = tmpqq$log
        ratios[[count]] = tmpqq$ratio
        means[[count]] = tmpqq$summary[['Median']]
        count = count + 1
    }
    multiplot(logs)
    log_plots = recordPlot()
    multiplot(ratios)
    ratio_plots = recordPlot()
    plots = list(logs=log_plots, ratios=ratio_plots, medians=means)
    return(plots)
}

#' Perform a qqplot between two columns of a matrix
#'
#' @param a data frame
#' @param or an expt!
#' @param x the first column
#' @param y the second column
#' @param labels whether to include the labels
#'
#' @return a list of the logs, ratios, and mean between the plots as ggplots.
#' @export
hpgl_qq_plot = function(df=NULL, expt=NULL, x=1, y=2, labels=TRUE) {
    hpgl_env = environment()    
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    if (is.null(df)) {
        hpgl_df = exprs(expt$expressionset)
    } else {
        hpgl_df = df
    }
    xlabel = colnames(hpgl_df)[x]
    ylabel = colnames(hpgl_df)[y]
    xvector = as.vector(hpgl_df[,x])
    yvector = as.vector(hpgl_df[,y])
    sorted_x = sort(xvector)
    sorted_y = sort(yvector)
    vector_ratio = sorted_x / sorted_y
    increment = as.vector(1:length(vector_ratio))
    ratio_df = data.frame(cbind(increment, sorted_x, sorted_y, vector_ratio))
    if (labels == "short") {
        y_string = paste(xlabel, " : ", ylabel, sep="")
    } else {
        y_string = paste("Ratio of sorted ", xlabel, " and ", ylabel, ".", sep="")
    }
    ratio_plot = ggplot2::ggplot(ratio_df, aes(x=increment, y=vector_ratio), environment=hpgl_env) +
        geom_point(colour=suppressWarnings(densCols(vector_ratio)), stat="identity", size=1, alpha=0.2, na.rm=TRUE) +
        scale_y_continuous(limits=c(0,2))
    if (isTRUE(labels)) {
        ratio_plot = ratio_plot + xlab("Sorted gene") + ylab(y_string) + theme(legend.position="none")
    } else if (labels == "short") {
        ratio_plot = ratio_plot + ylab(y_string) +
            theme(axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  legend.position="none",
                  panel.background=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank())                  
    } else {
        ratio_plot = ratio_plot + theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    }

    log_df = data.frame(cbind(log(sorted_x), log(sorted_y)))
    gg_max = max(log_df)
    colnames(log_df) = c(xlabel, ylabel)
    log_df$sub = log_df[,1] - log_df[,2]
    log_ratio_plot = ggplot2::ggplot(log_df, aes(x=get(xlabel), y=get(ylabel)), environment=hpgl_env) +
        geom_point(colour=suppressWarnings(densCols(sorted_x, sorted_y)), na.rm=TRUE) +
        scale_y_continuous(limits=c(0, gg_max)) + scale_x_continuous(limits=c(0, gg_max))
    if (isTRUE(labels)) {
        log_ratio_plot = log_ratio_plot + xlab(paste("log sorted ", xlabel)) + ylab(paste("log sorted ", ylabel)) + theme(legend.position="none")
    } else if (labels == "short") {
        log_ratio_plot = log_ratio_plot + xlab("gene") + ylab(y_string) +
            theme(axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  legend.position="none",
                  panel.background=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank())                          
    } else {
        log_ratio_plot = log_ratio_plot +
            theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  legend.position="none",
                  panel.background=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank())
    }
    log_summary = summary(log_df$sub)
    qq_plots = list(ratio=ratio_plot, log=log_ratio_plot, summary=log_summary)
    return(qq_plots)
}

#' Perform qq plots of every column against every other column of a dataset
#'
#' @param df the data
#' @param expt or an expt class
#'
#' @return a list containing the recordPlot() output of the ratios, logs, and means among samples
#'
#' @export
hpgl_qq_all_pairwise = function(df=NULL, expt=NULL, verbose=FALSE) {
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    hpgl_df = NULL
    if (is.null(df)) {
        hpgl_df = exprs(expt$expressionset)
    } else {
        hpgl_df = df
    }
    logs = list()
    ratios = list()
    rows = length(colnames(hpgl_df))
    means = matrix(nrow=rows, ncol=rows)
    count = 1
    for (i in 1:rows) {
        for (j in 1:rows) {
            ith = colnames(hpgl_df)[i]
            jth = colnames(hpgl_df)[j]
            if (verbose) {
                message(paste("Making plot of ", ith, "(", i, ") vs. ", jth, "(", j, ") as element: ", count, ".", sep=""))
            }
            tmp = hpgl_qq_plot(df=hpgl_df, x=i, y=j, labels=labels)
            logs[[count]] = tmp$log
            ratios[[count]] = tmp$ratio
            means[i,j] = tmp$summary[['Mean']]
            count = count + 1
        }
    }
    multiplot(logs)
    log_plots = recordPlot()
    multiplot(ratios)
    ratio_plots = recordPlot()
    heatmap.3(means, trace="none")
    means_heatmap = recordPlot()
    plots = list(logs=log_plots, ratios=ratio_plots, means=means_heatmap)
    return(plots)
}

## I thought multiplot() was a part of ggplot(), but no, weird:
## http://stackoverflow.com/questions/24387376/r-wired-error-could-not-find-function-multiplot
#' Make a grid of plots
#'
#' @param plots a list of plots
#' @param file a file to write to
#' @param cols the number of columns in the grid
#'
#' @return a multiplot!
#' @export
multiplot <- function(plots=NULL, file, cols=NULL, layout=NULL) {
  ## Make a list from the ... arguments and plotlist
  ##  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(cols)) {
      cols = ceiling(sqrt(length(plots)))
  }
  ## If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
      ## Make the panel
      ## ncol: Number of columns of plots
      ## nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
      print(plots[[1]])
  } else {
      ## Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

      ## Make each plot, in the correct location
      for (i in 1:numPlots) {
          ## Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                layout.pos.col = matchidx$col))
      }
  }
}
    
#' Make a heatmap.3 description of the correlation between samples.
#'
#' @param expt an expt set of samples
#' @param df alternately a data frame which must be accompanied by
#' @param design a design matrix and
#' @param colors a color scheme
#' @param method correlation statistic to use.  Defaults to pearson.
#' 
#' @return  corheat_plot a gplots heatmap describing how the samples
#' pairwise correlate with one another.
#' @seealso \code{\link{hpgl_cor}}, \code{\link{brewer.pal}},
#' \code{\link{heatmap.2}}, \code{\link{recordPlot}}
#' 
#' @export
#' @examples
#' ## corheat_plot = hpgl_corheat(expt=expt, method="robust")
#' ## corheat_plot
hpgl_corheat = function(df=NULL, colors=NULL, design=NULL, expt=NULL, method="pearson", names=NULL, row="batch", title=NULL, ...) {
    hpgl_heatmap(df=df, colors=colors, design=design, expt=expt, method=method, names=names, type="correlation", row=row, title=title, ...)
}

#' Make a heatmap.3 description of the similarity (euclildean distance) between samples.
#'
#' @param expt an expt set of samples
#' @param df alternately a data frame which must be accompanied by
#' @param design a design matrix and
#' @param colors a color scheme
#' @param method correlation statistic to use.  Defaults to euclidean.
#' 
#' @return a recordPlot() heatmap describing the distance between samples.
#' @seealso \code{\link{brewer.pal}},
#' \code{\link{heatmap.2}}, \code{\link{recordPlot}}
#' 
#' @export
#' @examples
#' ## disheat_plot = hpgl_disheat(expt=expt, method="euclidean")
#' ## disheat_plot
hpgl_disheat = function(df=NULL, colors=NULL, design=NULL, expt=NULL, method="euclidean", names=NULL, row="batch", title=NULL, ...) {
    hpgl_heatmap(df=df, colors=colors, design=design, expt=expt, method=method, names=names, type="distance", row=row, title=title, ...)
}

hpgl_heatmap = function(df=NULL, colors=NULL, design=NULL, expt=NULL, method="pearson", names=NULL, type="correlation", row="batch", title=NULL, ...) {
    hpgl_env = environment()
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    if (is.null(expt)) {
        hpgl_design = design
        hpgl_colors = colors
        hpgl_df = df
        hpgl_names = names
    } else if (is.null(df)) {
        hpgl_design = expt$design
        hpgl_colors = expt$colors
        hpgl_df = Biobase::exprs(expt$expressionset)
        hpgl_names = expt$names
    } else {
        stop("Both df and expt are defined, choose one.")
    }
    if (is.null(hpgl_colors)) {
        tt = ncol(hpgl_df)
        hpgl_colors = colorRampPalette(brewer.pal(tt,"Dark2"))(tt)
    }
    if (is.null(hpgl_names)) {
        hpgl_names = colnames(hpgl_df)
    }
    
    if (type == "correlation") {
        heatmap_data = hpgltools::hpgl_cor(hpgl_df, method=method)
        heatmap_colors = grDevices::colorRampPalette(brewer.pal(9, "OrRd"))(100)
    } else if (type == "distance") {
        heatmap_data = as.matrix(dist(t(hpgl_df)), method=method)
        heatmap_colors = grDevices::colorRampPalette(brewer.pal(9, "GnBu"))(100)
    }
    hpgl_colors = as.character(hpgl_colors)


    if (is.null(hpgl_design)) {
        row_colors = rep("white", length(hpgl_colors))
    } else if (length(as.integer(as.factor(as.data.frame(hpgl_design[ row ])[,1]))) >= 2) {
##        row_colors = brewer.pal(12, "Set3")[as.integer(as.list(hpgl_design[ row ]))]
        row_colors = RColorBrewer::brewer.pal(12, "Set3")[as.integer(as.factor(as.data.frame(hpgl_design[ row ])[,1]))]
    } else {
        row_colors = rep("green", length(hpgl_design[ row ]))
    }

    ## A temporary hack for the ribosome profiling paper
    ## Correlation heatmaps shouldn't set the colormap for the moment.
    if (type == "correlation") {
        hpgltools::heatmap.3(heatmap_data, keysize=2, labRow=hpgl_names,
                  labCol=hpgl_names, ColSideColors=hpgl_colors, RowSideColors=row_colors,
                  margins=c(8,8), scale="none", trace="none", linewidth=0.5, main=title)
    } else {
        hpgltools::heatmap.3(heatmap_data, keysize=2, labRow=hpgl_names, col=rev(heatmap_colors),
                  labCol=hpgl_names, ColSideColors=hpgl_colors, RowSideColors=row_colors,
                  margins=c(8,8), scale="none", trace="none", linewidth=0.5, main=title)
    }
    hpgl_heatmap_plot = recordPlot()    
    return(hpgl_heatmap_plot)
}

#' Make a heatmap.3 description of the similarity of the genes among samples.
#'
#' @param expt an expt set of samples
#' @param df alternately a data frame which must be accompanied by
#' @param design a design matrix and
#' @param colors a color scheme
#' @param method correlation statistic to use.  Defaults to euclidean.
#' 
#' @return a recordPlot() heatmap describing the samples.
#' @seealso \code{\link{brewer.pal}},
#' \code{\link{heatmap.3}}, \code{\link{recordPlot}}
#' 
#' @export
hpgl_sample_heatmap = function(df=NULL, colors=NULL, design=NULL, expt=NULL, method="pearson", names=NULL, type="correlation", row="batch", title=NULL, ...) {
    hpgl_env = environment()
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    if (is.null(expt)) {
        hpgl_df = df        
        hpgl_design = design
        if (is.null(colors)) {
            hpgl_colors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(dim(hpgl_df)[2])
        } else {
            hpgl_colors = colors
        }
        hpgl_names = names
    } else if (is.null(df)) {
        hpgl_design = expt$design
        hpgl_colors = expt$colors
        hpgl_df = Biobase::exprs(expt$expressionset)
        hpgl_names = expt$names
    } else {
        stop("Both df and expt are defined, choose one.")
    }

    heatmap_colors = redgreen(75)
    if (is.null(hpgl_names)) {
        hpgl_names = colnames(hpgl_df)
    }
    
    heatmap.3(hpgl_df, keysize=2, labRow=NA, col=heatmap_colors,
                   labCol=hpgl_names, margins=c(12,8), trace="none", linewidth=0.5, main=title)

    hpgl_heatmap_plot = recordPlot()    
    return(hpgl_heatmap_plot)
}

#' Make an R plot of the standard median correlation among samples
#'
#' @param expt an expt set of samples
#' @param df alternately a data frame which must be accompanied by
#' @param design a design matrix and
#' @param colors a color scheme
#' @param method a correlation method to use.  Defaults to pearson.
#' @param names use pretty names for the samples?
#' 
#' @return a recordPlot() of the standard median pairwise correlation
#' among the samples.  This will also write to an
#' open device.  The resulting plot measures the median correlation of
#' each sample among its peers.  It notes 1.5* the interquartile range
#' among the samples and makes a horizontal line at that correlation
#' coefficient.  Any sample which falls below this line is considered
#' for removal because it is much less similar to all of its peers.
#' 
#' @seealso \code{\link{hpgl_cor}}, \code{\link{matrixStats::rowMedians}},
#' \code{\link{quantile}}, \code{\link{diff}}, \code{\link{recordPlot}}
#' 
#' @export
#' @examples
#' ## smc_plot = hpgl_smc(expt=expt)
hpgl_smc = function(df=NULL, colors=NULL, expt=NULL, method="pearson", names=NULL, title=NULL, ...) {
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    if (is.null(expt)) {
        hpgl_colors = colors
        hpgl_df = df
    } else if (is.null(df)) {
        hpgl_colors = expt$colors
        hpgl_df = Biobase::exprs(expt$expressionset)
    } else {
        stop("Both df and expt are defined, choose one.")
    }
    if (is.null(expt$names)) {
        if (is.null(names)) {
            hpgl_names = colnames(hpgl_df)
        } else {
            hpgl_names = names
        }
    } else {
        hpgl_names = expt$names
    }
    if (is.null(hpgl_colors)) {
        hpgl_colors = colorRampPalette(brewer.pal(ncol(fun),"Dark2"))(ncol(fun))
    }
    hpgl_colors = as.character(hpgl_colors)
    correlations = hpgl_cor(hpgl_df, method=method)
    cor_median = matrixStats::rowMedians(correlations)
    cor_spread = stats::quantile(cor_median, p=c(1,3)/4)
    cor_iqr = diff(cor_spread)
    outer_limit = cor_spread[1] - 1.5 * cor_iqr
    ylimit = c(pmin(min(cor_median), outer_limit), max(cor_median))
    plot(cor_median, xaxt="n", ylim=ylimit,
         xlab="", main="", ylab="Median pairwise correlation",
         col=hpgl_colors, pch=16, cex=1.5)
    title(title)
    axis(side=1, at=seq(along=cor_median), labels=hpgl_names, las=2)
    abline(h=outer_limit, lty=2)
    abline(v=1:length(hpgl_names), lty=3, col="black")
    hpgl_smc_plot = recordPlot()    
    return(hpgl_smc_plot)
}

#' Make an R plot of the standard median distance among samples
#'
#' @param expt an expt set of samples
#' @param df alternately a data frame which must be accompanied by
#' @param design a design matrix and
#' @param colors a color scheme
#' @param method a distance metric to use.  Defaults to euclidean.
#' @param names use pretty names for the samples?
#' 
#' @return smd_plot a recordPlot of plot.  This will also write to an
#' open device.  This plot takes the median distance of each sample
#' with all of its peers.  It then calculates 1.5* the interquartile
#' range of distances.  Any sample which has a median distance greater
#' than this is considered for removal.
#' 
#' @seealso \code{\link{dist}}, \code{\link{quantile}},
#' \code{\link{diff}}, \code{\link{recordPlot}}
#' 
#' @export
#' @examples
#' ## smd_plot = hpgl_smd(expt=expt)
hpgl_smd = function(expt=NULL, df=NULL, colors=NULL, names=NULL, method="euclidean", title=NULL, ...) {
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    if (is.null(expt)) {
        hpgl_colors = colors
        hpgl_df = df
    } else if (is.null(df)) {
        hpgl_colors = expt$colors
        hpgl_df = Biobase::exprs(expt$expressionset)
    } else {
        stop("Both df and expt are defined, choose one.")
    }
    if (is.null(expt$names)) {
        if (is.null(names)) {
            hpgl_names = colnames(hpgl_df)
        } else {
            hpgl_names = names
        }
    } else {
        hpgl_names = expt$names
    }
    if (is.null(hpgl_colors)) {
        hpgl_colors = colorRampPalette(brewer.pal(ncol(fun),"Dark2"))(ncol(fun))
    }    
    hpgl_colors = as.character(hpgl_colors)
    dists = as.matrix(dist(t(hpgl_df)), method=method)
    dist_median = matrixStats::rowMedians(dists)
    dist_spread = stats::quantile(dist_median, p=c(1,3)/4)
    dist_iqr = diff(dist_spread)
    outer_limit = dist_spread[2] + (1.5 * dist_iqr)
##    ylimit = c(min(dist_median), max(dist_median))
    ylimit = c(min(dist_median), pmax(max(dist_median), outer_limit))
    plot(dist_median, xaxt="n", ylim=ylimit,
         xlab="", main="", ylab="Median pairwise distance",
         col=hpgl_colors, pch=16, cex=1.5)
    title(title)
    axis(side=1, at=seq(along=dist_median), labels=hpgl_names, las=2)
    abline(h=outer_limit, lty=2)
    abline(v=1:length(hpgl_names), lty=3, col="black")
    hpgl_smd_plot = recordPlot()        
    return(hpgl_smd_plot)
}

#' Make a ggplot PCA plot describing the samples' clustering
#'
#' @param expt an expt set of samples
#' @param df alternately a data frame which must be accompanied by
#' @param design a design matrix and
#' @param colors a color scheme
#' @param method a correlation method to use.  Defaults to pearson.
#' @param names use pretty names for the samples?
#' @param labels add labels?  Also, what type?  FALSE, "default", or "fancy"
#' 
#' @return a list containing the following:
#'   plot = ggplot2 pca_plot describing the principle component analysis of the samples.
#'   table = a table of the PCA plot data
#'   res = a table of the PCA res data
#'   variance = a table of the PCA plot variance
#' This makes use of cbcbSEQ and prints the table of variance by component.
#' 
#' @seealso \code{\link{makeSVD}}, \code{\link{pcRes}},
#' \code{\link{geom_dl}}
#' 
#' @export
#' @examples
#' ## pca_plot = hpgl_pca(expt=expt)
#' ## pca_plot
hpgl_pca = function(df=NULL, colors=NULL, design=NULL, expt=NULL, shapes="batch", title=NULL, labels=NULL, ...) {
    hpgl_env = environment()
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors.")
    }
    if (is.null(expt)) {
        hpgl_design = design
        hpgl_df = df
    } else if (is.null(df)) {
        hpgl_design = expt$design
        hpgl_df = Biobase::exprs(expt$expressionset)
    } else {
        stop("Both df and expt are defined, that is confusing.")
    }
    if (is.null(expt$names)) {
        hpgl_labels = colnames(hpgl_df)
    } else {
        hpgl_labels = expt$names
    }
    pca = cbcbSEQ::makeSVD(hpgl_df)  ## This is a part of cbcbSEQ
    if (length(levels(hpgl_design$batch)) == 1) {
        pca_res = cbcbSEQ::pcRes(pca$v, pca$d, hpgl_design$condition)
    } else {
        pca_res = cbcbSEQ::pcRes(pca$v, pca$d, hpgl_design$condition, hpgl_design$batch)
    }
    pca_variance = round((pca$d ^ 2) / sum(pca$d ^ 2) * 100, 2)
    xl = sprintf("PC1: %.2f%% variance", pca_variance[1])
    yl = sprintf("PC2: %.2f%% variance", pca_variance[2])
    pca_data = data.frame(SampleID=hpgl_labels,
        condition=hpgl_design$condition,
        batch=hpgl_design$batch,
        batch_int = as.integer(hpgl_design$batch),
        PC1=pca$v[,1],
        PC2=pca$v[,2])

    num_batches = length(hpgl_design$batch)
    pca_plot = ggplot(data=pca_data, environment=hpgl_env)
    if (num_batches > 6) { ## Then ggplot2 wants shapes specified manually...
        pca_plot = pca_plot + scale_shape_manual(values=1:num_batches)
    }
    pca_plot = pca_plot +
        geom_point(aes(x=PC1, y=PC2, color=hpgl_design$condition, shape=hpgl_design$batch), size=3) +
        scale_colour_discrete(name="Experimental\nCondition") +
        scale_shape_discrete(name="Experimental\nBatch") + 
        xlab(xl) + ylab(yl) + theme_bw()

    if (!is.null(labels)) {
        if (labels == "fancy") {
            pca_plot = pca_plot + directlabels::geom_dl(aes(label=hpgl_labels), method="smart.grid", colour=hpgl_design$condition)
        } else {
            pca_plot = pca_plot + geom_text(aes(x=PC1, y=PC2, label=labels), angle=45, size=4, vjust=2)
        }
    }
    if (!is.null(title)) {
        pca_plot = pca_plot + ggtitle(title)
    }
    pca_return = list(plot=pca_plot, table=pca_data, res=pca_res, variance=pca_variance)
    return(pca_return)
}

factor_rsquared = function(svd_v, factor) {
    svd_lm = lm(svd_v ~ factor)
    lm_summary = summary.lm(svd_lm)
    r_squared = lm_summary$r.squared
    result = round(r_squared * 100, 3)
    return(result)
}

#' Calculate some information useful for generating PCA plots
#'
#' @param df The data to analyze (usually exprs(somedataset))
#' @param design A dataframe describing the experimental design, containing columns with
#'   useful information like the conditions, batches, number of cells, whatever...
#' @param factors A character list of columns from the design matrix which will be queried
#'   for R^2 against the fast.svd of the data.  This defaults to c("condition","batch"), which
#'   is in all of my experimental designs thus far.
#' @param components A number of principle components to compare the design factors against.
#'   If left null, it will query the same number of components as factors asked for.
#' 
#' @return a list of fun pca information:
#'   svd_u/d/v: The u/d/v parameters from fast.svd
#'   rsquared_table: A table of the rsquared values between each factor and principle component
#'   pca_variance: A table of the pca variances
#'   pca_data: Coordinates for a pca plot
#'   pca_cor: A table of the correlations between the factors and principle components
#'   cor_heatmap: A heatmap from recordPlot() describing pca_cor.
#' @seealso \code{\link{fast.svd}}, \code{\link{lm}}
#' 
#' @export
#' @examples
#' ## pca_plot = hpgl_pca(expt=expt)
#' ## pca_plot
pca_information = function(df, design, factors=c("condition","batch"), num_components=NULL) {
    df = exprs(all_major_norm$expressionset)
    data = as.matrix(df)
    means = rowMeans(data)
    decomposed = fast.svd(data - means)
    positives = decomposed$d
    u = decomposed$u
    v = decomposed$v
    rownames(v) = colnames(data)
    component_variance = round((positives^2) / sum(positives^2) * 100, 3)
    cumulative_pc_variance = cumsum(component_variance)
    component_rsquared_table = data.frame(
        prop_var = component_variance,
        cumulative_prop_var = cumulative_pc_variance)

    if (is.null(factors)) {
        factors = colnames(design)
    } else if (factors == "all") {
        factors = colnames(design)
    }
    for (component in factors) {
        comp = as.factor(design[,component])
        column = apply(v, 2, factor_rsquared, factor=comp)
        component_rsquared_table[component] = column
    }
    pca_variance = round((positives ^ 2) / sum(positives ^2) * 100, 2)
    xl = sprintf("PC1: %.2f%% variance", pca_variance[1])
    yl = sprintf("PC2: %.2f%% variance", pca_variance[2])
    labels = rownames(design)

    pca_data = data.frame(SampleID=labels,
        condition=design$condition,
        batch=design$batch,
        batch_int=as.integer(design$batch))
    pc_df = data.frame(SampleID=labels)
    rownames(pc_df) = make.names(labels)
    
    if (is.null(num_components)) {
        num_components = length(factors)
    }
    for (pc in 1:num_components) {
        name = paste("PC", pc, sep="")
        pca_data[name] = v[,pc]
        pc_df[name] = v[,pc]
    }
    pc_df = pc_df[-1]

    factor_df = data.frame(SampleID=labels)
    rownames(factor_df) = make.names(labels)
    for (factor in factors) {
        factor_df[factor] = as.numeric(as.factor(design[,factor]))
    }
    factor_df = factor_df[-1]

    cor_df = data.frame()
    for (factor in factors) {
        for (pc in 1:num_components) {
            cor_df[factor,pc] = cor(factor_df[,factor], pc_df[,pc])
        }
    }
    rownames(cor_df) = colnames(factor_df)
    colnames(cor_df) = colnames(pc_df)
    cor_df = as.matrix(cor_df)
    silly_colors = grDevices::colorRampPalette(brewer.pal(9, "Purples"))(100)
    sillytime = heatmap.3(cor_df, scale="none", trace="none", linewidth=0.5, keysize=2, margins=c(8,8), col=silly_colors)
    pc_factor_corheat = recordPlot()

    pca_list = list(
        svd_d=positives,
        svd_u=u,
        svd_v=v,
        rsquared_table=component_rsquared_table,
        pca_variance=pca_variance,
        pca_data=pca_data,
        pca_cor=cor_df,
        cor_heatmap=pc_factor_corheat)
    return(pca_list)
}


#' Make a pretty MA plot from the output of voom/limma/eBayes/toptable
#'
#' @param counts df of linear-modelling, normalized counts by sample-type,
#' which is to say the output from voom/voomMod/hpgl_voom().
#' @param de_genes df from toptable or its friends containing p-values.
#' @param adjpval_cutoff a cutoff defining significant from not.
#' Defaults to 0.05.
#' @param alpha how transparent to make the dots.  Defaults to 0.6.
#' @param size how big are the dots?  Defaults to 2.
#' @param gvis_filename a filename to write a fancy html graph.
#' Defaults to NULL in which case the following parameter isn't needed.
#' @param tooltip_data a df of tooltip information for gvis
#' graphs. NULL by default.
#' 
#' @return a ggplot2 MA scatter plot.  This is defined as the rowmeans
#' of the normalized counts by type across all sample types on the
#' x-axis, and the log fold change between conditions on the y-axis.
#' Dots are colored depending on if they are 'significant.'  This will
#' make a fun clicky googleVis graph if requested.
#' 
#' @seealso \code{\link{hpgl_gvis_ma_plot}}, \code{\link{toptable}},
#' \code{\link{voom}}, \code{\link{voomMod}}, \code{\link{hpgl_voom}},
#' \code{\link{lmFit}}, \code{\link{makeContrasts}},
#' \code{\link{contrasts.fit}}
#' 
#' @export
#' @examples
#' ## hpgl_ma_plot(voomed_data, toptable_data, gvis_filename="html/fun_ma_plot.html")
#' ## Currently this assumes that a variant of toptable was used which
#' ## gives adjusted p-values.  This is not always the case and I should
#' ## check for that, but I have not yet.
hpgl_ma_plot = function(counts, de_genes, adjpval_cutoff=0.05, alpha=0.6, size=2, tooltip_data=NULL, gvis_filename=NULL, ...) {
    hpgl_env = environment()
    df = data.frame(AvgExp=rowMeans(counts[rownames(de_genes),]),
#        LogFC=de_genes$logFC, AdjPVal=de_genes$adj.P.Val)
        LogFC=de_genes$logFC, AdjPVal=de_genes$P.Value)
    plt = ggplot2::ggplot(df, aes(AvgExp, LogFC, color=(AdjPVal < adjpval_cutoff)), environment=hpgl_env) +
        geom_hline(yintercept=c(-1,1), color="Red", size=size) +
        geom_point(stat="identity", size=size, alpha=alpha) +
        theme(axis.text.x=element_text(angle=-90)) +
        xlab("Average Count (Millions of Reads)") +
        ylab("log fold change") +
        theme_bw()            
    if (!is.null(gvis_filename)) {
        hpgl_gvis_ma_plot(counts, de_genes, tooltip_data=tooltip_data, filename=gvis_filename, ...)
    }
    return(plt)
}

## Consider using these options for the kind of pretty graph Eva likes.
##ggplot(mydata) + aes(x=x, y=y) + scale_x_log10() + scale_y_log10() + 
##+   stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE) +
##+   scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))


#' Make a pretty scatter plot between two sets of numbers with a
#' cheesy distance metric and some statistics of the two sets.
#'
#' @param df a dataframe likely containing two columns
#' @param gvis_filename a filename to write a fancy html graph.
#' Defaults to NULL in which case the following parameter isn't needed.
#' @param tooltip_data a df of tooltip information for gvis
#' graphs. NULL by default.
#' 
#' @return a ggplot2 scatter plot.  This plot provides a "bird's eye"
#' view of two data sets.  This plot assumes the two data structures
#' are not correlated, and so it calculates the median/mad of each
#' axis and uses these to calculate a stupid, home-grown distance
#' metric away from both medians.  This distance metric is used to
#' color dots which are presumed the therefore be interesting because
#' they are far from 'normal.'  This will make a fun clicky googleVis
#' graph if requested.
#'
#' The distance metric should be codified and made more intelligent.
#' Currently it creates a dataframe of distances which are absolute
#' distances from each axis, multiplied by each other, summed by axis,
#' then normalized against the maximum.
#' 
#' @seealso \code{\link{hpgl_gvis_scatter}}, \code{\link{geom_scatter}},
#' \code{\link{hsv}}, \code{\link{hpgl_linear_scatter}}
#' 
#' @export
#' @examples
#' ## hpgl_dist_scatter(lotsofnumbers_intwo_columns, tooltip_data=tooltip_dataframe, gvis_filename="html/fun_scatterplot.html")
hpgl_dist_scatter = function(df, tooltip_data=NULL, gvis_filename=NULL, size=3) {
    hpgl_env = environment()
    df = data.frame(df[,c(1,2)])
    df = df[complete.cases(df),]
    df_columns = colnames(df)
    df_x_axis = df_columns[1]
    df_y_axis = df_columns[2]
    colnames(df) = c("first","second")
    first_median = summary(df[,1])["Median"]
    second_median = summary(df[,2])["Median"]
    first_mad = stats::mad(df[,1])
    second_mad = stats::mad(df[,2])
    mydist = sillydist(df[,1], df[,2], first_median, second_median)
    mydist$x = abs((mydist[,1] - first_median) / abs(first_median))
    mydist$y = abs((mydist[,2] - second_median) / abs(second_median))
    mydist$x = mydist$x / max(mydist$x)
    mydist$y = mydist$y / max(mydist$y)
    mydist$dist = mydist$x * mydist$y
    mydist$dist = mydist$dist / max(mydist$dist)
    line_size = size / 2
    first_vs_second = ggplot2::ggplot(df, aes(x=first, y=second), environment=hpgl_env) +
        xlab(paste("Expression of", df_x_axis)) +
        ylab(paste("Expression of", df_y_axis)) +
        geom_vline(color="grey", xintercept=(first_median - first_mad), size=line_size) +
        geom_vline(color="grey", xintercept=(first_median + first_mad), size=line_size) +
        geom_vline(color="darkgrey", xintercept=first_median, size=line_size) +
        geom_hline(color="grey", yintercept=(second_median - second_mad), size=line_size) +
        geom_hline(color="grey", yintercept=(second_median + second_mad), size=line_size) +
        geom_hline(color="darkgrey", yintercept=second_median, size=line_size) +
        geom_point(colour=hsv(mydist$dist, 1, mydist$dist), alpha=0.6, size=size) +
        theme(legend.position="none")
    if (!is.null(gvis_filename)) {
        hpgl_gvis_scatter(df, tooltip_data=tooltip_data, filename=gvis_filename)
    }
    return(first_vs_second)    
}

#' Make a pretty scatter plot between two sets of numbers
#'
#' @param df a dataframe likely containing two columns
#' @param gvis_filename a filename to write a fancy html graph.
#' Defaults to NULL in which case the following parameter isn't needed.
#' @param tooltip_data a df of tooltip information for gvis
#' graphs. NULL by default.
#' 
#' @return a ggplot2 scatter plot.
#' 
#' @seealso \code{\link{hpgl_gvis_scatter}}, \code{\link{geom_scatter}},
#' \code{\link{hpgl_linear_scatter}}
#' 
#' @export
#' @examples
#' ## hpgl_scatter(lotsofnumbers_intwo_columns, tooltip_data=tooltip_dataframe, gvis_filename="html/fun_scatterplot.html")
hpgl_scatter = function(df, tooltip_data=NULL, color="black", gvis_filename=NULL, size=3) {
    hpgl_env = environment()
    df = data.frame(df[,c(1,2)])
    df = df[complete.cases(df),]
    df_columns = colnames(df)
    df_x_axis = df_columns[1]
    df_y_axis = df_columns[2]
    colnames(df) = c("first","second")
    first_vs_second = ggplot2::ggplot(df, aes(x=first, y=second), environment=hpgl_env) +
        xlab(paste("Expression of", df_x_axis)) +
        ylab(paste("Expression of", df_y_axis)) +
        geom_point(colour=color, alpha=0.6, size=size) +
        theme(legend.position="none")
    if (!is.null(gvis_filename)) {
        hpgl_gvis_scatter(df, tooltip_data=tooltip_data, filename=gvis_filename)
    }
    return(first_vs_second)    
}

#' Make a pretty scatter plot between two sets of numbers with a
#' linear model superimposed and some supporting statistics.
#'
#' @param df a dataframe likely containing two columns
#' @param gvis_filename a filename to write a fancy html graph.
#' Defaults to NULL in which case the following parameter isn't needed.
#' @param tooltip_data a df of tooltip information for gvis
#' graphs. NULL by default.
#' @param cormethod what type of correlation to check?  Defaults to
#' 'pearson'
#' 
#' @return a list including a ggplot2 scatter plot and some
#' histograms.  This plot provides a "bird's eye"
#' view of two data sets.  This plot assumes a (potential) linear
#' correlation between the data, so it calculates the correlation
#' between them.  It then calculates and plots a robust linear model
#' of the data using an 'SMDM' estimator (which I don't remember how
#' to describe, just that the document I was reading said it is good).
#' The median/mad of each axis is calculated and plotted as well.  The
#' distance from the linear model is finally used to color the dots on
#' the plot.  Histograms of each axis are plotted separately and then
#' together under a single cdf to allow tests of distribution
#' similarity.  This will make a fun clicky googleVis graph if
#' requested.
#'
#' @seealso \code{\link{lmrob}}, \code{\link{weights}},
#' \code{\link{hsv}}, \code{\link{mad}}, \code{\link{hpgl_histogram}}
#'
#' @export
#' @examples
#' ## hpgl_linear_scatter(lotsofnumbers_intwo_columns, tooltip_data=tooltip_dataframe, gvis_filename="html/fun_scatterplot.html")
hpgl_linear_scatter = function(df, tooltip_data=NULL, gvis_filename=NULL, cormethod="pearson", size=2, verbose=FALSE, histargs=NULL, loess=FALSE, identity=FALSE, gvis_trendline=NULL, ...) {
    hpgl_env = environment()    
    df = data.frame(df[,c(1,2)])
    df = df[complete.cases(df),]
    correlation = cor.test(df[,1], df[,2], method=cormethod, exact=FALSE)
    df_columns = colnames(df)
    df_x_axis = df_columns[1]
    df_y_axis = df_columns[2]
    colnames(df) = c("first","second")
    linear_model = robustbase::lmrob(formula=second ~ first, data=df, method="SMDM")
    linear_model_summary = summary(linear_model)
    linear_model_weights = stats::weights(linear_model, type="robustness", na.action=NULL)
    linear_model_intercept = stats::coef(linear_model_summary)[1]
    linear_model_slope = stats::coef(linear_model_summary)[2]
    first_median = summary(df$first)["Median"]
    second_median = summary(df$second)["Median"]
    first_mad = stats::mad(df$first, na.rm=TRUE)
    second_mad = stats::mad(df$second, na.rm=TRUE)
    line_size = size / 2
    first_vs_second = ggplot2::ggplot(df, aes(x=first, y=second), environment=hpgl_env) +
        xlab(paste("Expression of", df_x_axis)) +
        ylab(paste("Expression of", df_y_axis)) +
        geom_vline(color="grey", xintercept=(first_median - first_mad), size=line_size) +
        geom_vline(color="grey", xintercept=(first_median + first_mad), size=line_size) +
        geom_hline(color="grey", yintercept=(second_median - second_mad), size=line_size) +
        geom_hline(color="grey", yintercept=(second_median + second_mad), size=line_size) +
        geom_hline(color="darkgrey", yintercept=second_median, size=line_size) +
        geom_vline(color="darkgrey", xintercept=first_median, size=line_size) +
        geom_abline(colour="grey", slope=linear_model_slope, intercept=linear_model_intercept, size=line_size) +
        geom_point(colour=hsv(linear_model_weights * 9/20,
                       linear_model_weights/20 + 19/20,
                       (1.0 - linear_model_weights)),
                       size=size, alpha=0.4)
    if (loess == TRUE) {
        first_vs_second = first_vs_second +
            geom_smooth(method="loess")
    }
    if (identity == TRUE) {
        first_vs_second = first_vs_second +
            geom_abline(colour="darkgreen", slope=1, intercept=0, size=1)
    }
    first_vs_second = first_vs_second +
        theme(legend.position="none")
    
    if (!is.null(gvis_filename)) {
        if (verbose) {
            message("Generating an interactive graph.")
        }
        hpgl_gvis_scatter(df, tooltip_data=tooltip_data, filename=gvis_filename, trendline=gvis_trendline)
    }
    x_histogram = hpgltools::hpgl_histogram(data.frame(df[,1]), verbose=verbose, fillcolor="lightblue", color="blue")
    y_histogram = hpgltools::hpgl_histogram(data.frame(df[,2]), verbose=verbose, fillcolor="pink", color="red")
    both_histogram = hpgltools::hpgl_multihistogram(df, verbose=verbose)
    plots = list(scatter=first_vs_second,
        x_histogram=x_histogram,
        y_histogram=y_histogram,
        both_histogram=both_histogram,
        correlation=correlation,
        lm_model=linear_model,
        lm_summary=linear_model_summary,
        lm_weights=linear_model_weights,
        first_median=first_median,
        first_mad=first_mad,
        second_median=second_median,
        second_mad=second_mad)
    if (verbose) {
        message(sprintf("Calculating correlation between the axes using:", cormethod))
        message(correlation)
        message("Calculating linear model between the axes")
        message(linear_model_summary)
        message("Generating histogram of the x axis.")
        message("Generating histogram of the y axis.")
        message("Generating a histogram comparing the axes.")
    }
    return(plots)
}

#' Make a pretty histogram of something
#'
#' @param df a dataframe of lots of pretty numbers
#' 
#' @return a ggplot histogram
#' @seealso \code{\link{geom_histogram}}, \code{\link{geom_density}},
#' 
#' @export
#' @examples
#' ## kittytime = hpgl_histogram(df)
hpgl_histogram = function(df, binwidth=NULL, log=FALSE, bins=500, verbose=FALSE, fillcolor="darkgrey", color="black") {
    hpgl_env = environment()
    if (class(df) == "data.frame") {
        colnames(df) = c("values")
    } else if (class(df) == "list") {
        df = data.frame(unlist(df))
        colnames(df) = c("values")
    } else if (class(df) == "numeric") {
        df = data.frame(unlist(df))
        colnames(df) = c("values")
    }
    if (is.null(binwidth)) {
        minval = min(df, na.rm=TRUE)
        maxval = max(df, na.rm=TRUE)
        binwidth = (maxval - minval) / bins
        if (verbose) {
            message(paste("No binwidth provided, setting it to ", binwidth, " in order to have ", bins, " bins.", sep=""))
        }
    }
    a_histogram = ggplot2::ggplot(df, aes(x=values), environment=hpgl_env) +
    geom_histogram(aes(y=..density..), stat="bin", binwidth=binwidth, colour=color, fill=fillcolor, position="identity") +
##    stat_bin(binwidth=binwidth, color=fillcolor) +
    geom_density(alpha=0.4, fill=fillcolor) +
    geom_vline(aes(xintercept=mean(values, na.rm=T)), color=color, linetype="dashed", size=1) +
    theme_bw()
    if (log) {
        log_histogram = try(a_histogram + scale_x_log10())
        if (log_histogram != 'try-error') {
            a_histogram = log_histogram
        }
    }
    return(a_histogram)
}

#' Make a pretty histogram of multiple datasets
#'
#' @param df a dataframe of lots of pretty numbers
#' 
#' @return a ggplot histogram comparing multiple data sets
#' Along the way this generates pairwise t tests of the columns of
#' data.
#' 
#' @seealso \code{\link{pairwise.t.test}}, \code{\link{ddply}},
#' \code{\link{rbind}}
#' 
#' @export
#' @examples
#' ## kittytime = hpgl_multihistogram(df)
hpgl_multihistogram = function(data, log=FALSE, binwidth=NULL, bins=NULL, verbose=FALSE) {
    if (is.data.frame(data)) {
        df = data
        columns = colnames(df)
        summary_df = summary(df)
        play_all = data.frame()
        for (col in 1:length(colnames(df))) {
            new_column = data.frame(expression=df[,col], cond=colnames(df)[col])
            play_all = BiocGenerics::rbind(play_all, new_column)
        }
    } else if (is.list(data)) {
        summary_df = summary(data)
        play_all = reshape::melt(data)
        colnames(play_all) = c("expression","cond")
    } else {
        stop("This can only work with a list or data frame.")
    }
    play_cdf = plyr::ddply(play_all, "cond", summarise, rating.mean=mean(expression, na.rm=TRUE))
    uncor_t = stats::pairwise.t.test(play_all$expression, play_all$cond, p.adjust="none")
    bon_t = try(stats::pairwise.t.test(play_all$expression, play_all$cond, p.adjust="bon", na.rm=TRUE))
    if (is.null(bins) & is.null(binwidth)) {
        minval = min(play_all$expression, na.rm=TRUE)
        maxval = max(play_all$expression, na.rm=TRUE)
        bins = 500
        binwidth = (maxval - minval) / bins
        message(paste("No binwidth nor bins provided, setting it to ", binwidth, " in order to have ", bins, " bins.", sep=""))        
    } else if  (is.null(binwidth)) {
        minval = min(play_all$expression, na.rm=TRUE)
        maxval = max(play_all$expression, na.rm=TRUE)
        binwidth = (maxval - minval) / bins
        message(paste("Setting binwidth to ", binwidth, " in order to have ", bins, " bins.", sep=""))        
    } else if (is.null(bins)) {
        message(paste("Setting binwidth to ", binwidth, sep=""))
    } else {
        message("Both bins and binwidth were provided, using binwidth: ", binwidth, sep="")
    }
    hpgl_multi = ggplot2::ggplot(play_all, aes(x=expression, fill=cond)) +
        geom_histogram(aes(y=..density..), binwidth=binwidth, alpha=0.4, position="identity") +
        xlab("Expression") +
        ylab("Number of observations") +
        geom_density(alpha=0.5) +                
        geom_vline(data=play_cdf, aes(xintercept=rating.mean,  colour=cond), linetype="dashed", size=0.75) +
        theme_bw()
    if (log) {
        logged = try(hpgl_multi + scale_x_log10())
        if (logged != 'try-error') {
            hpgl_multi = logged
        }
    }
    if (verbose) {
        message("Summarise the data.")
        message(summary_df)
        message("Uncorrected t test(s) between columns:")
        message(uncor_t)
        if (class(bon_t) == 'try-error') {
            message("Unable to perform corrected test.")
        } else {
            message("Bon Ferroni corrected t test(s) between columns:")
            message(bon_t)
        }
    }
    returns = list(plot=hpgl_multi,
        data_summary=summary_df,
        uncor_t=uncor_t,
        bon_t=bon_t)
    return(returns)
}

#' Density plots!
#'
#' @param df a data frame to look at
#' @param expt or an expt class!
#' @param names names of the samples
#' @param position stack or identity
#' @param fill NULL, but otherwise it fills in the colors for the rainbow coalition!
#'
#' @return a density plot!
#' @export
hpgl_density_plot = function(df=NULL, colors=NULL, expt=NULL, names=NULL, position="identity", fill=NULL, title=NULL) {
    hpgl_env = environment()
    message("This plot looks neat if you do position='fill' or position='stack'")
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df.")
    }
    if (is.null(expt)) {
        hpgl_colors = colors
        hpgl_df = df
        hpgl_names = names
    } else if (is.null(df)) {
        hpgl_colors = expt$colors
        hpgl_df = Biobase::exprs(expt$expressionset)
        hpgl_names = expt$names
    } else {
        stop("Both expt and df are defined, please choose one.")
    }

    if (!is.null(hpgl_names)) {
        colnames(hpgl_df) = make.names(hpgl_names, unique=TRUE)
    }
    ## If the columns lose the connectivity between the sample and values, then
    ## the ggplot below will fail with env missing.
    melted = reshape::melt(as.matrix(hpgl_df))
    if (dim(melted)[2] == 3) {
        colnames(melted) = c("id", "sample", "counts")        
    } else if (dim(melted)[2] == 2) {
        colnames(melted) = c("sample","counts")
    } else {
        stop("Could not properly melt the data.")
    }
    colors = factor(hpgl_colors)    
    if (is.null(hpgl_colors)) {
        hpgl_colors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(dim(hpgl_df)[2])
    }
    if (!is.null(fill)) {
        fill = "sample"
    }
    densityplot = ggplot2::ggplot(data=melted, aes(x=counts, colour=sample, fill=fill), environment=hpgl_env) +
        geom_density(aes(x=counts, y=..count..), position=position) +
        theme_bw()
    if (!is.null(title)) {
        densityplot = densityplot + ggplot2::ggtitle(title)
    }
    return(densityplot)
}

#' Make a pretty Volcano plot!
#'
#' @param toptable_data a dataframe from limma's toptable which
#' includes log(fold change) and an adjusted p-value. 
#' @param p_cutoff a cutoff defining significant from not.
#' Defaults to 0.05.
#' @param fc_cutoff a cutoff defining the minimum/maximum fold change
#' for interesting.  This is log, so I went with +/- 0.8 mostly
#' arbitrarily as the default.
#' @param alpha how transparent to make the dots.  Defaults to 0.6.
#' @param size how big are the dots?  Defaults to 2.
#' @param gvis_filename a filename to write a fancy html graph.
#' Defaults to NULL in which case the following parameter isn't needed.
#' @param tooltip_data a df of tooltip information for gvis
#' graphs. NULL by default.
#' 
#' @return a ggplot2 MA scatter plot.  This is defined as the
#' -log10(p-value) with respect to log(fold change).  The cutoff
#' values are delineated with lines and mark the boundaries between
#' 'significant' and not.  This will make a fun clicky googleVis graph
#' if requested.
#' 
#' @seealso \code{\link{hpgl_gvis_ma_plot}}, \code{\link{toptable}},
#' \code{\link{voom}}, \code{\link{voomMod}}, \code{\link{hpgl_voom}},
#' \code{\link{lmFit}}, \code{\link{makeContrasts}},
#' \code{\link{contrasts.fit}}
#' 
#' @export
#' @examples
#' ## hpgl_volcano_plot(toptable_data, gvis_filename="html/fun_ma_plot.html")
#' ## Currently this assumes that a variant of toptable was used which
#' ## gives adjusted p-values.  This is not always the case and I should
#' ## check for that, but I have not yet.
hpgl_volcano_plot = function(toptable_data, tooltip_data=NULL, gvis_filename=NULL, fc_cutoff=0.8, p_cutoff=0.05, size=2, alpha=0.6, ...) {
    hpgl_env = environment()
    low_vert_line = 0.0 - fc_cutoff
    horiz_line = -1 * log10(p_cutoff)
    toptable_data$modified_p = -1 * log10(toptable_data$P.Value)
    plt = ggplot2::ggplot(toptable_data, aes(x=logFC, y=modified_p, color=(P.Value <= p_cutoff)), environment=hpgl_env) +
        geom_hline(yintercept=horiz_line, color="black", size=size) +
        geom_vline(xintercept=fc_cutoff, color="black", size=size) +
        geom_vline(xintercept=low_vert_line, color="black", size=size) +
        geom_point(stat="identity", size=size, alpha=alpha) +
##        theme(axis.text.x=element_text(angle=-90)) +
        xlab("log fold change") +
        ylab("-log10(adjusted p value)") +
        theme(legend.position="none")
    if (!is.null(gvis_filename)) {
        hpgl_gvis_volcano_plot(toptable_data, fc_cutoff=fc_cutoff, p_cutoff=p_cutoff, tooltip_data=tooltip_data, filename=gvis_filename)
    }
    return(plt)
}

#' slight modification of heatmap.2
#'
#' @param holy crap so many parameters, same as heatmap.2
#'
#' @return a heatmap!
#' @seealso \code{\link{heatmap.2}}
#'
#' @export
heatmap.3 = function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
        w), symm = FALSE, scale = c("none", "row", "column"), 
    na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
    symbreaks = min(x < 0, na.rm = TRUE) || scale != "none", 
    col = "heat.colors", colsep, rowsep, sepcolor = "white", 
    sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
    na.color = par("bg"), trace = c("column", "row", "both", 
        "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
    linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors, 
    cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
        NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, key.title = NULL, 
    key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL, 
    key.par = list(), main = NULL, xlab = NULL, ylab = NULL, 
    lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, linewidth = 1.0, ...) 
{
    if (!is.null(main)) {
        if (main == FALSE) {
            main = NULL
        }
    }
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
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
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
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
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    if (is.null(srtCol)) 
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
            offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
            padj = adjCol[2])
    else {
        if (is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol)) 
                adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                strheight("M"), labels = labCol, adj = adjCol, 
                cex = cexCol, srt = srtCol)
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
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = linewidth, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = linewidth,
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

#' Steal EdgeR's plotBCV()
#'
#' @param expt
#'
#' @return a plot! of the BCV a la ggplot2
#' @export
hpgl_plot_bcv = function(expt) {
    data = exprs(expt$expressionset)^2
    data = DGEList(counts=data)
    edisp = estimateDisp(data)
    avg_log_cpm = edisp$AveLogCPM
    if (is.null(avg_log_cpm)) {
        avg_log_cpm = aveLogCPM(edisp$counts, offset=getOffset(edisp))
    }
    disper = getDispersion(edisp)
    if (is.null(disper)) {
        stop("No dispersions to plot")
    }
    if (attr(disper, "type") == "common") {
        disper = rep(disper, length = length(avg_log_cpm))
    }
    disp_df = data.frame(A=avg_log_cpm, disp=sqrt(disper))
    fitted_disp = gplots::lowess(disp_df$A, disp_df$disp, f=0.5)
    f = stats::approxfun(fitted_disp, rule=2)
    disp_plot = ggplot2::ggplot(disp_df, aes(x=A, y=disp)) +
        geom_point() +
        xlab("Average log(CPM)") +
        ylab("Dispersion of Biological Variance") +
        stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE, show_guide=FALSE) +
        scale_fill_gradientn(colours=colorRampPalette(c("white","black"))(256)) +
        geom_smooth(method="loess") +
        stat_function(fun=f, colour="red") +
        theme(legend.position="none")
    return(disp_plot)
}
