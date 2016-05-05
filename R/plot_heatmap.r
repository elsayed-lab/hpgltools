## Time-stamp: <Wed May  4 23:34:01 2016 Ashton Trey Belew (abelew@gmail.com)>

## plot_heatmap.r: Heatmaps separated by usage

#' Make a heatmap.3 description of the correlation between samples.
#'
#' Given a set of count tables and design, this will calculate the pairwise correlations and plot
#' them as a heatmap.  It attempts to standardize the inputs and eventual output.
#'
#' @param data Dataframe, expt, or expressionset to work with.
#' @param design Design matrix describing the experiment, not needed if this is an expt.
#' @param colors Color scheme for the samples, not needed if this is an expt.
#' @param method Correlation statistic to use. (pearson, spearman, kendall, robust).
#' @param names Alternate names to use for the samples.
#' @param title Title for the plot.
#' @param ... More options are wonderful!
#' @return Gplots heatmap describing describing how the samples are clustering vis a vis pairwise correlation.
#' @seealso \link{hpgl_cor} \link[RColorBrewer]{brewer.pal}
#' \link[grDevices]{recordPlot}
#' @examples
#' ## corheat_plot = hpgl_corheat(expt=expt, method="robust")
#' ## corheat_plot
#' @export
hpgl_corheat <- function(data, colors=NULL, design=NULL, method="pearson",
                         names=NULL, title=NULL, ...) {
    hpgl_heatmap(data, colors=colors, design=design, method=method,
                 names=names, type="correlation", title=title, ...)
}

#' Make a heatmap.3 description of the distances (euclidean by default) between samples.
#'
#' Given a set of count tables and design, this will calculate the pairwise distances and plot
#' them as a heatmap.  It attempts to standardize the inputs and eventual output.
#'
#' @param data Dataframe, expt, or expressionset to work with.
#' @param colors Color scheme (not needed if an expt is provided).
#' @param design Design matrix (not needed if an expt is provided).
#' @param method Distance metric to use.
#' @param names Alternate names to use for the samples.
#' @param title Title for the plot.
#' @param ... More parameters!
#' @return a recordPlot() heatmap describing the distance between samples.
#' @seealso \link[RColorBrewer]{brewer.pal} \link[gplots]{heatmap.2} \link[grDevices]{recordPlot}
#' @examples
#' \dontrun{
#'  disheat_plot = hpgl_disheat(expt=expt, method="euclidean")
#'  disheat_plot
#' }
#' @export
hpgl_disheat <- function(data, colors=NULL, design=NULL, method="euclidean",
                         names=NULL, title=NULL, ...) {
    hpgl_heatmap(data, colors=colors, design=design, method=method,
                 names=names, type="distance", row=row, title=title, ...)
}

#' Make a heatmap.3 plot, does the work for hpgl_disheat and hpgl_corheat.
#'
#' This does what is says on the tin.  Sets the colors for correlation or distance heatmaps, handles
#' the calculation of the relevant metrics, and plots the heatmap.
#'
#' @param data Dataframe, expt, or expressionset to work with.
#' @param colors Color scheme for the samples.
#' @param design Design matrix describing the experiment vis a vis conditions and batches.
#' @param method Distance or correlation metric to use.
#' @param names Alternate names to use for the samples.
#' @param type Defines the use of correlation, distance, or sample heatmap.
#' @param title Title for the plot.
#' @param ... I like elipses!
#' @return a recordPlot() heatmap describing the distance between samples.
#' @seealso \link[RColorBrewer]{brewer.pal} \link[grDevices]{recordPlot}
#' @export
hpgl_heatmap <- function(data, colors=NULL, design=NULL, method="pearson", names=NULL,
                         type="correlation", row="batch", title=NULL, ...) {
    arglist <- list(...)
    hpgl_env <- environment()
    data_class <- class(data)[1]
    if (data_class == "expt") {
        design <- data[["design"]]
        colors <- data[["colors"]]
        names <- data[["names"]]
        data <- Biobase::exprs(data[["expressionset"]])
    } else if (data_class == "ExpressionSet") {
        data <- Biobase::exprs(data)
    } else if (data_class == "matrix" | data_class == "data.frame") {
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    chosen_palette <- "Dark2"
    if (!is.null(arglist$palette)) {
        chosen_palette <- arglist$palette
    }

    if (is.null(colors)) {
        tt <- ncol(data)
        colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(tt, chosen_palette))(tt)
    }
    if (is.null(names)) {
        names <- colnames(data)
    }
    if (type == "correlation") {
        heatmap_data <- hpgl_cor(data, method=method)
        heatmap_colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(100)
    } else if (type == "distance") {
        heatmap_data <- as.matrix(dist(t(data)), method=method)
        heatmap_colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "GnBu"))(100)
    }
    colors <- as.character(colors)
    if (is.null(design)) {
        row_colors <- rep("white", length(colors))
    } else if (length(as.integer(as.factor(as.data.frame(design[row])[, 1]))) >= 2) {
        ## row_colors = brewer.pal(12, "Set3")[as.integer(as.list(hpgl_design[ row ]))]
        row_colors <- RColorBrewer::brewer.pal(12, "Set3")[as.integer(as.factor(as.data.frame(design[ row ])[, 1]))]
    } else {
        row_colors <- rep("green", length(design[row]))
    }
    if (type == "correlation") {
        heatmap.3(heatmap_data, keysize=2, labRow=names,
                  ##col=heatmap_colors,  ## OrRd is slightly different than what we have now
                  labCol=names, ColSideColors=colors, RowSideColors=row_colors,
                  margins=c(8,8), scale="none", trace="none", linewidth=0.5, main=title)
    } else {
        heatmap.3(heatmap_data, keysize=2, labRow=names, col=rev(heatmap_colors),
                  labCol=names, ColSideColors=colors, RowSideColors=row_colors,
                  margins=c(8,8), scale="none", trace="none", linewidth=0.5, main=title)
    }
    hpgl_heatmap_plot <- grDevices::recordPlot()
    return(hpgl_heatmap_plot)
}

hpgl_heatplus <- function(fundata) {
    heatmap_data <- hpgl_cor(fundata)
    heatmap_colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(100)
    funkytown <- Heatplus::annHeatmap2(heatmap_data)
    plot(funkytown)
    ret <- grDevices::recordPlot()
}

## Taken from https://plot.ly/ggplot2/ggdendro-dendrograms/
## Check out the following link for a neat dendrogram library.
## http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

ggplot2_heatmap <- function() {
    x <- as.matrix(scale(mtcars))
    dd.col <- as.dendrogram(hclust(dist(x)))
    dd.row <- as.dendrogram(hclust(dist(t(x))))
    dx <- dendro_data(dd.row)
    dy <- dendro_data(dd.col)
    ## helper function for creating dendograms
    ggdend <- function(df) {
        ggplot() +
            geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
            labs(x = "", y = "") + theme_minimal() +
            theme(axis.text = element_blank(), axis.ticks = element_blank(),
                  panel.grid = element_blank())
    }
    ## x/y dendograms
    px <- ggdend(dx$segments)
    py <- ggdend(dy$segments) + coord_flip()
    ## heatmap
    col.ord <- order.dendrogram(dd.col)
    row.ord <- order.dendrogram(dd.row)
    xx <- scale(mtcars)[col.ord, row.ord]
    xx_names <- attr(xx, "dimnames")
    df <- as.data.frame(xx)
    colnames(df) <- xx_names[[2]]
    df$car <- xx_names[[1]]
    df$car <- with(df, factor(car, levels=car, ordered=TRUE))
    mdf <- reshape2::melt(df, id.vars="car")
    p <- ggplot(mdf, aes(x = variable, y = car)) + geom_tile(aes(fill = value))
    ## hide axis ticks and grid lines
    eaxis <- list(
        showticklabels = FALSE,
        showgrid = FALSE,
        zeroline = FALSE
    )
}


#' Make a heatmap.3 description of the similarity of the genes among samples.
#'
#' Sometimes you just want to see how the genes of an experiment are related to each other.  This
#' can handle that.  These heatmap functions should probably be replaced with neatmaps or heatplus
#' or whatever it is, as the annotation dataframes in them are pretty awesome.
#'
#' @param data Expt/expressionset/dataframe set of samples.
#' @param colors Color scheme of the samples (not needed if input is an expt).
#' @param design Design matrix describing the experiment (gotten for free if an expt).
#' @param names Alternate samples names.
#' @param title Title of the plot!
#' @param Rowv Include the row names?
#' @param ... More parameters for a good time!
#' @return a recordPlot() heatmap describing the samples.
#' @seealso \link[RColorBrewer]{brewer.pal} \link[grDevices]{recordPlot}
#' @export
hpgl_sample_heatmap <- function(data, colors=NULL, design=NULL, names=NULL, title=NULL, Rowv=FALSE, ...) {
    hpgl_env <- environment()
    data_class <- class(data)[1]
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
    heatmap_colors <- gplots::redgreen(75)
    if (is.null(names)) {
        names <- colnames(data)
    }
    data <- as.matrix(data)
    heatmap.3(data, keysize=2, labRow=NA, col=heatmap_colors,
              labCol=names, margins=c(12,8), trace="none", linewidth=0.5, main=title, Rowv=Rowv)
    hpgl_heatmap_plot <- grDevices::recordPlot()
    return(hpgl_heatmap_plot)
}

#' a minor change to heatmap.2 makes heatmap.3
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
#' @export
heatmap.3 <- function (x, Rowv=TRUE, Colv=if (symm) "Rowv" else TRUE,
                       distfun=dist, hclustfun=hclust, dendrogram=c("both","row","column","none"),
                       reorderfun = function(d, w) reorder(d,w), symm = FALSE, scale = c("none", "row", "column"),
                       na.rm=TRUE, revC=identical(Colv,"Rowv"), add.expr, breaks,
                       symbreaks = min(x < 0, na.rm = TRUE) || scale != "none",
                       col = "heat.colors", colsep, rowsep, sepcolor = "white",
                       sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan",
                       na.color = par("bg"), trace = c("column", "row", "both", "none"),
                       tracecol = "cyan", hline = median(breaks), vline = median(breaks),
                       linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors,
                       cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
                       labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0,NA),
                       adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5,
                       key = TRUE, keysize = 1.5, density.info = c("histogram", "density", "none"),
                       denscol = tracecol, symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25,
                       key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL,
                       key.par = list(), main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL,
                       lwid = NULL, extrafun = NULL, linewidth = 1.0, ...) {
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
    if (!gtools::invalid(na.color) & any(is.na(x))) {
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

