## Time-stamp: <Fri May 13 14:59:10 2016 Ashton Trey Belew (abelew@gmail.com)>

## plot_bar.r: Some useful bar plots, currently only library size

#' Make a ggplot graph of library sizes.
#'
#' It is often useful to have a quick view of which samples have more/fewer reads.  This does that
#' and maintains one's favorite color scheme and tries to make it pretty!
#'
#' @param data Expt, dataframe, or expressionset of samples.
#' @param colors Color scheme if the data is not an expt.
#' @param names Alternate names for the x-axis.
#' @param text Add the numeric values inside the top of the bars of the plot?
#' @param title Title for the plot.
#' @param yscale Whether or not to log10 the y-axis.
#' @param ... More parameters for your good time!
#' @return a ggplot2 bar plot of every sample's size
#' @seealso \link[ggplot2]{geom_bar} \link[ggplot2]{geom_text}
#' \link{prettyNum} \link[ggplot2]{scale_y_log10}
#' @examples
#' \dontrun{
#'  libsize_plot = plot_libsize(expt=expt)
#'  libsize_plot  ## ooo pretty bargraph
#' }
#' @export
plot_libsize <- function(data, colors=NULL, names=NULL, text=TRUE, title=NULL,  yscale=NULL, ...) {
    hpgl_env <- environment()
    arglist <- list(...)
    if (is.null(text)) {
        text <- TRUE
    }

    ## In response to Keith's recent comment when there are more than 8 factors
    chosen_palette <- "Dark2"
    if (!is.null(arglist$palette)) {
        chosen_palette <- arglist$palette
    }

    data_class <- class(data)[1]
    if (data_class == "expt") {
        design <- data[["design"]]
        colors <- data[["colors"]]
        names <- data[["names"]]
        data <- Biobase::exprs(data[["expressionset"]])  ## Why does this need the simplifying
        ## method of extracting an element? (eg. data['expressionset'] does not work)
        ## that is _really_ weird!
    } else if (data_class == "ExpressionSet") {
        data <- Biobase::exprs(data)
    } else if (data_class == "matrix" | data_class == "data.frame") {
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }

    if (is.null(colors)) {
        colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(ncol(data), chosen_palette))(ncol(data))
    }
    tmp <- data.frame(id=colnames(data),
                      sum=colSums(data),
                      colors=colors)
    tmp$order <- factor(tmp$id, as.character(tmp$id))
    libsize_plot <- ggplot2::ggplot(data=tmp, ggplot2::aes_string(x="order", y="sum"),
                                    environment=hpgl_env, colour=colors) +
        ggplot2::geom_bar(ggplot2::aes_string(x="order"), stat="identity", colour="black", fill=tmp$colors) +
        ggplot2::xlab("Sample ID") +
        ggplot2::ylab("Library size in (pseudo)counts.") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1.5, vjust=0.5))
    if (isTRUE(text)) {
        ## ggplot2::geom_text(ggplot2::aes(data=tmp, label=prettyNum(sum, big.mark=",")),
        ##                                 angle=90, size=3, color="white", hjust=1.2)
        tmp$sum <- sprintf("%.2f", round(tmp$sum, 2))
        libsize_plot <- libsize_plot + ggplot2::geom_text(
                                           ggplot2::aes_string(x="order", label='prettyNum(as.character(tmp$sum), big.mark=",")'),
                                           angle=90, size=4, color="white", hjust=1.2)
    }
    if (!is.null(title)) {
        libsize_plot <- libsize_plot + ggplot2::ggtitle(title)
    }
    if (is.null(yscale)) {
        scale_difference <- max(as.numeric(tmp$sum)) / min(as.numeric(tmp$sum))
        if (scale_difference > 10.0) {
            message("The scale difference between the smallest and largest
   libraries is > 10. Assuming a log10 scale is better, set scale=FALSE if not.")
            scale = TRUE
        } else {
            scale = FALSE
        }
    }
    if (scale == TRUE) {
        libsize_plot <- libsize_plot +
            ggplot2::scale_y_log10()
    }
    if (!is.null(names)) {
        libsize_plot <- libsize_plot +
            ggplot2::scale_x_discrete(labels=names)
    }
    return(libsize_plot)
}

## EOF  Damners I don't have many bar plots, do I?
