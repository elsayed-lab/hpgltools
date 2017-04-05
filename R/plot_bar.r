## plot_bar.r: Some useful bar plots, currently only library size

#' Make a ggplot graph of library sizes.
#'
#' It is often useful to have a quick view of which samples have more/fewer reads.  This does that
#' and maintains one's favorite color scheme and tries to make it pretty!
#'
#' @param data Expt, dataframe, or expressionset of samples.
#' @param condition vector of sample condition names.
#' @param colors Color scheme if the data is not an expt.
#' @param names Alternate names for the x-axis.
#' @param text Add the numeric values inside the top of the bars of the plot?
#' @param title Title for the plot.
#' @param yscale Whether or not to log10 the y-axis.
#' @param ... More parameters for your good time!
#' @return a ggplot2 bar plot of every sample's size
#' @seealso \pkg{ggplot2}
#'  \code{\link[ggplot2]{geom_bar}} \code{\link[ggplot2]{geom_text}}
#'  \code{\link{prettyNum}} \code{\link[ggplot2]{scale_y_log10}}
#' @examples
#' \dontrun{
#'  libsize_plot <- plot_libsize(expt=expt)
#'  libsize_plot  ## ooo pretty bargraph
#' }
#' @export
plot_libsize <- function(data, condition=NULL, colors=NULL,
                         names=NULL, text=TRUE,
                         title=NULL,  yscale=NULL, ...) {
    hpgl_env <- environment()
    arglist <- list(...)
    if (is.null(text)) {
        text <- TRUE
    }

    ## In response to Keith's recent comment when there are more than 8 factors
    chosen_palette <- "Dark2"
    if (!is.null(arglist[["palette"]])) {
        chosen_palette <- arglist[["palette"]]
    }

    data_class <- class(data)[1]
    if (data_class == "expt") {
        condition <- data[["design"]][["condition"]]
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

    # Get conditions
    if (is.null(condition)) {
        stop("Missing condition label vector.")
    }

    colors <- as.character(colors)
    libsize_df <- data.frame("id" = colnames(data),
                             "sum" = colSums(data),
                             "condition" = condition,
                             "colors" = as.character(colors))
    libsize_df[["order"]] <- factor(libsize_df[["id"]], as.character(libsize_df[["id"]]))

    color_listing <- libsize_df[, c("condition", "colors")]
    color_listing <- unique(color_listing)
    color_list <- as.character(color_listing[["colors"]])
    names(color_list) <- as.character(color_listing[["condition"]])

    libsize_plot <- ggplot2::ggplot(data=libsize_df, environment=hpgl_env, colour=colors,
                                    aes_string(x="order",
                                               y="sum")) +
        ggplot2::geom_bar(stat="identity",
                          colour="black",
                          fill=libsize_df[["colors"]],
                          aes_string(x="order")) +
        ggplot2::xlab("Sample ID") +
        ggplot2::ylab("Library size in (pseudo)counts.") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1.5, vjust=0.5))

    if (isTRUE(text)) {
        libsize_df[["sum"]] <- sprintf("%.2f", round(as.numeric(libsize_df[["sum"]]), 2))
        ## newlabels <- prettyNum(as.character(libsize_df[["sum"]]), big.mark=",")
        libsize_plot <- libsize_plot +
            ggplot2::geom_text(parse=FALSE, angle=90, size=4, color="white", hjust=1.2,
                               ## ggplot2::aes_string(parse=FALSE,
                               ggplot2::aes_string(x="order",
                                                   label='prettyNum(as.character(libsize_df$sum), big.mark=",")'))
    }

    if (!is.null(title)) {
        libsize_plot <- libsize_plot + ggplot2::ggtitle(title)
    }
    if (is.null(yscale)) {
        scale_difference <- max(as.numeric(libsize_df[["sum"]])) / min(as.numeric(libsize_df[["sum"]]))
        if (scale_difference > 10.0) {
            message("The scale difference between the smallest and largest
   libraries is > 10. Assuming a log10 scale is better, set scale=FALSE if not.")
            scale <- TRUE
        } else {
            scale <- FALSE
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

#' Make relatively pretty bar plots of coverage in a genome.
#'
#' This was written for ribosome profiling coverage / gene.
#' It should however, work for any data with little or no modification.
#'
#' @param input  Coverage / position filename.
#' @param workdir  Where to put the resulting images.
#' @param output  Output image filename.
#' @param name  Gene name to print at the bottom of the plot.
#' @param start  Relative to 0, where is the gene's start codon.
#' @param end  Relative to 0, where is the gene's stop codon.
#' @param strand  Is this on the + or - strand? (+1/-1)
#' @param padding  How much space to provide on the sides?
#' @return coverage plot surrounging the ORF of interest
#' @seealso \pkg{ggplot2}
plot_rpm <- function(input, workdir="images", output="01.svg", name="LmjF.01.0010",
                    start=1000, end=2000, strand=1, padding=100) {

    mychr <- gsub("\\.\\d+$", "", name, perl=TRUE)
    plotted_start <- start - padding
    plotted_end <- end + padding
    my_start <- start
    my_end <- end
    ## These are good chances to use %>% I guess
    ## rpm_region = subset(input, chromosome==mychr & position >= plotted_start & position <= plotted_end)
    region_idx <- input[["chromosome"]] == mychr &
        input[["position"]] >= plotted_start &
        input[["position"]] <= plotted_end
    rpm_region <- rpm_region[region_idx, ]
    rpm_region <- rpm_region[, -1]
    rpm_region[["log"]] <- log2(rpm_region[["rpm"]] + 1)

    ## pre_start = subset(rpm_region, position < my_start)
    start_idx <- rpm_region[["position"]] < my_start
    pre_start <- rpm_region[start_idx, ]
    ## post_stop = subset(rpm_region, position > my_end)
    stop_idx <- rpm_region[["position"]] > my_end
    post_stop <- rpm_region[stop_idx, ]
    ## cds = subset(rpm_region, position >= my_start & position < my_end)
    cds_idx <- rpm_region[["position"]] >= my_start & rpm_region[["position"]] < my_end
    cds <- rpm_region[cds_idx, ]

    eval(substitute(
        expr = {
            stupid <- aes(y=0, yend=0, x=my_start, xend=my_end)
        },
        env <- list(my_start=my_start, my_end=my_end)))

    if (strand == "+") {
        gene_arrow <- grid::arrow(type="closed", ends="last")
    } else {
        gene_arrow <- grid::arrow(type="closed", ends="first")
    }
    xlabel_string <- paste(name, ": ", my_start, " to ", my_end)
    my_plot <- ggplot(rpm_region, aes_string(x="position", y="log")) +
        ggplot2::xlab(xlabel_string) +
        ggplot2::ylab("Log2(RPM) reads") +
        ggplot2::geom_bar(data=rpm_region, stat="identity", fill="black", colour="black") +
        ggplot2::geom_bar(data=pre_start, stat="identity", fill="red", colour="red") +
        ggplot2::geom_bar(data=post_stop, stat="identity", fill="red", colour="red") +
        ggplot2::geom_segment(data=rpm_region, mapping=stupid, arrow=gene_arrow, size=2, color="blue") +
        ggplot2::theme_bw()
    plot(my_plot)

}

#' Make a bar plot of the numbers of significant genes by contrast.
#' These plots are quite difficult to describe.
#'
#' @param ups  Set of up-regulated genes.
#' @param downs  Set of down-regulated genes.
#' @param maximum  Maximum/minimum number of genes to display.
#' @param text  Add text at the ends of the bars describing the number of genes >/< 0 fc.
#' @param color_list  Set of colors to use for the bars.
#' @param color_names  Categories associated with aforementioned colors.
#' @return weird significance bar plots
#' @seealso \pkg{ggplot2}
#'  \code{\link{extract_significant_genes}}
#' @export
plot_significant_bar <- function(ups, downs, maximum=NULL, text=TRUE,
                                 color_list=c("lightcyan", "lightskyblue", "dodgerblue",
                                              "plum1", "orchid", "purple4"),
                                 color_names=c("a_up_inner", "b_up_middle", "c_up_outer",
                                               "a_down_inner", "b_down_middle", "c_down_outer")) {
    choose_max <- function(u, d) {
        ## m is the maximum found in the ups/downs
        m <- 0
        ## which is extracted from ups and downs
        um <- max(as.numeric(u))
        dm <- max(as.numeric(d))
        ## choose the maximum by which is biggest!
        if (um >= dm) {
            m <- um
        } else {
            m <- dm
        }
        ## Figure out the number of digits in the number
        digits <- nchar(as.character(m))
        ## And the number of zeroes in it.
        num_zeroes <- digits - 1.0
        ## Add 1 to the first digit
        first_digit <- as.numeric(strsplit(x=as.character(m), split="")[[1]][[1]]) + 1.0
        ## And set maximum to that number * 10 to the number of zeroes.
        maximum <- first_digit * (10 ^ num_zeroes)
        return(maximum)
    }

    up_sums <- list()
    down_sums <- list()
    comp_names <- ups[ ups[["variable"]] == "a_up_inner", ][["comparisons"]]
    for (comp in 1:length(comp_names)) {
        comp_name <- comp_names[[comp]]
        up_sums[[comp_name]] <- sum(as.numeric(ups[ ups[["comparisons"]] == comp_name, ][["value"]]))
        down_sums[[comp_name]] <- sum(as.numeric(downs[ downs[["comparisons"]] == comp_name, ][["value"]])) * -1.0
    }

    if (is.null(maximum)) {
        maximum <- choose_max(up_sums, down_sums)
    }

    ## Try to ensure that ggplot orders my colors and bars in the specific order I want.
    ## holy ass crackers this is annoying and difficult to get correct,  as the ordering is (to my eyes) arbitrary.
    names(color_list) <- color_names
    levels(ups[["variable"]]) <- c("c_up_outer", "b_up_middle", "a_up_inner")
    levels(downs[["variable"]]) <- c("c_down_outer", "b_down_middle", "a_down_inner")

    sigbar_plot <- ggplot() +
        ggplot2::geom_col(data=ups, aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::geom_col(data=downs, aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::scale_fill_manual(values=c("a_up_inner"="lightcyan",
                                            "b_up_middle"="lightskyblue",
                                            "c_up_outer"="dodgerblue",
                                            "a_down_inner"="plum1",
                                            "b_down_middle"="orchid",
                                            "c_down_outer"="purple4")) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
        ggplot2::theme(legend.position="none")

    if (isTRUE(text)) {
        for (comp in 1:length(comp_names)) {
            comp_name <- comp_names[[comp]]
            upstring <- as.character(up_sums[[comp_name]])
            downstring <- as.character(down_sums[[comp_name]])
            sigbar_plot <- sigbar_plot +
                ggplot2::annotate("text", x=comp, y=maximum, label=upstring, angle=-90) +
                ggplot2::annotate("text", x=comp, y=maximum * -1, label=downstring, angle=90)
        }
    }
    return(sigbar_plot)
}

## EOF  Damners I don't have many bar plots, do I?
