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
plot_libsize <- function(data, colors=NULL,
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

    colors <- as.character(colors)
    libsize_df <- data.frame(id=colnames(data),
                             sum=colSums(data),
                             condition=design[["condition"]],
                             colors=as.character(colors))
    libsize_df[["order"]] <- factor(libsize_df[["id"]], as.character(libsize_df[["id"]]))

    color_listing <- libsize_df[, c("condition","colors")]
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
                               ggplot2::aes_string(parse=FALSE,
                                                   x="order",
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

plot_rpm = function(input,
                    output="~/riboseq/01.svg",
                    name="LmjF.01.0010",
                    start=1000,
                    end=2000,
                    strand=1,
                    padding=100) {
    head(genes)
    genes = genes[,c(2,3,5,11)]
    for(ch in 1:36) {
        mychr = paste("LmjF.", sprintf("%02d", ch), sep="")
        print(mychr)
        table_path = paste("~/riboseq/coverage/", mychr, "/testme.txt.cov.gz", sep="")

        rpms = read.table(table_path)
        colnames(rpms) = c("chromosome","position","rpm")
        chr_index = with(genes, grepl(mychr, Name))
        genes_on_chr = genes[chr_index, ]

        for(i in 1:nrow(genes_on_chr)) {
            row = genes_on_chr[i,]
            genename=row[["Name"]]
            output_file=paste("~/riboseq/coverage/", mychr, "/", genename, ".svg", sep="")
            svg(filename=output_file, height=2, width=8)
            plot_rpm(rpms, start=row[["start"]], end=row[["end"]], strand=row[["strand"]], output=output_file, name=row[["Name"]])
            dev.off()
        }
    }

    mychr = gsub("\\.\\d+$", "", name, perl=TRUE)
    plotted_start = start - padding
    plotted_end = end + padding
    my_start = start
    my_end = end
    ## These are good chances to use %>% I guess
    ## rpm_region = subset(input, chromosome==mychr & position >= plotted_start & position <= plotted_end)
    region_idx <- input[["chromosome"]] == mychr & input[["position"]] >= plotted_start & input[["position"]] <= plotted_end
    rpm_region <- rpm_region[region_idx, ]
    rpm_region = rpm_region[,-1]
    rpm_region[["log"]] = log2(rpm_region[["rpm"]] + 1)

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
            stupid = aes(y=0, yend=0, x=my_start, xend=my_end)
        },
        env = list(my_start=my_start, my_end=my_end)))

    if (strand == "+") {
        gene_arrow = grid::arrow(type="closed", ends="last")
    } else {
        gene_arrow = grid::arrow(type="closed", ends="first")
    }
    xlabel_string = paste(name, ": ", my_start, " to ", my_end)
    my_plot = ggplot(rpm_region, aes_string(x="position", y="log")) +
        ggplot2::xlab(xlabel_string) +
        ggplot2::ylab("Log2(RPM) reads") +
        ggplot2::geom_bar(data=rpm_region, stat="identity", fill="black", colour="black") +
        ggplot2::geom_bar(data=pre_start, stat="identity", fill="red", colour="red") +
        ggplot2::geom_bar(data=post_stop, stat="identity", fill="red", colour="red") +
        ggplot2::geom_segment(data=rpm_region, mapping=stupid, arrow=gene_arrow, size=2, color="blue") +
        ggplot2::theme_bw()
    plot(my_plot)

}

plot_significant_bar <- function(ups, downs, maximum=NULL, text=TRUE,
                                 color_list=c("lightcyan", "plum1", ## The light colors
                                              "lightskyblue", "orchid", ## The mid colors
                                              "dodgerblue", "purple4"), ## And the darks
                                 name_list=c("a_up_inner", "a_down_inner",
                                             "b_up_middle", "b_down_middle",
                                             "c_up_outer", "c_down_outer")) {
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

    names(color_list) <- name_list
    sigbar_plot <- ggplot() +
        ggplot2::geom_bar(data=ups, stat="identity",
                          aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::geom_bar(data=downs, stat="identity",
                          aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::scale_fill_manual(values=color_list) +
        ggplot2::scale_y_continuous(breaks=seq(maximum * -1, maximum, (maximum / 5))) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
        ggplot2::theme(legend.position="none")

    if (isTRUE(text)) {
        for (comp in 1:length(comp_names)) {
            comp_name <- comp_names[[comp]]
            upstring <- as.character(up_sums[[comp_name]])
            downstring <- as.character(down_sums[[comp_name]])
            sigbar_plot = sigbar_plot +
                ggplot2::annotate("text", x=comp, y=maximum, label=upstring, angle=-90) +
                ggplot2::annotate("text", x=comp, y=maximum * -1, label=downstring, angle=90)
        }
    }
    return(sigbar_plot)
}

## EOF  Damners I don't have many bar plots, do I?
