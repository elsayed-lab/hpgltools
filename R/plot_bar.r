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
            stupid = aes(y=0,yend=0,x=my_start,xend=my_end)
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

plot_significant_bar <- function(table_list, fc_cutoffs=c(0,1,2), maximum=7000,
                                 fc_column="limma_logfc", p_type="adj",
                                 p=0.05, z=NULL, order=NULL,
                                 color_list=c("purple4", "plum1", "orchid", "dodgerblue", "lightcyan", "lightskyblue")) {
    ## As currently written, this function is a piece of crap.
    ## It needs to be simplified so that the plotting is
    ## separated from the data collection.
    ## But as an initial implementation I just want it to work so...
    sig_lists_limma_up <- list()
    sig_lists_limma_down <- list()
    sig_lists_edger_up <- list()
    sig_lists_edger_down <- list()
    sig_lists_deseq_up <- list()
    sig_lists_deseq_down <- list()
    plots <- list()
    table_length <- 0
    fc_names <- c()
    for (fc in fc_cutoffs) {
        ## This is a bit weird and circuituous
        ## The most common caller of this function is in fact extract_significant_genes
        fc_sig <- sm(extract_significant_genes(table_list, fc=fc,
                                               p=p, z=z, n=NULL, excel=NULL,
                                               p_type=p_type, sig_bar=FALSE, ma=FALSE))

        table_length <- length(fc_sig[["limma"]][["ups"]])
        fc_name <- paste0("fc_", fc)
        fc_names <- append(fc_names, fc_name)

        for (tab in 1:table_length) { ## The table names are shared across methods and ups/downs
            table_name <- names(fc_sig[["limma"]][["ups"]])[tab]
            limma_table_up <- nrow(fc_sig[["limma"]][["ups"]][[table_name]])
            limma_table_down <- nrow(fc_sig[["limma"]][["downs"]][[table_name]])

            deseq_table_up <- nrow(fc_sig[["deseq"]][["ups"]][[table_name]])
            deseq_table_down <- nrow(fc_sig[["deseq"]][["downs"]][[table_name]])

            edger_table_up <- nrow(fc_sig[["edger"]][["ups"]][[table_name]])
            edger_table_down <- nrow(fc_sig[["edger"]][["downs"]][[table_name]])

            sig_lists_limma_up[[fc_name]][[table_name]] <- limma_table_up
            sig_lists_limma_down[[fc_name]][[table_name]] <- limma_table_down
            sig_lists_edger_up[[fc_name]][[table_name]] <- edger_table_up
            sig_lists_edger_down[[fc_name]][[table_name]] <- edger_table_down
            sig_lists_deseq_up[[fc_name]][[table_name]] <- deseq_table_up
            sig_lists_deseq_down[[fc_name]][[table_name]] <- deseq_table_down
        } ## End iterating through every table
    } ## End querying all fc cutoffs
    ## Now we need to collate the data and make the bars

    up_all_limma <- up_all_deseq <- up_all_edger <- numeric() ## The number of all genes FC > 0
    down_all_limma <- down_all_deseq <- down_all_edger <- numeric()  ## The number of all genes FC < 0
    up_mid_limma <- up_mid_deseq <- up_mid_edger <- numeric()  ## The number of genes 2<FC<4 (by default)
    down_mid_limma <- down_mid_deseq <- down_mid_edger <- numeric()  ## The number of genes -2>FC>-4
    up_max_limma <- up_max_deseq <- up_max_edger <- numeric()  ## The number of genes FC > 4
    down_max_limma <- down_max_deseq <- down_max_edger <- numeric()  ## The number of genes FC < -4
    ##  The bar graph looks like
    ######### #### #  <-- Total width is the number of all >1FC genes
    ##      ^    ^------- Total >0FC - the set between 4FC and 2FC
    ##      |------------ Total >0FC - the smallest set >4FC

    ups_limma <- ups_edger <- ups_deseq <- list()
    downs_limma <- downs_edger <- downs_deseq <- list()
    papa_bear <- fc_names[[1]]  ## Because it is the largest grouping
    mama_bear <- fc_names[[2]]  ## The middle grouping
    baby_bear <- fc_names[[3]]  ## And the smallest grouping
    for (t in 1:table_length) {
        table_name <- names(sig_lists_limma_up[[1]])[t]
        ## > 0 lfc
        everything_up <- sig_lists_limma_up[[papa_bear]][[table_name]]
        ## > 1 lfc
        mid_up <- sig_lists_limma_up[[mama_bear]][[table_name]]
        ## > 2 lfc
        exclusive_up <- sig_lists_limma_up[[baby_bear]][[table_name]]
        ## < 0 lfc
        everything_down <- sig_lists_limma_down[[papa_bear]][[table_name]]
        ## < 1 lfc
        mid_down <- sig_lists_limma_down[[mama_bear]][[table_name]]
        ## < 2 lfc
        exclusive_down <- sig_lists_limma_down[[baby_bear]][[table_name]]

        ## Ah, I think the problem is that by calculating the numbers a,b,c
        ## It is stacking them and so I am getting a final bar of the sum of a,b,c
        up_all_limma[[table_name]] <- everything_up
        down_all_limma[[table_name]] <- everything_down
        up_mid_limma[[table_name]] <- mid_up - exclusive_up
        down_mid_limma[[table_name]] <- mid_down - exclusive_down
        up_max_limma[[table_name]] <- exclusive_up
        down_max_limma[[table_name]] <- exclusive_down
        up_all_limma[[table_name]] <- up_all_limma[[table_name]] - up_mid_limma[[table_name]] - up_max_limma[[table_name]]
        down_all_limma[[table_name]] <- down_all_limma[[table_name]] - down_mid_limma[[table_name]] - down_max_limma[[table_name]]

        table_name <- names(sig_lists_deseq_up[[1]])[t]
        everything_up <- sig_lists_deseq_up[[papa_bear]][[table_name]]
        mid_up <- sig_lists_deseq_up[[mama_bear]][[table_name]]
        exclusive_up <- sig_lists_deseq_up[[baby_bear]][[table_name]]
        everything_down <- sig_lists_deseq_down[[papa_bear]][[table_name]]
        mid_down <- sig_lists_deseq_down[[mama_bear]][[table_name]]
        exclusive_down <- sig_lists_deseq_down[[baby_bear]][[table_name]]

        up_all_deseq[[table_name]] <- everything_up
        down_all_deseq[[table_name]] <- everything_down
        up_mid_deseq[[table_name]] <- mid_up - exclusive_up
        down_mid_deseq[[table_name]] <- mid_down - exclusive_down
        up_max_deseq[[table_name]] <- exclusive_up
        down_max_deseq[[table_name]] <- exclusive_down
        up_all_deseq[[table_name]] <- up_all_deseq[[table_name]] - up_mid_deseq[[table_name]] - up_max_deseq[[table_name]]
        down_all_deseq[[table_name]] <- down_all_deseq[[table_name]] - down_mid_deseq[[table_name]] - down_max_deseq[[table_name]]

        table_name <- names(sig_lists_edger_up[[1]])[t]
        everything_up <- sig_lists_edger_up[[papa_bear]][[table_name]]
        mid_up <- sig_lists_edger_up[[mama_bear]][[table_name]]
        exclusive_up <- sig_lists_edger_up[[baby_bear]][[table_name]]
        everything_down <- sig_lists_edger_down[[papa_bear]][[table_name]]
        mid_down <- sig_lists_edger_down[[mama_bear]][[table_name]]
        exclusive_down <- sig_lists_edger_down[[baby_bear]][[table_name]]

        up_all_edger[[table_name]] <- everything_up
        down_all_edger[[table_name]] <- everything_down
        up_mid_edger[[table_name]] <- mid_up - exclusive_up
        down_mid_edger[[table_name]] <- mid_down - exclusive_down
        up_max_edger[[table_name]] <- exclusive_up
        down_max_edger[[table_name]] <- exclusive_down
        up_all_edger[[table_name]] <- up_all_edger[[table_name]] - up_mid_edger[[table_name]] - up_max_edger[[table_name]]
        down_all_edger[[table_name]] <- down_all_edger[[table_name]] - down_mid_edger[[table_name]] - down_max_edger[[table_name]]

    }

    comparisons <- names(sig_lists_limma_up[[1]])
    up_limma <- cbind(comparisons, up_max_limma, up_mid_limma, up_all_limma)
    down_limma <- cbind(comparisons, down_max_limma, down_mid_limma, down_all_limma)
    up_limma <- as.data.frame(up_limma)
    down_limma <- as.data.frame(down_limma)
    colnames(up_limma) <- c("comparisons","up_max","up_mid","up_all")
    colnames(down_limma) <- c("comparisons","down_max","down_mid","down_all")
    up_limma <- reshape2::melt(up_limma, id.var="comparisons")
    down_limma <- reshape2::melt(down_limma, id.var="comparisons")
    up_limma[["comparisons"]] <- factor(up_limma[["comparisons"]], levels=comparisons)
    down_limma[["comparisons"]] <- factor(down_limma[["comparisons"]], levels=comparisons)
    up_limma[["value"]] <- as.numeric(up_limma[["value"]])
    down_limma[["value"]] <- as.numeric(down_limma[["value"]]) * -1

    limma_plot <- ggplot() +
        ggplot2::geom_bar(data=up_limma, stat="identity",
                          aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::geom_bar(data = down_limma, stat="identity",
                          aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::scale_fill_manual(values=color_list) +
        ggplot2::scale_y_continuous(breaks=seq(maximum * -1, maximum, (maximum / 5))) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
        ggplot2::theme(legend.position="none")

    comparisons <- names(sig_lists_deseq_up[[1]])
    up_deseq <- cbind(comparisons, up_max_deseq, up_mid_deseq, up_all_deseq)
    down_deseq <- cbind(comparisons, down_max_deseq, down_mid_deseq, down_all_deseq)
    up_deseq <- as.data.frame(up_deseq)
    down_deseq <- as.data.frame(down_deseq)
    colnames(up_deseq) <- c("comparisons","up_max","up_mid","up_all")
    colnames(down_deseq) <- c("comparisons","down_max","down_mid","down_all")
    up_deseq <- reshape2::melt(up_deseq, id.var="comparisons")
    down_deseq <- reshape2::melt(down_deseq, id.var="comparisons")
    up_deseq[["comparisons"]] <- factor(up_deseq[["comparisons"]], levels=comparisons)
    down_deseq[["comparisons"]] <- factor(down_deseq[["comparisons"]], levels=comparisons)
    up_deseq[["value"]] <- as.numeric(up_deseq[["value"]])
    down_deseq[["value"]] <- as.numeric(down_deseq[["value"]]) * -1

    deseq_plot <- ggplot() +
        ggplot2::geom_bar(data=up_deseq, stat="identity",
                          aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::geom_bar(data=down_deseq, stat="identity",
                          aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::scale_fill_manual(values=color_list) +
        ggplot2::scale_y_continuous(breaks=seq(maximum * -1, maximum, (maximum / 5))) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
        ggplot2::theme(legend.position="none")

    comparisons <- names(sig_lists_edger_up[[1]])
    up_edger <- cbind(comparisons, up_max_edger, up_mid_edger, up_all_edger)
    down_edger <- cbind(comparisons, down_max_edger, down_mid_edger, down_all_edger)
    up_edger <- as.data.frame(up_edger)
    down_edger <- as.data.frame(down_edger)
    colnames(up_edger) <- c("comparisons","up_max","up_mid","up_all")
    colnames(down_edger) <- c("comparisons","down_max","down_mid","down_all")
    up_edger <- reshape2::melt(up_edger, id.var="comparisons")
    down_edger <- reshape2::melt(down_edger, id.var="comparisons")
    up_edger[["comparisons"]] <- factor(up_edger[["comparisons"]], levels=comparisons)
    down_edger[["comparisons"]] <- factor(down_edger[["comparisons"]], levels=comparisons)
    up_edger[["value"]] <- as.numeric(up_edger[["value"]])
    down_edger[["value"]] <- as.numeric(down_edger[["value"]]) * -1

    edger_plot <- ggplot() +
        ggplot2::geom_bar(data=up_edger, stat="identity",
                          aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::geom_bar(data=down_edger, stat="identity",
                          aes_string(x="comparisons", y="value", fill="variable")) +
        ggplot2::scale_fill_manual(values=color_list) +
        ggplot2::scale_y_continuous(breaks=seq(maximum * -1, maximum, (maximum / 5))) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
        ggplot2::theme(legend.position="none")

    retlist <- list(
        "limma_up_table" = up_limma,
        "limma_down_table"= down_limma,
        "limma" = limma_plot,
        "deseq_up_table" = up_deseq,
        "deseq_down_table"= down_deseq,
        "deseq" = deseq_plot,
        "edger_up_table" = up_edger,
        "edger_down_table"= down_edger,
        "edger" = edger_plot)
    return(retlist)
}

## EOF  Damners I don't have many bar plots, do I?
