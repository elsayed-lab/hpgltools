## Time-stamp: <Wed May 18 13:22:11 2016 Ashton Trey Belew (abelew@gmail.com)>

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

plot_rpm = function(input, output="~/riboseq/01.svg", name="LmjF.01.0010", start=1000, end=2000, strand=1, padding=100) {

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
            genename=row$Name
            output_file=paste("~/riboseq/coverage/", mychr, "/", genename, ".svg", sep="")
            svg(file=output_file, height=2, width=8, units="in")
            plot_rpm(rpms, start=row$start, end=row$end, strand=row$strand, output=output_file, name=row$Name)
            dev.off()
        }
    }


    ##for(i in 1:nrow(genes)) {
    ##    row <- genes[i,]
    ##    genename=row$Name
    ##    chromosome_number =
    ##        output_file=paste("~/riboseq/rpm/", genename, ".svg", sep="")
    ##    svg(file=output_file, height=2, width=6)
    ##    plot_rpm(rpms, st=row$start, en=row$end, strand=row$strand, output=output_file, name=row$name)
    ##    dev.off()
    ##}
    ##dev.new(height=2, width=6)
    ##plot_rpm(rpms, st=2000, en=20000, strand="+", output=test, name=row$Name)
    ##dev.off()

    mychr = gsub("\\.\\d+$", "", name, perl=TRUE)
    plotted_start = start - padding
    plotted_end = end + padding
    my_start = start
    my_end = end
    rpm_region = subset(input, chromosome==mychr & position >= plotted_start & position <= plotted_end)
    rpm_region = rpm_region[,-1]
    rpm_region$log = log2(rpm_region$rpm + 1)

    pre_start = subset(rpm_region, position < my_start)
    post_stop = subset(rpm_region, position > my_end)
    cds = subset(rpm_region, position >= my_start & position < my_end)

    theme_set(theme_bw())

    ## Wow R is such a piece of shit.
    ## Apparently aes messes with its scope in arbitrary ways, so you have to do stupid workarounds like this...
    ## Ok, fuck you then.
    eval(substitute(
        expr = {
            stupid = aes(y=0,yend=0,x=my_start,xend=my_end)
        },
        env = list(my_start=my_start, my_end=my_end)))

    if (strand == "+") {
        gene_arrow = arrow(type="closed", ends="last")
    } else {
        gene_arrow = arrow(type="closed", ends="first")
    }
    xlabel_string = paste(name, ": ", my_start, " to ", my_end)
    my_plot = ggplot(rpm_region, aes(x=position, y=log)) +
        xlab(xlabel_string) +
            ylab("Log2(RPM) reads") +
                geom_bar(data=rpm_region, stat="identity", fill="black", colour="black") +
                    geom_bar(data=pre_start, stat="identity", fill="red", colour="red") +
                        geom_bar(data=post_stop, stat="identity", fill="red", colour="red") +
                            geom_segment(data=rpm_region, mapping=stupid, arrow=gene_arrow, size=2, color="blue")
    plot(my_plot)

}

## EOF  Damners I don't have many bar plots, do I?
