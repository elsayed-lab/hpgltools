## Time-stamp: <Mon Apr 25 15:04:17 2016 Ashton Trey Belew (abelew@gmail.com)>

## plot_hist.r: Histograms used in other functions

#'   Make a pretty histogram of something.
#'
#' @param df  a dataframe of lots of pretty numbers.
#' @param binwidth   width of the bins for the histogram.
#' @param log   replot on the log scale?
#' @param bins  bins for the histogram
#' @param verbose   be verbose?
#' @param fillcolor   change the fill colors of the plotted elements.
#' @param color   change the color of the lines of the plotted elements.
#' @return a ggplot histogram
#' @seealso \link[ggplot2]{geom_histogram} \link[ggplot2]{geom_density}
#' @examples
#' \dontrun{
#'  kittytime = hpgl_histogram(df)
#' }
#' @export
hpgl_histogram <- function(df, binwidth=NULL, log=FALSE, bins=500, verbose=FALSE,
                           fillcolor="darkgrey", color="black") {
    hpgl_env <- environment()
    if (class(df) == "data.frame") {
        colnames(df) <- c("values")
    } else if (class(df) == "list") {
        df <- data.frame(unlist(df))
        colnames(df) <- c("values")
    } else if (class(df) == "numeric") {
        df <- data.frame(unlist(df))
        colnames(df) <- c("values")
    }
    if (is.null(binwidth)) {
        minval <- min(df, na.rm=TRUE)
        maxval <- max(df, na.rm=TRUE)
        binwidth <- (maxval - minval) / bins
        if (verbose) {
            message(paste("No binwidth provided, setting it to ", binwidth, " in order to have ", bins, " bins.", sep=""))
        }
    }
    a_histogram <- ggplot2::ggplot(df, ggplot2::aes_string(x="values"), environment=hpgl_env) +
        ggplot2::geom_histogram(ggplot2::aes_string(y="..density.."), stat="bin", binwidth=binwidth,
                                colour=color, fill=fillcolor, position="identity") +
        ggplot2::geom_density(alpha=0.4, fill=fillcolor) +
        ggplot2::geom_vline(ggplot2::aes_string(xintercept="mean(values, na.rm=T)"), color=color, linetype="dashed", size=1) +
        ggplot2::theme_bw()
    if (log) {
        log_histogram <- try(a_histogram + ggplot2::scale_x_log10())
        if (log_histogram != "try-error") {
            a_histogram <- log_histogram
        }
    }
    return(a_histogram)
}


#' Make a pretty histogram of multiple datasets.
#'
#' @param data  a dataframe of lots of pretty numbers, this also accepts lists.
#' @param log   plot the data on the log scale?
#' @param bins   set a static # of bins of an unknown width?
#' @param binwidth   set a static bin width with an unknown # of bins?  If neither of these are provided, then bins is set to 500, if both are provided, then bins wins.
#' @param verbose   be verbose?
#' @return a ggplot histogram comparing multiple data sets
#' Along the way this generates pairwise t tests of the columns of
#' data.
#' @seealso \link[stats]{pairwise.t.test} \link[plyr]{ddply}
#' @examples
#' \dontrun{
#'  kittytime = hpgl_multihistogram(df)
#' }
#' @export
hpgl_multihistogram <- function(data, log=FALSE, binwidth=NULL, bins=NULL, verbose=FALSE) {
    if (is.data.frame(data)) {
        df <- data
        columns <- colnames(df)
        summary_df <- summary(df)
        play_all <- data.frame()
        for (col in 1:length(colnames(df))) {
            new_column <- data.frame(expression=df[,col], cond=colnames(df)[col])
            play_all <- BiocGenerics::rbind(play_all, new_column)
        }
    } else if (is.list(data)) {
        summary_df <- summary(data)
        play_all <- reshape2::melt(data)
        colnames(play_all) <- c("expression","cond")
    } else {
        stop("This can only work with a list or data frame.")
    }
    play_cdf <- plyr::ddply(play_all, "cond", plyr::summarise, rating.mean=mean(expression, na.rm=TRUE))
    uncor_t <- stats::pairwise.t.test(play_all$expression, play_all$cond, p.adjust="none")
    bon_t <- try(stats::pairwise.t.test(play_all$expression, play_all$cond, p.adjust="bon", na.rm=TRUE))
    if (is.null(bins) & is.null(binwidth)) {
        minval <- min(play_all$expression, na.rm=TRUE)
        maxval <- max(play_all$expression, na.rm=TRUE)
        bins <- 500
        binwidth <- (maxval - minval) / bins
        message(paste("Setting binwidth to ", binwidth, " in order to have ", bins, " bins.", sep=""))
    } else if  (is.null(binwidth)) {
        minval <- min(play_all$expression, na.rm=TRUE)
        maxval <- max(play_all$expression, na.rm=TRUE)
        binwidth <- (maxval - minval) / bins
        message(paste("Setting binwidth to ", binwidth, " in order to have ", bins, " bins.", sep=""))
    } else if (is.null(bins)) {
        message(paste("Setting binwidth to ", binwidth, sep=""))
    } else {
        message("Both bins and binwidth were provided, using binwidth: ", binwidth, sep="")
    }
    hpgl_multi <- ggplot2::ggplot(play_all, ggplot2::aes_string(x="expression", fill="cond")) +
        ggplot2::geom_histogram(ggplot2::aes_string(y="..density.."), binwidth=binwidth, alpha=0.4, position="identity") +
        ggplot2::xlab("Expression") +
        ggplot2::ylab("Observation likelihood") +
        ggplot2::geom_density(alpha=0.5) +
        ggplot2::geom_vline(data=play_cdf, ggplot2::aes_string(xintercept="rating.mean",  colour="cond"), linetype="dashed", size=0.75) +
        ggplot2::theme_bw()
    if (log) {
        logged <- try(hpgl_multi + ggplot2::scale_x_log10())
        if (class(logged) != "try-error") {
            hpgl_multi <- logged
        }
    }
    if (verbose) {
        message("Summarise the data.")
        message(summary_df)
        message("Uncorrected t test(s) between columns:")
        message(uncor_t)
        if (class(bon_t) == "try-error") {
            message("Unable to perform corrected test.")
        } else {
            message("Bon Ferroni corrected t test(s) between columns:")
            message(bon_t)
        }
    }
    returns <- list(plot=hpgl_multi, data_summary=summary_df,
                    uncor_t=uncor_t, bon_t=bon_t)
    return(returns)
}

