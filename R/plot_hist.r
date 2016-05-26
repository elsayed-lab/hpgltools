## plot_hist.r: Histograms used in other functions

#' Make a pretty histogram of something.
#'
#' A shortcut to make a ggplot2 histogram which makes an attempt to set reasonable bin widths and
#' set the scale to log if that seems a good idea.
#'
#' @param df Dataframe of lots of pretty numbers.
#' @param binwidth Width of the bins for the histogram.
#' @param log Replot on the log scale?
#' @param bins Number of bins for the histogram.
#' @param fillcolor Change the fill colors of the plotted elements?
#' @param color Change the color of the lines of the plotted elements?
#' @return Ggplot histogram.
#' @seealso \link[ggplot2]{geom_histogram} \link[ggplot2]{geom_density}
#' @examples
#' \dontrun{
#'  kittytime = plot_histogram(df)
#' }
#' @export
plot_histogram <- function(df, binwidth=NULL, log=FALSE, bins=500,
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
        message(paste("No binwidth provided, setting it to ", binwidth, " in order to have ", bins, " bins.", sep=""))
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
#' If there are multiple data sets, it might be useful to plot them on a histogram together and look
#' at the t.test results between distributions.
#'
#' @param data Dataframe of lots of pretty numbers, this also accepts lists.
#' @param log Plot the data on the log scale?
#' @param bins Set a static # of bins of an unknown width?
#' @param binwidth Set a static bin width with an unknown # of bins?  If neither of these are
#'     provided, then bins is set to 500, if both are provided, then bins wins.
#' @return List of the ggplot histogram and some statistics describing the distributions.
#' @seealso \link[stats]{pairwise.t.test} \link[plyr]{ddply}
#' @examples
#' \dontrun{
#'  kittytime = plot_multihistogram(df)
#' }
#' @export
plot_multihistogram <- function(data, log=FALSE, binwidth=NULL, bins=NULL) {
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

    returns <- list(
        "plot" = hpgl_multi,
        "data_summary" = summary_df,
        "uncor_t" = uncor_t,
        "bon_t" = bon_t)
    return(returns)
}

## EOF
