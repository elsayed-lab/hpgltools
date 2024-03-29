## plot_hist.r: Histograms used in other functions

#' Make a pretty histogram of something.
#'
#' A shortcut to make a ggplot2 histogram which makes an attempt to set
#' reasonable bin widths and set the scale to log if that seems a good idea.
#'
#' @param df Dataframe of lots of pretty numbers.
#' @param binwidth Width of the bins for the histogram.
#' @param log Replot on the log scale?
#' @param bins Number of bins for the histogram.
#' @param adjust The prettification parameter in the ggplot2 density.
#' @param fillcolor Change the fill colors of the plotted elements?
#' @param color Change the color of the lines of the plotted elements?
#' @return Ggplot histogram.
#' @seealso [ggplot2]
#' @examples
#' \dontrun{
#'  kittytime = plot_histogram(df)
#' }
#' @export
plot_histogram <- function(df, binwidth = NULL, log = FALSE, bins = 500, adjust = 1,
                           fillcolor = "darkgrey", color = "black") {
  if (class(df) == "data.frame") {
    colnames(df) <- c("values")
  } else if (class(df) == "list") {
    df <- data.frame(unlist(df))
    colnames(df) <- c("values")
  } else if (class(df) == "numeric" || class(df) == "integer") {
    df <- data.frame(unlist(df))
    colnames(df) <- c("values")
  }
  if (is.null(binwidth)) {
    minval <- min(df, na.rm = TRUE)
    maxval <- max(df, na.rm = TRUE)
    binwidth <- (maxval - minval) / bins
  }
  density <- NULL
  a_histogram <- ggplot(df, aes(x = .data[["values"]])) +
    ggplot2::geom_histogram(aes(y = ggplot2::after_stat(density)),
                            binwidth = binwidth,
                            colour = color, fill = fillcolor, position = "identity") +
    ggplot2::geom_density(alpha = 0.4, fill = fillcolor, adjust = adjust) +
    ggplot2::geom_vline(aes(xintercept = mean(.data[["values"]], na.rm = TRUE)),
                        color = color, linetype = "dashed", size = 1) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"))
  if (log) {
    log_histogram <- try(a_histogram +
                           ggplot2::scale_x_continuous(trans = "log10"))
    if (log_histogram != "try-error") {
      a_histogram <- log_histogram
    }
  }
  return(a_histogram)
}

#' Make a pretty histogram of multiple datasets.
#'
#' If there are multiple data sets, it might be useful to plot them on a
#' histogram together and look at the t.test results between distributions.
#'
#' @param data Dataframe of lots of pretty numbers, this also accepts lists.
#' @param log Plot the data on the log scale?
#' @param binwidth Set a static bin width with an unknown # of bins?  If neither of these are
#'  provided, then bins is set to 500, if both are provided, then bins wins.
#' @param bins Set a static # of bins of an unknown width?
#' @param colors Change the default colors of the densities?
#' @return List of the ggplot histogram and some statistics describing the distributions.
#' @seealso [stats::pairwise.t.test()] [ggplot2]
#' @examples
#' \dontrun{
#'  kittytime = plot_multihistogram(df)
#' }
#' @export
plot_multihistogram <- function(data, log = FALSE, binwidth = NULL, bins = NULL, colors = NULL) {
  if (is.data.frame(data)) {
    df <- data
    columns <- colnames(df)
    summary_df <- summary(df)
    play_all <- data.frame()
    for (col in seq_along(colnames(df))) {
      new_column <- data.frame(expression = df[, col], cond = colnames(df)[col])
      play_all <- BiocGenerics::rbind(play_all, new_column)
    }
  } else if (is.list(data)) {
    summary_df <- summary(data)
    play_all <- reshape2::melt(data)
    colnames(play_all) <- c("expression", "cond")
  } else {
    stop("This can only work with a list or data frame.")
  }
  play_all[["expression"]] <- as.numeric(play_all[["expression"]])
  play_all[["cond"]] <- as.factor(play_all[["cond"]])
  play_cdf <- plyr::ddply(play_all, "cond",
                          plyr::summarise, rating.mean = mean(expression, na.rm = TRUE))
  uncor_t <- stats::pairwise.t.test(play_all[["expression"]],
                                    play_all[["cond"]], p.adjust = "none")
  bon_t <- try(stats::pairwise.t.test(play_all[["expression"]], play_all[["cond"]],
                                      p.adjust = "bon", na.rm = TRUE))
  if (is.null(bins) && is.null(binwidth)) {
    ##minval <- min(play_all[["expression"]], na.rm = TRUE)
    ##maxval <- max(play_all[["expression"]], na.rm = TRUE)
    ##bins <- 500
    ##binwidth <- (maxval - minval) / bins
    binwidth <- stats::density(play_all[["expression"]], na.rm = TRUE)[["bw"]]
  } else if  (is.null(binwidth)) {
    minval <- min(play_all[["expression"]], na.rm = TRUE)
    maxval <- max(play_all[["expression"]], na.rm = TRUE)
    binwidth <- (maxval - minval) / bins
  } else {
    message("Both bins and binwidth were provided, using binwidth: ", binwidth, sep = "")
  }
  multi <- ggplot(play_all,
                  aes(x = .data[["expression"]], fill = .data[["cond"]])) +
    ggplot2::geom_histogram(aes(y = ggplot2::after_stat(density)), binwidth = binwidth,
                            alpha = 0.4, position = "identity") +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(data = play_cdf,
                        aes(xintercept = .data[["rating.mean"]],  colour = .data[["cond"]]),
                        linetype = "dashed", size = 0.75) +
    ggplot2::xlab("Expression") +
    ggplot2::ylab("Observation likelihood") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"))
  if (!is.null(colors)) {
    multi <- multi +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_color_manual(values = colors)
  }
  if (log) {
    logged <- try(multi + ggplot2::scale_x_continuous(trans = "log10"))
    if (class(logged) != "try-error") {
      multi <- logged
    }
  }
  ##  if (class(bon_t) == "try-error") {
  ##    message("Unable to perform corrected test.")
  ##  } else {
  ##    message("Used Bonferroni corrected t test(s) between columns.")
  ##  }

  returns <- list(
    "plot" = multi,
    "data_summary" = summary_df,
    "uncor_t" = uncor_t,
    "bon_t" = bon_t)
  return(returns)
}

## EOF
