#' Validation function when creating a circos class.
#'
#' This is the one of the first steps taken to make the circos plot
#' builder into an object oriented set of functions.  Thank you,
#' Theresa!
#'
#' @param object The object to check for validity.
#' @return TRUE or FALSE
check_circos <- function(object) {
  ret <- c()
  base_dir <- dirname(object@data_dir)
  conf_dir <- dirname(object@cfg_file)
  data_dir <- object@data_dir

  if (!file.exists(data_dir)) {
    msg <- message(data_dir, " does not exist. Creating the data directory now.")
    created <- dir.create(data_dir, recursive = TRUE)
    if (isFALSE(created)) {
      ret <- c(ret, data_dir)
    }
  }

  if (!file.exists(conf_dir)) {
    msg <- message("The circos directory does not exist, creating: ", conf_dir)
    created <- dir.create(conf_dir, recursive = TRUE)
    if (isFALSE(created)) {
      ret <- c(ret, conf_dir)
    }
  }

  if (length(ret) == 0) {
    ret <- TRUE
  }

  return(ret)
}

#' Create a class for circos data
setClass("circos",
         representation(
             name = "character",
             data_dir = "character",
             cfg_file = "character",
             karyotype_cfg_file = "character",
             ideogram_cfg_file = "character",
             tick_cfg_file = "character",
             plus_cfg_file = "character",
             plus_data_file = "character",
             minus_cfg_file = "character",
             minus_data_file = "character",
             annotation = "data.frame",
             annot = "data.frame",
             plus_df = "data.frame",
             minus_df = "data.frame"),
         validity = check_circos)

#' Generic method to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setGeneric("iDA", signature = c("object"),
           function(object, ...) {
             standardGeneric("iDA")
           })

#' Set method for matrix to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setMethod("iDA", "matrix",
          function(object, ...) {
            iDAoutput <- iDA_core(object, ...)
            return(iDAoutput)
          })

## Methods for plotting functions.
##setGeneric("plot_libsize",
##           valueClass = "matrix",
##           function(data, ...) {
##             standardGeneric("plot_libsize")
##             ## plot_libsize(data, ...)
##           })

#' Send a SummarizedExperiment to plot_libsize().
#'
#' @param data SummarizedExperiment presumably created by create_se().
#' @param condition Set of conditions observed in the metadata, overriding
#'  the metadata in the SE.
#' @param colors Set of colors for the plot, overriding the SE metadata.
#' @param text Print text with the counts/sample observed at the top of the bars?
#' @param order Optionally redefine the order of the bars of the plot.
#' @param title Plot title!
#' @param yscale Explicitly set the scale on the log or base10 scale.
#' @param expt_names Optionally change the names of the bars.
#' @param label_chars If the names of the bars are larger than this, abbreviate them.
#' @param ... Additonal arbitrary arguments.
#' @return Plot of library sizes and a couple tables describing the data.
#' @export
setMethod("plot_libsize",
          signature = signature(data = "SummarizedExperiment"),
          definition = function(data, condition = NULL, colors = NULL, text = TRUE,
                                order = NULL, plot_title = NULL, yscale = NULL,
                                expt_names = NULL, label_chars = 10, ...) {
            mtrx <- as.matrix(assay(data))
            condition <- metadata(data)[["conditions"]]
            colors <- metadata(data)[["colors"]]
            plot_libsize(mtrx, condition = condition, colors = colors, text = text,
                         order = order, plot_title = plot_title, yscale = yscale,
                         expt_names = expt_names, label_chars = label_chars,
                         ...)
          })
setMethod("plot_libsize",
          signature = signature(data = "expt"),
          definition = function(data, condition = NULL, colors = NULL, text = TRUE,
                                order = NULL, plot_title = NULL, yscale = NULL,
                                expt_names = NULL, label_chars = 10, ...) {
            mtrx <- exprs(data)
            condition <- pData(data)[["condition"]]
            colors = data[["colors"]]
            plot_libsize(mtrx, condition = condition, colors = colors, text = text,
                         order = order, plot_title = plot_title, yscale = yscale,
                         expt_names = expt_names, label_chars = label_chars, ...)
          })
setMethod("plot_libsize",
          signature = signature(data = "ExpressionSet"),
          definition = function(data, condition = NULL, colors = NULL, text = TRUE,
                                order = NULL, plot_title = NULL, yscale = NULL,
                                expt_names = NULL, label_chars = 10, ...) {
            mtrx <- exprs(data)
            condition <- pData(data)[["conditions"]]
            plot_libsize(mtrx, condition = condition, colors = colors,
                         text = text, order = order, plot_title = plot_title,
                         yscale = yscale, expt_names = expt_names, label_chars = label_chars,
                         ...)
          })
setMethod("plot_libsize",
          signature = signature(data = "data.frame", condition = "factor", colors = "character"),
          definition = function(data, condition, colors, text = TRUE,
                                order = NULL, plot_title = NULL, yscale = NULL,
                                expt_names = NULL, label_chars = 10, ...) {
            data <- as.matrix(data)
            plot_libsize(data, condition = condition, colors = colors,
                         text = text, order = order, plot_title = plot_title, yscale = yscale,
                         expt_names = expt_names, label_chars = label_chars, ...) # , ...)
          })

## The following seems reasonable to me, but is a good example of me not understanding
## R's method-based OO style, because it causes an infinite loop in function dispatch.

##setMethod("plot_libsize",
##          signature = signature(data = "matrix", condition = "factor", colors = "character"),
##          definition = function(data, condition, colors, text = TRUE,
##                                order = NULL, plot_title = NULL, yscale = NULL,
##                                expt_names = NULL, label_chars = 10, ...) {
##            data <- as.matrix(data)
##            plot_libsize(data, condition = condition, colors = colors,
##                         text = text, order = order, plot_title = plot_title, yscale = yscale,
##                         expt_names = expt_names, label_chars = label_chars, ...) # , ...)
##          })
