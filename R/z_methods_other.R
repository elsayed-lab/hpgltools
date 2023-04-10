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
            iDAoutput <- iDA::iDA_core(object, ...)
            return(iDAoutput)
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
