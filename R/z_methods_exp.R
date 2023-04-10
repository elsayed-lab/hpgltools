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


#' A series of setMethods for expts, ExpressionSets, and SummarizedExperiments.
#'
#' @param object One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param value New value to add to the object.
#' @param x I absoluately do not understand why R explicitly requires
#'  this.
#' @param i Nor this.
#' @importFrom SummarizedExperiment assay assay<- colData colData<- rowData rowData<-
#' @importFrom Biobase annotation annotation<-
setMethod("annotation", signature = "expt",
          function(object) {
            Biobase::annotation(object[["expressionset"]])
          })
setMethod("annotation<-", signature = "expt",
          function(object, value) {
            fData(object[["expressionset"]]) <- value
            return(object)
          })
setMethod("assay", signature = "expt",
          function(x, withDimnames = TRUE, ...) {
            mtrx <- Biobase::exprs(x[["expressionset"]])
            return(mtrx)
          })
setMethod("assay", "ExpressionSet",
          function(x, withDimnames = TRUE, ...) {
            Biobase::exprs(x)
          })
setMethod("assay<-", signature = "expt",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::exprs(x[["expressionset"]]) <- value
            return(x)
          })
setMethod("assay<-", signature = "ExpressionSet",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::exprs(x) <- value
            return(x)
          })
setMethod("colData", "expt",
          function(x, withDimnames = TRUE, ...) {
            Biobase::pData(x[["expressionset"]])
          })
setMethod("colData", "ExpressionSet",
          function(x, withDimnames = TRUE, ...) {
            Biobase::pData(x)
          })
setMethod("colData<-", signature = "expt",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::pData(x[["expressionset"]]) <- value
            return(x)
          })
setMethod("colData<-", signature = "ExpressionSet",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::pData(x) <- value
            return(x)
          })
setMethod("exprs", signature = "expt",
          function(object) {
            Biobase::exprs(object[["expressionset"]])
          })
setMethod("exprs<-", signature = "expt",
          function(object, value) {
            exprs(object[["expressionset"]]) <- value
            return(object)
          })
setMethod("fData", signature = "expt",
          function(object) {
            Biobase::fData(object[["expressionset"]])
          })
setMethod("fData<-", signature = "expt",
          function(object, value) {
            fData(object[["expressionset"]]) <- value
            return(object)
          })
setMethod("notes", signature = "expt",
          function(object) {
            Biobase::notes(object[["expressionset"]])
          })
setMethod("pData", signature = "expt",
          function(object) {
            Biobase::pData(object[["expressionset"]])
          })
setMethod("pData<-", signature = "expt",
          function(object, value) {
            pData(object[["expressionset"]]) <- value
            return(object)
          })
setMethod("rowData", "expt",
          function(x, withDimnames = TRUE, ...) {
            Biobase::fData(x[["expressionset"]])
          })
setMethod("rowData", "ExpressionSet",
          function(x, withDimnames = TRUE, ...) {
            Biobase::fData(x)
          })
setMethod("rowData<-", signature = "expt",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::fData(x[["expressionset"]]) <- value
            return(x)
          })
setMethod("sampleNames", signature = "expt",
          function(object) {
            Biobase::sampleNames(object[["expressionset"]])
          })
setMethod("sampleNames<-", signature = "expt",
          function(object, value) {
            set_expt_samplenames(object, value)
          })
