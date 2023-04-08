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
