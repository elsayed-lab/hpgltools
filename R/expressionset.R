setMethod("assay",
          signature = "ExpressionSet",
          definition = function(x, withDimnames = TRUE, ...) {
            Biobase::exprs(x)
          })
setMethod("assay<-",
          signature = "ExpressionSet",
          definition = function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::exprs(x) <- value
            return(x)
          })
setMethod(
  "backup_expression_data", signature = signature(expt = "ExpressionSet"),
  definition = function(expt) {
    message("I do not know a good way to backup the expressionset data, returning.")
    return(expt)
  })
setMethod("colData",
          signature = "ExpressionSet",
          definition = function(x, withDimnames = TRUE, ...) {
            Biobase::pData(x)
          })
setMethod("colData<-",
          signature = "ExpressionSet",
          definition = function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::pData(x) <- value
            return(x)
          })
setMethod("get_backup_expression_data",
          signature = "ExpressionSet",
          definition = function(expt) {
            message("The expressionset currently does not keep a backup.")
            return(expt)
          })
setMethod("rowData",
          signature = "ExpressionSet",
          definition = function(x, withDimnames = TRUE, ...) {
            Biobase::fData(x)
          })
