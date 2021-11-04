
check_circos_prefix <- function(object, data_dir = "data", conf_dir = "conf",
                                karyotype_dir = "conf/karyotypes") {
  errors <- character()
  if (!file.exists(data_dir)) {
    msg <- message(data_dir, " does not exist. Creating the data directory now.")
    dir.create(data_dir, recursive = TRUE)
    errors <- c(errors, msg)
  }

  if (!file.exists(conf_dir)) {
    msg <- message("The circos directory does not exist, creating: ", conf_dir)
    dir.create(conf_dir, recursive = TRUE)
    errors <- c(errors, msg)
  }

  if (!file.exists(karyotype_dir)) {
    msg <- message("The karyotype directory does not exist, creating: ", karyotype_dir)
    dir.create(karyotype_dir, recursive = TRUE)
    errors <- c(errors, msg)
  }
  ideogram_dir <- file.path(conf_dir, "ideograms")
  if (!file.exists(ideogram_dir)) {
    msg <- message("The ideogram directory does not exist, creating: ", ideogram_dir)
    dir.create(ideogram_dir, recursive = TRUE)
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) TRUE else errors
}

## This causes roxygen to explode, and I am not sure why.
## Here is the error:
## Error in substituteFunctionArgs(validity, "object", functionName = sprintf("validity method for class '%s'",  :
##   trying to change the argument list of for validity method for class 'circos_prefix' with 4 arguments to have arguments (object)
## Calls: suppressPackageStartupMessages ... makeClassRepresentation -> .makeValidityMethod -> substituteFunctionArgs

##setClass("circos_prefix",
##         representation( name = "character",
##                         basedir = "character",
##                         cfg_file = "character",
##                         conf_dir = "character",
##                         data_dir = "character",
##                         annotation = "character",
##                         annot = "character",
##                         karyotype_cfg_file = "character",
##                         ideogram_cfg_file = "character",
##                         tick_cfg_file = "character",
##                         plus_df = "character",
##                         plus_cfg_file = "character",
##                         plus_data_file = "character",
##                         minus_df = "character",
##                         minus_cfg_file = "character",
##                         minus_data_file = "character"),
##         validity = check_circos_prefix)

#' Generic method to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setGeneric("iDA", signature=c("object"),
           function(object, ...) standardGeneric("iDA"))

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
