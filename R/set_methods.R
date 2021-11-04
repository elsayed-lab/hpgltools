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
