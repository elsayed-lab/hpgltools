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
