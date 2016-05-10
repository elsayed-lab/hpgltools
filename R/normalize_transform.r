## Time-stamp: <Tue May 10 12:20:04 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Perform a simple transformation of a count table (log2)
#'
#' the add argument is only important if the data was previously cpm'd because that does a +1, thus
#' this will avoid a double+1 on the data.
#'
#' @param count_table  A matrix of count data
#' @param transform   A type of transformation to perform: log2/log10/log
#' @param base   for other log scales
#' @return dataframe of logx(counts)
#' @examples
#' \dontrun{
#' filtered_table = transform_counts(count_table, transform='log2', converted='cpm')
#' }
#' @export
transform_counts <- function(count_table, transform="raw",
                             base=NULL, ...) {
    arglist <- list(...)
    ## Short circuit this if we are going with raw data.
    if (transform == "raw") {
        libsize <- colSums(count_table)
        counts <- list(count_table=count_table, libsize=libsize)
        return(counts)
    }

    ## If we are performing a transformation, then the minimum value I want is 1 before performing the logn
    less_zero <- sum(count_table < 0)
    if (less_zero > 0) {
        message(paste0("transform_counts: Found ", less_zero, " values less than 0."))
    }

    num_zero <- sum(count_table == 0)
    if (num_zero > 0) {
        message(paste0("transform_counts: Found ", num_zero, " values equal to 0, adding 1 to the matrix."))
        count_table <- count_table + 1
    }

    if (!is.null(base)) {
        count_table <- (log(count_table) / log(base))
    } else if (transform == "log2") {
        count_table <- log2(count_table)
    } else if (transform == "log10") {
        count_table <- log10(count_table)
    } else if (transform == "log") {  ## Natural log
        count_table <- log(count_table)  ## Apparently log1p does this.
    } else {
        message("Did not recognize the transformation, leaving the table.
 Recognized transformations include: 'log2', 'log10', 'log'
")
    }
    libsize <- colSums(count_table)
    counts <- list(
        "count_table" = count_table,
        "libsize" = libsize)
    return(counts)
}

