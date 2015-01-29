## Time-stamp: "Wed Jan 28 16:40:35 2015 Ashton Trey Belew (abelew@gmail.com)"
## misc_functions.R contains a host of (hopefully) useful functions when
## dealing with genomic data in R.

### This file is intended to keep a set of functions which will prove
### elsewhere in the analysis of ribosome profiling data.  These are
### mostly stolen from other people's code, notably Ramzi Temanni,
### Keith Hughitt, Laura Dillon, Yuan Li, and Kwame.

#' Wrap cor() to include robust correlations
#'
#' @param df a data frame to test
#' @param method Correlation method to use.  Defaults to pearson.
#' Includes pearson, spearman, kendal, robust.
#' 
#' @return  correlation some fun correlation statistics
#' @seealso \code{\link{cor}}, \code{\link{cov}}, \code{\link{covRob}}
#' 
#' @export
#' @examples
#' ## hpgl_cor(df=df)
#' ## hpgl_cor(df=df, method="robust")
hpgl_cor = function(df=NULL, method="pearson", ...) {
    if (method == "robust") {
        robust_cov = robust::covRob(df, corr=TRUE)
        correlation = robust_cov$cov
    } else {
        correlation = stats::cor(df, method=method, ...)
    }
    return(correlation)
}


    
#' Write a dataframe to an excel spreadsheet sheet.
#'
#' @param data a dataframe of information
#' @param sheet the name of an excel sheet in a workbook.
#' @param file an excel workbook to which to write.  Defaults to "excel/workbook.xls"
#' @param rowname include rownames?  Defalts to no.
#' 
#' @return NULL, on the say it creates a workbook if necessary,
#' creates a sheet, and writes the data to it.
#' 
#' @seealso \code{\link{loadWorkbook}}, \code{\link{createSheet}},
#' \code{\link{writeWorksheet}}, \code{\link{saveWorkbook}}
#' 
#' @export
#' @examples
#' ## write_xls(dataframe, "hpgl_data")
#' ## Sometimes it is a good idea to go in and delete the workbook and
#' ## re-create it if this is used heavily, because it will get crufty.
write_xls = function(data=NULL, sheet="first", file="excel/workbook.xls", rowname="rownames", overwrite=FALSE) {
    excel_dir = dirname(file)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }
    if (is.na(file.info(file)$size)) {
        xls = XLConnect::loadWorkbook(file, create=TRUE)
    } else if (file.info(file)$size == 0) {
        file.remove(file)
        xls = XLConnect::loadWorkbook(file, create=TRUE)        
    } else {
        xls = XLConnect::loadWorkbook(file)
    }

    write_sheet = 0
    if (isTRUE(existsSheet(xls, sheet))) {
        if (isTRUE(overwrite)) {
            write_sheet = 1
        } else {
            write_sheet = 0
        }
    } else {
        write_sheet = 1
    }

    if (write_sheet == 1) {
        XLConnect::createSheet(xls, name=sheet)
        if (is.na(rowname)) {
            XLConnect::writeWorksheet(xls, data, sheet=sheet)
        } else {
            XLConnect::writeWorksheet(xls, data, sheet=sheet, rowname=rowname)
        }
        XLConnect::saveWorkbook(xls)
    }
}


#' A stupid distance function of a point against two axes
#'
#' @param firstterm the x-values of the points
#' @param secondterm the y-values of the points
#' @param firstaxis the x-value of the vertical axis
#' @param secondaxis the y-value of the second axis
#'
#' @return dataframe of the distances
#' @export
sillydist = function(firstterm, secondterm, firstaxis, secondaxis) {
    dataframe = data.frame(firstterm, secondterm)
    dataframe$x = (abs(dataframe[,1]) - abs(firstaxis)) / abs(firstaxis)
    dataframe$y = abs((dataframe[,2] - secondaxis) / secondaxis)
    dataframe$x = abs(dataframe[,1] / max(dataframe$x))
    dataframe$y = abs(dataframe[,2] / max(dataframe$y))
    dataframe$dist = abs(dataframe$x * dataframe$y)
    dataframe$dist = dataframe$dist / max(dataframe$dist)
    return(dataframe)
}

Beta.NA = function(y,X) {
    des = X[!is.na(y),]
    y1 = y[!is.na(y)]
    B = solve(t(des)%*%des)%*%t(des)%*%y1
    B
}
