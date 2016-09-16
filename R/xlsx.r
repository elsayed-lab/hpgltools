#' Write a dataframe to an excel spreadsheet sheet.
#'
#' I like to give folks data in any format they prefer, even though I sort
#' of hate excel.  Most people I work with use it, so therefore I do too.
#' This function has been through many iterations, first using XLConnect,
#' then xlsx, and now openxlsx.  Hopefully this will not change again.
#'
#' @param data Data frame to print.
#' @param wb Workbook to which to write.
#' @param sheet Name of the sheet to write.
#' @param rownames  Include row names in the output?
#' @param start_row First row of the sheet to write. Useful if writing multiple tables.
#' @param start_col First column to write.
#' @param ...  Set of extra arguments given to openxlsx.
#' @return List containing the sheet and workbook written as well as the bottom-right coordinates of
#'     the last row/column written to the worksheet.
#' @seealso \pkg{openxlsx} \link[openxlsx]{writeDataTable}
#' @examples
#' \dontrun{
#'  xls_coords <- write_xls(dataframe, sheet="hpgl_data")
#'  xls_coords <- write_xls(another_df, sheet="hpgl_data", start_row=xls_coords$end_col)
#' }
#' @export
write_xls <- function(data, wb=NULL, sheet="first", rownames=TRUE,
                      start_row=1, start_col=1, ...) {
    arglist <- list(...)
    if (class(data) == 'matrix') {
        data <- as.data.frame(data)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="hpgltools")
    } else if (class(wb)[[1]] == "list") { ## In case the return from write_xls() was passed to write_xls()
        wb <- wb[["workbook"]]
    } else if (class(wb)[[1]] != "Workbook") {
        stop("A workbook was passed to this, but the format is not understood.")
    }

    newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
    if (class(newsheet) == 'try-error') {
        if (grepl(pattern="too long", x=newsheet[1])) {
            message("The sheet name was too long for Excel, truncating it by removing vowels.")
            sheet <- gsub(pattern="a|e|i|o|u", replacement="", x=sheet)
            newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
            if (class(newsheet) == 'try-error') {
                message("hmph, still too long, truncating to 30 characters.")
                sheet <- substr(sheet, start=1, stop=30)
                newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet))
            }
        } else {
            message(paste0("The sheet ", sheet, " already exists, it will get overwritten"))
        }
    }
    hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")
    new_row <- start_row
    new_col <- start_col
    ##print(paste0("GOT HERE openxlswrite, title? ", arglist$title))
    if (!is.null(arglist[["title"]])) {
        openxlsx::writeData(wb, sheet, x=arglist[["title"]], startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
    }

    ## I might have run into a bug in openxlsx, in WorkbookClass.R there is a call to is.nan() for a data.frame
    ## and it appears to me to be called oddly and causing problems
    ## I hacked the writeDataTable() function in openxlsx and sent a bug report.
    ## Another way to trip this up is for a column in the table to be of class 'list'
    for (col in colnames(data)) {
        ## data[[col]] <- as.character(data[[col]])
        ## print(paste0("TESTME: ", class(data[[col]])))
        if (class(data[[col]]) == 'list' | class(data[[col]]) == 'vector' |
            class(data[[col]]) == 'factor' | class(data[[col]]) == 'AsIs') {
            ## message(paste0("Converted ", col, " to characters."))
            data[[col]] <- as.character(data[[col]])
        }
    }
    openxlsx::writeDataTable(wb, sheet, x=data, tableStyle="TableStyleMedium9",
                             startRow=new_row, rowNames=rownames, startCol=new_col)
    new_row <- new_row + nrow(data) + 2
    ## Set the column lengths, hard set the first to 20,
    ## then try to set it to auto if the length is not too long.
    for (col in 1:ncol(data)) {
        if (col == 1) {
            openxlsx::setColWidths(wb, sheet=sheet, widths=20, cols=col)
        } else if (max(nchar(data[[col]]), na.rm=TRUE) > 30) {
            openxlsx::setColWidths(wb, sheet=sheet, widths=30, cols=col)
        } else {
            openxlsx::setColWidths(wb, sheet=sheet, widths="auto", cols=col)
        }
    }
    end_col <- ncol(data) + 1
    ret <- list(
        "workbook" = wb,
        "end_row" = new_row,
        "end_col" = end_col)
    return(ret)
}

