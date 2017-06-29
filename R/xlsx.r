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
#'  the last row/column written to the worksheet.
#' @seealso \pkg{openxlsx}
#' @examples
#'  \dontrun{
#'   xls_coords <- write_xls(dataframe, sheet="hpgl_data")
#'   xls_coords <- write_xls(another_df, sheet="hpgl_data", start_row=xls_coords$end_col)
#'  }
#' @export
write_xls <- function(data="undef", wb=NULL, sheet="first", rownames=TRUE,
                      start_row=1, start_col=1, ...) {
    arglist <- list(...)
    if (class(data) == "matrix" | class(data) == "character") {
        data <- as.data.frame(data)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="hpgltools")
    } else if (class(wb)[[1]] == "list") {
        ## In case the return from write_xls() was passed to write_xls()
        wb <- wb[["workbook"]]
    } else if (class(wb)[[1]] != "Workbook") {
        stop("A workbook was passed to this, but the format is not understood.")
    }

    newsheet <- NULL
    ##current_sheets <- names(wb@.xData$worksheets)
    current_sheets <- wb@.xData[[".->sheet_names"]]
    if (sheet %in% current_sheets) {
        message(paste0("The sheet: ", sheet, " is in ", toString(current_sheets), "."))
    } else {
        newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
        if (class(newsheet) == "try-error") {
            if (grepl(pattern="already exists", x=newsheet[1])) {
                message("The sheet already exists.")
                tt <- openxlsx::removeWorksheet(wb, sheet)
                newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
            } else if (grepl(pattern="too long", x=newsheet[1])) {
                message("The sheet name was too long for Excel, truncating it by removing vowels.")
                sheet <- gsub(pattern="a|e|i|o|u", replacement="", x=sheet)
                newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
                if (class(newsheet) == "try-error") {
                    if (grepl(pattern="already exists", x=newsheet[1])) {
                        tt <- openxlsx::removeWorksheet(wb, sheet)
                        newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
                    } else if (grepl(pattern="too long", x=newsheet[1])) {
                        message("hmph, still too long, truncating to 30 characters.")
                        sheet <- substr(sheet, start=1, stop=30)
                        newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet))
                        if (grepl(pattern="already exists", x=newsheet[1])) {
                            tt <- openxlsx::removeWorksheet(wb, sheet)
                            newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
                        }
                    }
                }
            } else {
                message(paste0("Unknown error: ", newsheet))
            }
        }
    }
    hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold",
                                 border="Bottom", fontSize="30")
    new_row <- start_row
    new_col <- start_col
    if (!is.null(arglist[["title"]])) {
        xl_result <- openxlsx::writeData(wb, sheet, arglist[["title"]],
                                         startRow=new_row, startCol=new_col)
        openxlsx::addStyle(wb, sheet, hs1, new_row, new_col)
        new_row <- new_row + 1
    }

    ## I might have run into a bug in openxlsx, in WorkbookClass.R there is a call to is.nan()
    ## for a data.frame and it appears to me to be called oddly and causing problems
    ## I hacked the writeDataTable() function in openxlsx and sent a bug report.
    ## Another way to trip this up is for a column in the table to be of class 'list'
    test_column <- 0
    for (col in colnames(data)) {
        test_column <- test_column + 1
        colnames(data)[test_column] <- paste0(colnames(data)[test_column], "_", test_column)
        ## Originally, this was a single test condition, but I fear I need to do separate tasks for each data type.
        ## If that proves to be the case, I am ready, but until then it remains a series of as.character() castings.
        if (class(data[, test_column]) == "list") {
          ##  message(paste0("The column: ", col, " is a list."))
            ## The above did not work, trying what I found in:
            ## https://stackoverflow.com/questions/15930880/unlist-all-list-elements-in-a-dataframe
            ##list_entries <- is.list(data[, test_column])
            ##ListCols <- sapply(data, is.list)
            ##cbind(data[!ListCols], t(apply(data[ListCols], 1, unlist)))
            data[, test_column] <- as.character(data[, test_column])
        } else if (class(data[, test_column]) == "vector") {
          ##  message(paste0("The column: ", col, " is a vector."))
            data[, test_column] <- as.character(data[, test_column])
        } else if (class(data[, test_column]) == "factor") {
          ##  message(paste0("The column: ", col, " is a factor."))
            data[, test_column] <- as.character(data[, test_column])
        } else if (class(data[, test_column]) == "AsIs") {
          ##  message(paste0("The column: ", col, " is an AsIs."))
            data[, test_column] <- as.character(data[, test_column])
        }
    }  ## Finished adjusting stupid column types.

    final_colnames <- colnames(data)
    final_colnames <- tolower(final_colnames)
    final_colnames <- make.names(final_colnames, unique=TRUE)
    colnames(data) <- final_colnames
    wtf <- try(openxlsx::writeDataTable(wb, sheet, data, startCol=new_col,
                                        startRow=new_row, tableStyle="TableStyleMedium9",
                                        rowNames=rownames, colNames=TRUE))
    new_row <- new_row + nrow(data) + 2
    ## Set the column lengths, hard set the first to 20,
    ## then try to set it to auto if the length is not too long.
    for (data_col in 1:ncol(data)) {
        ## Make an explicit check that the data is not null, which comes out here as character(0)
        test_column <- as.character(data[[data_col]])
        test_column[is.na(test_column)] <- ""
        test_null <- identical(as.character(test_column), character(0))
        test_max <- 4

        if (isTRUE(test_null)) {
            test_max <- 1
        } else {
            test_max <- max(nchar(as.character(test_column)), na.rm=TRUE)
        }

        ## Keep in mind that if we are going to set the column widths
        ## and we set a start_col, then the actual column we will be changing is start_col + data_col.
        current_col <- start_col + data_col - 1  ## start_col is 1 indexed.
        if (data_col == 1) {
            openxlsx::setColWidths(wb, sheet, current_col, 20)
        } else if (test_max > 30) {
            openxlsx::setColWidths(wb, sheet, current_col, 30)
        } else {
            openxlsx::setColWidths(wb, sheet, current_col, "auto")
        }
    }
    end_col <- ncol(data) + 1
    ret <- list(
        "workbook" = wb,
        "sheet" = sheet,
        "end_row" = new_row,
        "end_col" = end_col)
    return(ret)
}

#' An attempt to improve the behaivor of openxlsx's plot inserter.
#'
#' The functions provided by openxlsx for adding plots to xlsx files are quite
#' nice, but they can be a little annoying.  This attempt to catch some corner cases
#' and potentially save an extra svg-version of each plot inserted.
#'
#' @param a_plot  The plot provided
#' @param wb  Workbook to which to write.
#' @param sheet  Name or number of the sheet to which to add the plot.
#' @param width  Plot width in the sheet.
#' @param height  Plot height in the sheet.
#' @param res  Resolution of the png image inserted into the sheet.
#' @param plotname  Prefix of the pdf file created.
#' @param savedir  Directory to which to save pdf copies of the plots.
#' @param fancy_type  Plot publication quality images in this format.
#' @param start_row  Row on which to place the plot in the sheet.
#' @param start_col  Column on which to place the plot in the sheet.
#' @param file_type  Currently this only does pngs, but perhaps I will parameterize this.
#' @param units  Units for the png plotter.
#' @param ...  Extra arguments are passed to arglist (Primarily for vennerable plots which are odd)
#' @return  A list containing the result of the tryCatch{} used to invoke the plot prints.
#' @seealso \pkg{openxlsx}
#' @examples
#'  \dontrun{
#'   fun_plot <- plot_pca(stuff)$plot
#'   try_results <- xlsx_plot_png(fun_plot)
#' }
#' @export
xlsx_plot_png <- function(a_plot, wb=NULL, sheet=1, width=6, height=6, res=90,
                          plotname="plot", savedir="saved_plots", fancy_type="pdf",
                          start_row=1, start_col=1, file_type="png", units="in", ...) {
    arglist <- list(...)
    if (is.null(a_plot)) {
        return(NULL)
    }
    if (!is.null(arglist[["doWeights"]])) {
        requireNamespace(package="Vennerable")
        library("Vennerable")
    }

    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="hpgltools")
    } else if (class(wb)[[1]] == "list") {
        ## In case the return from write_xls() was passed to write_xls()
        wb <- wb[["workbook"]]
    } else if (class(wb)[[1]] != "Workbook") {
        stop("A workbook was passed to this, but the format is not understood.")
    }
    high_quality <- paste0(savedir, "/", plotname, ".", fancy_type)
    fancy_print_ret <- png_print_ret <- NULL
    if (!is.null(savedir)) {
        if (!file.exists(savedir)) {
            dir.create(savedir, recursive=TRUE)
        }
        high_quality <- paste0(savedir, "/", plotname, ".", fancy_type)
        if (fancy_type == "pdf") {
            fancy_ret <- try(pdf(file=high_quality))
        } else if (fancy_type == "ps") {
            fancy_ret <- try(postscript(file=high_quality))
        } else if (fancy_type == "svg") {
            fancy_ret <- try(svg(filename=high_quality))
        } else if (fancy_type == "emf") {
            fancy_ret <- try(devEMF::emf(file=high_quality))
        } else {
            ## Default to pdf
            high_quality_renamed <- gsub(pattern="\\..*$", replacement="\\.pdf", x=high_quality)
            fancy_ret <- try(pdf(file=high_quality_renamed))
        }

        ## I do not understand why some images are plot()ed while others
        ## seem to need to be print()ed.  Adding a try to attempt
        ## to work around this concern.
        if (class(a_plot)[[1]] == "Venn") {
            pdf_print_ret <- try(Vennerable::plot(a_plot, doWeights=FALSE))
        } else {
            pdf_print_ret <- try(print(a_plot))
        }
        if (class(pdf_print_ret)[[1]] == "try-error") {
            pdf_print_ret <- try(plot(a_plot, ...))
        }
        dev.off()
    }
    fileName <- tempfile(pattern = "figureImage", fileext = paste0(".", file_type))
    png_ret <- try(png(filename=fileName,
                       width=width,
                       height=height,
                       units=units,
                       res=res))

    if (class(a_plot)[[1]] == "Venn") {
        pdf_print_ret <- try(Vennerable::plot(a_plot, doWeights=FALSE))
    } else {
        png_print_ret <- try(print(a_plot))
    }
    if (class(png_print_ret)[[1]] == "try-error") {
        png_print_ret <- try(plot(a_plot, ...))
    }
    dev.off()
    insert_ret <- try(openxlsx::insertImage(wb=wb, sheet=sheet, file=fileName, width=width,
                                     height=height, startRow=start_row, startCol=start_col,
                                     units=units, dpi=res))
    if (class(insert_ret) == "try-error") {
        message(paste0("There was an error inserting the image at: ", fileName))
    }
    ret <- list(
        "png_print" = png_print_ret,
        "pdf_print" = pdf_print_ret,
        "openxlsx" = insert_ret)
    return(ret)
}

## EOF
