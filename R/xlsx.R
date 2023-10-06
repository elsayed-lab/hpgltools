## xlsx.r: Functions to simplify working with the xlsx format.  Most of the
## people with whom we work prefer Excel files. I am not particularly a fan and
## so wanted a way to create reasonably nice workbooks without intervention and
## therefore hopefully without significant chance of shenanigans.

#' Create the named worksheet in a workbook, this function was not well named.
#'
#' This tries to make sure that some of the problems of creating new worksheets
#' do not occur.  E.g. Names must be less than something and must be unique.
#'
#' @param wb Workbook to modify
#' @param sheet Sheet to check/create.
#' @return The workbook object hopefully with a new worksheet.
#' @seealso [openxlsx::addWorksheet()]
check_xlsx_worksheet <- function(wb, sheet) {
  newsheet <- NULL
  if (is.null(wb)) {
    return(NULL)
  }
  current_sheets <- wb@.xData[[".->sheet_names"]]
  found_sheets <- 0
  if (sheet %in% current_sheets) {
    found_sheets <- found_sheets + 1
    retlist <- list("wb" = wb, "sheet" = sheet)
    return(retlist)
  }

  newsheet <- try(openxlsx::addWorksheet(wb, sheetName = sheet), silent = TRUE)
  if (class(newsheet)[1] == "try-error") {
    if (grepl(pattern = "already exists", x = newsheet[1])) {
      message("The sheet already exists.")
      tt <- openxlsx::removeWorksheet(wb, sheet)
      newsheet <- try(openxlsx::addWorksheet(wb, sheetName = sheet), silent = TRUE)
    } else if (grepl(pattern = "too long", x = newsheet[1])) {
      sheet <- abbreviate(sheet, minlength = 28)
      newsheet <- try(openxlsx::addWorksheet(wb, sheetName = sheet), silent = TRUE)
    } else {
      message("Unknown error: ", newsheet)
    }
  }
  retlist <- list("wb" = wb, sheet = sheet)
  return(retlist)
}

#' Initialize an xlsx file with a little bit of logic to make sure there are no
#' annoying downstream errors.
#'
#' @param excel Excel file to create.
#' @return List containing the basename of the excel file along with the
#'  openxlsx workbook data structure.
#' @seealso [openxlsx::createWorkbook()]
#' @export
init_xlsx <- function(excel = "excel/something.xlsx") {
  if (isFALSE(excel)) {
    excel <- NULL
  }
  ## Thank you, Najib for this new and more robust regular expression.
  excel_basename <- basename(gsub(pattern = "\\.xlsx$|\\.xl.+$", replacement = "", x = excel))

  if (is.null(excel)) {
    return(NULL)
  }

  excel_dir <- dirname(as.character(excel))
  if (!file.exists(excel_dir)) {
    created <- dir.create(excel_dir, recursive = TRUE)
  }
  write_permission <- as.numeric(file.access(excel_dir, 2))
  if (write_permission < 0) {
    warning("The directory: ", excel_dir, " does not have write permission, this will fail.")
  }

  if (file.exists(excel)) {
    message("Deleting the file ", excel, " before writing the tables.")
    removed <- file.remove(excel)
  }
  wb <- openxlsx::createWorkbook(creator = "hpgltools")
  retlist <- list(
      "dirname" = excel_dir,
      "basename" = excel_basename,
      "wb" = wb)
  return(retlist)
}

#' Sanitize unreliable presentation of percent values from excel.
#'
#' In a recent sample sheet, we had some percentage values which were
#' '0.5', '5%', ' 6%' and a few other weirdo things. This function
#' should sanitize such shenanigans.
#'
#' @param numbers Either a vector of excel crap, or a column
#'  name/number.
#' @param df When provided, a data frame from which to extract the
#'  numbers.
#' @return Either the numbers or dataframe with the sanitized information.
#' @export
sanitize_percent <- function(numbers, df = NULL) {
  number_column <- NULL
  if (!is.null(df)) {
    number_column <- numbers
    numbers <- df[[numbers]]
  }

  ## Count up the things we changed
  num_pct_removed <- 0
  num_gt_one <- 0
  num_na <- 0
  num_other <- 0
  num_numeric <- 0

  numbers <- gsub(pattern = "\\s+", replacement = "", x = numbers)

  encoded <- grepl(x = numbers, pattern = "\\%")
  numericable <- suppressWarnings(as.numeric(numbers))
  num_encoded <- sum(encoded)
  for (n in seq_along(numbers)) {
    start <- numbers[n]
    encode <- encoded[n]

    if (isTRUE(encode)) {
      re_encoded <- gsub(x = start, pattern = "^.*?(\\d+)\\%.*$",
                         replacement = "\\1", perl = TRUE)
      re_encoded
      re_encoded <- as.numeric(re_encoded) / 100.0
      num_pct_removed <- num_pct_removed + 1
    } else if (is.na(start))  {
      re_encoded <- NA
      num_na <- num_na + 1
    } else if (is.na(numericable[n])) {
      re_encoded <- NA
      num_other <- num_other + 1
    } else if (as.numeric(start) > 1) {
      re_encoded <- as.numeric(start) / 100.0
      num_gt_one <- num_gt_one + 1
    } else {
      re_encoded <- as.numeric(start)
      num_numeric <- num_numeric + 1
    }
    numbers[n] <- re_encoded
  }
  numbers <- as.numeric(numbers)

  mesg("Re-encoded the following: ")
  mesg("Contained a '%': ", num_pct_removed, ".")
  mesg("Written as greater than 1: ", num_gt_one, ".")
  mesg("Written as a non-number: ", num_other, ".")
  mesg("Written as NA: ", num_na, ".")
  mesg("Written as a normal number from 0-1: ", num_numeric, ".")

  ## If a df was provided, return that instead of the number vector.
  if (!is.null(df)) {
    df[[number_column]] <- numbers
    numbers <- df
  }
  return(numbers)
}

#' Re-encode South American numbers so they fit my preconceptions.
#'
#' Ensure that we handle numbers encoded as '4.012.321,10' are
#' properly (from the perspective of R using my encoding system)
#' interpreted as 'four million twelve thousand three hundred
#' twenty-one and one tenth.'
#'
#' @param numbers Column of numbers.
#' @param df optional df rather than just a vector.
sanitize_number_encoding <- function(numbers, df = NULL) {
  number_column <- NULL
  if (!is.null(df)) {
    number_column <- numbers
    numbers <- df[[numbers]]
  }
  numbers <- gsub(pattern = "\\s+", replacement = "", x = numbers)

  ## If any numbers in the set have a . followed by a ,
  ## Then let us first get rid of any .'s and then replace
  ## the comma with a .
  encoded <- grepl(x = numbers, pattern = "\\..*\\,")
  num_encoded <- sum(encoded)
  if (num_encoded > 0) {
    removed_dot <- gsub(x = numbers, pattern = "\\.",
                        replacement = "")
    replaced_comma <- gsub(x = removed_dot, pattern = "\\,",
                           replacement = "\\.")
    numbers <- as.numeric(replaced_comma)
  }

  if (!is.null(df)) {
    df[[number_column]] <- numbers
    numbers <- df
  }
  return(numbers)
}

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
#' @param excel Filename of final excel workbook to write
#' @param rownames Include row names in the output?
#' @param start_row First row of the sheet to write. Useful if writing multiple tables.
#' @param start_col First column to write.
#' @param title Title for this xlsx table.
#' @param ... Set of extra arguments given to openxlsx.
#' @return List containing the sheet and workbook written as well as the
#'  bottom-right coordinates of the last row/column written to the worksheet.
#' @seealso [openxlsx] [openxlsx::createWorkbook()] [openxlsx::writeData()]
#'  [openxlsx::writeDataTable()] [openxlsx::saveWorkbook()]
#' @examples
#'  \dontrun{
#'   xls_coords <- write_xlsx(dataframe, sheet = "hpgl_data", excel = "testing.xlsx")
#'   xls_coords <- write_xlsx(another_df, wb = xls_coords$workbook,
#'                           sheet = "hpgl_data", start_row = xls_coords$end_col)
#'  }
#' @export
write_xlsx <- function(data = NULL, wb = NULL, sheet = "first", excel = NULL,
                       rownames = TRUE, start_row = 1, start_col = 1,
                       title = NULL, number_format = "0.000", data_table = TRUE,
                       freeze_first_row = TRUE, freeze_first_column = TRUE, ...) {
  arglist <- list(...)
  if (is.null(data)) {
    return(NULL)
  }
  if (class(data)[[1]] == "list") {
    written <- NULL
    for (element in seq_along(data)) {
      sheet_name <- names(data)[element]
      one_df <- data[[sheet_name]]
      if (!is.null(written[["workbook"]])) {
        wb <- written[["workbook"]]
        print(names(wb))
      }
      written <- write_xlsx(
        data = one_df, wb = wb, sheet = sheet_name, excel = excel, rownames = rownames,
        start_row = start_row, start_col = start_col,
        freeze_first_column = freeze_first_column, title = title,
        freeze_first_row = freeze_first_row)
      print(names(written[["workbook"]]))
    }

    return(written)
  }
  if ("matrix" %in% class(data) || "character" %in% class(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  if ("data.table" %in% class(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  if (!is.null(excel)) {
    xlsx <- init_xlsx(excel = excel)
    if (is.null(wb)) {
      wb <- xlsx[["wb"]]
    }
  }

  if (!is.null(number_format)) {
    old_options <- options("openxlsx.numFmt" = number_format)
  }
  ## Heading style 1 (For titles)
  hs1 <- openxlsx::createStyle(fontColour = "#000000", halign = "LEFT",
                               textDecoration = "bold", border = "Bottom",
                               fontSize = "30")

  ## Create the new worksheet.
  wb_sheet <- check_xlsx_worksheet(wb, sheet)
  wb <- wb_sheet[["wb"]]
  sheet <- wb_sheet[["sheet"]]

  new_row <- start_row
  new_col <- start_col
  if (!is.null(title)) {
    openxlsx::addStyle(wb = wb, sheet = sheet, style = hs1, rows = new_row,
                       cols = new_col, gridExpand = FALSE, stack = FALSE)
    ## An important caveat:
    ## glue'd() strings are not just class character, but their own thing,
    ## which means that just dumping them in this situation leads to unexpected results.
    xl_result <- openxlsx::writeData(
                               wb = wb, sheet = sheet,
                               x = as.character(title),
                               startCol = new_col, startRow = new_row)
    new_row <- new_row + 1
  }

  ## I might have run into a bug in openxlsx, in WorkbookClass.R there is a call to is.nan()
  ## for a data.frame and it appears to me to be called oddly and causing problems
  ## I hacked the writeDataTable() function in openxlsx and sent a bug report.
  ## Another way to trip this up is for a column in the table to be of class 'list'
  test_column <- 0
  final_colnames <- colnames(data)
  final_colnames <- tolower(final_colnames)
  final_colnames <- make.unique(final_colnames, sep = "_")
  final_colnames <- make.names(final_colnames)
  colnames(data) <- final_colnames

  final_data <- as.data.frame(data)
  for (col in seq_len(ncol(final_data))) {
    column_name <- colnames(final_data)[col]
    if ("list" %in% class(final_data[[col]])) {
      ## The above did not work, trying what I found in:
      ## https://stackoverflow.com/questions/15930880/unlist-all-list-elements-in-a-dataframe
      ##list_entries <- is.list(data[, test_column])
      ##ListCols <- sapply(data, is.list)
      ##cbind(data[!ListCols], t(apply(data[ListCols], 1, unlist)))
      final_data[[col]] <- as.character(final_data[[col]])
    } else if ("vector" %in% class(final_data[[col]])) {
      final_data[[col]] <- as.character(final_data[[col]])
    } else if ("factor" %in% class(final_data[[col]])) {
      final_data[[col]] <- as.character(final_data[[col]])
    } else if ("AsIs" %in% class(final_data[[col]])) {
      final_data[[col]] <- as.character(final_data[[col]])
    }
  }  ## Finished adjusting stupid column types.
  ## One more check that there are no '.'s in the column names,
  ## apparently excel doesn't like that sometimes (but others it doesn't care)
  colnames(final_data) = gsub(x = colnames(final_data),
                              pattern = "\\.", replacement = "_")
  written <- NULL
  if (isTRUE(data_table)) {
  written <- try(openxlsx::writeDataTable(
    wb = wb, sheet = sheet, x = final_data,
    startCol = new_col, startRow = new_row, tableStyle = table_style,
    rowNames = rownames, colNames = TRUE))
  } else {
  written <- try(openxlsx::writeData(
    wb = wb, sheet = sheet, x = final_data,
    startCol = new_col, startRow = new_row,
    rowNames = rownames, colNames = TRUE))
  }
  new_row <- new_row + nrow(final_data) + 2

  ## Set the column lengths, hard set the first to 20,
  ## then try to set it to auto if the length is not too long.
  for (data_col in seq_len(ncol(final_data))) {
    ## Make an explicit check that the data is not null, which comes out here as character(0)
    test_column <- as.character(final_data[[data_col]])
    test_column[is.na(test_column)] <- ""
    test_null <- identical(as.character(test_column), character(0))
    test_max <- 4

    if (isTRUE(test_null)) {
      test_max <- 1
      final_data[[data_col]] <- NULL  ## Drop the offending column.
    } else {
      test_max <- max(nchar(as.character(test_column)), na.rm = TRUE)
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
  end_col <- new_col + ncol(final_data) + 1

  new_options <- options(old_options)

  frozen <- NULL
  if (isTRUE(freeze_first_row) && isTRUE(freeze_first_column)) {
    frozen <- openxlsx::freezePane(wb, sheet, firstCol = TRUE, firstRow = TRUE)
  } else if (isTRUE(freeze_first_column)) {
    frozen <- openxlsx::freezePane(wb, sheet, firstCol = TRUE)
  } else if (isTRUE(freeze_first_row)) {
    frozen <- openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  }

  ret <- list(
      "workbook" = wb,
      "sheet" = sheet,
      "frozen" = frozen,
      "end_row" = new_row,
      "end_col" = end_col)
  if (!is.null(excel)) {
    mesg("Saving to: ", excel)
    save_result <- openxlsx::saveWorkbook(wb, excel, overwrite = TRUE)
    ret[["save_result"]] <- save_result
  }
  return(ret)
}

#' An attempt to improve the behaivor of openxlsx's plot inserter.
#'
#' The functions provided by openxlsx for adding plots to xlsx files are quite
#' nice, but they can be a little annoying.  This attempt to catch some corner cases
#' and potentially save an extra svg-version of each plot inserted.
#'
#' @param a_plot The plot provided
#' @param wb Workbook to which to write.
#' @param sheet Name or number of the sheet to which to add the plot.
#' @param width Plot width in the sheet.
#' @param height Plot height in the sheet.
#' @param res Resolution of the png image inserted into the sheet.
#' @param plotname Prefix of the pdf file created.
#' @param savedir Directory to which to save pdf copies of the plots.
#' @param fancy Plot fancy plots with the xlsx file?
#' @param fancy_type Plot publication quality images in this format.
#' @param start_row Row on which to place the plot in the sheet.
#' @param start_col Column on which to place the plot in the sheet.
#' @param file_type Currently this only does pngs, but perhaps I will parameterize this.
#' @param units Units for the png plotter.
#' @param ... Extra arguments are passed to arglist (Primarily for vennerable
#'  plots which are odd)
#' @return List containing the result of the tryCatch{} used to invoke the plot prints.
#' @seealso [openxlsx::insertImage()]
#' @examples
#'  \dontrun{
#'   fun_plot <- plot_pca(stuff)$plot
#'   df <- some_data_frame
#'   wb <- write_xlsx(df, excel = "funkytown.xlsx")$workbook
#'   try_results <- xlsx_insert_png(fun_plot, wb = wb)
#' }
#' @export
xlsx_insert_png <- function(a_plot, wb = NULL, sheet = 1, width = 6, height = 6, res = 90,
                          plotname = "plot", savedir = "saved_plots", fancy = FALSE,
                          fancy_type = "pdf", start_row = 1, start_col = 1,
                          file_type = "png", units = "in", ...) {
  arglist <- list(...)
  if (is.null(a_plot)) {
    return(NULL)
  }
  if (!is.null(arglist[["doWeights"]])) {
    requireNamespace(package = "Vennerable")
    ## library("Vennerable")
  }

  if (is.null(wb)) {
    wb <- openxlsx::createWorkbook(creator = "hpgltools")
  } else if (class(wb)[1] == "list") {
    ## In case the return from write_xlsx() was passed to write_xlsx()
    wb <- wb[["workbook"]]
  } else if (class(wb)[1] != "Workbook") {
    stop("A workbook was passed to this, but the format is not understood.")
  }

  insert_ret <- fancy_ret <- png_ret <- print_ret <- NULL
  if (isTRUE(fancy)) {
    high_quality <- file.path(savedir, glue("{plotname}.{fancy_type}"))
    if (!is.null(savedir)) {
      if (!file.exists(savedir)) {
        created <- dir.create(savedir, recursive = TRUE)
      }
    }
    if (fancy_type == "pdf") {
      fancy_ret <- try(pdf(file = high_quality))
    } else if (fancy_type == "ps") {
      fancy_ret <- try(postscript(file = high_quality))
    } else if (fancy_type == "svg") {
      fancy_ret <- try(svg(filename = high_quality))
    } else if (fancy_type == "emf") {
      fancy_ret <- try(devEMF::emf(file = high_quality))
    } else {
      ## Default to pdf
      high_quality_renamed <- gsub(pattern = "\\..*$", replacement = "\\.pdf", x = high_quality)
      fancy_ret <- try(pdf(file = high_quality_renamed))
    }
    ## I do not understand why some images are plot()ed while others
    ## seem to need to be print()ed.  Adding a try to attempt
    ## to work around this concern.
    if (class(a_plot)[1] == "Venn") {
      fancy_ret <- try(Vennerable::plot(a_plot, doWeights = FALSE))
    } else {
      fancy_ret <- try(suppressWarnings(print(a_plot)))
    }
    if (class(fancy_ret)[1] == "try-error") {
      fancy_ret <- try(suppressWarnings(plot(a_plot, ...)))
    }
    dev.off()
  }
  png_name <- try(tmpmd5file(pattern = "figureImage", fileext = glue(".{file_type}")))
  if ("try-error" %in% class(png_name)) {
    warning("There are too many tmp files in your current Rtmp directory.")
    warning("You need to clean it out ASAP.")
    png_name <- try(tmpmd5file(pattern = "figureImage2", fileext = glue(".{file_type}")))
  }
  png_ret <- try(png(filename = png_name,
                     width = width,
                     height = height,
                     units = units,
                     res = res))

  if (class(a_plot)[1] == "Venn") {
    print_ret <- try(Vennerable::plot(a_plot, doWeights = FALSE))
  } else {
    print_ret <- try(suppressWarnings(print(a_plot)))
  }
  if (class(print_ret)[1] == "try-error") {
    print_ret <- try(suppressWarnings(plot(a_plot, ...)))
  }
  dev.off()

  ## Check that the worksheet exists and add the plot.
  wb_sheet <- check_xlsx_worksheet(wb, sheet)
  wb <- wb_sheet[["wb"]]
  sheet <- wb_sheet[["sheet"]]

  if (file.exists(png_name)) {
    insert_ret <- try(openxlsx::insertImage(wb = wb, sheet = sheet, file = png_name,
                                            width = width, height = height,
                                            startRow = start_row, startCol = start_col,
                                            units = units, dpi = res))
    if (class(insert_ret)[1] == "try-error") {
      message("There was an error inserting the image at: ", png_name)
    }
  } else {
    message("The png file name did not exist: ", png_name)
  }
  ret <- list(
      "filename" = png_name,
      "png_fh" = png_ret,
      "png_print" = print_ret,
      "openxlsx" = insert_ret)
  if (isTRUE(fancy)) {
    ret[["fancy_print"]] <- fancy_ret
  }
  return(ret)
}

## EOF
