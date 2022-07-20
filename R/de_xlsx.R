## de_xlsx.r: Attempt to export DE results to consistent xlsx reports. This is,
## in effect, a dumping ground for the various metrics which people have asked
## me to include when attempting to decide if a given DE result is
## 'significant'.  I therefore try in this file to find an appropriate place to
## dump each metric into an excel workbook.

#' Combine portions of deseq/limma/edger table output.
#'
#' This hopefully makes it easy to compare the outputs from
#' limma/DESeq2/EdgeR on a table-by-table basis.
#'
#' @param apr Output from all_pairwise().
#' @param extra_annot Add some annotation information?
#' @param excel Filename for the excel workbook, or null if not printed.
#' @param excel_title Title for the excel sheet(s).  If it has the
#'  string 'YYY', that will be replaced by the contrast name.
#' @param keepers List of reformatted table names to explicitly keep
#'  certain contrasts in specific orders and orientations.
#' @param excludes List of columns and patterns to use for excluding genes.
#' @param adjp Perhaps you do not want the adjusted p-values for plotting?
#' @param include_limma Include limma analyses in the table?
#' @param include_deseq Include deseq analyses in the table?
#' @param include_edger Include edger analyses in the table?
#' @param include_ebseq Include ebseq analyses in the table?
#' @param include_basic Include my stupid basic logFC tables?
#' @param rownames Add rownames to the xlsx printed table?
#' @param add_plots Add plots to the end of the sheets with expression values?
#' @param loess Add time intensive loess estimation to plots?
#' @param plot_dim Number of inches squared for the plot if added.
#' @param compare_plots Add some plots comparing the results.
#' @param padj_type Add a consistent p adjustment of this type.
#' @param lfc_cutoff In this context, only used for plotting volcano/MA plots.
#' @param p_cutoff In this context, used for volcano/MA plots.
#' @param de_types Used for plotting pvalue/logFC cutoffs.
#' @param rda Write a rda file of the results.
#' @param ... Arguments passed to significance and abundance tables.
#' @return Table combining limma/edger/deseq outputs.
#' @seealso [all_pairwise()] [extract_significant_genes()]
#' @examples
#' \dontrun{
#'  expt <- create_expt(metadata="some_metadata.xlsx", gene_info=funkytown)
#'  big_result <- all_pairwise(expt, model_batch=FALSE)
#'  pretty <- combine_de_tables(big_result, table='t12_vs_t0')
#'  pretty <- combine_de_tables(big_result, table='t12_vs_t0',
#'                              keepers=list("avsb"=c("a","b")))
#'  pretty <- combine_de_tables(big_result, table='t12_vs_t0',
#'                              keepers=list("avsb"=c("a","b")),
#'                              excludes=list("description"=c("sno","rRNA")))
#' }
#' @export
combine_de_tables <- function(apr, extra_annot = NULL,
                              excel = NULL, excel_title = "Table SXXX: Combined Differential Expression of YYY",
                              keepers = "all", excludes = NULL, adjp = TRUE, include_limma = TRUE,
                              include_deseq = TRUE, include_edger = TRUE, include_ebseq = TRUE,
                              include_basic = TRUE, rownames = TRUE, add_plots = TRUE, loess = FALSE,
                              plot_dim = 6, compare_plots = TRUE, padj_type = "ihw",
                              lfc_cutoff = 1, p_cutoff = 0.05, de_types = c("limma", "deseq", "edger"),
                              rda = NULL, ...) {
  arglist <- list(...)
  retlist <- NULL
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  do_excel <- TRUE
  if (is.null(wb)) {
    do_excel <- FALSE
  }

  ## Create a list of image files so that they may be properly cleaned up
  ## after writing the xlsx file.
  image_files <- c()

  ## First pull out the data for each tool
  limma <- apr[["limma"]]
  deseq <- apr[["deseq"]]
  edger <- apr[["edger"]]
  ebseq <- apr[["ebseq"]]
  basic <- apr[["basic"]]

  ## If any of the tools failed, then we cannot plot stuff with confidence.
  if (!isTRUE(include_limma) | !isTRUE(include_deseq) |
      !isTRUE(include_edger) | !isTRUE(include_basic)) {
    add_plots <- FALSE
    mesg("One or more methods was excluded.  Not adding the plots.")
  }
  if ("try-error" %in% class(limma) | is.null(limma)) {
    include_limma <- FALSE
    add_plots <- FALSE
    message("Not adding plots, limma had an error.")
  }
  if ("try-error" %in% class(deseq) | is.null(deseq)) {
    include_deseq <- FALSE
    add_plots <- FALSE
    message("Not adding plots, deseq had an error.")
  }
  if ("try-error" %in% class(edger) | is.null(edger)) {
    include_edger <- FALSE
    add_plots <- FALSE
    message("Not adding plots, edger had an error.")
  }
  if ("try-error" %in% class(ebseq) | is.null(ebseq)) {
    include_ebseq <- FALSE
  }
  if ("try-error" %in% class(basic) | is.null(basic)) {
    include_basic <- FALSE
    add_plots <- FALSE
    message("Not adding plots, basic had an error.")
  }

  ## A common request is to have the annotation data added to the table.  Do that here.
  annot_df <- fData(apr[["input"]])
  if (!is.null(extra_annot)) {
    annot_df <- merge(annot_df, extra_annot, by = "row.names", all.x = TRUE)
    rownames(annot_df) <- annot_df[["Row.names"]]
    annot_df <- annot_df[, -1, drop = FALSE]
  }

  ## Write the legend.
  legend <- write_combined_legend(wb, excel_basename, plot_dim, apr,
                                  limma, include_limma,
                                  deseq, include_deseq,
                                  edger, include_edger,
                                  ebseq, include_ebseq,
                                  basic, include_basic,
                                  padj_type)
  image_files <- c(legend[["image_files"]], image_files)
  table_names <- legend[["table_names"]]

  ## Because we sometimes use different names for tables in the keepers
  ## This is probably not the correct way to handle the list of coefficients.
  ## Instead, we should just grab the conditions list from deseq/edger/etc
  all_coefficients <- unlist(strsplit(x = table_names, split = "_vs_"))
  ## But I want to see what happend before I change this.
  ## all_coefficients <- apr[["deseq"]][["conditions"]]

  ## Extract and write the data tables.
  extracted <- list(
      "data" = list(),
      "table_names" = list(),
      "plots" = list(),
      "summaries" = data.frame())
  if ("list" %in% class(keepers)) {
    ## Here, we will look for only those elements in the keepers list.
    ## In addition, if someone wanted a_vs_b, but we did b_vs_a, then this will
    ## flip the logFCs.
    extracted <- extract_keepers_lst(extracted, keepers,
                                     legend[["table_names"]],
                                     all_coefficients,
                                     limma, edger, ebseq, deseq, basic,
                                     adjp, annot_df,
                                     include_deseq, include_edger,
                                     include_ebseq, include_limma,
                                     include_basic, excludes, padj_type,
                                     loess =  loess, lfc_cutoff = lfc_cutoff,
                                     p_cutoff =  p_cutoff)
  } else if (class(keepers)[1] == "character" & keepers[1] == "all") {
    ## If you want all the tables in a dump
    ## The logic here is the same as above without worrying about a_vs_b, but
    ## instead just iterating through every returned table, combining them, and
    ## printing them to the excel.
    extracted <- extract_keepers_all(apr, extracted, keepers,
                                     legend[["table_names"]],
                                     all_coefficients,
                                     limma, edger, ebseq, deseq, basic,
                                     adjp, annot_df,
                                     include_deseq, include_edger,
                                     include_ebseq, include_limma,
                                     include_basic, excludes, padj_type,
                                     loess = loess, lfc_cutoff = lfc_cutoff,
                                     p_cutoff =  p_cutoff)
  } else if (class(keepers)[1] == "character") {
    ## Finally, the simplest case, just print a single table.  Otherwise the logic
    ## should be identical to the first case above.
    extracted <- extract_keepers_single(apr, extracted, keepers,
                                        legend[["table_names"]],
                                        all_coefficients,
                                        limma, edger, ebseq, deseq, basic,
                                        adjp, annot_df,
                                        include_deseq, include_edger,
                                        include_ebseq, include_limma,
                                        include_basic, excludes, padj_type,
                                        loess =  loess, lfc_cutoff = lfc_cutoff,
                                        p_cutoff = p_cutoff)
  } else {
    stop("I don't know what to do with your specification of tables to keep.")
  } ## End different types of things to keep.
  numerators <- extracted[["numerators"]]
  denominators <- extracted[["denominators"]]

  ## At this point, we have done everything we can to combine the requested
  ## tables. So lets dump the tables to the excel file and compare how the
  ## various tools performed with some venn diagrams, and finally dump the plots
  ## from above into the sheet.
  venns <- list()
  venns_sig <- list()
  comp <- list()

  ## The following if() is too long and should be split into its own function.
  if (isTRUE(do_excel)) {
    ## Starting a new counter of sheets.
    ## I am considering adding some logic to collect the linear models
    ## Then check to see if the slopes/intercepts are duplicated across any
    ## of the contrasts, if this is true, then it is highly likely a mistake was made
    ## when setting up the contrasts such that something got duplicated.
    tnames <- names(extracted[["table_names"]])
    tsources <- as.character(extracted[["table_names"]])
    for (x in 1:length(tnames)) {
      tab <- tnames[x]
      written_table <- extracted[["data"]][[tab]]
      if (! "data.frame" %in% class(written_table)) {
        message("There is no data for ", tab, ", skipping it.")
        next
      }
      final_excel_title <- gsub(pattern = "YYY", replacement = tab, x = excel_title)
      final_excel_title <- glue("{final_excel_title}; Contrast numerator: {numerators[x]}.  Contrast denominator: {denominators[x]}.")
      ## Dump each table to the appropriate excel sheet
      xls_result <- write_xlsx(data = written_table, wb = wb, sheet = tab,
                               title = final_excel_title, rownames = rownames)
      ## The function write_xlsx has some logic in it to get around excel name
      ## limitations (30 characters), therefore set the sheetname to what was
      ## returned in case it had to change the sheet's name.
      sheetname <- xls_result[["sheet"]]
      current_row <- 1
      venn_rows <- 16
      plot_rows <- 31
      current_column <- xls_result[["end_col"]] + 2
      venn_columns <- 4
      plot_columns <- 10
      if (isTRUE(add_plots)) {
        ## Text on row 1, plots from 2-17 (15 rows)
        mesg("Adding venn plots for ", tnames[x], ".")
        ## Make some venn diagrams comparing deseq/limma/edger!
        venn_nop_lfc0 <- try(de_venn(written_table, lfc = 0, adjp = FALSE, p = 1.0))
        venn_nop <- try(de_venn(written_table, lfc = 1, adjp = FALSE, p = 1.0))
        venn_list <- try(de_venn(written_table, lfc = 0, adjp = adjp))
        venn_sig_list <- try(de_venn(written_table, lfc = 1, adjp = adjp))
        venns[[tab]] <- list(venn_nop_lfc0, venn_nop, venn_list, venn_sig_list)
        names(venns[[tab]]) <- c("nop_lfc0", "nop_lfc1", "p_lfc0", "p_lfc1")
        ## If they worked, add them to the excel sheets after the data,
        ## but make them smaller than other graphs.
        if (class(venn_list)[1] != "try-error") {
          ## First row of plots all going up
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     x = "Venn of all genes, lfc > 0.",
                                     startRow = current_row, startCol = current_column)
          up_plot <- venn_nop_lfc0[["up_venn"]]
          try_result <- xlsx_plot_png(
              up_plot, wb = wb, sheet = sheetname, width = (plot_dim / 2), height = (plot_dim / 2),
              start_col = current_column, plotname = "lfc0upvennnop", savedir = excel_basename,
              start_row = current_row + 1, doWeights = FALSE)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          current_column <- current_column + venn_columns
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     x = "Venn of all genes, lfc > 1.",
                                     startRow = 1, startCol = current_column)
          up_plot <- venn_nop[["up_venn"]]
          try_result <- xlsx_plot_png(
              up_plot, wb = wb, sheet = sheetname, width = (plot_dim / 2), height = (plot_dim / 2),
              start_col = current_column, plotname = "upvennnop", savedir = excel_basename,
              start_row = current_row + 1, doWeights = FALSE)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          current_column <- current_column + venn_columns
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     x = "Venn of p-value up genes, lfc > 0.",
                                     startRow = 1, startCol = current_column)
          up_plot <- venn_list[["up_venn"]]
          try_result <- xlsx_plot_png(
              up_plot, wb = wb, sheet = sheetname, width = (plot_dim / 2), height = (plot_dim / 2),
              start_col = current_column, plotname = "upvenn", savedir = excel_basename,
              start_row = current_row + 1, doWeights = FALSE)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          current_column <- current_column + venn_columns
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     x = "Venn of p-value up genes, lfc > 1.",
                                     startRow = current_row, startCol = current_column)
          sig_up_plot <- venn_sig_list[["up_venn"]]
          try_result <- xlsx_plot_png(
              sig_up_plot, wb = wb, sheet = sheetname, width = (plot_dim / 2), height =  (plot_dim / 2),
              start_col = current_column, plotname = "upvenn", savedir = excel_basename,
              start_row = current_row + 1, doWeights = FALSE)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          sig_methods <- c()
          if (isTRUE(include_limma)) {
            sig_methods <- c("limma", sig_methods)
          }
          if (isTRUE(include_edger)) {
            sig_methods <- c("edger", sig_methods)
          }
          if (isTRUE(include_deseq)) {
            sig_methods <- c("deseq", sig_methods)
          }
          siggene_lst <- try(plot_num_siggenes(written_table, methods =  sig_methods))
          current_column <- current_column + venn_columns
          if (class(siggene_lst)[1] != "try-error") {
            xl_result <- openxlsx::writeData(
                                       wb = wb, sheet = sheetname,
                                       x = "Significant genes by fc going up.",
                                       startRow = current_row, startCol = current_column)
            try_result <- xlsx_plot_png(
                siggene_lst[["up"]], wb = wb, sheet =  sheetname,
                width = plot_dim, height = (plot_dim / 2),
                start_col =  current_column, plotname =  "siggenesup",
                savedir = excel_basename, start_row = current_row + 1, doWeights = FALSE)
            if (! "try-error" %in% class(try_result)) {
              image_files <- c(image_files, try_result[["filename"]])
            }
            current_column <- current_column + plot_columns
            xl_result <- openxlsx::writeData(
                                       wb = wb, sheet = sheetname,
                                       x = "Significant genes by p going up.",
                                       startRow = 1, startCol = current_column)
            try_result <- xlsx_plot_png(
                siggene_lst[["pup"]], wb = wb, sheet = sheetname,
                width = plot_dim, height = (plot_dim / 2),
                start_col = current_column, plotname = "siggenespup",
                savedir = excel_basename, start_row = current_row + 1, doWeights = FALSE)
            if (! "try-error" %in% class(try_result)) {
              image_files <- c(image_files, try_result[["filename"]])
            }
          }
          ## Plot down venns etc, so reset the column but not the rows!
          current_column <- xls_result[["end_col"]] + 2
          current_row <- current_row + venn_rows
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     x = "Venn of all genes, lfc < 0.",
                                     startRow = current_row, startCol = current_column)
          down_plot <- venn_nop_lfc0[["down_venn"]]
          try_result <- xlsx_plot_png(
              down_plot, wb = wb, sheet = sheetname, width = plot_dim / 2, height = plot_dim / 2,
              start_col = current_column, plotname = "lfc0downvennnop", savedir = excel_basename,
              start_row = current_row + 1, doWeights = FALSE)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          current_column <- current_column + venn_columns
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     x = "Venn of all genes, lfc < -1.",
                                     startRow = current_row, startCol = current_column)
          down_plot <- venn_nop[["down_venn"]]
          try_result <- xlsx_plot_png(
              down_plot, wb = wb, sheet = sheetname, width = plot_dim / 2, height = plot_dim / 2,
              start_col = current_column, plotname = "downvennnop", savedir = excel_basename,
              start_row = current_row + 1, doWeights = FALSE)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          current_column <- current_column + venn_columns
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     x = "Venn of p-value down genes, lfc < 0.",
                                     startRow = current_row, startCol = current_column)
          down_plot <- venn_list[["down_venn"]]
          try_result <- xlsx_plot_png(
              down_plot, wb = wb, sheet = sheetname, width = plot_dim / 2, height = plot_dim / 2,
              start_col = current_column, plotname = "downvenn", savedir = excel_basename,
              start_row = current_row + 1, doWeights = FALSE)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          current_column <- current_column + venn_columns
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     x = "Venn of p-value down genes, lfc < -1.",
                                     startRow = current_row, startCol = current_column)
          down_plot <- venn_sig_list[["down_venn"]]
          try_result <- xlsx_plot_png(
              down_plot, wb = wb, sheet = sheetname, width = plot_dim / 2, height = plot_dim / 2,
              start_col = current_column, plotname = "downvenn", savedir = excel_basename,
              start_row = current_row + 1, doWeights = FALSE)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          current_column <- current_column + venn_columns
          if (class(siggene_lst)[1] != "try-error") {
            xl_result <- openxlsx::writeData(
                                       wb = wb, sheet = sheetname,
                                       x = "Significant genes by fc going down.",
                                       startRow = current_row, startCol = current_column)
            try_result <- xlsx_plot_png(
                siggene_lst[["down"]], wb = wb, sheet = sheetname,
                width = plot_dim, height = (plot_dim / 2),
                start_col = current_column, plotname = "siggenesup",
                savedir = excel_basename, start_row = current_row + 1, doWeights = FALSE)
            if (! "try-error" %in% class(try_result)) {
              image_files <- c(image_files, try_result[["filename"]])
            }
            current_column <- current_column + plot_columns
            xl_result <- openxlsx::writeData(
                                       wb = wb, sheet = sheetname,
                                       x = "Significant genes by p going down.",
                                       startRow = current_row, startCol = current_column)
            try_result <- xlsx_plot_png(
                siggene_lst[["pdown"]], wb = wb, sheet = sheetname,
                width = plot_dim, height = (plot_dim / 2),
                start_col = current_column, plotname = "siggenespdown",
                savedir = excel_basename, start_row = current_row + 1, doWeights = FALSE)
            if (! "try-error" %in% class(try_result)) {
              image_files <- c(image_files, try_result[["filename"]])
            }
          }
          current_row <- current_row + venn_rows
        } ## End checking on venns
        ## Now add the coefficients, ma, and volcanoes below the venns.
        ## Text on row 18, plots from 19-49 (30 rows)
        for (t in 1:length(de_types)) {
          type <- de_types[t]
          sc <- paste0(type, "_scatter_plots")
          ma <- paste0(type, "_ma_plots")
          vo <- paste0(type, "_vol_plots")
          pp <- paste0(type, "_p_plots")
          short <- substr(type, 0, 2)
          cap <- R.utils::capitalize(type)
          plt <- extracted[["plots"]][[sheetname]][[sc]]
          ma_plt <- extracted[["plots"]][[sheetname]][[ma]]
          vol_plt <- extracted[["plots"]][[sheetname]][[vo]]
          p_plt <- extracted[["plots"]][[sheetname]][[pp]]
          current_row <- current_row + 2
          current_column <- xls_result[["end_col"]] + 2
          if (class(plt)[1] != "try-error" & !is.null(plt)) {
            printme <- as.character(
                glue::glue("{cap} expression coefficients for {tnames[x]}; R^2: \\
                          {signif(x=plt[['lm_rsq']], digits=3)}; equation: \\
                          {ymxb_print(plt[['lm_model']])}"))
            xl_result <- openxlsx::writeData(
                                       wb = wb, sheet = sheetname, x = printme,
                                       startRow = current_row, startCol = current_column)
            plotname <- paste0(short, "scatter")
            try_result <- xlsx_plot_png(
                plt[["scatter"]], wb = wb, sheet = sheetname,
                width = plot_dim, height = plot_dim, start_col = current_column,
                plotname = plotname, savedir = excel_basename, start_row = current_row + 1)
            if (! "try-error" %in% class(try_result)) {
              image_files <- c(image_files, try_result[["filename"]])
            }
            current_column <- current_column + plot_columns
            xl_result <- openxlsx::writeData(
                                       wb = wb, sheet = sheetname, x = paste0(type, " MA plot"),
                                       startRow = current_row, startCol = current_column)
            plotname <- paste0(short, "ma")
            try_ma_result <- xlsx_plot_png(
                ma_plt[["plot"]], wb = wb, sheet = sheetname, width = plot_dim,
                height = plot_dim, start_col = current_column, plotname = plotname,
                savedir = excel_basename, start_row = current_row + 1)
            if (! "try-error" %in% class(try_ma_result)) {
              image_files <- c(image_files, try_ma_result[["filename"]])
            }
            current_column <- current_column + plot_columns
            xl_result <- openxlsx::writeData(
                                       wb = wb, sheet = sheetname, x = paste0(type, " volcano plot"),
                                       startRow = current_row, startCol = current_column)
            plotname <- paste0(short, "vol")
            try_vol_result <- xlsx_plot_png(
                vol_plt[["plot"]], wb = wb, sheet = sheetname, width = plot_dim,
                height = plot_dim, start_col = current_column, pltname = plotname,
                savedir = excel_basename, start_row = current_row + 1)
            if (! "try-error" %in% class(try_vol_result)) {
              image_files <- c(image_files, try_vol_result[["filename"]])
            }
            current_column <- current_column + plot_columns
            xl_result <- openxlsx::writeData(
                                       wb = wb, sheet = sheetname, x = paste0(type, " p-value plot"),
                                       startRow = current_row, startCol = current_column)
            plotname <- paste0(short, "p")
            try_p_result <- xlsx_plot_png(
                p_plt, wb = wb, sheet = sheetname, width = plot_dim,
                height = plot_dim, start_col = current_column, pltname = plotname,
                savedir = excel_basename, start_row = current_row + 1)
            if (! "try-error" %in% class(try_p_result)) {
              image_files <- c(image_files, try_p_result[["filename"]])
            }
            current_row <- current_row + plot_rows
          }
        } ## End adding limma, deseq, and edger plots.
      } ## End checking whether to add plots
    }  ## End for loop iterating over every kept table.

    ## Now add some summary data and some plots comparing the tools.
    xls_result <- write_combined_summary(wb, excel_basename, apr, extracted, compare_plots,
                                         lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff)
    image_files <- c(image_files, xls_result[["image_files"]])

    ## Finish up.
    if (!is.null(apr[["original_pvalues"]])) {
      mesg("Appending a data frame of the original pvalues before sva messed with them.")
      xls_result <- write_xlsx(
          wb, data = apr[["original_pvalues"]], sheet = "original_pvalues",
          title = "Original pvalues for all contrasts before sva adjustment.",
          start_row = 1, rownames = rownames)
    }

    design_result <- write_sample_design(wb, apr)

    mesg("Performing save of ", excel, ".")
    save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
    if (class(save_result)[1] == "try-error") {
      warning("Saving xlsx failed.")
    }
  } ## End if !is.null(excel)

  ## We have finished!  Dump the important stuff into a return list.
  ret <- NULL
  if (is.null(retlist)) {
    ret <- list(
        "data" = extracted[["data"]],
        "table_names" = extracted[["table_names"]],
        "input" = apr,
        "plots" = extracted[["plots"]],
        "comp_plot" = comp,
        "venns" = venns,
        "keepers" = keepers,
        ## Kept is currently broken.
        "kept" = extracted[["kept"]],
        "de_summary" = extracted[["summaries"]])
  } else {
    ret <- retlist
  }
  class(ret) <- c("combined_de", "list")

  if (!is.null(rda)) {
    saved <- save(list = "ret", file = rda)
  }
  ## Cleanup the saved image files.
  for (img in image_files) {
    removed <- try(suppressWarnings(file.remove(img)), silent = TRUE)
  }
  return(ret)
}

#' Gather data required to make MA/Volcano plots for pairwise comparisons.
#'
#' It should be possible to significantly simplify the arguments passed to this
#' function, but I have thus far focused only on getting it to work with the
#' newly split-apart combine_de_tables() functions.
#'
#' @param name Name of the table to plot.
#' @param combined Modified pairwise result, containing the various DE methods.
#' @param denominator Name of the denominator coefficient.
#' @param numerator Name of the numerator coefficient.
#' @param plot_inputs The individual outputs from limma etc.
#' @param include_basic Add basic data?
#' @param include_deseq Add deseq data?
#' @param include_edger Add edger data?
#' @param include_limma Add limma data?
#' @param include_ebseq Add ebseq data?
#' @param loess Add a loess estimation?
#' @param logfc For Volcano/MA plot lines.
#' @param p For Volcano/MA plot lines.
#' @param do_inverse Flip the numerator/denominator?
#' @param found_table The table name actually used.
combine_extracted_plots <- function(name, combined, denominator, numerator, plot_inputs,
                                    include_basic = TRUE, include_deseq = TRUE,
                                    include_edger = TRUE, include_limma = TRUE,
                                    include_ebseq = FALSE, loess = FALSE, logfc = 1, p = 0.05,
                                    do_inverse = FALSE, found_table = NULL, p_type = "all") {
  combined_data <- combined[["data"]]
  plots <- list()
  if (isTRUE(include_deseq)) {
    plots[["deseq_scatter_plots"]] <- list()
    plots[["deseq_ma_plots"]] <- list()
    plots[["deseq_vol_plots"]] <- list()
    plots[["deseq_p_plots"]] <- list()
  }
  if (isTRUE(include_edger)) {
    plots[["edger_scatter_plots"]] <- list()
    plots[["edger_ma_plots"]] <- list()
    plots[["edger_vol_plots"]] <- list()
    plots[["edger_p_plots"]] <- list()
  }
  if (isTRUE(include_limma)) {
    plots[["limma_scatter_plots"]] <- list()
    plots[["limma_ma_plots"]] <- list()
    plots[["limma_vol_plots"]] <- list()
    plots[["limma_p_plots"]] <- list()
  }
  types <- c()
  if (isTRUE(include_limma)) {
    types <- c("limma", types)
  }
  if (isTRUE(include_edger)) {
    types <- c("edger", types)
  }
  if (isTRUE(include_deseq)) {
    types <- c("deseq", types)
  }
  for (t in 1:length(types)) {
    type <- types[t]
    sc_name <- paste0(type, "_scatter_plots")
    ma_name <- paste0(type, "_ma_plots")
    vol_name <- paste0(type, "_vol_plots")
    p_name <- paste0(type, "_p_plots")
    if (is.null(sc_name)) {
      mesg("Skipping scatter plot for ", type, ".")
    } else {
      x_y <- strsplit(x = found_table, split = "_vs_")[[1]]
      coef_scatter <- sm(try(extract_coefficient_scatter(
          plot_inputs[[type]], type = type, invert = do_inverse,
          x = x_y[1], y = x_y[2])))
      plots[[sc_name]] <- coef_scatter
    }
    if (is.null(ma_name)) {
      mesg("Skipping volcano/MA plot for ", type, ".")
    } else {
      ma_vol <- try(extract_de_plots(
          plot_inputs[[type]], type = type, invert = do_inverse,
          logfc = logfc, p = p, table = found_table), silent = TRUE)
      if ("try-error" %in% class(ma_vol)) {
        plots[[ma_name]] <- NULL
        plots[[vol_name]] <- NULL
      } else {
        plots[[ma_name]] <- ma_vol[["ma"]]
        plots[[vol_name]] <- ma_vol[["volcano"]]
      }
    }
    if (is.null(p_name)) {
      mesg("Skipping p-value plot for ", t, ".")
    } else {
      ## If one sets the padj_type, then one will need to pull the correct columns
      ## from the data at this point.  In my vignette, I set padj_type to 'BH'
      ## and as a result I have a series of columns: 'limma_adj_bh' etc.
      ## Therefore we need to get that information to this function call.
      pval_plot <- plot_de_pvals(combined[["data"]], type = type,
                                 p_type = p_type)
      plots[[p_name]] <- pval_plot
    }
  }
  return(plots)
}

#' Given a limma, edger, and deseq table, combine them into one.
#'
#' This combines the outputs from the various differential expression
#' tools and formalizes some column names to make them a little more
#' consistent.
#'
#' @param li Limma output table.
#' @param ed Edger output table.
#' @param eb EBSeq output table
#' @param de DESeq2 output table.
#' @param ba Basic output table.
#' @param table_name Name of the table to merge.
#' @param final_table_names Vector of the final table names.
#' @param annot_df Add some annotation information?
#' @param do_inverse Invert the fold changes?
#' @param adjp Use adjusted p-values?
#' @param padj_type Add this consistent p-adjustment.
#' @param include_deseq Include tables from deseq?
#' @param include_edger Include tables from edger?
#' @param include_ebseq Include tables from ebseq?
#' @param include_limma Include tables from limma?
#' @param include_basic Include the basic table?
#' @param lfc_cutoff Preferred logfoldchange cutoff.
#' @param p_cutoff Preferred pvalue cutoff.
#' @param excludes Set of genes to exclude from the output.
#' @param sheet_count What sheet is being written?
#' @return List containing a) Dataframe containing the merged
#'  limma/edger/deseq/basic tables, and b) A summary of how many
#'  genes were observed as up/down by output table.
#' @seealso [data.table] [hpgl_padjust()] [extract_keepers_all()] [extract_keepers_lst()]
#'  [extract_keepers_single()]
combine_single_de_table <- function(li = NULL, ed = NULL, eb = NULL, de = NULL, ba = NULL,
                                    table_name = "", final_table_names = c(),
                                    annot_df = NULL, do_inverse = FALSE,
                                    adjp = TRUE, padj_type = "fdr",
                                    include_deseq = TRUE, include_edger = TRUE,
                                    include_ebseq = TRUE, include_limma = TRUE,
                                    include_basic = TRUE, lfc_cutoff = 1,
                                    p_cutoff = 0.05, excludes = NULL, sheet_count = 0) {
  if (padj_type[1] != "ihw" & !padj_type %in% p.adjust.methods) {
    warning("The p adjustment ", padj_type, " is not in the set of p.adjust.methods.
Defaulting to fdr.")
    padj_type <- "fdr"
  }

  num_den_names <- strsplit(x = table_name, split = "_vs_")[[1]]
  extra <- FALSE
  if (length(num_den_names) == 1) {
    ## Then this should be coming from extra_contrasts
    num_name = NULL
    den_name = NULL
    inverse_name = NULL
    extra <- TRUE
  } else {
    num_name <- num_den_names[[1]]
    den_name <- num_den_names[[2]]
    inverse_name <- glue::glue("{den_name}_vs_{num_name}")
  }

  lidf <- data.frame("limma_logfc" = 0, "limma_ave" = 0, "limma_t" = 0,
                     "limma_p" = 0, "limma_adjp" = 0, "limma_b" = 0)
  rownames(lidf) <- "dropme"
  dedf <- data.frame("deseq_basemean" = 0, "deseq_logfc" = 0, "deseq_lfcse" = 0,
                     "deseq_stat" = 0, "deseq_p" = 0, "deseq_adjp" = 0)
  rownames(dedf) <- "dropme"
  eddf <- data.frame("edger_logfc" = 0, "edger_logcpm" = 0, "edger_lr" = 0,
                     "edger_p" = 0, "edger_adjp" = 0)
  rownames(eddf) <- "dropme"
  ebdf <- data.frame("ebseq_fc" = 0, "ebseq_logfc" = 0, "ebseq_c1mean" = 0,
                     "ebseq_c2mean" = 0, "ebseq_mean" = 0, "ebseq_var" = 0,
                     "ebseq_postfc" = 0, "ebseq_ppee" = 0, "ebseq_ppde" = 0,
                     "ebseq_adjp" = 0)
  rownames(ebdf) <- "dropme"
  badf <- data.frame("numerator_median" = 0, "denominator_median" = 0, "numerator_var" = 0,
                     "denominator_var" = 0, "logFC" = 0, "t" = 0, "p" = 0, "adjp" = 0)
  rownames(badf) <- "dropme"

  ## Check that the limma result is valid.
  if (is.null(li) | class(li)[1] == "try-error") {
    mesg("The limma table is null.")
    li <- NULL
    include_limma <- FALSE
  } else if (isTRUE(extra)) {
    test_lidf <- li[["all_tables"]][[table_name]]
    if (!is.null(test_lidf)) {
      lidf <- test_lidf
    }
  } else {
    test_lidf <- li[["all_tables"]][[table_name]]
    if (is.null(test_lidf)) {
      test_lidf <- li[["all_tables"]][[inverse_name]]
      mesg("Used the inverse table, might need to -1 the logFC.")
      if (is.null(test_lidf)) {
        warning("The limma table seems to be missing.")
      } else {
        lidf <- test_lidf
      }
    } else {
      lidf <- test_lidf
    }
  }

  ## Check that the deseq result is valid.
  if (is.null(de) | class(de)[1] == "try-error") {
    mesg("The deseq table is null.")
    de <- NULL
    include_deseq <- FALSE
  } else if (isTRUE(extra)) {
    test_dedf <- de[["all_tables"]][[table_name]]
    if (!is.null(test_dedf)) {
      dedf <- test_dedf
    }
  } else {
    test_dedf <- de[["all_tables"]][[table_name]]
    if (is.null(test_dedf)) {
      test_dedf <- de[["all_tables"]][[inverse_name]]
      mesg("Used the inverse table, might need to -1 the logFC and stat.")
      if (is.null(test_dedf)) {
        warning("The deseq table seems to be missing.")
      } else {
        dedf <- test_dedf
      }
    } else {
      dedf <- test_dedf
    }
  }

  ## Check that the edger result is valid.
  if (is.null(ed) | class(ed)[1] == "try-error") {
    mesg("The edger table is null.")
    ed <- NULL
    include_edger <- FALSE
  } else if (isTRUE(extra)) {
    test_eddf <- ed[["all_tables"]][[table_name]]
    if (!is.null(test_eddf)) {
      eddf <- test_eddf
    }
  } else {
    test_eddf <- ed[["all_tables"]][[table_name]]
    if (is.null(test_eddf)) {
      test_eddf <- ed[["all_tables"]][[inverse_name]]
      mesg("Used the inverse table, might need to -1 the logFC.")
      if (is.null(test_eddf)) {
        warning("The edger table seems to be missing.")
      } else {
        eddf <- test_eddf
      }
    } else {
      eddf <- test_eddf
    }
  }

  ## Check that the ebseq result is valid.
  if (is.null(eb) | class(eb)[1] == "try-error") {
    mesg("The ebseq table is null.")
    eb <- NULL
    include_ebseq <- FALSE
  } else if (isTRUE(extra)) {
    test_ebdf <- eb[["all_tables"]][[table_name]]
    if (!is.null(test_ebdf)) {
      ebdf <- test_ebdf
    }
  } else {
    test_ebdf <- eb[["all_tables"]][[table_name]]
    if (is.null(test_ebdf)) {
      test_ebdf <- eb[["all_tables"]][[inverse_name]]
      mesg("Used the inverse table, might need to -1 the logFC.")
      if (is.null(test_ebdf)) {
        warning("The ebseq table seems to be missing.")
      } else {
        ebdf <- test_ebdf
      }
    } else {
      ebdf <- test_ebdf
    }
  }

  ## And finally, check that my stupid basic result is valid.
  if (is.null(ba) | class(ba)[1] == "try-error") {
    mesg("The basic table is null.")
    ba <- NULL
    include_basic <- FALSE
  } else if (isTRUE(extra)) {
    test_badf <- ba[["all_tables"]][[table_name]]
    if (!is.null(test_badf)) {
      badf <- test_badf
    }
  } else {
    test_badf <- ba[["all_tables"]][[table_name]]
    if (is.null(test_badf)) {
      test_badf <- ba[["all_tables"]][[inverse_name]]
      mesg("Used the inverse table, might need to -1 the logFC.")
      if (is.null(test_badf)) {
        warning("The basic table seems to be missing.")
      } else {
        badf <- test_badf
      }
    } else {
      badf <- test_badf
    }
  }

  colnames(lidf) <- c("limma_logfc", "limma_ave", "limma_t", "limma_p",
                      "limma_adjp", "limma_b")
  li_stats <- lidf[, c("limma_ave", "limma_t", "limma_b", "limma_p")]
  li_lfc_adjp <- lidf[, c("limma_logfc", "limma_adjp")]

  colnames(dedf) <- c("deseq_basemean", "deseq_logfc", "deseq_lfcse",
                      "deseq_stat", "deseq_p", "deseq_adjp")
  de_stats <- dedf[, c("deseq_basemean", "deseq_lfcse", "deseq_stat", "deseq_p")]
  de_lfc_adjp <- dedf[, c("deseq_logfc", "deseq_adjp")]

  colnames(eddf) <- c("edger_logfc", "edger_logcpm", "edger_lr", "edger_p", "edger_adjp")
  ed_stats <- eddf[, c("edger_logcpm", "edger_lr", "edger_p")]
  ed_lfc_adjp <- eddf[, c("edger_logfc", "edger_adjp")]
  colnames(ebdf) <- c("ebseq_fc", "ebseq_logfc", "ebseq_c1mean",
                      "ebseq_c2mean", "ebseq_mean", "ebseq_var",
                      "ebseq_postfc", "ebseq_ppee", "ebseq_ppde",
                      "ebseq_adjp")

  ## I recently changed basic to optionally do means or medians.  I need to take that into
  ## account when working with these tables.  For the moment, I think I will simply rename
  ## the column to _median to avoid confusion.
  mean_idx <- "numerator_mean" == colnames(badf)
  if (sum(mean_idx) > 0) {
    colnames(badf)[mean_idx] <- "numerator_median"
    mean_idx <- "denominator_mean" == colnames(badf)
    colnames(badf)[mean_idx] <- "denominator_median"
  }
  ba_stats <- badf[, c("numerator_median", "denominator_median", "numerator_var",
                       "denominator_var", "logFC", "t", "p", "adjp")]
  colnames(ba_stats) <- c("basic_nummed", "basic_denmed", "basic_numvar", "basic_denvar",
                          "basic_logfc", "basic_t", "basic_p", "basic_adjp")

  datalst <- list()
  statslst <- list()
  if (isTRUE(include_basic)) {
    statslst[["basic"]] <- data.table::as.data.table(ba_stats)
    statslst[["basic"]][["rownames"]] <- rownames(ba_stats)
  }
  if (isTRUE(include_deseq)) {
    datalst[["deseq"]] <- data.table::as.data.table(de_lfc_adjp)
    datalst[["deseq"]][["rownames"]] <- rownames(de_lfc_adjp)
    statslst[["deseq"]] <- data.table::as.data.table(de_stats)
    statslst[["deseq"]][["rownames"]] <- rownames(de_stats)
  }
  if (isTRUE(include_ebseq)) {
    statslst[["ebseq"]] <- data.table::as.data.table(ebdf)
    statslst[["ebseq"]][["rownames"]] <- rownames(ebdf)
  }
  if (isTRUE(include_edger)) {
    datalst[["edger"]] <- data.table::as.data.table(ed_lfc_adjp)
    datalst[["edger"]][["rownames"]] <- rownames(ed_lfc_adjp)
    statslst[["edger"]] <- data.table::as.data.table(ed_stats)
    statslst[["edger"]][["rownames"]] <- rownames(ed_stats)
  }
  if (isTRUE(include_limma)) {
    datalst[["limma"]] <- data.table::as.data.table(li_lfc_adjp)
    datalst[["limma"]][["rownames"]] <- rownames(li_lfc_adjp)
    statslst[["limma"]] <- data.table::as.data.table(li_stats)
    statslst[["limma"]][["rownames"]] <- rownames(li_stats)
  }

  ## Make the initial data structure
  wanted <- names(datalst)[[1]]
  num_stats <- length(statslst)
  num_data <- length(datalst)
  if (num_data == 1) {
    if (num_stats == 1) {
      ## Then this is a chunk of data and associated stats.
      comb <- merge(datalst[[1]], statslst[[1]], by = "rownames", all = TRUE)
    } else {
      ## Then there must only be a chunk of data.
      comb <- datalst[[1]]
    }
  } else {
    ## There is more than one set of data to merge.
    comb <- datalst[[1]]
    for (i in 2:length(datalst)) {
      comb <- merge(comb, datalst[[i]], by = "rownames", all = TRUE)
    }
    if (length(statslst) > 0) {
      for (j in 1:length(statslst)) {
        comb <- merge(comb, statslst[[j]], by = "rownames", all = TRUE)
      }
    }
  }
  ## Doing the merge in the way above will lead to a single row which is essentially blank
  ## It is probably the first row.

  ## The next lines are intended to drop that blank row.
  comb <- as.data.frame(comb)
  rownames(comb) <- comb[["rownames"]]
  dropme <- rownames(comb) == "dropme"
  comb <- comb[!dropme, ]

  keepers <- colnames(comb) != "rownames"
  comb <- comb[, keepers, drop = FALSE]
  comb[is.na(comb)] <- 0
  if (isTRUE(do_inverse)) {
    if (isTRUE(include_basic)) {
      comb[["basic_logfc"]] <- comb[["basic_logfc"]] * -1.0
    }
    if (isTRUE(include_limma)) {
      comb[["limma_logfc"]] <- comb[["limma_logfc"]] * -1.0
    }
    if (isTRUE(include_deseq)) {
      comb[["deseq_logfc"]] <- comb[["deseq_logfc"]] * -1.0
      comb[["deseq_stat"]] <- comb[["deseq_stat"]] * -1.0
    }
    if (isTRUE(include_edger)) {
      comb[["edger_logfc"]] <- comb[["edger_logfc"]] * -1.0
    }
    if (isTRUE(include_ebseq)) {
      comb[["ebseq_logfc"]] <- comb[["ebseq_logfc"]] * -1.0
    }
  }

  ## Add one final p-adjustment to ensure a consistent and user defined value.
  if (!is.null(comb[["limma_p"]])) {
    colname <- glue::glue("limma_adjp_{padj_type}")
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "limma_p", mean_column = "limma_ave",
                                    method = padj_type, significance = 0.05)
    comb[[colname]] <- format(x = comb[[colname]], digits = 4, scientific = TRUE, trim = TRUE)
  }
  if (!is.null(comb[["deseq_p"]])) {
    colname <- glue::glue("deseq_adjp_{padj_type}")
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "deseq_p", mean_column = "deseq_basemean",
                                    method = padj_type, significance = 0.05)
    comb[[colname]] <- format(x=comb[[colname]], digits = 4, scientific = TRUE, trim = TRUE)
  }
  if (!is.null(comb[["edger_p"]])) {
    colname <- glue::glue("edger_adjp_{padj_type}")
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "edger_p", mean_column = "edger_logcpm",
                                    method = padj_type, significance = 0.05)
    comb[[colname]] <- format(x = comb[[colname]], digits = 4, scientific = TRUE, trim = TRUE)
  }
  if (!is.null(comb[["ebseq_ppde"]])) {
    colname <- glue::glue("ebseq_adjp_{padj_type}")
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "ebseq_ppde", mean_column = "ebseq_mean",
                                    method = padj_type, significance = 0.05)
    comb[[colname]] <- format(x = comb[[colname]], digits = 4, scientific = TRUE, trim = TRUE)
  }
  if (!is.null(comb[["basic_p"]])) {
    colname <- glue::glue("basic_adjp_{padj_type}")
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "basic_p", mean_column = "basic_nummed",
                                    method = padj_type, significance = 0.05)
    comb[[colname]] <- format(x = comb[[colname]], digits = 4, scientific = TRUE, trim = TRUE)
  }

  ## I made an odd choice in a moment to normalize.quantiles the combined fold changes
  ## This should be reevaluated
  temp_fc <- data.frame()
  if (isTRUE(include_limma) & isTRUE(include_deseq) & isTRUE(include_edger)) {
    temp_fc <- cbind(as.numeric(comb[["limma_logfc"]]),
                     as.numeric(comb[["edger_logfc"]]),
                     as.numeric(comb[["deseq_logfc"]]))
    temp_fc <- preprocessCore::normalize.quantiles(as.matrix(temp_fc))
    comb[["lfc_meta"]] <- rowMeans(temp_fc, na.rm = TRUE)
    comb[["lfc_var"]] <- genefilter::rowVars(temp_fc, na.rm = TRUE)
    comb[["lfc_varbymed"]] <- comb[["lfc_var"]] / comb[["lfc_meta"]]
    temp_p <- cbind(as.numeric(comb[["limma_p"]]),
                    as.numeric(comb[["edger_p"]]),
                    as.numeric(comb[["deseq_p"]]))
    comb[["p_meta"]] <- rowMeans(temp_p, na.rm = TRUE)
    comb[["p_var"]] <- genefilter::rowVars(temp_p, na.rm = TRUE)
    comb[["lfc_meta"]] <- signif(x=comb[["lfc_meta"]], digits = 4)
    comb[["lfc_var"]] <- format(x=comb[["lfc_var"]], digits = 4, scientific = TRUE, trim = TRUE)
    comb[["lfc_varbymed"]] <- format(x=comb[["lfc_varbymed"]], digits = 4, scientific = TRUE, trim = TRUE)
    comb[["p_var"]] <- format(x=comb[["p_var"]], digits = 4, scientific = TRUE, trim = TRUE)
    comb[["p_meta"]] <- format(x=comb[["p_meta"]], digits = 4, scientific = TRUE, trim = TRUE)
  }
  if (!is.null(annot_df)) {
    colnames(annot_df) <- gsub(pattern = "[[:punct:]]",
                               replacement = "",
                               x = colnames(annot_df))
    comb <- merge(annot_df, comb, by = "row.names", all.y = TRUE)
    rownames(comb) <- comb[["Row.names"]]
    comb <- comb[, -1, drop = FALSE]
    colnames(comb) <- make.names(tolower(colnames(comb)), unique = TRUE)
  }

  ## Exclude rows based on a list of unwanted columns/strings
  if (!is.null(excludes)) {
    for (colnum in 1:length(excludes)) {
      col <- names(excludes)[colnum]
      for (exclude_num in 1:length(excludes[[col]])) {
        exclude <- excludes[[col]][exclude_num]
        remove_column <- comb[[col]]
        remove_idx <- grep(pattern = exclude, x = remove_column, perl = TRUE, invert = TRUE)
        removed_num <- sum(as.numeric(remove_idx))
        message("Removed ", removed_num, " genes using ",
                exclude, " as a string against column ", col, ".")
        comb <- comb[remove_idx, ]
      }  ## End iterating through every string to exclude
    }  ## End iterating through every element of the exclude list
  }

  up_fc <- lfc_cutoff
  down_fc <- -1.0 * lfc_cutoff
  summary_table_name <- table_name
  if (isTRUE(do_inverse)) {
    summary_table_name <- glue::glue("{summary_table_name}-inverted")
  }
  limma_p_column <- "limma_adjp"
  deseq_p_column <- "deseq_adjp"
  edger_p_column <- "edger_adjp"
  if (!isTRUE(adjp)) {
    limma_p_column <- "limma_p"
    deseq_p_column <- "deseq_p"
    edger_p_column <- "edger_p"
  }
  summary_lst <- list(
      "table" = summary_table_name,
      "total" = nrow(comb),
      "limma_up" = sum(comb[["limma_logfc"]] >= up_fc),
      "limma_sigup" = sum(
          comb[["limma_logfc"]] >= up_fc & as.numeric(comb[[limma_p_column]]) <= p_cutoff),
      "deseq_up" = sum(comb[["deseq_logfc"]] >= up_fc),
      "deseq_sigup" = sum(
          comb[["deseq_logfc"]] >= up_fc & as.numeric(comb[[deseq_p_column]]) <= p_cutoff),
      "edger_up" = sum(comb[["edger_logfc"]] >= up_fc),
      "edger_sigup" = sum(
          comb[["edger_logfc"]] >= up_fc & as.numeric(comb[[edger_p_column]]) <= p_cutoff),
      "basic_up" = sum(comb[["basic_logfc"]] >= up_fc),
      "basic_sigup" = sum(
          comb[["basic_logfc"]] >= up_fc & as.numeric(comb[["basic_p"]]) <= p_cutoff),
      "limma_down" = sum(comb[["limma_logfc"]] <= down_fc),
      "limma_sigdown" = sum(
          comb[["limma_logfc"]] <= down_fc & as.numeric(comb[[limma_p_column]]) <= p_cutoff),
      "deseq_down" = sum(comb[["deseq_logfc"]] <= down_fc),
      "deseq_sigdown" = sum(
          comb[["deseq_logfc"]] <= down_fc & as.numeric(comb[[deseq_p_column]]) <= p_cutoff),
      "edger_down" = sum(comb[["edger_logfc"]] <= down_fc),
      "edger_sigdown" = sum(
          comb[["edger_logfc"]] <= down_fc & as.numeric(comb[[edger_p_column]]) <= p_cutoff),
      "basic_down" = sum(comb[["basic_logfc"]] <= down_fc),
      "basic_sigdown" = sum(
          comb[["basic_logfc"]] <= down_fc & as.numeric(comb[["basic_p"]]) <= p_cutoff),
      "meta_up" = sum(comb[["fc_meta"]] >= up_fc),
      "meta_sigup" = sum(
          comb[["lfc_meta"]] >= up_fc & as.numeric(comb[["p_meta"]]) <= p_cutoff),
      "meta_down" = sum(comb[["lfc_meta"]] <= down_fc),
      "meta_sigdown" = sum(
          comb[["lfc_meta"]] <= down_fc & as.numeric(comb[["p_meta"]]) <= p_cutoff)
  )

  ret <- list(
      "data" = comb,
      "summary" = summary_lst)
  class(ret) <- c("combined_table", "list")
  return(ret)
}

#' When no set of 'keeper' contrasts is specified, grab them all.
#'
#' This has a couple of cousin functions, extract_keepers_list and _single.
#' These handle extracting one or more contrasts out of the various tables
#' produced by all_pairwise().
#'
#' @param apr Result from all_pairwise().
#' @param extracted Table of extracted data.
#' @param keepers In this case, one may assume either NULL or 'all'.
#' @param table_names The set of tables produced by all_pairwise().
#' @param all_coefficients The set of all experimental conditions in the
#'  experimental metadata.
#' @param limma The limma data from all_pairwise().
#' @param edger The edger data from all_pairwise().
#' @param ebseq The ebseq data from all_pairwise().
#' @param deseq The deseq data from all_pairwise().
#' @param basic The basic data from all_pairwise().
#' @param adjp Pull out the adjusted p-values from the data?
#' @param annot_df What annotations should be added to the table?
#' @param include_deseq Whether or not to include the deseq data.
#' @param include_edger Whether or not to include the edger data.
#' @param include_ebseq Whether or not to include the ebseq data.
#' @param include_limma Whether or not to include the limma data.
#' @param include_basic Whether or not to include the basic data.
#' @param excludes Set of genes to exclude.
#' @param padj_type Choose a specific p adjustment.
#' @param loess Include a loess estimator in the plots?
#' @param lfc_cutoff Passed for plotting volcano/MA plots.
#' @param p_cutoff Passed for volcano/MA plots.
extract_keepers_all <- function(apr, extracted, keepers, table_names,
                                all_coefficients,
                                limma, edger, ebseq, deseq, basic,
                                adjp, annot_df,
                                include_deseq, include_edger,
                                include_ebseq, include_limma,
                                include_basic, excludes, padj_type, loess = FALSE,
                                lfc_cutoff = 1, p_cutoff = 0.05) {
  names_length <- length(table_names)
  kept_tables <- list()
  numerators <- denominators <- c()
  for (a in 1:names_length) {
    name <- table_names[a]
    kept_tables[a] <- name
    mesg("Working on table ", a, "/", names_length, ": ", name)
    splitted <- strsplit(x=name, split="_vs_")
    denominator <- splitted[[1]][2]
    numerator <- splitted[[1]][1]
    numerators[a] <- numerator
    denominators[a] <- denominator
    combined <- combine_single_de_table(
        li = limma, ed = edger, eb = ebseq, de = deseq, ba = basic,
        include_basic = include_basic, include_deseq = include_deseq,
        include_edger = include_edger, include_ebseq = include_ebseq,
        include_limma = include_limma,
        table_name = name, adjp = adjp, annot_df = annot_df,
        excludes = excludes, padj_type = padj_type,
        lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff)
    extracted[["data"]][[name]] <- combined[["data"]]
    extracted[["table_names"]][[name]] <- combined[["summary"]][["table"]]
    extracted[["kept"]] <- kept_tables
    extracted[["keepers"]] <- keepers
    plot_inputs <- list(
        "limma" = limma,
        "deseq" = deseq,
        "edger" = edger,
        "ebseq" = ebseq,
        "basic" = basic)
    extracted[["plots"]][[name]] <- combine_extracted_plots(
        name, combined, denominator, numerator, plot_inputs,
        include_basic = include_basic, include_deseq = include_deseq,
        include_edger = include_edger, include_limma = include_limma,
        include_ebseq = include_ebseq, loess = loess, logfc = lfc_cutoff, p = p_cutoff,
        found_table = name, p_type = padj_type)
    extracted[["summaries"]] <- rbind(extracted[["summaries"]],
                                      as.data.frame(combined[["summary"]]))
    extracted[["numerators"]] <- numerators
    extracted[["denominators"]] <- denominators
  } ## End for list
  return(extracted)
}

#' When a list of 'keeper' contrasts is specified, extract it from the data.
#'
#' This is the most interesting of the extract_keeper functions.  It must check
#' that the numerators and denominators match the desired contrast and flip the
#' signs in the logFCs when appropriate.
#'
#' @param extracted Tables extracted from the all_pairwise data.
#' @param keepers In this case, one may assume either NULL or 'all'.
#' @param table_names The set of tables produced by all_pairwise().
#' @param all_coefficients The set of all experimental conditions in the
#'  experimental metadata.
#' @param limma The limma data from all_pairwise().
#' @param edger The edger data from all_pairwise().
#' @param ebseq The ebseq data from all_pairwise().
#' @param deseq The deseq data from all_pairwise().
#' @param basic The basic data from all_pairwise().
#' @param adjp Pull out the adjusted p-values from the data?
#' @param annot_df What annotations should be added to the table?
#' @param include_deseq Whether or not to include the deseq data.
#' @param include_edger Whether or not to include the edger data.
#' @param include_ebseq Whether or not to include the ebseq data.
#' @param include_limma Whether or not to include the limma data.
#' @param include_basic Whether or not to include the basic data.
#' @param excludes Set of genes to exclude.
#' @param padj_type Choose a specific p adjustment.
#' @param loess Add a loess to plots?
#' @param lfc_cutoff Passed for volcano/MA plots.
#' @param p_cutoff Passed for volcano/MA plots.
extract_keepers_lst <- function(extracted, keepers, table_names,
                                all_coefficients,
                                limma, edger, ebseq, deseq, basic,
                                adjp, annot_df,
                                include_deseq, include_edger,
                                include_ebseq, include_limma,
                                include_basic, excludes, padj_type,
                                loess = FALSE, lfc_cutoff = 1, p_cutoff = 0.05) {
  ## First check that your set of kepers is in the data
  all_keepers <- as.character(unlist(keepers))
  keeper_names <- names(keepers)
  found_keepers <- sum(all_keepers %in% all_coefficients)

  ## Add a little logic in case we have weirdo contrasts, in which case
  ## the name of the table should remain constant and the numerator/denominator will be ignored
  ## because they should have been defined when setting up the weirdo contrast.
  found_names_edger_idx <- names(edger[["all_tables"]]) %in% keeper_names
  found_names_deseq_idx <- names(deseq[["all_tables"]]) %in% keeper_names
  found_names_limma_idx <- names(limma[["all_tables"]]) %in% keeper_names
  extra_names <- c()
  if (sum(found_names_edger_idx) > 0) {
    extra_names <- names(edger[["all_tables"]])[found_names_edger_idx]
  }
  if (sum(found_names_deseq_idx) > 0) {
    extra_names <- c(extra_names, names(deseq[["all_tables"]])[found_names_deseq_idx])
  }
  if (sum(found_names_deseq_idx) > 0) {
    extra_names <- c(extra_names, names(deseq[["all_tables"]])[found_names_deseq_idx])
  }
  extra_names <- unique(extra_names)

  ## Just make sure we have something to work with.
  if (sum(found_keepers, length(extra_names)) == 0) {
    message("The keepers has no elements in the coefficients.")
    message("Here are the keepers: ", toString(all_keepers))
    message("Here are the coefficients: ", toString(all_coefficients))
    stop("Unable to find the set of contrasts to keep, fix this and try again.")
  }

  numerators <- denominators <- c()
  keeper_len <- length(names(keepers))
  contrast_list <- names(keepers)
  kept_tables <- list()
  for (a in 1:length(names(keepers))) {
    name <- names(keepers)[a]
    extracted[["data"]][[name]] <- list()  ## Provided by combine_de_table() shortly.

    ## Initially, set same_string to the name of the table, then if there is a
    ## table with that, just use it directly rather than hunt for numerator/denominator.
    same_string <- name
    ## The numerators and denominators will be used to check that we are a_vs_b or b_vs_a
    numerator <- keepers[[name]][1]
    if (is.null(numerator)) {
      ## Then this should be the name of a special table, e.g. from extra_contrasts.
      message("Assuming ", name, " is a table from extra_contrasts.")
      inverse_string <- same_string
    } else {
      numerators[name] <- numerator
      denominator <- keepers[[name]][2]
      denominators[name] <- denominator
      mesg("Working on ", a, "/", keeper_len, ": ",  name,
           " which is: ", numerator, "/", denominator, ".")
      same_string <- numerator
      inverse_string <- numerator
      if (!is.na(denominator)) {
        same_string <- glue::glue("{numerator}_vs_{denominator}")
        inverse_string <- glue::glue("{denominator}_vs_{numerator}")
      }
    }

    ## Do the actual table search, checking for the same_string (a_vs_b) and
    ## inverse (b_vs_a) Set a flag do_inverse appropriately, this will be used
    ## later to flip some numbers.
    found <- 0
    found_table <- NULL
    do_inverse <- FALSE
    for (t in 1:length(table_names)) {
      tab <- table_names[t]
      if (tab == same_string) {
        do_inverse <- FALSE
        found <- found + 1
        found_table <- same_string
        mesg("Found table with ", same_string)
        kept_tables[tab] <- same_string
      } else if (tab == inverse_string) {
        do_inverse <- TRUE
        found <- found + 1
        found_table <- inverse_string
        kept_tables[tab] <- inverse_string
        mesg("Found inverse table with ", inverse_string)
      } else {
        kept_tables[tab] <- "Not found"
      }
    }
    if (found == 0) {
      warning("FOUND NEITHER ", same_string, " NOR ", inverse_string, "!")
      break
    }

    ## Now make a single table from the limma etc results.
    if (found > 0) {
      combined <- combine_single_de_table(
          li = limma, ed = edger, eb = ebseq, de = deseq, ba = basic,
          include_deseq = include_deseq, include_edger = include_edger,
          include_ebseq = include_ebseq, include_limma = include_limma,
          include_basic = include_basic,
          table_name = found_table, do_inverse = do_inverse,
          adjp = adjp, annot_df = annot_df,
          excludes = excludes, padj_type = padj_type,
          lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff)
      extracted[["data"]][[name]] <- combined[["data"]]
      extracted[["table_names"]][[name]] <- combined[["summary"]][["table"]]
      extracted[["kept"]] <- kept_tables
      extracted[["keepers"]] <- keepers
      plot_inputs <- list(
          "limma" = limma,
          "deseq" = deseq,
          "edger" = edger,
          "ebseq" = ebseq,
          "basic" = basic)
      ## Changing this to a try() for when we have weirdo extra_contrasts.
      extracted_plots <- try(combine_extracted_plots(
          name, combined, denominator, numerator, plot_inputs,
          include_basic = include_basic, include_deseq = include_deseq,
          include_edger = include_edger, include_limma = include_limma,
          include_ebseq = include_ebseq, loess = loess, logfc = lfc_cutoff, p = p_cutoff,
          do_inverse = do_inverse, found_table = found_table, p_type = padj_type), silent = TRUE)
      if ("try-error" %in% class(extracted_plots)) {
        extracted[["plots"]][[name]] <- NULL
      } else {
        extracted[["plots"]][[name]] <- extracted_plots
      }
      extracted[["summaries"]] <- rbind(extracted[["summaries"]],
                                        as.data.frame(combined[["summary"]]))
      extracted[["numerators"]] <- numerators
      extracted[["denominators"]] <- denominators
    } else {
      ## End checking that we found the numerator/denominator
      warning("Did not find either ", same_string, " nor ", inverse_string, ".")
      message("Did not find either ", same_string, " nor ", inverse_string, ".")
      break
    }
  } ## Ending the for loop of elements in the keepers list.
  return(extracted)
}

#' When a single 'keeper' contrast is specified, find and extract it.
#'
#' @param apr Data from all_pairwise().
#' @param extracted Tables extracted in combine_de_tables().
#' @param keepers In this case, one may assume either NULL or 'all'.
#' @param table_names The set of tables produced by all_pairwise().
#' @param all_coefficients The set of all experimental conditions in the
#'  experimental metadata.
#' @param limma The limma data from all_pairwise().
#' @param edger The edger data from all_pairwise().
#' @param ebseq The ebseq data from all_pairwise().
#' @param deseq The deseq data from all_pairwise().
#' @param basic The basic data from all_pairwise().
#' @param adjp Pull out the adjusted p-values from the data?
#' @param annot_df What annotations should be added to the table?
#' @param include_deseq Whether or not to include the deseq data.
#' @param include_edger Whether or not to include the edger data.
#' @param include_ebseq Whether or not to include the ebseq data.
#' @param include_limma Whether or not to include the limma data.
#' @param include_basic Whether or not to include the basic data.
#' @param excludes Set of genes to exclude.
#' @param padj_type Choose a specific p adjustment.
#' @param loess Add a loess to plots?
#' @param lfc_cutoff Passed for volcano/MA plots.
#' @param p_cutoff Passed for volcano/MA plots.
extract_keepers_single <- function(apr, extracted, keepers, table_names,
                                   all_coefficients,
                                   limma, edger, ebseq, deseq, basic,
                                   adjp, annot_df,
                                   include_deseq, include_edger,
                                   include_ebseq, include_limma,
                                   include_basic, excludes, padj_type,
                                   loess = FALSE, lfc_cutoff = 1, p_cutoff = 0.05) {
  splitted <- strsplit(x = keepers, split = "_vs_")
  numerator <- splitted[[1]][1]
  denominator <- splitted[[1]][2]
  inverse_keeper <- paste0(denominator, "_vs_", numerator)
  table <- keepers
  do_inverse <- FALSE
  if (keepers[1] %in% table_names) {
    mesg("I found ", keepers, " in the available contrasts.")
  } else if (inverse_keeper %in% table_names) {
    mesg("I found ", inverse_keeper, " the inverse keeper in the contrasts.")
    table <- inverse_keeper
    do_inverse <- TRUE
  } else {
    mesg("I did not find ", keepers, " in the available contrasts.", warn = TRUE)
    mesg("The available tables are: ", table_names, ".", warn = TRUE)
    table <- table_names[1]
    mesg("Choosing the first table: ", table)
  }
  combined <- combine_single_de_table(
      li = limma, ed = edger, eb = ebseq, de = deseq, ba = basic,
      include_deseq = include_deseq, include_edger = include_edger,
      include_ebseq = include_ebseq, include_limma = include_limma,
      include_basic = include_basic,
      table_name = table,
      adjp = adjp, annot_df = annot_df,
      excludes = excludes, padj_type = padj_type,
      lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff)
  extracted[["data"]][[table]] <- combined[["data"]]
  extracted[["table_names"]][[table]] <- combined[["summary"]][["table"]]
  extracted[["kept"]] <- table
  extracted[["keepers"]] <- keepers
  plot_inputs <- list(
      "limma" = limma,
      "deseq" = deseq,
      "edger" = edger,
      "ebseq" = ebseq,
      "basic" = basic)
  extracted[["plots"]][[table]] <- combine_extracted_plots(
      ## vv I changed this from 'name', I think that is correct but am uncertain.
      table, combined, denominator, numerator, plot_inputs,
      include_basic = include_basic, include_deseq = include_deseq,
      include_edger = include_edger, include_limma = include_limma,
      include_ebseq = include_ebseq, loess = loess, found_table = table,
      logfc = lfc_cutoff, p = p_cutoff, do_inverse = do_inverse, p_type = padj_type)
  extracted[["summaries"]] <- rbind(extracted[["summaries"]],
                                    as.data.frame(combined[["summary"]]))
  extracted[["numerators"]] <- numerator
  extracted[["denominators"]] <- denominator
  return(extracted)
}

#' Extract the sets of genes which are significantly more abundant than the rest.
#'
#' Given the output of something_pairwise(), pull out the genes for each contrast
#' which are the most/least abundant.  This is in contrast to extract_significant_genes().
#' That function seeks out the most changed, statistically significant genes.
#'
#' @param pairwise Output from _pairwise()().
#' @param according_to What tool(s) define 'most?'  One may use deseq, edger,
#'  limma, basic, all.
#' @param n How many genes to pull?
#' @param z Instead take the distribution of abundances and pull those past the
#'  given z score.
#' @param unique One might want the subset of unique genes in the top-n which
#'  are unique in the set of available conditions.  This will attempt to
#'  provide that.
#' @param excel Excel file to write.
#' @param ... Arguments passed into arglist.
#' @return The set of most/least abundant genes by contrast/tool.
#' @seealso \pkg{openxlsx}
#' @export
extract_abundant_genes <- function(pairwise, according_to = "all", n = 200,
                                   z = NULL, unique = FALSE,
                                   excel = "excel/abundant_genes.xlsx", ...) {
  arglist <- list(...)
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  abundant_lists <- list()
  final_list <- list()

  data <- NULL
  if (according_to[[1]] == "all") {
    according_to <- c("limma", "deseq", "edger", "basic")
  }

  for (type in according_to) {
    datum <- pairwise[[type]]
    abundant_lists[[type]] <- get_abundant_genes(datum, type = type, n = n, z = z,
                                                 unique = unique)
  }

  if (class(excel)[1] == "character") {
    mesg("Writing a legend of columns.")
    legend <- data.frame(rbind(
        c("The first ~3-10 columns of each sheet:",
          "are annotations provided by our chosen annotation source for this experiment."),
        c("Next column", "The most/least abundant genes.")),
        stringsAsFactors = FALSE)
    colnames(legend) <- c("column name", "column definition")
    xls_result <- write_xlsx(wb, data = legend, sheet = "legend", rownames = FALSE,
                             title = "Columns used in the following tables.")
  }

  ## Now make the excel sheet for each method/coefficient
  for (according in names(abundant_lists)) {
    for (coef in names(abundant_lists[[according]][["high"]])) {
      sheetname <- glue::glue("{according}_high_{coef}")
      annotations <- fData(pairwise[["input"]])
      high_abundances <- abundant_lists[[according]][["high"]][[coef]]
      kept_annotations <- names(high_abundances)
      kept_idx <- rownames(annotations) %in% kept_annotations
      kept_annotations <- annotations[kept_idx, ]
      high_data <- data.frame()
      if (nrow(annotations) > 0 & ncol(annotations) > 0) {
        high_data <- merge(data.frame(high_abundances), annotations, by = "row.names", all.x = TRUE)
        rownames(high_data) <- high_data[["Row.names"]]
        high_data[["Row.names"]] <- NULL
      } else {
        high_data <- as.data.frame(high_abundances)
      }
      start_row <- 1
      if (class(excel)[1] == "character") {
        title <- glue::glue("Table SXXX: High abundance genes in {coef} according to {according}.")
        xls_result <- write_xlsx(data = high_data, wb = wb, sheet = sheetname, title = title)
        start_row <- start_row + xls_result[["end_row"]] + 2
      }

      sheetname <- glue::glue("{according}_low_{coef}")
      low_abundances <- abundant_lists[[according]][["low"]][[coef]]
      kept_annotations <- names(low_abundances)
      kept_idx <- rownames(annotations) %in% kept_annotations
      kept_annotations <- annotations[kept_idx, ]
      low_data <- data.frame()
      if (nrow(annotations) > 0 & ncol(annotations) > 0) {
        low_data <- merge(data.frame(low_abundances), annotations,
                          by = "row.names", all.x = TRUE)
        rownames(low_data) <- low_data[["Row.names"]]
        low_data[["Row.names"]] <- NULL
      } else {
        low_data <- as.data.frame(low_abundances)
      }
      if (class(excel)[1] == "character") {
        title <- glue::glue("Table SXXX: Low abundance genes in {coef} according to {according}.")
        xls_result <- write_xlsx(data = low_data, wb = wb, sheet = sheetname, title = title)
      }
    } ## End the for loop
  }

  if (class(excel)[1] == "character") {
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }
  ret <- list(
      "with_annotations" = final_list,
      "abundances" = abundant_lists)
  return(ret)
}

#' Alias for extract_significant_genes because I am dumb.
#'
#' @param ... The parameters for extract_significant_genes()
#' @return It should return a reminder for me to remember my function names or
#'   change them to something not stupid.
#' @export
extract_siggenes <- function(...) {
  extract_significant_genes(...)
}
#' Extract the sets of genes which are significantly up/down regulated
#' from the combined tables.
#'
#' Given the output from combine_de_tables(), extract the genes in
#' which we have the greatest likely interest, either because they
#' have the largest fold changes, lowest p-values, fall outside a
#' z-score, or are at the top/bottom of the ranked list.
#'
#' @param combined Output from combine_de_tables().
#' @param according_to What tool(s) decide 'significant?'  One may use
#'  the deseq, edger, limma, basic, meta, or all.
#' @param lfc Log fold change to define 'significant'.
#' @param p (Adjusted)p-value to define 'significant'.
#' @param sig_bar Add bar plots describing various cutoffs of 'significant'?
#' @param z Z-score to define 'significant'.
#' @param n Take the top/bottom-n genes.
#' @param top_percent Use a percentage to get the top-n genes.
#' @param ma Add ma plots to the sheets of 'up' genes?
#' @param p_type use an adjusted p-value?
#' @param invert_barplots Invert the significance barplots as per Najib's request?
#' @param excel Write the results to this excel file, or NULL.
#' @param siglfc_cutoffs Set of cutoffs used to define levels of 'significant.'
#' @param ... Arguments passed into arglist.
#' @return The set of up-genes, down-genes, and numbers therein.
#' @seealso \code{\link{combine_de_tables}}
#' @export
extract_significant_genes <- function(combined, according_to = "all", lfc = 1.0,
                                      p = 0.05, sig_bar = TRUE, z = NULL, n = NULL, top_percent = NULL,
                                      ma = TRUE, p_type = "adj", invert_barplots = FALSE,
                                      excel = NULL, fc_column = NULL, p_column = NULL,
                                      siglfc_cutoffs = c(0, 1, 2), column_suffix = TRUE,
                                      gmt = NULL, category = "category",
                                      phenotype_name = "phenotype", set_name = "set",
                                      current_id = "ENSEMBL", required_id = "ENTREZID", ...) {
  arglist <- list(...)
  image_files <- c()  ## For cleaning up tmp image files after saving the xlsx file.

  ## The following two if statements are a little silly, but will be used later
  ## To check for the existence of a column in the input data via p_column %in% colnames()
  if (is.null(fc_column)) {
    fc_column <- ""
  }
  if (is.null(p_column)) {
    p_column <- ""
  }

  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  if (is.null(excel)) {
    ma <- FALSE
  }
  num_tables <- 0
  table_names <- NULL
  all_tables <- NULL
  table_mappings <- NULL
  if (class(combined)[1] == "data.frame") {
    ## Then this is just a data frame.
    all_tables[["all"]] <- combined
    table_names <- "all"
    num_tables <- 1
    table_mappings <- table_names
  } else if (class(combined)[1] == "combined_de") {
    ## Then this is the result of combine_de_tables()
    num_tables <- length(names(combined[["data"]]))
    ##table_names <- names(combined[["data"]])
    table_names <- combined[["table_names"]]
    all_tables <- combined[["data"]]
    table_mappings <- table_names
  } else if (!is.null(combined[["contrast_list"]])) {
    ## Then this is the result of all_pairwise()
    num_tables <- length(combined[["contrast_list"]])
    ## Extract the names of the tables which filled combined
    table_names <- names(combined[["data"]])
    ## Pull the table list
    all_tables <- combined[["data"]]
    ## Get the mappings of contrast_name -> table_name
    table_mappings <- combined[["keepers"]]
  } else {
    ## Then this is just a data frame.
    all_tables[["all"]] <- combined
    table_names <- "all"
    num_tables <- 1
    table_mappings <- table_names
  }

  if (!is.null(top_percent)) {
    n <- floor(nrow(all_tables[[1]]) * (top_percent / 100))
    mesg("Setting n to ", n)
  }
  if (!is.null(n)) {
    lfc <- NULL
    p <- NULL
  }

  logfc_suffix <- "_logfc"
  p_suffix <- "_p"
  adjp_suffix <- "_adjp"
  if (!isTRUE(column_suffix)) {
    logfc_suffix <- ""
    p_suffix <- ""
    adjp_suffix <- ""
  }

  trimmed_up <- list()
  trimmed_down <- list()
  up_titles <- list()
  down_titles <- list()
  sig_list <- list()
  title_append <- ""
  if (!is.null(lfc)) {
    title_append <- glue::glue("{title_append} |log2fc| >= {lfc}")
  }
  if (!is.null(p)) {
    title_append <- glue::glue("{title_append} p <= {p}")
  }
  if (!is.null(z)) {
    title_append <- glue::glue("{title_append} |z| >= {z}")
  }
  if (!is.null(n)) {
    title_append <- glue::glue("{title_append} top|bottom n={n}")
  }

  if (according_to[[1]] == "all") {
    according_to <- c("limma", "edger", "deseq", "ebseq", "basic")
  }

  if ("character" %in% class(excel)) {
    written <- write_sig_legend(wb)
    xls_result <- written[["xls_result"]]
  }

  ret <- list()
  sheet_count <- 0
  according_kept <- according_to
  chosen_columns <- list()
  ## Iterate over the according_to entries to see if there are valid p-value and logfc columns
  for (summary_count in 1:length(according_to)) {
    according <- according_to[summary_count]
    test_fc_param <- fc_column  ## Check if a column was provided by the user
    ## Otherwise make a column name from the method employed followed by the suffix.
    test_fc_column <- glue::glue("{according}{logfc_suffix}")
    test_p_param <- p_column
    test_p_column <- glue::glue("{according}{p_suffix}")
    test_adjp_column <- glue::glue("{according}{adjp_suffix}")
    skip <- TRUE
    if (test_fc_param %in% colnames(combined[["data"]][[1]])) {
      chosen_columns[[according]][["fc"]] <- test_fc_param
      if (test_p_param %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_p_param
      } else if (test_p_column %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_p_column
      } else if (test_adjp_column %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_adjp_column
      }
    } else if (test_fc_column %in% colnames(combined[["data"]][[1]])) {
      skip <- FALSE
      chosen_columns[[according]][["fc"]] <- test_fc_column
      if (test_p_param %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_p_param
      } else if (test_p_column %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_p_column
      } else if (test_adjp_column %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_adjp_column
      }
    }
    ## We have given some chances to not skip this group, if skip never got set to FALSE
    ## then this is probably not a useful criterion for searching.
    if (isTRUE(skip)) {
      mesg("Did not find the ", test_fc_column, ", skipping ", according, ".")
      according_kept <- according_to[!according == according_to]
      next
    }
  } ## End checking set of according_tos for valid entries.
  ## Our list of chosen_columns should have two levels, the first by method
  ## The second with logfc and p, which contain the columns of interest.

  according <- NULL
  according_to <- according_kept
  for (summary_count in 1:length(according_to)) {
    according <- according_to[summary_count]
    ret[[according]] <- list()
    ma_plots <- list()
    change_counts_up <- list()
    change_counts_down <- list()
    this_fc_column <- chosen_columns[[according]][["fc"]]
    this_p_column <- chosen_columns[[according]][["p"]]
    for (table_count in 1:length(table_names)) {
      table_name <- names(table_names)[table_count]
      plot_name <- as.character(table_names)[table_count]
      plot_name <- gsub(pattern = "-inverted", replacement = "", x = plot_name)
      ## Extract the MA data if requested.
      if (isTRUE(ma)) {
        ma_lfc <- lfc
        ma_p <- p
        if (is.null(ma_lfc)) {
          ma_lfc <- 1
        }
        if (is.null(ma_p)) {
          ma_p <- 0.05
        }
        single_ma <- try(extract_de_plots(
            combined[["input"]], type = according, table = plot_name, lfc = ma_lfc, p = ma_p), silent = TRUE)
        if ("try-error" %in% single_ma) {
          ma_plots[[table_name]] <- NULL
        } else {
          ma_plots[[table_name]] <- single_ma[["ma"]][["plot"]]
        }
      }

      this_table <- all_tables[[table_name]]
      trimming <- get_sig_genes(
          this_table, lfc = lfc, p = p, z = z, n = n, column = this_fc_column,
          p_column = this_p_column)

      trimmed_up[[table_name]] <- trimming[["up_genes"]]
      change_counts_up[[table_name]] <- nrow(trimmed_up[[table_name]])
      trimmed_down[[table_name]] <- trimming[["down_genes"]]
      change_counts_down[[table_name]] <- nrow(trimmed_down[[table_name]])
      up_title <- glue::glue("Table SXXX: Genes deemed significantly up in \\
                       {table_name} with {title_append} according to {according}.")
      up_titles[[table_name]] <- up_title
      down_title <- glue::glue("Table SXXX: Genes deemed significantly down in \\
                         {table_name} with {title_append} according to {according}.")
      down_titles[[table_name]] <- down_title
    } ## End extracting significant genes for loop

    change_counts <- as.data.frame(cbind(as.numeric(change_counts_up),
                                         as.numeric(change_counts_down)))
    colnames(change_counts) <- c("up", "down")
    rownames(change_counts)[table_count] <- table_name

    summary_title <- glue::glue("Counting the number of changed genes by contrast according to \\
                          {according} with {title_append}.")

    ret[[according]] <- list(
        "ups" = trimmed_up,
        "downs" = trimmed_down,
        "counts" = change_counts,
        "up_titles" = up_titles,
        "down_titles" = down_titles,
        "counts_title" = summary_title,
        "ma_plots" = ma_plots)

    ## I want to start writing out msigdb compatible gmt files and therefore
    ## want to start creating gene set collections from our data.

    do_excel <- TRUE
    if (is.null(excel)) {
      do_excel <- FALSE
    }
    if (isFALSE(excel)) {
      do_excel <- FALSE
    } else {
      mesg("Printing significant genes to the file: ", excel)
      xlsx_ret <- print_ups_downs(ret[[according]], wb, excel_basename, according = according,
                                  summary_count = summary_count, ma = ma)
      image_files <- c(xlsx_ret[["images_files"]], image_files)
      ## This is in case writing the sheet resulted in it being shortened.
      ## wb <- xlsx_ret[["workbook"]]
    } ## End of an if whether to print the data to excel
  } ## End list of according_to's

  ## the extraneous message() statements and instead fill that information into
  ## this data frame.
  name_element <- according_to[1]
  summary_df <- data.frame(row.names = names(ret[[name_element]][["ups"]]))

  sig_bar_plots <- NULL
  if (isTRUE(do_excel) & isTRUE(sig_bar)) {
    ## This needs to be changed to get_sig_genes()
    sig_bar_plots <- significant_barplots(
        combined, lfc_cutoffs = siglfc_cutoffs, invert = invert_barplots,
        p = p, z = z, p_type = p_type,
        according_to = according_to)
    plot_row <- 1
    plot_col <- 1
    mesg("Adding significance bar plots.")

    num_tables <- length(according_to)
    plot_row <- plot_row + ((nrow(change_counts) + 1) * num_tables) + 4
    ## The +4 is for the number of tools.
    ## I know it is silly to set the row in this very explicit fashion, but I
    ## want to make clear the fact that the table has a title, a set of
    ## headings, a length corresponding to the number of contrasts,  and then
    ## the new stuff should be added. Now add in a table summarizing the numbers
    ## in the plot. The information required to make this table is in
    ## sig_bar_plots[["ups"]][["limma"]] and sig_bar_plots[["downs"]][["limma"]]
    plot_row <- plot_row + 3
    ## I messed up something here.  The plots and tables
    ## at this point should start:
    ## 5(blank spaces and titles) + 4(table headings) + 4 * the number of contrasts.
    checked <- check_xlsx_worksheet(wb, "number_changed")

    for (according in according_to) {
      tmp_df <- ret[[according]][["counts"]]
      rownames(tmp_df) <- names(ret[[according]][["ups"]])
      colnames(tmp_df) <- paste0(according, "_", colnames(tmp_df))
      summary_df <- cbind(summary_df, tmp_df)
      sig_message <- as.character(glue::glue("Significant {according} genes."))
      xls_result <- openxlsx::writeData(
                                  wb = wb, sheet = "number_changed", x = sig_message,
                                  startRow = plot_row, startCol = plot_col)
      plot_row <- plot_row + 1
      plotname <- glue::glue("sigbar_{according}")
      try_result <- xlsx_plot_png(
          a_plot = sig_bar_plots[[according]], wb = wb, sheet = "number_changed", plotname = plotname,
          savedir = excel_basename, width = 9, height = 6, start_row = plot_row, start_col = plot_col)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      summary_row <- plot_row
      summary_col <- plot_col + 11
      de_summary <- summarize_ups_downs(sig_bar_plots[["ups"]][[according]],
                                        sig_bar_plots[["downs"]][[according]])
      xls_summary <- write_xlsx(
          data = de_summary, wb = wb, sheet = "number_changed", rownames = TRUE,
          start_row = summary_row, start_col = summary_col)
      plot_row <- plot_row + 30
    } ## End for loop writing out significance bar plots
  } ## End if we want significance bar plots
  ret[["sig_bar_plots"]] <- sig_bar_plots
  summary_df[["rownames"]] <- NULL
  ret[["summary_df"]] <- summary_df

  if (isTRUE(do_excel)) {
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }

  for (img in image_files) {
    removed <- try(suppressWarnings(file.remove(img)), silent = TRUE)
  }
  class(ret) <- c("sig_genes", "list")

  ## If there is no gmt fielname provided, but an excel filename is provided,
  ## save a .gmt with the xlsx file's basename
  if (is.null(gmt)) {
    if (!is.null(excel_basename)) {
      output_dir <- xlsx[["dirname"]]
      output_base <- paste0(xlsx[["basename"]], ".gmt")
      gmt <- file.path(output_dir, output_base)
    }
  }

  if (!is.null(gmt)) {
    message("Going to attempt to create gmt files from these results.")
    annotation_name <- annotation(combined[["input"]][["input"]])
    gsc <- make_gsc_from_pairwise(ret, according_to = according_to, orgdb = annotation_name,
                                  pair_names = c("ups", "downs"), category_name = category,
                                  phenotype_name = phenotype_name, set_name = set_name,
                                  current_id = current_id, required_id = required_id)
    types <- c("colored", "up", "down")
    for (t in types) {
      if (!is.null(gsc[[t]])) {
        for (g in 1:length(gsc[[t]])) {
          datum <- gsc[[t]][[g]]
          contrast_name <- names(gsc[[t]])[g]
          dname <- dirname(gmt)
          bname <- gsub(x = basename(gmt), pattern = "\\.gmt", replacement = "")
          newname <- paste0(bname, "_", t, "_", contrast_name, ".gmt")
          write_to <- file.path(dname, newname)
          written <- GSEABase::toGmt(datum, write_to)
        }
      }
    }
  }
  return(ret)
}

summarize_ups_downs <- function(ups, downs) {
  ## The ups and downs tables have 1 row for each contrast, 3 columns of numbers named
  ## 'a_up_inner', 'b_up_middle', 'c_up_outer'.
  ups <- ups[, -1]
  downs <- downs[, -1]
  ups[[1]] <- as.numeric(ups[[1]])
  ups[[2]] <- as.numeric(ups[[2]])
  ups[[3]] <- as.numeric(ups[[3]])
  ups[["up_sum"]] <- rowSums(ups)
  downs[[1]] <- as.numeric(downs[[1]])
  downs[[2]] <- as.numeric(downs[[2]])
  downs[[3]] <- as.numeric(downs[[3]])
  downs[["down_sum"]] <- rowSums(downs)
  summary_table <- as.data.frame(cbind(ups, downs))
  summary_table <- summary_table[, c(1, 2, 3, 5, 6, 7, 4, 8)]
  colnames(summary_table) <- c("up_from_0_to_2", "up_from_2_to_4", "up_gt_4",
                               "down_from_0_to_2", "down_from_2_to_4", "down_gt_4",
                               "sum_up", "sum_down")
  summary_table[["up_gt_2"]] <- summary_table[["up_from_2_to_4"]] +
    summary_table[["up_gt_4"]]
  summary_table[["down_gt_2"]] <- summary_table[["down_from_2_to_4"]] +
    summary_table[["down_gt_4"]]
  summary_table_idx <- rev(rownames(summary_table))
  summary_table <- summary_table[summary_table_idx, ]
  return(summary_table)
}

#' Find the sets of intersecting significant genes
#'
#' Use extract_significant_genes() to find the points of agreement between
#' limma/deseq/edger.
#'
#' @param combined Result from combine_de_tables().
#' @param lfc Define significant via fold-change.
#' @param p Or p-value.
#' @param padding_rows How much space to put between groups of data?
#' @param z Use a z-score filter?
#' @param p_type Use normal or adjusted p-values.
#' @param selectors List of methods to intersect.
#' @param order When set to the default 'inverse', go from the set with the most
#' least intersection to the most. E.g. Start with abc,bc,ac,c,ab,b,a as
#'  opposed to a,b,ab,c,ac,bc,abc.
#' @param excel An optional excel workbook to which to write.
#' @param ... Extra arguments for extract_significant_genes() and friends.
#' @return List containing the intersections between the various DE methods for
#'  both the up and down sets of genes.  It should also provide some venn
#'  diagrams showing the degree of similarity between the methods.
#' @examples
#'  \dontrun{
#'   expt <- create_expt(metadata="some_metadata.xlsx", gene_info=funkytown)
#'   big_result <- all_pairwise(expt, model_batch=FALSE)
#'   pretty <- combine_de_tables(big_result, excel="excel/combined_expt.xlsx")
#'   intersect <- intersect_significant(pretty, excel="excel/intersecting_genes.xlsx")
#'  }
#' @export
intersect_significant <- function(combined, lfc = 1.0, p = 0.05, padding_rows = 2,
                                  z = NULL, p_type = "adj", selectors = c("limma", "deseq", "edger"),
                                  order = "inverse", excel = "excel/intersect_significant.xlsx",
                                  ...) {
  ## Check the set of first->sixth and see that they exist in the combined table.
  arglist <- list(...)
  image_files <- c()
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  chosen_selectors <- c()
  extract_selectors <- c()
  alternate_selectors <- c()
  ## Use the column names from the first table to make this decision.
  possible_columns <- colnames(combined[["data"]][[1]])
  for (c in 1:length(selectors)) {
    rawname <- selectors[c]
    conjugate <- glue::glue("{rawname}_logfc")
    if (rawname %in% possible_columns) {
      chosen_selectors <- c(chosen_selectors, rawname)
      alternate_selectors <- c(alternate_selectors, rawname)
    } else if (conjugate %in% possible_columns) {
      chosen_selectors <- c(chosen_selectors, rawname)
      extract_selectors <- c(extract_selectors, rawname)
    } else if (arglist[["fc_column"]] %in% possible_columns &
               arglist[["p_column"]] %in% possible_columns) {
      chosen_selectors <- c(chosen_selectors, "alternate")
      alternate_selectors <- c(alternate_selectors, "alternate")
    } else {
      mesg("Skipping ", rawname)
    }
  }
  if (length(chosen_selectors) < 2) {
    stop("This requires more selectors in order to create an intersection.")
  }

  sig_genes <- sm(extract_significant_genes(combined, lfc = lfc, p = p,
                                            z = z, p_type = p_type, excel = NULL))
  if (length(alternate_selectors) > 0) {
    alt_genes <- extract_significant_genes(combined, according_to = "alternate", lfc = lfc,
                                           fc_column = arglist[["fc_column"]],
                                           p_column = arglist[["p_column"]], ma = FALSE,
                                           p = p, z = z, p_type = p_type, excel = FALSE)
    sig_genes[["alternate"]] <- alt_genes[["alternate"]]
  }

  xls_result <- NULL
  ## Set up the base data structure, a list of ups and a list of downs.
  lst <- list("ups" = list(), "downs" = list())
  set_names <- c()
  venn_row <- 1
  venn_col <- 1
  ## Arbitrarily pull the set of ups from the first 'chooser'
  table_lst <- names(sig_genes[[1]][["ups"]])
  for (t in 1:length(table_lst)) {
    table <- table_lst[t]
    lst[["ups"]][[table]] <- list()
    lst[["downs"]][[table]] <- list()
    for (dir in c("ups", "downs")) {
      for (i in 1:length(chosen_selectors)) {
        name <- chosen_selectors[i]
        lst[[dir]][[table]][[name]] <- rownames(sig_genes[[name]][[dir]][[table]])
      } ## End pulling the significants by selectors.
      sets <- Vennerable::Venn(Sets = lst[[dir]][[table]])
      intersections <- sets@IntersectionSets
      tmp_file <- tempfile(pattern = "venn", fileext = ".png")
      this_plot <- png(filename = tmp_file)
      controlled <- dev.control("enable")
      plt <- Vennerable::plot(sets, doWeights = FALSE)
      rec <- grDevices::recordPlot()
      dev.off()
      removed <- file.remove(tmp_file)
      lst[[dir]][[table]][["sets"]] <- sets
      lst[[dir]][[table]][["intersections"]] <- intersections
      lst[[dir]][[table]][["plot"]] <- rec
      print_order <- names(lst[[dir]][[table]][["intersections"]])
      if (order == "inverse") {
        print_order <- rev(print_order)
      }
      set_names <- list()
      inner_count <- 0
      for (symbolic in print_order) {
        inner_count <- inner_count + 1
        symbols <- strsplit(as.character(symbolic), "")[[1]]
        name <- c()
        for (s in 1:length(symbols)) {
          symbol <- symbols[s]
          if (symbol == 1) {
            name <- c(name, chosen_selectors[s])
          }
        }
        name <- toString(name)
        set_names[[inner_count]] <- name
      } ## End for symbolic in names(elements_per_set)
      set_names[length(set_names)] <- "none"
      set_names[1] <- "all"

      ## This has been the source of some confusion for me.
      ## I think that the
      if (order == "inverse") {
        names(set_names) <- rev(names(lst[[dir]][[table]][["intersections"]]))
      } else {
        set_names <- rev(set_names)
        names(invert_names) <- names(lst[[dir]][[table]][["intersections"]])
      }

      lst[[dir]][[table]][["set_names"]] <- set_names
      table_rows <- combined[["data"]][[table]]
      xlsx_row <- 1
      xlsx_table <- ""
      lst[[dir]][[table]][["data"]] <- list()
      for (s in 1:length(set_names)) {
        sname <- set_names[s]
        set_selection_name <- names(sname)
        clean_sname <- sname
        clean_sname <- gsub(pattern = "\\, ", replacement = "_", x = sname)

        lst[[dir]][[table]][["data"]][[clean_sname]] <- data.frame()
        ## This next operation was previously being done on the list index
        ## number, which is exactly wrong.
        if (length(lst[[dir]][[table]][["intersections"]][[set_selection_name]]) > 0) {
          idx <- lst[[dir]][[table]][["intersections"]][[set_selection_name]]
          table_subset <- table_rows[idx, ]
          lst[[dir]][[table]][["data"]][[clean_sname]] <- table_subset
          if (dir == "ups") {
            text_dir <- "up"
          } else {
            text_dir <- "down"
          }
          xlsx_table <- glue::glue("{text_dir}_{table}")
          xlsx_title <- glue::glue("Genes deemed {text_dir} significant via logFC: {lfc}\\
                             , p-value: {p}; by {clean_sname}.")
          if (!is.null(excel)) {
            xl_result <- write_xlsx(data = table_subset, wb = wb, sheet = xlsx_table,
                                    start_row = xlsx_row, title = xlsx_title)
            xlsx_row <- xlsx_row + nrow(table_subset) + padding_rows + 2
          } ## End checking to write the excel file.
        }
      } ## End iterating over writing the set names.
    } ## End ups vs. downs

    up_plot <- lst[["ups"]][[table]][["plot"]]
    down_plot <- lst[["downs"]][[table]][["plot"]]
    venn_title <- glue::glue("Summary of intersections among {toString(chosen_selectors)} \\
                       for {table}.")

    summary_df <- rbind(t(Vennerable::Weights(lst[["ups"]][[table]][["sets"]])),
                        t(Vennerable::Weights(lst[["downs"]][[table]][["sets"]])))
    rownames(summary_df) <- c("up", "down")
    tmp_colnames <- rev(lst[["ups"]][[table]][["set_names"]])
    tmp_colnames[length(tmp_colnames)] <- "all"
    colnames(summary_df) <- tmp_colnames
    summary_df <- summary_df[, -1]
    lst[["summary"]] <- summary_df
    if (!is.null(excel)) {
      xl_result <- write_xlsx(wb = wb,
                              data = summary_df,
                              sheet = "summary", title = venn_title,
                              start_row = venn_row, start_col = venn_col)
      venn_row <- venn_row + 4
      try_result <- xlsx_plot_png(
          up_plot, wb = wb, sheet = "summary", width = 6, height = 6,
          start_col = venn_col, start_row = venn_row)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      venn_col <- venn_col + 12
      try_result <- xlsx_plot_png(
          down_plot, wb = wb, sheet = "summary", width = 6, height = 6,
          start_col = venn_col, start_row = venn_row)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      venn_row <- venn_row + 26
    }
  }  ## End iterating over the tables

  if (!is.null(excel)) {
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }
  for (img in image_files) {
    removed <- try(suppressWarnings(file.remove(img)), silent = TRUE)
  }

  class(lst) <- c("sig_intersect", "list")
  return(lst)
}

#' Reprint the output from extract_significant_genes().
#'
#' I found myself needing to reprint these excel sheets because I
#' added some new information. This shortcuts that process for me.
#'
#' @param upsdowns Output from extract_significant_genes().
#' @param wb Workbook object to use for writing, or start a new one.
#' @param excel_basename Used when including plots in the xlsx sheet.
#' @param according Use limma, deseq, or edger for defining 'significant'.
#' @param summary_count For spacing sequential tables one after another.
#' @param ma Include ma plots?
#' @return Return from write_xlsx.
#' @seealso \code{\link{combine_de_tables}}
#' @export
print_ups_downs <- function(upsdowns, wb, excel_basename, according = "limma",
                            summary_count = 1, ma = FALSE) {
  image_files <- c()
  ups <- upsdowns[["ups"]]
  downs <- upsdowns[["downs"]]
  up_titles <- upsdowns[["up_titles"]]
  down_titles <- upsdowns[["down_titles"]]
  summary <- as.data.frame(upsdowns[["counts"]])
  summary_title <- upsdowns[["counts_title"]]
  ma_plots <- upsdowns[["ma_plots"]]
  table_count <- 0
  summary_count <- summary_count - 1
  num_tables <- length(names(ups))
  summary_start <- ((num_tables + 2) * summary_count) + 1
  xls_summary_result <- write_xlsx(wb = wb, data = summary, start_col = 1, start_row = summary_start,
                                   sheet = "number_changed", title = summary_title)
  xls_result <- NULL
  for (table_count in 1:length(names(ups))) {
    base_name <- names(ups)[table_count]
    up_name <- glue::glue("up_{according}_{base_name}")
    down_name <- glue::glue("down_{according}_{base_name}")
    up_table <- ups[[table_count]]
    up_title <- up_titles[[table_count]]
    if (nrow(up_table) > 0) {
      mesg(table_count, "/", num_tables, ": Creating significant table ", up_name)
      xls_result <- write_xlsx(data = up_table, wb = wb, sheet = up_name, title = up_title)
      ## This is in case the sheet name is past the 30 character limit.
      sheet_name <- xls_result[["sheet"]]
      if (isTRUE(ma)) {
        ma_row <- 1
        ma_col <- xls_result[["end_col"]] + 1
        if (!is.null(ma_plots[[base_name]])) {
          plot_name <- glue::glue("ma_{according}_{base_name}")
          try_result <- xlsx_plot_png(ma_plots[[base_name]], wb = wb, sheet = sheet_name,
                                      plotname = plot_name, savedir = excel_basename,
                                      start_row = ma_row, start_col = ma_col)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
        }
      }
    } else {
      mesg("The up table ", base_name, " is empty.")
    }


    down_table <- downs[[table_count]]
    down_title <- down_titles[[table_count]]
    if (nrow(down_table) > 0) {
      xls_result <- write_xlsx(data = down_table, wb = wb, sheet = down_name, title = down_title)
    } else {
      mesg("The down table ", base_name, " is empty.")
    }
  } ## End for each name in ups
  retlist <- list(
      "result" = xls_result,
      "image_files" = image_files)
  return(retlist)
}

#' Write the legend of an excel file for combine_de_tables()
#'
#' @param wb Workbook to write
#' @param excel_basename Where to write it
#' @param plot_dim Default plot size.
#' @param apr The all_pairwise() result.
#' @param limma The limma result, which is redundant.
#' @param include_limma Include the limma result?
#' @param deseq The deseq result, which is redundant.
#' @param include_deseq Include the deseq result?
#' @param edger The edger result, which is redundant.
#' @param include_edger Include the edger result?
#' @param ebseq The ebseq result, which is redundant.
#' @param include_ebseq Include the ebseq result?
#' @param basic Basic data
#' @param include_basic Include the basic result?
#' @param padj_type P-adjustment employed.
write_combined_legend <- function(wb, excel_basename, plot_dim, apr,
                                  limma, include_limma,
                                  deseq, include_deseq,
                                  edger, include_edger,
                                  ebseq, include_ebseq,
                                  basic, include_basic,
                                  padj_type) {
  ## I want to print a string reminding the user what kind of model was used in
  ## the analysis. Do that here.  Noting that if 'batch' is actually from a
  ## surrogate variable, then we will not have TRUE/FALSE but instead a matrix.
  do_excel <- TRUE
  if (length(excel_basename) == 0) {
    do_excel <- FALSE
  }
  reminder_model_cond <- apr[["model_cond"]]
  reminder_model_batch <- apr[["model_batch"]]
  reminder_extra <- apr[["extra_contrasts"]]
  reminder_string <- NULL
  if (class(reminder_model_batch)[1] == "matrix") {
    reminder_string <- "The contrasts were performed using surrogates from sva/ruv/etc."
  } else if (isTRUE(reminder_model_batch) & isTRUE(reminder_model_cond)) {
    reminder_string <- "The contrasts were performed with condition and batch in the model."
  } else if (isTRUE(reminder_model_cond)) {
    reminder_string <- "The contrasts were performed with only condition in the model."
  } else {
    reminder_string <- "The contrasts were performed in a strange way, beware!"
  }

  ## The next large set of data.frame() calls create the first sheet, containing a legend.
  mesg("Writing a legend of columns.")
  legend <- data.frame(rbind(
      c("", reminder_string),
      c("The first ~3-10 columns of each sheet:",
        "are annotations provided by our chosen annotation source for this experiment."),
      c("Next 6 columns", "The logFC and p-values reported by limma, edger, and deseq2.")
  ),
  stringsAsFactors = FALSE)
  deseq_legend <- data.frame(rbind(
      c("The next 7 columns", "Statistics generated by DESeq2."),
      c("deseq_logfc_rep", "The log2 fold change reported by DESeq2, again."),
      c("deseq_adjp_rep", "The adjusted-p value reported by DESeq2, again."),
      c("deseq_basemean", "The base mean of all samples according to DESeq2."),
      c("deseq_lfcse", "The standard error observed given the log2 fold change."),
      c("deseq_stat", "T-statistic reported by DESeq2 given the log2FC and observed variances."),
      c("deseq_p", "Resulting p-value."),
      c(glue::glue("deseq_adjp_{padj_type}"), glue::glue("p-value adjusted with {padj_type}")),
      c("deseq_q", "False-positive corrected p-value.")
  ),
  stringsAsFactors = FALSE)
  edger_legend <- data.frame(rbind(
      c("The next 7 columns", "Statistics generated by edgeR."),
      c("edger_logfc_rep", "The log2 fold change reported by edgeR, again."),
      c("edger_adjp_rep", "The adjusted-p value reported by edgeR, again."),
      c("edger_logcpm",
        "Similar DESeq2's basemean, except only including the samples in the comparison."),
      c("edger_lr", "Undocumented, I think it is the T-statistic calculated by edgeR."),
      c("edger_p", "The observed p-value from edgeR."),
      c(glue::glue("edger_adjp_{padj_type}"), glue::glue("p-value adjusted with {padj_type}")),
      c("edger_q", "The observed corrected p-value from edgeR.")
  ),
  stringsAsFactors = FALSE)
  ebseq_legend <- data.frame(rbind(
      c("The next 6 columns", "Statistics generated by EBSeq."),
      c("ebseq_FC", "The fold change reported by EBSeq."),
      c("ebseq_logfc", "The log2 fold change from EBSeq."),
      c("ebseq_postfc", "The post-probability fold change."),
      c("ebseq_mean", "Mean of the EBSeq values."),
      c("PPEE", "Post-probability that the numerator/denominator are equivalent."),
      c("PPDE",  "Post-probability that the numerator/denominator are different."),
      c("ebseq_adjp",  "Attempt at FDR correction of the PPEE, this is just copied from PPEE
and is in _no_ way statistically valid, but added as a plotting conveinence.")
),
stringsAsFactors = FALSE)
  limma_legend <- data.frame(rbind(
      c("The next 7 columns", "Statistics generated by limma."),
      c("limma_logfc_rep", "The log2 fold change reported by limma, again."),
      c("limma_adjp_rep", "The adjusted-p value reported by limma, again."),
      c("limma_ave", "Average log2 expression observed by limma across all samples."),
      c("limma_t", "T-statistic reported by limma given the log2FC and variances."),
      c("limma_p", "Derived from limma_t, the p-value asking 'is this logfc significant?'"),
      c(glue::glue("limma_adjp_{padj_type}"), glue::glue("p-value adjusted with {padj_type}")),
      c("limma_b", "Bayesian estimate of the log-odds significance."),
      c("limma_q", "A q-value FDR adjustment of the p-value above.")
  ),
  stringsAsFactors = FALSE)
  basic_legend <- data.frame(rbind(
      c("The next 8 columns", "Statistics generated by the basic analysis."),
      c("basic_nummed", "log2 median values of the numerator for this comparison."),
      c("basic_denmed", "log2 median values of the denominator for this comparison."),
      c("basic_numvar", "Variance observed in the numerator values."),
      c("basic_denvar", "Variance observed in the denominator values."),
      c("basic_logfc", "The log2 fold change observed by the basic analysis."),
      c("basic_t", "T-statistic from basic."),
      c("basic_p", "Resulting p-value."),
      c(glue::glue("basic_adjp_{padj_type}"), glue::glue("p-value adjusted with {padj_type}")),
      c("basic_adjp", "BH correction of the p-value.")
  ),
  stringsAsFactors = FALSE)
  summary_legend <- data.frame(rbind(
      c("The next 5 columns", "Summaries of the limma/deseq/edger results."),
      c("lfc_meta", "The mean fold-change value of limma/deseq/edger."),
      c("lfc_var", "The variance between limma/deseq/edger."),
      c("lfc_varbymed", "The ratio of the variance/median (lower means better agreement.)"),
      c("p_meta", "A meta-p-value of the mean p-values."),
      c("p_var", "Variance among the 3 p-values."),
      c("The last columns: top plot left",
        "Venn of genes with logFC > 0 and p-value <= 0.05 for limma/DESeq/Edger."),
      c("The last columns: top plot right",
        "Venn of genes with logFC < 0 and p-value <= 0.05 for limma/DESeq/Edger."),
      c("The last columns: second plot",
        "Scatter plot of the voom-adjusted/normalized counts for each coefficient."),
      c("The last columns: third plot",
        "Scatter plot of the adjusted/normalized counts for each coefficient from edgeR."),
      c("The last columns: fourth plot",
        "Scatter plot of the adjusted/normalized counts for each coefficient from DESeq."),
      c("", "If this data was adjusted with sva, check for a sheet 'original_pvalues' at the end.")
  ),
  stringsAsFactors = FALSE)

  ## Here we including only those columns which are relevant to the analysis performed.
  if (isTRUE(include_limma)) {
    legend <- rbind(legend,
                    c("limma_logfc", "The log2 fold change reported by limma."),
                    c("limma_adjp", "The adjusted-p value reported by limma."))
  }
  if (isTRUE(include_deseq)) {
    legend <- rbind(legend,
                    c("deseq_logfc", "The log2 fold change reported by DESeq2."),
                    c("deseq_adjp", "The adjusted-p value reported by DESeq2."))
  }
  if (isTRUE(include_edger)) {
    legend <- rbind(legend,
                    c("edger_logfc", "The log2 fold change reported by edgeR."),
                    c("edger_adjp", "The adjusted-p value reported by edgeR."))
  }

  ## Simplify the code a little by bringing the table definitions up here.
  ## Figure out the set of possible contrasts by name.
  ## I am changing this logic to have it include the union of table names
  ## rather than just pick up the first set that worked.
  table_names <- c()
  if (isTRUE(include_limma)) {
    legend <- rbind(legend, limma_legend)
    table_names <- c(table_names, limma[["contrasts_performed"]])
  } else {
    limma <- NULL
  }
  if (isTRUE(include_deseq)) {
    legend <- rbind(legend, deseq_legend)
    table_names <- c(table_names, deseq[["contrasts_performed"]])
  } else {
    deseq <- NULL
  }
  if (isTRUE(include_edger)) {
    legend <- rbind(legend, edger_legend)
    table_names <- c(table_names, names(edger[["contrast_list"]]))
  } else {
    edger <- NULL
  }
  if (isTRUE(include_ebseq)) {
    legend <- rbind(legend, ebseq_legend)
    table_names <- c(table_names, names(ebseq[["all_tables"]]))
  } else {
    ebseq <- NULL
  }
  if (isTRUE(include_basic)) {
    legend <- rbind(legend, basic_legend)
    table_names <- c(table_names, names(basic[["all_tables"]]))
  } else {
    basic <- NULL
  }
  ## Drop the set of tables names to the unique set.
  table_names <- unique(table_names)

  ## Make sure there were no errors and die if things went catastrophically wrong.
  if (!isTRUE(include_limma) & !isTRUE(include_deseq) & !isTRUE(include_edger) &
      !isTRUE(include_basic) & !isTRUE(include_ebseq)) {
    stop("None of the DE tools appear to have worked.")
  }
  if (length(table_names) == 0) {
    stop("Could not find the set of table names.")
  }

  if (isTRUE(include_limma) & isTRUE(include_deseq) &
      isTRUE(include_edger) & isTRUE(include_basic)) {
    legend <- rbind(legend, summary_legend)
  }
  colnames(legend) <- c("column name", "column definition")
  xls_result <- write_xlsx(
      wb, data = legend, sheet = "legend", rownames = FALSE,
      title = "Columns used in the following tables.")

  ## Some folks have asked for some PCA showing the before/after surrogates.
  ## Put that on the first sheet, then.
  ## This if (isTRUE()) is a little odd, perhaps it should be removed or moved up.
  xls_result <- NULL
  image_files <- c()
  if (isTRUE(do_excel)) {
    mesg("Printing pca plots before and after surrogate|batch estimation.")
    ## Add PCA before/after
    chosen_estimate <- apr[["batch_type"]]
    xl_result <- openxlsx::writeData(
                               wb = wb, sheet = "legend",
                               x = "PCA plot before surrogate estimation.",
                               startRow = 1, startCol = 10)
    try_result <- xlsx_plot_png(
        apr[["pre_batch"]][["plot"]], wb = wb, sheet = "legend", start_row = 2,
        width = (plot_dim * 3/2), height = plot_dim, start_col = 10,
        plotname = "pre_pca", savedir = excel_basename)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    xl_result <- openxlsx::writeData(
                               wb = wb, sheet = "legend",
                               x = as.character(glue::glue("PCA after surrogate estimation with: {chosen_estimate}")),
                               startRow = 36, startCol = 10)
    try_result <- xlsx_plot_png(
        apr[["post_batch"]][["plot"]], wb = wb, sheet = "legend", start_row = 37,
        width = (plot_dim * 3/2), height = plot_dim, start_col = 10,
        plotname = "pre_pca", savedir = excel_basename)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    pre_table <- write_xlsx(
        wb, data = apr[["pre_batch"]][["table"]],
        sheet = "legend", title = "Pre-Batch PCA table.",
        start_row = 66, start_col = 10)
    xls_result <- write_xlsx(
        wb, data = apr[["post_batch"]][["table"]],
        sheet = "legend", title = "Pre-Batch PCA table.",
        start_row = pre_table[["end_row"]] + 2, start_col = 10)
  }
  retlist <- list(
      "image_files" = image_files,
      "table_names" = table_names,
      "xls_result" = xls_result)
  return(retlist)
}

#' Internal function to write a summary of some combined data
#'
#' @param wb xlsx workbook to which to write.
#' @param excel_basename basename for printing plots.
#' @param apr a pairwise result
#' @param extracted table extracted from the pairwise result
#' @param compare_plots series of plots to print out.
#' @param lfc_cutoff Used for volcano/MA plots.
#' @param p_cutoff Used for volcano/MA plots.
write_combined_summary <- function(wb, excel_basename, apr, extracted, compare_plots,
                                   lfc_cutoff = 1, p_cutoff = 0.05) {
  image_files <- c()
  if (length(apr[["comparison"]]) == 0) {
    compare_plots <- FALSE
  }
  if (isTRUE(compare_plots)) {
    sheetname <- "pairwise_summary"
    xls_result <- write_xlsx(
        wb, data = extracted[["summaries"]], sheet = sheetname,
        title = "Summary of contrasts (lfc cutoff:{lfc_cutoff} p cutoff: {p_cutoff}).")
    new_row <- xls_result[["end_row"]] + 2
    xls_result <- write_xlsx(
        wb, data = apr[["comparison"]][["comp"]], sheet = sheetname, start_row = new_row,
        title = "Pairwise correlation coefficients among differential expression tools.")

    new_row <- xls_result[["end_row"]] + 2
    if (class(apr[["comparison"]][["heat"]])[1] == "recordedplot") {
      try_result <- xlsx_plot_png(
          apr[["comparison"]][["heat"]], wb = wb, sheet = sheetname, plotname = "pairwise_summary",
          savedir = excel_basename, start_row = new_row + 1, start_col = 1)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
    }
    logfc_comparisons <- try(compare_logfc_plots(extracted), silent = TRUE)
    if (class(logfc_comparisons) != "try-error") {
      logfc_names <- names(logfc_comparisons)
      new_row <- new_row + 2
      for (c in 1:length(logfc_names)) {
        lname <- logfc_names[c]
        new_row <- new_row + 32
        le <- logfc_comparisons[[c]][["le"]]
        ld <- logfc_comparisons[[c]][["ld"]]
        de <- logfc_comparisons[[c]][["de"]]
        tmpcol <- 1
        if (!is.null(le)) {
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     startRow = new_row - 2, startCol = tmpcol,
                                     x = glue::glue("Comparing DE tools for the \\
                                          comparison of: {logfc_names[c]}"))
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     startRow = new_row - 1, startCol = tmpcol,
                                     x="Log2FC(Limma vs. EdgeR)")
          try_result <- xlsx_plot_png(
              le, wb = wb, sheet = "pairwise_summary", plotname = "compare_le",
              savedir = excel_basename, start_row = new_row, start_col = tmpcol)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          tmpcol <- 8
        }
        if (!is.null(ld)) {
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     startRow = new_row - 1, startCol = tmpcol,
                                     x = "Log2FC(Limma vs. DESeq2)")
          try_result <- xlsx_plot_png(
              ld, wb = wb, sheet = sheetname, plotname = "compare_ld", savedir = excel_basename,
              start_row = new_row, start_col = tmpcol)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          tmpcol <- 15
        }
        if (!is.null(de)) {
          xl_result <- openxlsx::writeData(
                                     wb = wb, sheet = sheetname,
                                     startRow = new_row - 1, startCol = tmpcol,
                                     x = "Log2FC(DESeq2 vs. EdgeR)")
          try_result <- xlsx_plot_png(
              de, wb = wb, sheet = sheetname, plotname = "compare_ld", savedir = excel_basename,
              start_row = new_row, start_col = tmpcol)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
        }
      } ## End iterating over the comparison of logfc plots.
    } ## End checking if printing the logfc comparison plots worked.
  } ## End if compare_plots is TRUE
  retlist  <- list(
      "image_files" = image_files,
      "result" = xl_result)
  return(retlist)
}

#' Writes out the results of a single pairwise comparison.
#'
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the
#'     contrast matrix.
#' 2.  Write out the results() output for them in separate sheets in excel.
#' 3.  Since I have been using qvalues a lot for other stuff, add a column.
#'
#' Tested in test_24deseq.R
#' Rewritten in 2016-12 looking to simplify combine_de_tables().  That function
#' is far too big, this should become a template for that.
#'
#' @param data Output from results().
#' @param type Which DE tool to write.
#' @param excel Filename into which to save the xlsx data.
#' @param ... Parameters passed downstream, dumped into arglist and passed,
#'  notably the number of genes (n), the coefficient column (coef)
#' @return List of data frames comprising the toptable output for each
#'  coefficient, I also added a qvalue entry to these toptable() outputs.
#' @seealso \code{\link{write_xlsx}}
#' @examples
#' \dontrun{
#'  finished_comparison <- eBayes(deseq_output)
#'  data_list <- write_deseq(finished_comparison, workbook="excel/deseq_output.xls")
#' }
#' @export
write_de_table <- function(data, type = "limma", excel = "de_table.xlsx", ...) {
  arglist <- list(...)
  if (!is.null(data[[type]])) {
    data <- data[[type]]
  }
  n <- arglist[["n"]]
  if (is.null(n)) {
    n <- 0
  }
  coef <- arglist[["coef"]]
  if (is.null(coef)) {
    coef <- data[["contrasts_performed"]]
  } else {
    coef <- as.character(coef)
  }

  ## Figure out the number of genes if not provided
  if (n == 0) {
    n <- nrow(data[["coefficients"]])
  }

  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]

  return_data <- list()
  end <- length(coef)
  for (c in 1:end) {
    comparison <- coef[c]
    mesg("Writing ", c, "/", end, ": table: ", comparison, ".")
    table <- data[["all_tables"]][[c]]

    written <- try(write_xlsx(
        data = table, wb = wb, sheet = comparison,
        title = glue::glue("{type} results for: {comparison}.")))
  }

  save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  return(save_result)
}

#' Internal function to write a legend for significant gene tables.
#'
#' @param wb xlsx workbook object from openxlsx.
write_sig_legend <- function(wb) {
  legend <- data.frame(rbind(
      c("The first ~3-10 columns of each sheet:",
        "are annotations provided by our chosen annotation source for this experiment."),
      c("Next 6 columns", "The logFC and p-values reported by limma, edger, and deseq2."),
      c("limma_logfc", "The log2 fold change reported by limma."),
      c("deseq_logfc", "The log2 fold change reported by DESeq2."),
      c("edger_logfc", "The log2 fold change reported by edgeR."),
      c("limma_adjp", "The adjusted-p value reported by limma."),
      c("deseq_adjp", "The adjusted-p value reported by DESeq2."),
      c("edger_adjp", "The adjusted-p value reported by edgeR."),
      c("The next 5 columns", "Statistics generated by limma."),
      c("limma_ave", "Average log2 expression observed by limma across all samples."),
      c("limma_t", "T-statistic reported by limma given the log2FC and variances."),
      c("limma_p", "Derived from limma_t, the p-value asking 'is this logfc significant?'"),
      c("limma_b", "Use a Bayesian estimate to calculate log-odds significance."),
      c("limma_q", "A q-value FDR adjustment of the p-value above."),
      c("The next 5 columns", "Statistics generated by DESeq2."),
      c("deseq_basemean", "Analagous to limma's ave column, the base mean of all samples."),
      c("deseq_lfcse", "The standard error observed given the log2 fold change."),
      c("deseq_stat", "T-statistic reported by DESeq2 given the log2FC and observed variances."),
      c("deseq_p", "Resulting p-value."),
      c("deseq_q", "False-positive corrected p-value."),
      c("The next 4 columns", "Statistics generated by edgeR."),
      c("edger_logcpm",
        "Similar DESeq2's basemean, only including the samples in the comparison."),
      c("edger_lr", "Undocumented, I think it is the T-statistic calculated by edgeR."),
      c("edger_p", "The observed p-value from edgeR."),
      c("edger_q", "The observed corrected p-value from edgeR."),
      c("The next 7 columns", "Statistics generated by ebseq."),
      c("ebseq_fc", "Fold-change reported by ebseq."),
      c("ebseq_logfc", "Oddly, ebseq also reports a log2fc."),
      c("ebseq_postfc", "Post-analysis fold change from ebseq."),
      c("ebseq_mean", "Mean values in the ebseq analysis."),
      c("ebseq_ppee", "Prior probability that the numerators and denominators are the same."),
      c("ebseq_ppde", "... that they are different (eg. 1-ppee)."),
      c("ebseq_adjp", "Currently just a copy of ppee until I figure them out."),
      c("The next 8 columns", "Statistics generated by the basic analysis."),
      c("basic_nummed", "log2 median values of the numerator (like edgeR's basemean)."),
      c("basic_denmed", "log2 median values of the denominator for this comparison."),
      c("basic_numvar", "Variance observed in the numerator values."),
      c("basic_denvar", "Variance observed in the denominator values."),
      c("basic_logfc", "The log2 fold change observed by the basic analysis."),
      c("basic_t", "T-statistic from basic."),
      c("basic_p", "Resulting p-value."),
      c("basic_adjp", "BH correction of the p-value."),
      c("The next 5 columns", "Summaries of the limma/deseq/edger results."),
      c("lfc_meta", "The mean fold-change value of limma/deseq/edger."),
      c("lfc_var", "The variance between limma/deseq/edger."),
      c("lfc_varbymed", "The ratio of the variance/median (lower means better agreement.)"),
      c("p_meta", "A meta-p-value of the mean p-values."),
      c("p_var", "Variance among the 3 p-values."),
      c("The last columns: top plot left",
        "Venn diagram of the genes with logFC > 0 and p-value <= 0.05 for limma/DESeq/Edger."),
      c("The last columns: top plot right",
        "Venn diagram of the genes with logFC < 0 and p-value <= 0.05 for limma/DESeq/Edger."),
      c("The last columns: second plot",
        "Scatter plot of the voom-adjusted/normalized counts for each coefficient."),
      c("The last columns: third plot",
        "Scatter plot of the adjusted/normalized counts for each coefficient from edgeR."),
      c("The last columns: fourth plot",
        "Scatter plot of the adjusted/normalized counts for each coefficient from DESeq."),
      c("", "If this data was adjusted with sva, look 'original_pvalues' at the end.")
  ),
  stringsAsFactors = FALSE)

  colnames(legend) <- c("column name", "column definition")
  xls_result <- write_xlsx(wb, data = legend, sheet = "legend", rownames = FALSE,
                           title = "Columns used in the following tables.")
  retlist <- list("wb" = wb, xls_result = xls_result)
  return(retlist)
}

#' Put the metadata at the end of combined_de_tables()
#'
#' For the moment this is a stupidly short function.  I am betting we will
#' elaborate on this over time.
#'
#' @param wb workbook object.
#' @param apr Pairwise result.
write_sample_design <- function(wb, apr) {
  meta_df <- pData(apr[["input"]])
  xls_meta_result <- write_xlsx(wb = wb, data = meta_df,
                                sheet = "metadata", title = "Experiment metadata.")
  return(xls_meta_result)
}

## EOF
