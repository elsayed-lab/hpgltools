gather_ontology_genes <- function(table, mappings, new="ORF") {
    strings <- table[["geneID"]]
    separate <- strsplit(x=strings, split="/")
    rownames(mappings) <- make.names(mappings[["ENTREZID"]], unique=TRUE)

    ## Taken from https://stackoverflow.com/questions/23420331/replacing-elements-in-a-list-of-lists
    a <- as.relistable(separate)
    u <- paste0("X", unlist(a))
    u <- mappings[u, new]
    mapped <- relist(u, a)
    return(mapped)
}

#' Make a pretty table of clusterprofiler data in excel.
#'
#' It is my intention to make a function like this for each ontology tool in my repetoire
#'
#' @param cp_result A set of results from simple_clusterprofiler().
#' @param excel An excel file to which to write some pretty results.
#' @param wb  Workbook object to write to.
#' @param add_trees  Include topgoish ontology trees?
#' @param pval Choose a cutoff for reporting by p-value.
#' @param add_plots Include some pvalue plots in the excel output?
#' @param height  Height of included plots.
#' @param width and their width.
#' @param ... Extra arguments are passed to arglist.
#' @return The result from openxlsx in a prettyified xlsx file.
#' @seealso \pkg{openxlsx} \pkg{goseq}
#' @export
write_cp_data <- function(cp_result, excel="excel/clusterprofiler.xlsx", wb=NULL, add_trees=TRUE,
                             pval=0.1, add_plots=TRUE, height=15, width=10, ...) {
    arglist <- list(...)
    table_style <- "TableStyleMedium9"
    if (!is.null(arglist[["table_style"]])) {
        table_style <- arglist[["TableStyleMedium9"]]
    }
    excel_basename <- gsub(pattern="\\.xlsx", replacement="", x=excel)
    excel_dir <- dirname(excel)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold",
                                     border="Bottom", fontSize="30")
    }
    pval_column <- "limma_adjp"
    if (!is.null(arglist[["pval_column"]])) {
        pval_column <- arglist[["pval_column"]]
    }

    if (class(excel) == "character") {
        message("Writing a sheet containing the legend.")
        wb <- openxlsx::createWorkbook(creator="hpgltools")
        legend <- data.frame(rbind(
            c("Ontology", "Molecular Function, Biological Process, or Cellular Component."),
            c("Category", "Gene ontology Identifier."),
            c("Term", "Short definition of the category."),
            c("Over p-value", "Estimate of cp over-representation in the row's ontology category."),
            c("Q-value", "False discovery rate correction of the p-value."),
            c("DE genes in cat", "What genes provided are in this specific category?"),
            c("All genes in cat", "The full set of gene annotations included in this ontology category."),
            c("Num. de", "The number of genes in column 'F'."),
            c("Num. in cat", "The number of genes in column 'G'.")
        ))
        colnames(legend) <- c("column name", "column definition")
        xls_result <- write_xls(wb, data=legend, sheet="legend", rownames=FALSE,
                                title="Columns used in the following tables.")
        summary_row <- nrow(legend) + 5
        summary_df <- data.frame(rbind(
            c("Queried BP ontologies", nrow(cp_result[["enrich_go"]][["BP_all"]])),
            c("Significant BP ontologies", nrow(cp_result[["enrich_go"]][["BP_sig"]])),
            c("Queried MF ontologies", nrow(cp_result[["enrich_go"]][["MF_all"]])),
            c("Significant MF ontologies", nrow(cp_result[["enrich_go"]][["MF_sig"]])),
            c("Queried CC ontologies", nrow(cp_result[["enrich_go"]][["CC_all"]])),
            c("Significant CC ontologies", nrow(cp_result[["enrich_go"]][["CC_sig"]]))))
        colnames(summary_df) <- c("Ontology type", "Number found")
        xls_result <- write_xls(wb, data=summary_df, sheet="legend", rownames=FALSE,
                                title="Summary of the cp search.", start_row=1, start_col=4)
    }  ## End making sure that an excel is desired.

    ## Pull out the relevant portions of the cp data
    ## For this I am using the same (arbitrary) rules as in gather_ontology_genes()
    cp_mf <- cp_result[["enrich_go"]][["MF_sig"]]
    cp_mf <- cp_mf[ cp_mf[["pvalue"]] <= pval, ]
    cp_mf_genes <- gather_ontology_genes(cp_result[["enrich_go"]][["MF_sig"]], cp_result[["all_mappings"]])
    cp_mf[["named_genes"]] <- cp_mf_genes
    mf_idx <- order(cp_mf[["qvalue"]])
    cp_mf <- cp_mf[mf_idx, ]
    cp_mf[["Ontology"]] <- "MF"

    cp_bp <- cp_result[["enrich_go"]][["BP_sig"]]
    cp_bp <- cp_bp[ cp_bp[["pvalue"]] <= pval, ]
    cp_bp_genes <- gather_ontology_genes(cp_result[["enrich_go"]][["BP_sig"]], cp_result[["all_mappings"]])
    cp_bp[["named_genes"]] <- cp_bp_genes
    bp_idx <- order(cp_bp[["qvalue"]])
    cp_bp <- cp_bp[bp_idx, ]
    cp_bp[["Ontology"]] <- "BP"

    cp_cc <- cp_result[["enrich_go"]][["CC_sig"]]
    cp_cc <- cp_cc[ cp_cc[["pvalue"]] <= pval, ]
    cp_cc_genes <- gather_ontology_genes(cp_result[["enrich_go"]][["CC_sig"]], cp_result[["all_mappings"]])
    cp_cc[["named_genes"]] <- cp_cc_genes
    cc_idx <- order(cp_cc[["qvalue"]])
    cp_cc <- cp_cc[cc_idx, ]
    cp_cc[["Ontology"]] <- "CC"

    kept_columns <- c("ID", "Ontology", "Description", "GeneRatio", "BgRatio", "pvalue",
                      "p.adjust", "qvalue", "Count", "geneID", "named_genes")
    cp_mf <- cp_mf[, kept_columns]
    cp_bp <- cp_bp[, kept_columns]
    cp_cc <- cp_cc[, kept_columns]
    new_columns <- c("ID", "Ontology", "Description", "Ratio", "BgRatio", "P value",
                     "Adjusted P", "Q value", "Count", "gene ID", "Named ID")
    colnames(cp_mf) <- new_columns
    colnames(cp_bp) <- new_columns
    colnames(cp_cc) <- new_columns

    cp_kegg <- cp_result[["kegg_data"]][["kegg_sig"]]
    kegg_idx <- order(cp_kegg[["qvalue"]])
    cp_kegg <- cp_kegg[kegg_idx, ]

    cp_david <- cp_result[["david_data"]]
    david_idx <- order(cp_kegg[["qvalue"]])
    cp_david <- cp_david[david_idx, ]

    new_row <- 1
    message("Writing the BP data.")
    sheet <- "BP"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "BP Results from cp.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=cp_bp, tableStyle=table_style, startRow=new_row)
    ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
    if (isTRUE(add_plots)) {
        a_plot <- cp_result[["plots"]][["ego_sig_bp"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(cp_bp) + 2, start_row=new_row,
                                  plotname="bp_plot", savedir=excel_basename, doWeights=FALSE)
        b_plot <- cp_result[["plots"]][["tree_sig_bp"]]
        if (!is.null(b_plot)) {
            plot_try <- xlsx_plot_png(b_plot, wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(cp_bp) + 2, start_row=80, res=210,
                                      plotname="bp_trees", savedir=excel_basename)
        }
    }
    new_row <- new_row + nrow(cp_bp) + 2
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)

    new_row <- 1
    message("Writing the MF data.")
    sheet <- "MF"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "MF Results from cp.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=cp_mf, tableStyle=table_style, startRow=new_row)
    ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
    if (isTRUE(add_plots)) {
        a_plot <- cp_result[["plots"]][["ego_sig_mf"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(cp_mf) + 2, start_row=new_row,
                                  plotname="mf_plot", savedir=excel_basename, doWeights=FALSE)
        b_plot <- cp_result[["plots"]][["tree_sig_mf"]]
        if (!is.null(b_plot)) {
            plot_try <- xlsx_plot_png(b_plot, wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(cp_mf) + 2, start_row=80, res=210,
                                      plotname="mf_trees", savedir=excel_basename)
        }
    }
    new_row <- new_row + nrow(cp_mf) + 2
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)

    new_row <- 1
    message("Writing the CC data.")
    sheet <- "CC"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "CC Results from cp.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=cp_cc, tableStyle=table_style, startRow=new_row)
    ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
    if (isTRUE(add_plots)) {
        a_plot <- cp_result[["plots"]][["ego_sig_cc"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(cp_cc) + 2, start_row=new_row,
                                  plotname="cc_plot", savedir=excel_basename, doWeights=FALSE)
        b_plot <- cp_result[["plots"]][["tree_sig_cc"]]
        if (!is.null(b_plot)) {
            plot_try <- xlsx_plot_png(b_plot, wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(cp_cc) + 2, start_row=80, res=210,
                                      plotname="cc_trees", savedir=excel_basename)
        }
    }
    new_row <- new_row + nrow(cp_cc) + 2
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)

    new_row <- 1
    message("Writing the KEGG data.")
    sheet <- "KEGG"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "KEGG Results from cp.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=cp_kegg, tableStyle=table_style, startRow=new_row)
    ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology

    new_row <- 1
    message("Writing the DAVID data.")
    sheet <- "DAVID"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "DAVID Results from cp.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=cp_david, tableStyle=table_style, startRow=new_row)

    res <- openxlsx::saveWorkbook(wb, excel, overwrite=TRUE)
    message("Finished writing excel file.")
    return(res)
}

#' Make a pretty table of goseq data in excel.
#'
#' It is my intention to make a function like this for each ontology tool in my repetoire
#'
#' @param goseq_result  A set of results from simple_goseq().
#' @param excel  An excel file to which to write some pretty results.
#' @param wb   Workbook object to write to.
#' @param add_trees  Include topgoish ontology trees?
#' @param pval  Choose a cutoff for reporting by p-value.
#' @param add_plots  Include some pvalue plots in the excel output?
#' @param height  Height of included plots.
#' @param width  and their width.
#' @param ...  Extra arguments are passed to arglist.
#' @return  The result from openxlsx in a prettyified xlsx file.
#' @seealso \pkg{openxlsx} \pkg{goseq}
#' @export
write_goseq_data <- function(goseq_result, excel="excel/goseq.xlsx", wb=NULL, add_trees=TRUE,
                             pval=0.1, add_plots=TRUE, height=15, width=10, ...) {
    arglist <- list(...)
    table_style <- "TableStyleMedium9"
    if (!is.null(arglist[["table_style"]])) {
        table_style <- arglist[["TableStyleMedium9"]]
    }
    excel_basename <- gsub(pattern="\\.xlsx", replacement="", x=excel)
    excel_dir <- dirname(excel)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold",
                                     border="Bottom", fontSize="30")
    }
    pval_column <- "limma_adjp"
    if (!is.null(arglist[["pval_column"]])) {
        pval_column <- arglist[["pval_column"]]
    }

    if (class(excel) == "character") {
        message("Writing a sheet containing the legend.")
        wb <- openxlsx::createWorkbook(creator="hpgltools")
        legend <- data.frame(rbind(
            c("Ontology", "Molecular Function, Biological Process, or Cellular Component."),
            c("Category", "Gene ontology Identifier."),
            c("Term", "Short definition of the category."),
            c("Over p-value", "Estimate of goseq over-representation in the row's ontology category."),
            c("Q-value", "False discovery rate correction of the p-value."),
            c("DE genes in cat", "What genes provided are in this specific category?"),
            c("All genes in cat", "The full set of gene annotations included in this ontology category."),
            c("Num. de", "The number of genes in column 'F'."),
            c("Num. in cat", "The number of genes in column 'G'.")
        ))
        colnames(legend) <- c("column name", "column definition")
        xls_result <- write_xls(wb, data=legend, sheet="legend", rownames=FALSE,
                                title="Columns used in the following tables.")
        summary_row <- nrow(legend) + 5
        summary_df <- data.frame(rbind(
            c("Queried BP ontologies", nrow(goseq_result[["bp_subset"]])),
            c("Significant BP ontologies", nrow(goseq_result[["bp_interesting"]])),
            c("Queried MF ontologies", nrow(goseq_result[["mf_subset"]])),
            c("Significant MF ontologies", nrow(goseq_result[["mf_interesting"]])),
            c("Queried CC ontologies", nrow(goseq_result[["cc_subset"]])),
            c("Significant CC ontologies", nrow(goseq_result[["cc_interesting"]]))))
        colnames(summary_df) <- c("Ontology type", "Number found")
        xls_result <- write_xls(wb, data=summary_df, sheet="legend", rownames=FALSE,
                                title="Summary of the goseq search.", start_row=1, start_col=4)
        if (isTRUE(add_plots)) {
            printme <- "Histogram of observed ontology (adjusted) p-values by goseq."
            xl_result <- openxlsx::writeData(wb, "legend", x=printme,
                                             startRow=summary_row - 1, startCol=1)
            plot_try <- xlsx_plot_png(goseq_result[["pvalue_histogram"]], wb=wb, sheet="legend",
                                      start_col=1, start_row=summary_row, plotname="p_histogram",
                                      savedir=excel_basename)
        }
    }  ## End making sure that an excel is desired.

    trees <- NULL
    if (isTRUE(add_trees)) {
        trees <- goseq_trees(goseq_result, pval_column=pval_column)
    }

    ## Pull out the relevant portions of the goseq data
    ## For this I am using the same (arbitrary) rules as in gather_ontology_genes()
    goseq_mf <- goseq_result[["mf_subset"]]
    goseq_mf <- goseq_mf[ goseq_mf[["over_represented_pvalue"]] <= pval, ]
    goseq_mf_genes <- gather_ontology_genes(goseq_result, ontology="MF", pval=pval)
    mf_genes <- as.data.frame(goseq_mf_genes)
    rownames(mf_genes) <- rownames(goseq_mf_genes)
    goseq_mf <- merge(goseq_mf, mf_genes, by="row.names")
    rownames(goseq_mf) <- goseq_mf[["Row.names"]]
    goseq_mf <- goseq_mf[-1]
    mf_idx <- order(goseq_mf[["qvalue"]])
    goseq_mf <- goseq_mf[mf_idx, ]

    goseq_bp <- goseq_result[["bp_subset"]]
    goseq_bp <- goseq_bp[ goseq_bp[["over_represented_pvalue"]] <= pval, ]
    goseq_bp_genes <- gather_ontology_genes(goseq_result, ontology="BP", pval=pval)
    bp_genes <- as.data.frame(goseq_bp_genes)
    rownames(bp_genes) <- rownames(goseq_bp_genes)
    goseq_bp <- merge(goseq_bp, bp_genes, by="row.names")
    rownames(goseq_bp) <- goseq_bp[["Row.names"]]
    goseq_bp <- goseq_bp[-1]
    bp_idx <- order(goseq_bp[["qvalue"]])
    goseq_bp <- goseq_bp[bp_idx, ]

    goseq_cc <- goseq_result[["cc_subset"]]
    goseq_cc <- goseq_cc[ goseq_cc[["over_represented_pvalue"]] <= pval, ]
    goseq_cc_genes <- gather_ontology_genes(goseq_result, ontology="CC", pval=pval)
    cc_genes <- as.data.frame(goseq_cc_genes)
    rownames(cc_genes) <- rownames(goseq_cc_genes)
    goseq_cc <- merge(goseq_cc, cc_genes, by="row.names")
    rownames(goseq_cc) <- goseq_cc[["Row.names"]]
    goseq_cc <- goseq_cc[-1]
    cc_idx <- order(goseq_cc[["qvalue"]])
    goseq_cc <- goseq_cc[cc_idx, ]

    kept_columns <- c("ontology", "category", "term", "over_represented_pvalue",
                      "qvalue", "sig", "all", "numDEInCat", "numInCat",
                      "limma_sigfc", "deseq_sigfc", "edger_sigfc")
    goseq_mf <- goseq_mf[, kept_columns]
    goseq_bp <- goseq_bp[, kept_columns]
    goseq_cc <- goseq_cc[, kept_columns]
    new_columns <- c("Ontology", "Category", "Term", "Over p-value", "Q-value",
                     "DE genes in cat", "All genes in cat", "Num. DE", "Num. in cat.",
                     "FC from limma", "FC from DESeq", "FC from edgeR")
    colnames(goseq_mf) <- new_columns
    colnames(goseq_bp) <- new_columns
    colnames(goseq_cc) <- new_columns

    new_row <- 1
    message("Writing the BP data.")
    sheet <- "BP"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "BP Results from goseq.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=goseq_bp, tableStyle=table_style, startRow=new_row)
    ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
    if (isTRUE(add_plots)) {
        a_plot <- goseq_result[["pvalue_plots"]][["bpp_plot_over"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(goseq_bp) + 2, start_row=new_row,
                                  plotname="bp_plot", savedir=excel_basename, doWeights=FALSE)
        if (!is.null(trees[["BP_over"]])) {
            plot_try <- xlsx_plot_png(trees[["BP_over"]], wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(goseq_bp) + 2, start_row=80, res=210,
                                      plotname="bp_trees", savedir=excel_basename)
        }
    }
    new_row <- new_row + nrow(goseq_bp) + 2
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)

    new_row <- 1
    message("Writing the MF data.")
    sheet <- "MF"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "MF Results from goseq.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=goseq_mf, tableStyle=table_style, startRow=new_row)
    if (isTRUE(add_plots)) {
        a_plot <- goseq_result[["pvalue_plots"]][["mfp_plot_over"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(goseq_mf) + 2, start_row=new_row,
                                  plotname="mf_plot", savedir=excel_basename, doWeights=FALSE)
        if (!is.null(trees[["MF_over"]])) {
            plot_try <- xlsx_plot_png(trees[["MF_over"]], wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(goseq_bp) + 2, start_row=80, res=210,
                                      plotname="mf_trees", savedir=excel_basename)
        }
    }
    new_row <- new_row + nrow(goseq_mf) + 2
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)

    new_row <- 1
    message("Writing the CC data.")
    sheet <- "CC"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "CC Results from goseq.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=goseq_cc, tableStyle=table_style, startRow=new_row)
    if (isTRUE(add_plots)) {
        a_plot <- goseq_result[["pvalue_plots"]][["ccp_plot_over"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(goseq_cc) + 2, start_row=new_row,
                                  plotname="cc_plot", savedir=excel_basename, doWeights=FALSE)
        if (!is.null(trees[["CC_over"]])) {
            plot_try <- xlsx_plot_png(trees[["CC_over"]], wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(goseq_bp) + 2, start_row=80, res=210,
                                      plotname="cc_trees", savedir=excel_basename)
        }
    }
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)
    new_row <- new_row + nrow(goseq_cc) + 2

    res <- openxlsx::saveWorkbook(wb, excel, overwrite=TRUE)
    message("Finished writing excel file.")
    return(res)
}

#' Make a pretty table of gostats data in excel.
#'
#' It is my intention to make a function like this for each ontology tool in my repetoire
#'
#' @param gostats_result  A set of results from simple_gostats().
#' @param excel  An excel file to which to write some pretty results.
#' @param wb  Workbook object to write to.
#' @param add_trees  Include topgoish ontology trees?
#' @param pval  Choose a cutoff for reporting by p-value.
#' @param add_plots  Include some pvalue plots in the excel output?
#' @param height  Height of included plots.
#' @param width  and their width.
#' @param ...  Extra arguments are passed to arglist.
#' @return  The result from openxlsx in a prettyified xlsx file.
#' @seealso \pkg{openxlsx} \pkg{gostats}
#' @export
write_gostats_data <- function(gostats_result, excel="excel/gostats.xlsx", wb=NULL, add_trees=TRUE,
                             pval=0.1, add_plots=TRUE, height=15, width=10, ...) {
    arglist <- list(...)
    table_style <- "TableStyleMedium9"
    if (!is.null(arglist[["table_style"]])) {
        table_style <- arglist[["TableStyleMedium9"]]
    }
    excel_basename <- gsub(pattern="\\.xlsx", replacement="", x=excel)
    excel_dir <- dirname(excel)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold",
                                     border="Bottom", fontSize="30")
    }
    pval_column <- "limma_adjp"
    if (!is.null(arglist[["pval_column"]])) {
        pval_column <- arglist[["pval_column"]]
    }

    if (class(excel) == "character") {
        message("Writing a sheet containing the legend.")
        wb <- openxlsx::createWorkbook(creator="hpgltools")
        legend <- data.frame(rbind(
            c("Ontology", "Molecular Function, Biological Process, or Cellular Component."),
            c("Category", "Gene ontology Identifier."),
            c("Term", "Short definition of the category."),
            c("Over p-value", "Estimate of gostats over-representation in the row's ontology category."),
            c("Q-value", "False discovery rate correction of the p-value."),
            c("DE genes in cat", "What genes provided are in this specific category?"),
            c("All genes in cat", "The full set of gene annotations included in this ontology category."),
            c("Num. de", "The number of genes in column 'F'."),
            c("Num. in cat", "The number of genes in column 'G'.")
        ))
        colnames(legend) <- c("column name", "column definition")
        xls_result <- write_xls(wb, data=legend, sheet="legend", rownames=FALSE,
                                title="Columns used in the following tables.")
        summary_row <- nrow(legend) + 5
        summary_df <- data.frame(rbind(
            c("Queried BP ontologies", nrow(gostats_result[["bp_subset"]])),
            c("Significant BP ontologies", nrow(gostats_result[["bp_interesting"]])),
            c("Queried MF ontologies", nrow(gostats_result[["mf_subset"]])),
            c("Significant MF ontologies", nrow(gostats_result[["mf_interesting"]])),
            c("Queried CC ontologies", nrow(gostats_result[["cc_subset"]])),
            c("Significant CC ontologies", nrow(gostats_result[["cc_interesting"]]))))
        colnames(summary_df) <- c("Ontology type", "Number found")
        xls_result <- write_xls(wb, data=summary_df, sheet="legend", rownames=FALSE,
                                title="Summary of the gostats search.", start_row=1, start_col=4)
        if (isTRUE(add_plots)) {
            printme <- "Histogram of observed ontology (adjusted) p-values by gostats."
            xl_result <- openxlsx::writeData(wb, "legend", x=printme,
                                             startRow=summary_row - 1, startCol=1)
            plot_try <- xlsx_plot_png(gostats_result[["pvalue_histogram"]], wb=wb, sheet="legend",
                                      start_col=1, start_row=summary_row, plotname="p_histogram",
                                      savedir=excel_basename)
        }
    }  ## End making sure that an excel is desired.

    trees <- NULL
    if (isTRUE(add_trees)) {
        trees <- gostats_trees(gostats_result, pval_column=pval_column)
    }

    ## Pull out the relevant portions of the gostats data
    ## For this I am using the same (arbitrary) rules as in gather_ontology_genes()
    gostats_mf <- gostats_result[["mf_subset"]]
    gostats_mf <- gostats_mf[ gostats_mf[["over_represented_pvalue"]] <= pval, ]
    gostats_mf_genes <- gather_ontology_genes(gostats_result, ontology="MF", pval=pval)
    mf_genes <- as.data.frame(gostats_mf_genes)
    rownames(mf_genes) <- rownames(gostats_mf_genes)
    gostats_mf <- merge(gostats_mf, mf_genes, by="row.names")
    rownames(gostats_mf) <- gostats_mf[["Row.names"]]
    gostats_mf <- gostats_mf[-1]
    mf_idx <- order(gostats_mf[["qvalue"]])
    gostats_mf <- gostats_mf[mf_idx, ]

    gostats_bp <- gostats_result[["bp_subset"]]
    gostats_bp <- gostats_bp[ gostats_bp[["over_represented_pvalue"]] <= pval, ]
    gostats_bp_genes <- gather_ontology_genes(gostats_result, ontology="BP", pval=pval)
    bp_genes <- as.data.frame(gostats_bp_genes)
    rownames(bp_genes) <- rownames(gostats_bp_genes)
    gostats_bp <- merge(gostats_bp, bp_genes, by="row.names")
    rownames(gostats_bp) <- gostats_bp[["Row.names"]]
    gostats_bp <- gostats_bp[-1]
    bp_idx <- order(gostats_bp[["qvalue"]])
    gostats_bp <- gostats_bp[bp_idx, ]

    gostats_cc <- gostats_result[["cc_subset"]]
    gostats_cc <- gostats_cc[ gostats_cc[["over_represented_pvalue"]] <= pval, ]
    gostats_cc_genes <- gather_ontology_genes(gostats_result, ontology="CC", pval=pval)
    cc_genes <- as.data.frame(gostats_cc_genes)
    rownames(cc_genes) <- rownames(gostats_cc_genes)
    gostats_cc <- merge(gostats_cc, cc_genes, by="row.names")
    rownames(gostats_cc) <- gostats_cc[["Row.names"]]
    gostats_cc <- gostats_cc[-1]
    cc_idx <- order(gostats_cc[["qvalue"]])
    gostats_cc <- gostats_cc[cc_idx, ]

    kept_columns <- c("ontology", "category", "term", "over_represented_pvalue",
                      "qvalue", "sig", "all", "numDEInCat", "numInCat",
                      "limma_sigfc", "deseq_sigfc", "edger_sigfc")
    gostats_mf <- gostats_mf[, kept_columns]
    gostats_bp <- gostats_bp[, kept_columns]
    gostats_cc <- gostats_cc[, kept_columns]
    new_columns <- c("Ontology", "Category", "Term", "Over p-value", "Q-value",
                     "DE genes in cat", "All genes in cat", "Num. DE", "Num. in cat.",
                     "FC from limma", "FC from DESeq", "FC from edgeR")
    colnames(gostats_mf) <- new_columns
    colnames(gostats_bp) <- new_columns
    colnames(gostats_cc) <- new_columns

    new_row <- 1
    message("Writing the BP data.")
    sheet <- "BP"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "BP Results from gostats.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=gostats_bp, tableStyle=table_style, startRow=new_row)
    ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
    if (isTRUE(add_plots)) {
        a_plot <- gostats_result[["pvalue_plots"]][["bpp_plot_over"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(gostats_bp) + 2, start_row=new_row,
                                  plotname="bp_plot", savedir=excel_basename, doWeights=FALSE)
        if (!is.null(trees[["BP_over"]])) {
            plot_try <- xlsx_plot_png(trees[["BP_over"]], wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(gostats_bp) + 2, start_row=80, res=210,
                                      plotname="bp_trees", savedir=excel_basename)
        }
    }
    new_row <- new_row + nrow(gostats_bp) + 2
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)

    new_row <- 1
    message("Writing the MF data.")
    sheet <- "MF"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "MF Results from gostats.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=gostats_mf, tableStyle=table_style, startRow=new_row)
    if (isTRUE(add_plots)) {
        a_plot <- gostats_result[["pvalue_plots"]][["mfp_plot_over"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(gostats_mf) + 2, start_row=new_row,
                                  plotname="mf_plot", savedir=excel_basename, doWeights=FALSE)
        if (!is.null(trees[["MF_over"]])) {
            plot_try <- xlsx_plot_png(trees[["MF_over"]], wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(gostats_bp) + 2, start_row=80, res=210,
                                      plotname="mf_trees", savedir=excel_basename)
        }
    }
    new_row <- new_row + nrow(gostats_mf) + 2
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)

    new_row <- 1
    message("Writing the CC data.")
    sheet <- "CC"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "CC Results from gostats.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=gostats_cc, tableStyle=table_style, startRow=new_row)
    if (isTRUE(add_plots)) {
        a_plot <- gostats_result[["pvalue_plots"]][["ccp_plot_over"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(gostats_cc) + 2, start_row=new_row,
                                  plotname="cc_plot", savedir=excel_basename, doWeights=FALSE)
        if (!is.null(trees[["CC_over"]])) {
            plot_try <- xlsx_plot_png(trees[["CC_over"]], wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(gostats_bp) + 2, start_row=80, res=210,
                                      plotname="cc_trees", savedir=excel_basename)
        }
    }
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)
    new_row <- new_row + nrow(gostats_cc) + 2

    res <- openxlsx::saveWorkbook(wb, excel, overwrite=TRUE)
    message("Finished writing excel file.")
    return(res)
}

#' Write some excel results from a gprofiler search.
#'
#' Gprofiler is pretty awesome.  This function will attempt to write its results to an excel file.
#'
#' @param gprofiler_result  The result from simple_gprofiler().
#' @param wb  Optional workbook object, if you wish to append to an existing workbook.
#' @param excel  Excel file to which to write.
#' @param add_plots  Add some pvalue plots?
#' @param height  Height of included plots?
#' @param width  And their width.
#' @param ... More options, not currently used I think.
#' @return A prettyified table in an xlsx document.
#' @seealso \pkg{openxlsx} \pkg{gProfiler}
#' @export
write_gprofiler_data <- function(gprofiler_result, wb=NULL, excel="excel/gprofiler_result.xlsx",
                                 add_plots=TRUE, height=15, width=10, ...) {
    arglist <- list(...)
    table_style <- "TableStyleMedium9"
    if (!is.null(arglist[["table_style"]])) {
        table_style <- arglist[["TableStyleMedium9"]]
    }

    excel_basename <- gsub(pattern="\\.xlsx", replacement="", x=excel)
    excel_dir <- dirname(excel)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT",
                                     textDecoration="bold", border="Bottom", fontSize="30")
    }

    bp_data <- NULL
    mf_data <- NULL
    cc_data <- NULL
    kegg_data <- NULL
    tf_data <- NULL
    react_data <- NULL
    mi_data <- NULL
    hp_data <- NULL
    corum_data <- NULL

    do_go <- TRUE
    if (is.null(gprofiler_result[["go"]])) {
        do_go <- FALSE
    }
    if (nrow(gprofiler_result[["go"]]) == 0) {
        do_go <- FALSE
    }
    if (isTRUE(do_go)) {
        new_row <- 1
        sheet <- "GO_BP"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        go_data <- gprofiler_result[["go"]]
        bp_data <- go_data[go_data[["domain"]]=="BP", ]
        mf_data <- go_data[go_data[["domain"]]=="MF", ]
        cc_data <- go_data[go_data[["domain"]]=="CC", ]

        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=bp_data,
                                                tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["pvalue_plots"]][["bpp_plot_over"]]
            plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                      start_col=ncol(bp_data) + 2, start_row=new_row,
                                      plotname="bp_plot", savedir=excel_basename, doWeights=FALSE)
        }
        new_row <- new_row + nrow(bp_data) + 2

        new_row <- 1
        sheet <- "GO_MF"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=mf_data,
                                                tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["pvalue_plots"]][["mfp_plot_over"]]
            plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                      start_col=ncol(mf_data) + 2, start_row=new_row,
                                      plotname="mf_plot", savedir=excel_basename, doWeights=FALSE)
        }
        new_row <- new_row + nrow(mf_data) + 2

        new_row <- 1
        sheet <- "GO_CC"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=cc_data,
                                                tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["pvalue_plots"]][["ccp_plot_over"]]
            plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                      start_col=ncol(cc_data) + 2, start_row=new_row,
                                      plotname="cc_plot", savedir=excel_basename, doWeights=FALSE)
        }
        new_row <- new_row + nrow(cc_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking if go data is null

    do_kegg <- TRUE
    if (is.null(gprofiler_result[["kegg"]])) {
        do_kegg <- FALSE
    }
    if (nrow(gprofiler_result[["kegg"]]) == 0) {
        do_kegg <- FALSE
    }
    if (isTRUE(do_kegg)) {
        new_row <- 1
        sheet <- "KEGG"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        kegg_data <- gprofiler_result[["kegg"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=kegg_data,
                                                tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["pvalue_plots"]][["kegg_plot_over"]]
            plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                      start_col=ncol(kegg_data) + 2, start_row=new_row,
                                      plotname="kegg_plot", savedir=excel_basename, doWeights=FALSE)
        }
        new_row <- new_row + nrow(kegg_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking KEGG data


    do_tf <- TRUE
    if (is.null(gprofiler_result[["tf"]])) {
        do_tf <- FALSE
    }
    if (nrow(gprofiler_result[["tf"]]) == 0) {
        do_tf <- FALSE
    }
    if (isTRUE(do_tf)) {
        new_row <- 1
        sheet <- "transcription_factor"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        tf_data <- gprofiler_result[["tf"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=tf_data,
                                                tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["pvalue_plots"]][["tf_plot_over"]]
            plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                      start_col=ncol(tf_data) + 2, start_row=new_row,
                                      plotname="tf_plot", savedir=excel_basename, doWeights=FALSE)
        }
        new_row <- new_row + nrow(tf_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking tf data

    do_reactome <- TRUE
    if (is.null(gprofiler_result[["reactome"]])) {
        do_reactome <- FALSE
    }
    if (nrow(gprofiler_result[["reactome"]]) == 0) {
        do_reactome <- FALSE
    }
    if (isTRUE(do_reactome)) {
        new_row <- 1
        sheet <- "reactome"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        react_data <- gprofiler_result[["reactome"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=react_data,
                                                tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["pvalue_plots"]][["reactome_plot_over"]]
            plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                      start_col=ncol(react_data) + 2, start_row=new_row,
                                      plotname="react_plot", savedir=excel_basename, doWeights=FALSE)
        }
        new_row <- new_row + nrow(react_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking reactome data

    do_mi <- TRUE
    if (is.null(gprofiler_result[["mi"]])) {
        do_mi <- FALSE
    }
    if (nrow(gprofiler_result[["mi"]]) == 0) {
        do_mi <- FALSE
    }
    if (isTRUE(do_mi)) {
        new_row <- 1
        sheet <- "mirna"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        mi_data <- gprofiler_result[["mi"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=mi_data,
                                                tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["pvalue_plots"]][["mi_plot_over"]]
            plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                      start_col=ncol(mi_data) + 2, start_row=new_row,
                                      plotname="mi_plot", savedir=excel_basename, doWeights=FALSE)
        }
        new_row <- new_row + nrow(mi_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking mirna data

    do_hp <- TRUE
    if (is.null(gprofiler_result[["hp"]])) {
        do_hp <- FALSE
    }
    if (nrow(gprofiler_result[["hp"]]) == 0) {
        do_hp <- FALSE
    }
    if (isTRUE(do_hp)) {
        new_row <- 1
        sheet <- "hp"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        hp_data <- gprofiler_result[["hp"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=hp_data,
                                                tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["pvalue_plots"]][["hp_plot_over"]]
            plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                      start_col=ncol(hp_data) + 2, start_row=new_row,
                                      plotname="hp_plot", savedir=excel_basename, doWeights=FALSE)
        }
        new_row <- new_row + nrow(hp_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking corum data

    do_corum <- TRUE
    if (is.null(gprofiler_result[["corum"]])) {
        do_corum <- FALSE
    }
    if (nrow(gprofiler_result[["corum"]]) == 0) {
        do_corum <- FALSE
    }
    if (isTRUE(do_corum)) {
        new_row <- 1
        sheet <- "corum"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        corum_data <- gprofiler_result[["corum"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=corum_data,
                                                tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["pvalue_plots"]][["corum_plot_over"]]
            plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                      start_col=ncol(corum_data) + 2, start_row=new_row,
                                      plotname="corum_plot", savedir=excel_basename, doWeights=FALSE)
        }
        new_row <- new_row + nrow(corum_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking corum data

    excel_ret <- NULL
    if (!is.null(excel)) {
        excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite=TRUE))
    }
    if (class(excel_ret) == "try-error") {
        ## In case of an error with zip(1) in making an excel file.
        csvfile <- paste0(excel, "_gomf.csv")
        write.csv(mf_data, file=csvfile)
        csvfile <- paste0(excel, "_gobp.csv")
        write.csv(bp_data, file=csvfile)
        csvfile <- paste0(excel, "_gocc.csv")
        write.csv(cc_data, file=csvfile)
        csvfile <- paste0(excel, "_kegg.csv")
        write.csv(kegg_data, file=csvfile)
        csvfile <- paste0(excel, "_tf.csv")
        write.csv(tf_data, file=csvfile)
        csvfile <- paste0(excel, "_reactome.csv")
        write.csv(react_data, file=csvfile)
        csvfile <- paste0(excel, "_mi.csv")
        write.csv(mi_data, file=csvfile)
        csvfile <- paste0(excel, "_hp.csv")
        write.csv(hp_data, file=csvfile)
        csvfile <- paste0(excel, "_corum.csv")
        write.csv(corum_data, file=csvfile)
        excel_ret <- "csv"
    }
    message("Finished writing excel file.")
    return(excel_ret)
}

#' Make a pretty table of topgo data in excel.
#'
#' It is my intention to make a function like this for each ontology tool in my repetoire
#'
#' @param topgo_result  A set of results from simple_topgo().
#' @param excel  An excel file to which to write some pretty results.
#' @param wb  Workbook object to write to.
#' @param add_trees  Include topgoish ontology trees?
#' @param pval  Choose a cutoff for reporting by p-value.
#' @param add_plots  Include some pvalue plots in the excel output?
#' @param height  Height of included plots.
#' @param width  and their width.
#' @param ...  Extra arguments are passed to arglist.
#' @return The result from openxlsx in a prettyified xlsx file.
#' @seealso \pkg{openxlsx} \pkg{topgo}
#' @export
write_topgo_data <- function(topgo_result, excel="excel/topgo.xlsx", wb=NULL, add_trees=TRUE,
                             pval=0.1, add_plots=TRUE, height=15, width=10, ...) {
    arglist <- list(...)
    table_style <- "TableStyleMedium9"
    if (!is.null(arglist[["table_style"]])) {
        table_style <- arglist[["TableStyleMedium9"]]
    }
    excel_basename <- gsub(pattern="\\.xlsx", replacement="", x=excel)
    excel_dir <- dirname(excel)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold",
                                     border="Bottom", fontSize="30")
    }
    pval_column <- "limma_adjp"
    if (!is.null(arglist[["pval_column"]])) {
        pval_column <- arglist[["pval_column"]]
    }

    if (class(excel) == "character") {
        message("Writing a sheet containing the legend.")
        wb <- openxlsx::createWorkbook(creator="hpgltools")
        legend <- data.frame(rbind(
            c("Ontology", "Molecular Function, Biological Process, or Cellular Component."),
            c("Category", "Gene ontology Identifier."),
            c("Term", "Short definition of the category."),
            c("Over p-value", "Estimate of topgo over-representation in the row's ontology category."),
            c("Q-value", "False discovery rate correction of the p-value."),
            c("DE genes in cat", "What genes provided are in this specific category?"),
            c("All genes in cat", "The full set of gene annotations included in this ontology category."),
            c("Num. de", "The number of genes in column 'F'."),
            c("Num. in cat", "The number of genes in column 'G'.")
        ))
        colnames(legend) <- c("column name", "column definition")
        xls_result <- write_xls(wb, data=legend, sheet="legend", rownames=FALSE,
                                title="Columns used in the following tables.")
        summary_row <- nrow(legend) + 5
        summary_df <- data.frame(rbind(
            c("Queried BP ontologies", nrow(topgo_result[["bp_subset"]])),
            c("Significant BP ontologies", nrow(topgo_result[["bp_interesting"]])),
            c("Queried MF ontologies", nrow(topgo_result[["mf_subset"]])),
            c("Significant MF ontologies", nrow(topgo_result[["mf_interesting"]])),
            c("Queried CC ontologies", nrow(topgo_result[["cc_subset"]])),
            c("Significant CC ontologies", nrow(topgo_result[["cc_interesting"]]))))
        colnames(summary_df) <- c("Ontology type", "Number found")
        xls_result <- write_xls(wb, data=summary_df, sheet="legend", rownames=FALSE,
                                title="Summary of the topgo search.", start_row=1, start_col=4)
        if (isTRUE(add_plots)) {
            printme <- "Histogram of observed ontology (adjusted) p-values by topgo."
            xl_result <- openxlsx::writeData(wb, "legend", x=printme,
                                             startRow=summary_row - 1, startCol=1)
            plot_try <- xlsx_plot_png(topgo_result[["pvalue_histogram"]], wb=wb, sheet="legend",
                                      start_col=1, start_row=summary_row, plotname="p_histogram",
                                      savedir=excel_basename)
        }
    }  ## End making sure that an excel is desired.

    trees <- NULL
    if (isTRUE(add_trees)) {
        trees <- topgo_trees(topgo_result)
    }

    ## Pull out the relevant portions of the topgo data
    ## For this I am using the same (arbitrary) rules as in gather_ontology_genes()
    topgo_mf <- topgo_result[["mf_subset"]]
    topgo_mf <- topgo_mf[ topgo_mf[["over_represented_pvalue"]] <= pval, ]
    topgo_mf_genes <- gather_ontology_genes(topgo_result, ontology="MF", pval=pval)
    mf_genes <- as.data.frame(topgo_mf_genes)
    rownames(mf_genes) <- rownames(topgo_mf_genes)
    topgo_mf <- merge(topgo_mf, mf_genes, by="row.names")
    rownames(topgo_mf) <- topgo_mf[["Row.names"]]
    topgo_mf <- topgo_mf[-1]
    mf_idx <- order(topgo_mf[["qvalue"]])
    topgo_mf <- topgo_mf[mf_idx, ]

    topgo_bp <- topgo_result[["bp_subset"]]
    topgo_bp <- topgo_bp[ topgo_bp[["over_represented_pvalue"]] <= pval, ]
    topgo_bp_genes <- gather_ontology_genes(topgo_result, ontology="BP", pval=pval)
    bp_genes <- as.data.frame(topgo_bp_genes)
    rownames(bp_genes) <- rownames(topgo_bp_genes)
    topgo_bp <- merge(topgo_bp, bp_genes, by="row.names")
    rownames(topgo_bp) <- topgo_bp[["Row.names"]]
    topgo_bp <- topgo_bp[-1]
    bp_idx <- order(topgo_bp[["qvalue"]])
    topgo_bp <- topgo_bp[bp_idx, ]

    topgo_cc <- topgo_result[["cc_subset"]]
    topgo_cc <- topgo_cc[ topgo_cc[["over_represented_pvalue"]] <= pval, ]
    topgo_cc_genes <- gather_ontology_genes(topgo_result, ontology="CC", pval=pval)
    cc_genes <- as.data.frame(topgo_cc_genes)
    rownames(cc_genes) <- rownames(topgo_cc_genes)
    topgo_cc <- merge(topgo_cc, cc_genes, by="row.names")
    rownames(topgo_cc) <- topgo_cc[["Row.names"]]
    topgo_cc <- topgo_cc[-1]
    cc_idx <- order(topgo_cc[["qvalue"]])
    topgo_cc <- topgo_cc[cc_idx, ]

    kept_columns <- c("ontology", "category", "term", "over_represented_pvalue",
                      "qvalue", "sig", "all", "numDEInCat", "numInCat",
                      "limma_sigfc", "deseq_sigfc", "edger_sigfc")
    topgo_mf <- topgo_mf[, kept_columns]
    topgo_bp <- topgo_bp[, kept_columns]
    topgo_cc <- topgo_cc[, kept_columns]
    new_columns <- c("Ontology", "Category", "Term", "Over p-value", "Q-value",
                     "DE genes in cat", "All genes in cat", "Num. DE", "Num. in cat.",
                     "FC from limma", "FC from DESeq", "FC from edgeR")
    colnames(topgo_mf) <- new_columns
    colnames(topgo_bp) <- new_columns
    colnames(topgo_cc) <- new_columns

    new_row <- 1
    message("Writing the BP data.")
    sheet <- "BP"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "BP Results from topgo.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=topgo_bp, tableStyle=table_style, startRow=new_row)
    ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
    if (isTRUE(add_plots)) {
        a_plot <- topgo_result[["pvalue_plots"]][["bpp_plot_over"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(topgo_bp) + 2, start_row=new_row,
                                  plotname="bp_plot", savedir=excel_basename, doWeights=FALSE)
        if (!is.null(trees[["BP_over"]])) {
            plot_try <- xlsx_plot_png(trees[["BP_over"]], wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(topgo_bp) + 2, start_row=80, res=210,
                                      plotname="bp_trees", savedir=excel_basename)
        }
    }
    new_row <- new_row + nrow(topgo_bp) + 2
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)

    new_row <- 1
    message("Writing the MF data.")
    sheet <- "MF"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "MF Results from topgo.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=topgo_mf, tableStyle=table_style, startRow=new_row)
    if (isTRUE(add_plots)) {
        a_plot <- topgo_result[["pvalue_plots"]][["mfp_plot_over"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(topgo_mf) + 2, start_row=new_row,
                                  plotname="mf_plot", savedir=excel_basename, doWeights=FALSE)
        if (!is.null(trees[["MF_over"]])) {
            plot_try <- xlsx_plot_png(trees[["MF_over"]], wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(topgo_bp) + 2, start_row=80, res=210,
                                      plotname="mf_trees", savedir=excel_basename)
        }
    }
    new_row <- new_row + nrow(topgo_mf) + 2
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)

    new_row <- 1
    message("Writing the CC data.")
    sheet <- "CC"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "CC Results from topgo.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=topgo_cc, tableStyle=table_style, startRow=new_row)
    if (isTRUE(add_plots)) {
        a_plot <- topgo_result[["pvalue_plots"]][["ccp_plot_over"]]
        plot_try <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=width, height=height,
                                  start_col=ncol(topgo_cc) + 2, start_row=new_row,
                                  plotname="cc_plot", savedir=excel_basename, doWeights=FALSE)
        if (!is.null(trees[["CC_over"]])) {
            plot_try <- xlsx_plot_png(trees[["CC_over"]], wb=wb, sheet=sheet, width=12, height=12,
                                      start_col=ncol(topgo_bp) + 2, start_row=80, res=210,
                                      plotname="cc_trees", savedir=excel_basename)
        }
    }
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)
    new_row <- new_row + nrow(topgo_cc) + 2

    res <- openxlsx::saveWorkbook(wb, excel, overwrite=TRUE)
    message("Finished writing excel file.")
    return(res)
}

#'   Write gene ontology tables for data subsets
#'
#' Given a set of ontology results, this attempts to write them to an excel
#' workbook in a consistent and relatively easy-to-read fashion.
#'
#' @param kept_ontology  A result from subset_ontology_search()
#' @param outfile   Workbook to which to write.
#' @param dated    Append the year-month-day-hour to the workbook.
#' @param n   How many ontology categories to write for each search
#' @param overwritefile   Overwrite an existing workbook?
#' @param add_plots   Add the various p-value plots to the end of each sheet?
#' @param table_style   The chosen table style for excel
#' @param ...  some extra parameters
#' @return a set of excel sheet/coordinates
#' @seealso \pkg{openxlsx}
#' @examples
#' \dontrun{
#'  all_contrasts <- all_pairwise(expt, model_batch=TRUE)
#'  keepers <- list(bob = ('numerator','denominator'))
#'  kept <- combine_de_tables(all_contrasts, keepers=keepers)
#'  changed <- extract_significant_genes(kept)
#'  kept_ontologies <- subset_ontology_search(changed, lengths=gene_lengths,
#'                                            goids=goids, gff=gff, gff_type='gene')
#'  go_writer <- write_subset_ontologies(kept_ontologies)
#' }
#' @export
write_subset_ontologies <- function(kept_ontology, outfile="excel/subset_go", dated=TRUE,
                                    n=NULL, overwritefile=TRUE,
                                    add_plots=TRUE, table_style="TableStyleMedium9", ...) {
    arglist <- list(...)
    if (is.null(table_style)) {
        table_style <- "TableStyleMedium9"
    }
    if (is.null(outfile)) {
        outfile <- "excel/subset_go"
    }
    excel_dir <- dirname(outfile)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }

    suffix <- ".xlsx"
    outfile <- gsub(pattern="\\.xlsx", replacement="", outfile, perl=TRUE)
    outfile <- gsub(pattern="\\.xls", replacement="", outfile, perl=TRUE)
    excel_basename <- outfile  ## for making pdf plots

    types_list <- c("up_goseq", "down_goseq", "up_cluster", "down_cluster",
                    "up_topgo", "down_topgo", "up_gostats", "down_gostats",
                    "up_gprofiler", "down_gprofiler")
    ## names_list doesn't exist at this point, I losted it
    ## It is buried not very deep in kept_ontology I think
    names_list <- names(kept_ontology[["up_goseq"]])
    count <- 0
    for (name in names_list) {
        count <- count + 1
        up_filename <- paste0(outfile, "_up-", name)
        down_filename <- paste0(outfile, "_down-", name)
        if (isTRUE(dated)) {
            timestamp <- format(Sys.time(), "%Y%m%d%H")
            up_filename <- paste0(up_filename, "-", timestamp, suffix)
            down_filename <- paste0(down_filename, "-", timestamp, suffix)
        } else {
            up_filename <- paste0(up_filename, suffix)
            down_filename <- paste0(down_filename, suffix)
        }

        onts <- c("bp", "mf", "cc")
        up_stuff <- list()
        down_stuff <- list()
        for (ont in onts) {
            ONT <- toupper(ont)
            ## The goseq columns are probably wrong because I dropped one, remember that.
            varname <- paste0(ont, "_subset")

            if (!identical(list(), kept_ontology[["up_goseq"]])) {
                goseq_up <- kept_ontology[["up_goseq"]][[count]]
                goseq_up_ont <- goseq_up[[varname]]
                if (!is.null(n)) {
                    goseq_up_ont <- head(goseq_up_ont, n=n)
                }
                goseq_up_ont <- goseq_up_ont[, c(7, 1, 6, 2, 4, 5, 8)]
                colnames(goseq_up_ont) <- c("Ontology", "Category", "Term", "Over p-value",
                                            "Num. DE", "Num. in cat.", "Q-value")
                goseq_down <- kept_ontology[["down_goseq"]][[count]]
                goseq_down_ont <- goseq_down[[varname]]
                if (!is.null(n)) {
                    goseq_down_ont <- head(goseq_down_ont, n=n)
                }
                goseq_down_ont <- goseq_down_ont[, c(7, 1, 6, 2, 4, 5, 8)]
                colnames(goseq_down_ont) <- c("Ontology", "Category", "Term", "Over p-value",
                                              "Num. DE", "Num. in cat.", "Q-value")
                element_name <- paste0("goseq_", ont)
                up_stuff[[element_name]] <- goseq_up_ont
                down_stuff[[element_name]] <- goseq_down_ont
            }

            if (!identical(list(), kept_ontology[["up_cluster"]])) {
                varname <- paste0(ont, "_all")
                cluster_up <- kept_ontology[["up_cluster"]][[count]]
                cluster_up_ont <- as.data.frame(cluster_up[[varname]]@result)
                if (!is.null(n)) {
                    cluster_up_ont <- head(cluster_up_ont, n=n)
                }
                cluster_up_ont[["geneID"]] <- gsub(cluster_up_ont[["geneID"]],
                                                   pattern="/", replacement=" ")
                cluster_up_ont[["ontology"]] <- ONT
                cluster_up_ont <- cluster_up_ont[, c(10, 1, 2, 5, 3, 4, 6, 7, 9, 8)]
                colnames(cluster_up_ont) <- c("Ontology", "Category", "Term", "Over p-value",
                                              "Gene ratio", "BG ratio", "Adj. p-value", "Q-value",
                                              "Count", "Genes")
                cluster_down <- kept_ontology[["down_cluster"]][[count]]
                cluster_down_ont <- as.data.frame(cluster_down[[varname]]@result)
                if (!is.null(n)) {
                    cluster_down_ont <- head(cluster_down_ont, n=n)
                }
                cluster_down_ont[["geneID"]] <- gsub(cluster_down_ont[["geneID"]],
                                                     pattern="/", replacement=" ")
                cluster_down_ont[["ontology"]] <- ONT
                cluster_down_ont <- cluster_down_ont[, c(10, 1, 2, 5, 3, 4, 6, 7, 9, 8)]
                colnames(cluster_down_ont) <- c("Ontology", "Category", "Term", "Over p-value",
                                                "Gene ratio", "BG ratio", "Adj. p-value", "Q-value",
                                                "Count", "Genes")
                element_name <- paste0("cluster_", ont)
                up_stuff[[element_name]] <- cluster_up_ont
                down_stuff[[element_name]] <- cluster_down_ont
            }

            if (!identical(list(), kept_ontology[["up_topgo"]])) {
                varname <- paste0(ont, "_interesting")
                topgo_up <- kept_ontology[["up_topgo"]][[count]]
                topgo_up_ont <- topgo_up[["tables"]][[varname]]
                if (!is.null(n)) {
                    topgo_up_ont <- head(topgo_up_ont, n=n)
                }
                topgo_up_ont <- topgo_up_ont[, c(2, 1, 11, 6, 7, 8, 9, 10, 4, 3, 5)]
                colnames(topgo_up_ont) <- c("Ontology", "Category", "Term", "Fisher p-value",
                                            "Q-value", "KS score", "EL score", "Weight score",
                                            "Num. DE", "Num. in cat.", "Exp. in cat.")
                topgo_down <- kept_ontology[["down_topgo"]][[count]]
                topgo_down_ont <- topgo_down[["tables"]][[varname]]
                if (!is.null(n)) {
                    topgo_down_ont <- head(topgo_down_ont, n=n)
                }
                topgo_down_ont <- topgo_down_ont[, c(2, 1, 11, 6, 7, 8, 9, 10, 4, 3, 5)]
                colnames(topgo_down_ont) <- c("Ontology", "Category", "Term", "Fisher p-value",
                                              "Q-value", "KS score", "EL score", "Weight score",
                                              "Num. DE", "Num. in cat.", "Exp. in cat.")
                element_name <- paste0("topgo_", ont)
                up_stuff[[element_name]] <- topgo_up_ont
                down_stuff[[element_name]] <- topgo_down_ont
            }

            if (!identical(list(), kept_ontology[["up_gostats"]])) {
                varname <- paste0(ont, "_over_all")
                gostats_up <- kept_ontology[["up_gostats"]][[count]]
                gostats_up_ont <- gostats_up[[varname]]
                if (!is.null(n)) {
                    gostats_up_ont <- head(gostats_up_ont, n=n)
                }
                gostats_up_ont[["t"]] <- gsub(gostats_up_ont[["Term"]],
                                              pattern=".*\">(.*)</a>", replacement="\\1")
                gostats_up_ont[["Term"]] <- gsub(gostats_up_ont[["Term"]],
                                                 pattern="<a href=\"(.*)\">.*", replacement="\\1")
                gostats_up_ont[["ont"]] <- ONT
                gostats_up_ont <- gostats_up_ont[, c(10, 1, 9, 2, 5, 6, 3, 4, 8, 7)]
                colnames(gostats_up_ont) <- c("Ontology", "Category", "Term", "Fisher p-value",
                                              "Num. DE", "Num. in cat.", "Odds ratio", "Exp. in cat.",
                                              "Q-value", "Link")
                gostats_down <- kept_ontology[["down_gostats"]][[count]]
                gostats_down_ont <- gostats_down[[varname]]
                if (!is.null(n)) {
                    gostats_down_ont <- head(gostats_down_ont, n=n)
                }
                gostats_down_ont[["t"]] <- gsub(gostats_down_ont[["Term"]],
                                                pattern=".*\">(.*)</a>", replacement="\\1")
                gostats_down_ont[["Term"]] <- gsub(gostats_down_ont[["Term"]],
                                                   pattern="<a href=\"(.*)\">.*", replacement="\\1")
                gostats_down_ont[["ont"]] <- ONT
                gostats_down_ont <- gostats_down_ont[, c(10, 1, 9, 2, 5, 6, 3, 4, 8, 7)]
                colnames(gostats_down_ont) <- c("Ontology", "Category", "Term", "Fisher p-value",
                                                "Num. DE", "Num. in cat.", "Odds ratio", "Exp. in cat.",
                                                "Q-value", "Link")
                element_name <- paste0("gostats_", ont)
                up_stuff[[element_name]] <- gostats_up_ont
                down_stuff[[element_name]] <- gostats_down_ont
            }

            if (!identical(list(), kept_ontology[["up_gprofiler"]])) {
                gprofiler_up <- kept_ontology[["up_gprofiler"]][[count]]
                gprofiler_up_ont <- gprofiler_up[[varname]]
                if (!is.null(n)) {
                    gprofiler_up_ont <- head(gprofiler_up_ont, n=n)
                }
                gprofiler_up_ont <- gprofiler_up_ont[, c(7, 1, 6, 2, 4, 5, 8)]
                colnames(gprofiler_up_ont) <- c("Ontology", "Category", "Term", "Over p-value",
                                                "Num. DE", "Num. in cat.", "Q-value")
                gprofiler_down <- kept_ontology[["down_gprofiler"]][[count]]
                gprofiler_down_ont <- gprofiler_down[[varname]]
                if (!is.null(n)) {
                    gprofiler_down_ont <- head(gprofiler_down_ont, n=n)
                }
                gprofiler_down_ont <- gprofiler_down_ont[, c(7, 1, 6, 2, 4, 5, 8)]
                colnames(gprofiler_down_ont) <- c("Ontology", "Category", "Term", "Over p-value",
                                                  "Num. DE", "Num. in cat.", "Q-value")
                element_name <- paste0("gprofiler_", ont)
                up_stuff[[element_name]] <- gprofiler_up_ont
                down_stuff[[element_name]] <- gprofiler_down_ont
            }
        } ## End MF/BP/CC loop

        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold",
                                     border="Bottom", fontSize="30")
        ## This stanza will be repeated so I am just incrementing the new_row

        ## Write goseq data
        new_row <- 1
        sheet <- "goseq"
        openxlsx::addWorksheet(wb, sheetName=sheet)

        ## Write goseq BP data
        if (!is.null(up_stuff[["goseq_bp"]])) {
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["goseq_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_goseq"]][[name]][["pvalue_plots"]][["bpp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["goseq_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["goseq_bp"]]) + 2,
                                            start_row=new_row, plotname="goseq_bp",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["goseq_bp"]]) + 2
        }

        ## write goseq MF data
        if (!is.null(up_stuff[["goseq_mf"]])) {
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["goseq_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_goseq"]][[name]][["pvalue_plots"]][["mfp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["goseq_mf"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["goseq_mf"]]) + 2,
                                            start_row=new_row, plotname="goseq_mf",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["goseq_mf"]]) + 2
        }

        ## write goseq CC data
        if (!is.null(up_stuff[["goseq_cc"]])) {
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["goseq_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_goseq"]][[name]][["pvalue_plots"]][["ccp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["goseq_cc"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["goseq_cc"]]) + 2,
                                            start_row=new_row, plotname="goseq_cc",
                                            savedir=excel_basename)
            }
            openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
        }

        ## Move to  cluster profiler
        new_row <- 1
        sheet <- "clusterProfiler"
        openxlsx::addWorksheet(wb, sheetName=sheet)

        ## Write clusterprofiler BP data
        if (!is.null(up_stuff[["cluster_bp"]])) {
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["cluster_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_cluster"]][[name]][["pvalue_plots"]][["bp_all_barplot"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["cluster_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["cluster_bp"]]) + 2,
                                            start_row=new_row, plotname="cluster_bp",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["cluster_bp"]]) + 2
        }

        ## Write clusterprofiler MF data
        if (!is.null(up_stuff[["cluster_mf"]])) {
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["cluster_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_cluster"]][[name]][["pvalue_plots"]][["mf_all_barplot"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["cluster_mf"]]) + 2,
                ##                     startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["cluster_mf"]]) + 2,
                                            start_row=new_row, plotname="cluster_mf",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["cluster_mf"]]) + 2
        }

        ## Write clusterprofiler CC data
        if (!is.null(up_stuff[["cluster_cc"]])) {
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["cluster_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_cluster"]][[name]][["pvalue_plots"]][["cc_all_barplot"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["cluster_cc"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["cluster_cc"]]) + 2,
                                            start_row=new_row, plotname="cluster_cc",
                                            savedir=excel_basename)
            }
        }
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")

        ## Move to topgo
        new_row <- 1
        sheet <- "topgo"
        openxlsx::addWorksheet(wb, sheetName=sheet)

        ## Write topgo BP results
        if (!is.null(up_stuff[["topgo_bp"]])) {
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["topgo_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_topgo"]][[name]][["pvalue_plots"]][["bpp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["topgo_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["topgo_bp"]]) + 2,
                                            start_row=new_row, plotname="topgo_bp",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["topgo_bp"]]) + 2
        }

        ## write topgo MF results
        if (!is.null(up_stuff[["topgo_mf"]])) {
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["topgo_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_topgo"]][[name]][["pvalue_plots"]][["mfp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["topgo_mf"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["topgo_mf"]]) + 2,
                                            start_row=new_row, plotname="topgo_mf",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["topgo_mf"]]) + 2
        }

        ## and cc
        if (!is.null(up_stuff[["topgo_cc"]])) {
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["topgo_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_topgo"]][[name]][["pvalue_plots"]][["ccp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["topgo_cc"]]) + 2,
                ##                     startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["topgo_cc"]]) + 2,
                                            start_row=new_row, plotname="topgo_cc",
                                            savedir=excel_basename)
            }
            openxlsx::setColWidths(wb, sheet=sheet, cols=2:11, widths="auto")
        }

        ## move to gostats
        new_row <- 1
        sheet <- "gostats"
        openxlsx::addWorksheet(wb, sheetName=sheet)

        ## Write gostats BP stuff
        if (!is.null(up_stuff[["gostats_bp"]])) {
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["gostats_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            links <- up_stuff[["gostats_bp"]][["Link"]]
            class(links) <- "hyperlink"
            names(links) <- up_stuff[["gostats_bp"]][["Category"]]
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_gostats"]][[name]][["pvalue_plots"]][["bp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["gostats_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["gostats_bp"]]) + 2,
                                            start_row=new_row, plotname="gostats_bp",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["gostats_bp"]]) + 2
        }

        ## and mf
        if (!is.null(up_stuff[["gostats_mf"]])) {
            openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
            message("The previous line was a warning about overwriting existing data because of a link.")
            new_row <- new_row + nrow(up_stuff[["gostats_bp"]]) + 2
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["gostats_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            links <- up_stuff[["gostats_mf"]][["Link"]]
            class(links) <- "hyperlink"
            names(links) <- up_stuff[["gostats_mf"]][["Category"]]
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_gostats"]][[name]][["pvalue_plots"]][["mf_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["gostats_mf"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["gostats_mf"]]) + 2,
                                            start_row=new_row, plotname="gostats_mf",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["gostats_mf"]]) + 2
        }

        ## and cc
        if (!is.null(up_stuff[["gostats_cc"]])) {
            openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
            new_row <- new_row + nrow(up_stuff[["gostats_mf"]]) + 2
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["gostats_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            links <- up_stuff[["gostats_cc"]][["Link"]]
            class(links) <- 'hyperlink'
            names(links) <- up_stuff[["gostats_cc"]][["Category"]]
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_gostats"]][[name]][["pvalue_plots"]][["cc_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["gostats_cc"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["gostats_cc"]]) + 2,
                                            start_row=new_row, plotname="gostats_cc",
                                            savedir=excel_basename)
            }
            openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
        }

        ## Write gprofiler data
        new_row <- 1
        sheet <- "gprofiler"
        openxlsx::addWorksheet(wb, sheetName=sheet)

        ## Write gprofiler BP data
        if (!is.null(up_stuff[["gprofiler_bp"]])) {
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["gprofiler_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_gprofiler"]][[name]][["pvalue_plots"]][["bpp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["gprofiler_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["gprofiler_bp"]]) + 2,
                                            start_row=new_row, plotname="gprofiler_bp",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["gprofiler_bp"]]) + 2
        }

        ## write gprofiler MF data
        if (!is.null(up_stuff[["gprofiler_mf"]])) {
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["gprofiler_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_gprofiler"]][[name]][["pvalue_plots"]][["mfp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["gprofiler_mf"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["gprofiler_mf"]]) + 2,
                                            start_row=new_row, plotname="gprofiler_mf",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(up_stuff[["gprofiler_mf"]]) + 2
        }

        ## write gprofiler CC data
        if (!is.null(up_stuff[["gprofiler_cc"]])) {
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            cnew_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=up_stuff[["gprofiler_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["up_gprofiler"]][[name]][["pvalue_plots"]][["ccp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(up_stuff[["gprofiler_cc"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(up_stuff[["gprofiler_cc"]]) + 2,
                                            start_row=new_row, plotname="gprofiler_cc",
                                            savedir=excel_basename)
            }
            openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
        }

        ## Now the down data.

        openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
        res <- openxlsx::saveWorkbook(wb, up_filename, overwrite=TRUE)

        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold",
                                     border="Bottom", fontSize="30")
        ## This stanza will be repeated so I am just incrementing the new_row

        ## Starting with goseq
        new_row <- 1
        sheet <- "goseq"
        openxlsx::addWorksheet(wb, sheetName=sheet)

        ## Goseq down BP data
        if (!is.null(down_stuff[["goseq_bp"]])) {
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["goseq_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_goseq"]][[name]][["pvalue_plots"]][["bpp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["goseq_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["goseq_bp"]]) + 2,
                                            start_row=new_row, plotname="goseq_bp",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(down_stuff[["goseq_bp"]]) + 2
        }

        ## Goseq down MF data
        if (!is.null(down_stuff[["goseq_mf"]])) {
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["goseq_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_goseq"]][[name]][["pvalue_plots"]][["mfp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["goseq_mf"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["goseq_mf"]]) + 2,
                                            start_row=new_row, plotname="goseq_mf",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(down_stuff[["goseq_mf"]]) + 2
        }

        ## Goseq down CC data
        if (!is.null(down_stuff[["goseq_cc"]])) {
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["goseq_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_goseq"]][[name]][["pvalue_plots"]][["ccp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["goseq_cc"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["goseq_cc"]]) + 2,
                                            start_row=new_row, plotname="goseq_cc",
                                            savedir=excel_basename)
            }
        }
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")

        ## Now clusterprofiler data
        new_row <- 1
        sheet <- "clusterProfiler"

        ## cp down bp
        if (!is.null(down_stuff[["cluster_bp"]])) {
            openxlsx::addWorksheet(wb, sheetName=sheet)
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["cluster_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_cluster"]][[name]][["pvalue_plots"]][["bp_all_barplot"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["cluster_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["cluster_bp"]]) + 2,
                                            start_row=new_row, plotname="cluster_bp",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(down_stuff[["cluster_bp"]]) + 2
        }

        ## cp down mf
        if (!is.null(down_stuff[["cluster_mf"]])) {
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["cluster_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_cluster"]][[name]][["pvalue_plots"]][["mf_all_barplot"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["cluster_mf"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["cluster_mf"]]) + 2,
                                            start_row=new_row, plotname="cluster_mf",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(down_stuff[["cluster_mf"]]) + 2
        }

        ## cp down cc
        if (!is.null(down_stuff[["cluster_cc"]])) {
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["cluster_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_cluster"]][[name]][["cc_all_barplot"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["cluster_cc"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["cluster_cc"]]) + 2,
                                            start_row=new_row, plotname="cluster_cc",
                                            savedir=excel_basename)
            }
        }

        ## Move to topgo
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
        new_row <- 1
        sheet <- "topgo"
        openxlsx::addWorksheet(wb, sheetName=sheet)

        ## tp down bp
        if (!is.null(down_stuff[["topgo_bp"]])) {
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["topgo_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_topgo"]][[name]][["pvalue_plots"]][["bpp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["topgo_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["topgo_bp"]]) + 2,
                                            start_row=new_row, plotname="topgo_bp",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(down_stuff[["topgo_bp"]]) + 2
        }

        ## tp down mf
        if (!is.null(down_stuff[["topgo_mf"]])) {
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["topgo_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_topgo"]][[name]][["pvalue_plots"]][["mfp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["topgo_mf"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["topgo_mf"]]) + 2,
                                            start_row=new_row, plotname="topgo_mf",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(down_stuff[["topgo_mf"]]) + 2
        }

        ## tp down cc
        if (!is.null(down_stuff[["topgo_cc"]])) {
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["topgo_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_topgo"]][[name]][["pvalue_plots"]][["ccp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["topgo_cc"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["topgo_cc"]]) + 2,
                                            start_row=new_row, plotname="topgo_cc",
                                            savedir=excel_basename)
            }
        }

        ## Move to gostats
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:11, widths="auto")
        new_row <- 1
        sheet <- "gostats"
        openxlsx::addWorksheet(wb, sheetName=sheet)

        ## gs down bp
        if (!is.null(down_stuff[["gostats_bp"]])) {
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["gostats_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_gostats"]][[name]][["pvalue_plots"]][["bp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["gostats_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["gostats_bp"]]) + 2,
                                            start_row=new_row, plotname="gotats_bp",
                                            savedir=excel_basename)
            }
            links <- down_stuff[["gostats_bp"]][["Link"]]
            class(links) <- 'hyperlink'
            names(links) <- down_stuff[["gostats_bp"]][["Category"]]
            openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
            new_row <- new_row + nrow(down_stuff[["gostats_bp"]]) + 2
        }

        ## gs down mf
        if (!is.null(down_stuff[["gostats_mf"]])) {
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["gostats_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_gostats"]][[name]][["pvalue_plots"]][["mf_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["gostats_mf"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["gostats_mf"]]) + 2,
                                            start_row=new_row, plotname="gostats_mf",
                                            savedir=excel_basename)
            }
            links <- down_stuff[["gostats_mf"]][["Link"]]
            class(links) <- 'hyperlink'
            names(links) <- down_stuff[["gostats_mf"]][["Category"]]
            openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
            new_row <- new_row + nrow(down_stuff[["gostats_mf"]]) + 2
        }

        ## gs down cc
        if (!is.null(down_stuff[["gostats_cc"]])) {
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["gostats_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_gostats"]][[name]][["pvalue_plots"]][["cc_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["gostats_cc"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["gostats_cc"]]) + 2,
                                            start_row=new_row, plotname="gostats_cc",
                                            savedir=excel_basename)
            }
            links <- down_stuff[["gostats_cc"]][["Link"]]
            class(links) <- 'hyperlink'
            names(links) <- down_stuff[["gostats_cc"]][["Category"]]
            openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
        }
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")

        ## Write gprofiler data
        new_row <- 1
        sheet <- "gprofiler"
        openxlsx::addWorksheet(wb, sheetName=sheet)

        ## Write gprofiler BP data
        if (!is.null(down_stuff[["gprofiler_bp"]])) {
            openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["gprofiler_bp"]],
                                     tableStyle=table_style, startRow=new_row)
            ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_gprofiler"]][[name]][["pvalue_plots"]][["bpp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["gprofiler_bp"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["gprofiler_bp"]]) + 2,
                                            start_row=new_row, plotname="gprofiler_bp",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(down_stuff[["gprofiler_bp"]]) + 2
        }

        ## write gprofiler MF data
        if (!is.null(down_stuff[["gprofiler_mf"]])) {
            openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["gprofiler_mf"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_gprofiler"]][[name]][["pvalue_plots"]][["mfp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["gprofiler_mf"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["gprofiler_mf"]]) + 2,
                                            start_row=new_row, plotname="gprofiler_mf",
                                            savedir=excel_basename)
            }
            new_row <- new_row + nrow(down_stuff[["gprofiler_mf"]]) + 2
        }

        ## write gprofiler CC data
        if (!is.null(down_stuff[["gprofiler_cc"]])) {
            openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
            openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
            new_row <- new_row + 1
            openxlsx::writeDataTable(wb, sheet, x=down_stuff[["gprofiler_cc"]],
                                     tableStyle=table_style, startRow=new_row)
            if (isTRUE(add_plots)) {
                a_plot <- kept_ontology[["down_gprofiler"]][[name]][["pvalue_plots"]][["ccp_plot_over"]]
                ## print(a_plot)
                ## openxlsx::insertPlot(wb, sheet, width=6, height=6,
                ##                      startCol=ncol(down_stuff[["gprofiler_cc"]]) + 2,
                ##                      startRow=new_row, fileType="png", units="in")
                try_result <- xlsx_plot_png(a_plot, wb=wb, sheet=sheet, width=6, height=6,
                                            start_col=ncol(down_stuff[["gprofiler_cc"]]) + 2,
                                            start_row=new_row, plotname="gprofiler_cc",
                                            savedir=excel_basename)
            }
            openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
        }


        res <- openxlsx::saveWorkbook(wb, down_filename, overwrite=TRUE)
    }  ## End of name_list
}

#' Write gene ontology tables for excel
#'
#' Combine the results from goseq, cluster profiler, topgo, and gostats and drop them into excel
#' Hopefully with a relatively consistent look.
#'
#' @param goseq  The goseq result from simple_goseq()
#' @param cluster The result from simple_clusterprofiler()
#' @param topgo  Guess
#' @param gostats  Yep, ditto
#' @param gprofiler  woo hoo!
#' @param file   the file to save the results.
#' @param dated   date the excel file
#' @param n   the number of ontology categories to include in each table.
#' @param overwritefile   overwrite an existing excel file
#' @return the list of ontology information
#' @seealso \pkg{openxlsx} \pkg{goseq} \pkg{clusterProfiler} \pkg{goStats} \pkg{topGO} \pkg{gProfiler}
#' @export
write_go_xls <- function(goseq, cluster, topgo, gostats, gprofiler, file="excel/merged_go",
                         dated=TRUE, n=30, overwritefile=TRUE) {
    n <- get0('n')
    if (is.null(n)) {
        n <- 30
    }
    file <- get0('file')
    if (is.null(file)) {
        file <- "excel/merged_go"
    }
    excel_dir <- dirname(file)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }

    suffix <- ".xlsx"
    file <- gsub(pattern="\\.xlsx", replacement="", file, perl=TRUE)
    file <- gsub(pattern="\\.xls", replacement="", file, perl=TRUE)
    filename <- NULL
    if (isTRUE(dated)) {
        timestamp <- format(Sys.time(), "%Y%m%d%H")
        filename <- paste0(file, "-", timestamp, suffix)
    } else {
        filename <- paste0(file, suffix)
    }

    if (file.exists(filename)) {
        if (isTRUE(overwritefile)) {
            backup_file(filename)
        }
    }

    ## Massage the goseq tables to match Najib's request
    goseq_mf <- head(goseq[["mf_subset"]], n=n)
    goseq_bp <- head(goseq[["bp_subset"]], n=n)
    goseq_cc <- head(goseq[["cc_subset"]], n=n)
    goseq_mf <- goseq_mf[, c(7,1,6,2,4,5,8)]
    goseq_bp <- goseq_bp[, c(7,1,6,2,4,5,8)]
    goseq_cc <- goseq_cc[, c(7,1,6,2,4,5,8)]
    colnames(goseq_mf) <- c("Ontology", "Category", "Term", "Over p-value",
                            "Num. DE", "Num. in cat.", "Q-value")
    colnames(goseq_bp) <- c("Ontology", "Category", "Term", "Over p-value",
                            "Num. DE", "Num. in cat.", "Q-value")
    colnames(goseq_cc) <- c("Ontology", "Category", "Term", "Over p-value",
                            "Num. DE", "Num. in cat.", "Q-value")

    ## Massage the clusterProfiler tables similarly
    cluster_mf <- head(as.data.frame(cluster[["mf_all"]]@result), n=n)
    cluster_bp <- head(as.data.frame(cluster[["bp_all"]]@result), n=n)
    cluster_cc <- head(as.data.frame(cluster[["cc_all"]]@result), n=n)
    cluster_mf[["geneID"]] <- gsub(cluster_mf[["geneID"]], pattern="/", replacement=" ")
    cluster_bp[["geneID"]] <- gsub(cluster_bp[["geneID"]], pattern="/", replacement=" ")
    cluster_cc[["geneID"]] <- gsub(cluster_cc[["geneID"]], pattern="/", replacement=" ")
    cluster_mf[["ontology"]] <- "MF"
    cluster_bp[["ontology"]] <- "BP"
    cluster_cc[["ontology"]] <- "CC"
    cluster_mf <- cluster_mf[,c(10, 1, 2, 5, 3, 4, 6, 7, 9, 8)]
    cluster_bp <- cluster_bp[,c(10, 1, 2, 5, 3, 4, 6, 7, 9, 8)]
    cluster_cc <- cluster_cc[,c(10, 1, 2, 5, 3, 4, 6, 7, 9, 8)]
    colnames(cluster_mf) <- c("Ontology", "Category", "Term", "Over p-value", "Gene ratio",
                              "BG ratio", "Adj. p-value", "Q-value", "Count", "Genes")
    colnames(cluster_bp) <- c("Ontology", "Category", "Term", "Over p-value", "Gene ratio",
                              "BG ratio", "Adj. p-value", "Q-value", "Count", "Genes")
    colnames(cluster_cc) <- c("Ontology", "Category", "Term", "Over p-value", "Gene ratio",
                              "BG ratio", "Adj. p-value", "Q-value", "Count", "Genes")

    ## Now do the topgo data
    topgo_mf <- head(topgo[["tables"]][["mf_interesting"]], n=n)
    topgo_bp <- head(topgo[["tables"]][["bp_interesting"]], n=n)
    topgo_cc <- head(topgo[["tables"]][["cc_interesting"]], n=n)
    topgo_mf <- topgo_mf[,c(2, 1, 11, 6, 7, 8, 9, 10, 4, 3, 5)]
    topgo_bp <- topgo_bp[,c(2, 1, 11, 6, 7, 8, 9, 10, 4, 3, 5)]
    topgo_cc <- topgo_cc[,c(2, 1, 11, 6, 7, 8, 9, 10, 4, 3, 5)]
    colnames(topgo_mf) <- c("Ontology", "Category", "Term", "Fisher p-value", "Q-value", "KS score",
                            "EL score", "Weight score", "Num. DE", "Num. in cat.", "Exp. in cat.")
    colnames(topgo_bp) <- c("Ontology", "Category", "Term", "Fisher p-value", "Q-value", "KS score",
                            "EL score", "Weight score", "Num. DE", "Num. in cat.", "Exp. in cat.")
    colnames(topgo_cc) <- c("Ontology", "Category", "Term", "Fisher p-value", "Q-value", "KS score",
                            "EL score", "Weight score", "Num. DE", "Num. in cat.", "Exp. in cat.")

    ## And the gostats data
    gostats_mf <- head(gostats[["mf_over_all"]], n=n)
    gostats_bp <- head(gostats[["bp_over_all"]], n=n)
    gostats_cc <- head(gostats[["cc_over_all"]], n=n)
    gostats_mf[["t"]] <- gsub(gostats_mf[["Term"]], pattern=".*\">(.*)</a>", replacement="\\1")
    gostats_bp[["t"]] <- gsub(gostats_bp[["Term"]], pattern=".*\">(.*)</a>", replacement="\\1")
    gostats_cc[["t"]] <- gsub(gostats_cc[["Term"]], pattern=".*\">(.*)</a>", replacement="\\1")
    gostats_mf[["Term"]] <- gsub(gostats_mf[["Term"]], pattern="<a href=\"(.*)\">.*", replacement="\\1")
    gostats_bp[["Term"]] <- gsub(gostats_bp[["Term"]], pattern="<a href=\"(.*)\">.*", replacement="\\1")
    gostats_cc[["Term"]] <- gsub(gostats_cc[["Term"]], pattern="<a href=\"(.*)\">.*", replacement="\\1")
    gostats_mf[["ont"]] <- "MF"
    gostats_bp[["ont"]] <- "BP"
    gostats_cc[["ont"]] <- "CC"
    gostats_mf <- gostats_mf[,c(10, 1, 9, 2, 5, 6, 3, 4, 8, 7)]
    gostats_bp <- gostats_bp[,c(10, 1, 9, 2, 5, 6, 3, 4, 8, 7)]
    gostats_cc <- gostats_cc[,c(10, 1, 9, 2, 5, 6, 3, 4, 8, 7)]
    colnames(gostats_mf) <- c("Ontology", "Category", "Term", "Fisher p-value", "Num. DE",
                              "Num. in cat.", "Odds ratio", "Exp. in cat.", "Q-value", "Link")
    colnames(gostats_bp) <- c("Ontology", "Category", "Term", "Fisher p-value", "Num. DE",
                              "Num. in cat.", "Odds ratio", "Exp. in cat.", "Q-value", "Link")
    colnames(gostats_cc) <- c("Ontology", "Category", "Term", "Fisher p-value", "Num. DE",
                              "Num. in cat.", "Odds ratio", "Exp. in cat.", "Q-value", "Link")

    lst <- list("goseq_mf" = goseq_mf,
                "goseq_bp" = goseq_bp,
                "goseq_cc" = goseq_cc,
                "cluster_mf" = cluster_mf,
                "cluster_bp" = cluster_bp,
                "cluster_cc" = cluster_cc,
                "topgo_mf" = topgo_mf,
                "topgo_bp" = topgo_bp,
                "topgo_cc" = topgo_cc,
                "gostats_mf" = gostats_mf,
                "gostats_bp" = gostats_bp,
                "gostats_cc" = gostats_cc)
    ## require.auto("awalker89/openxlsx")
    wb <- openxlsx::createWorkbook(creator="atb")
    hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT",
                                 textDecoration="bold", border="Bottom", fontSize="30")

    ## This stanza will be repeated so I am just incrementing the new_row
    new_row <- 1
    sheet <- "goseq"
    openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["goseq_mf"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst[["goseq_mf"]]) + 2
    openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["goseq_cc"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")

    new_row <- 1
    sheet <- "clusterProfiler"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["cluster_bp"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst[["cluster_bp"]]) + 2
    openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["cluster_mf"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst[["cluster_mf"]]) + 2
    openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["cluster_cc"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")

    new_row <- 1
    sheet <- "topgo"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["topgo_bp"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst[["topgo_bp"]]) + 2
    openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["topgo_mf"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst[["topgo_mf"]]) + 2
    openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["topgo_cc"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:11, widths="auto")

    new_row <- 1
    sheet <- "gostats"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["gostats_bp"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    links <- lst[["gostats_bp"]][["Link"]]
    class(links) <- 'hyperlink'
    names(links) <- lst[["gostats_bp"]][["Category"]]
    openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
    new_row <- new_row + nrow(lst[["gostats_bp"]]) + 2
    openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["gostats_mf"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    links <- lst[["gostats_mf"]][["Link"]]
    class(links) <- 'hyperlink'
    names(links) <- lst[["gostats_mf"]][["Category"]]
    openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
    new_row <- new_row + nrow(lst[["gostats_mf"]]) + 2
    openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst[["gostats_cc"]],
                             tableStyle="TableStyleMedium9", startRow=new_row)
    links <- lst[["gostats_cc"]][["Link"]]
    class(links) <- 'hyperlink'
    names(links) <- lst[["gostats_cc"]][["Category"]]
    openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")

    res <- openxlsx::saveWorkbook(wb, file, overwrite=TRUE)
    return(res)
}

## EOF
