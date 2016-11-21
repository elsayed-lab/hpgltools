#' Make a pretty table of goseq data in excel.
#'
#' It is my intention to make a function like this for each ontology tool in my repetoire
#'
#' @param goseq A set of results from simple_goseq().
#' @param file An excel file to which to write some pretty results.
#' @param pval Choose a cutoff for reporting by p-value.
#' @param add_plots Include some pvalue plots in the excel output?
#' @return The result from openxlsx
#' @export
write_goseq_data <- function(goseq, excel="excel/goseq.xlsx", wb=NULL,
                             pval=0.1, add_plots=TRUE, height=15, width=10, ...) {
    arglist <- list(...)
    table_style <- "TableStyleMedium9"
    if (!is.null(arglist[["table_style"]])) {
        table_style <- arglist[["TableStyleMedium9"]]
    }
    excel_dir <- dirname(excel)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")
    }

    if (class(excel) == "character") {
        message("Writing a legend of columns.")
        wb <- openxlsx::createWorkbook(creator="hpgltools")
        legend <- data.frame(rbind(
            c("Ontology", "Which portion of the ontology tree is being examined?  Molecular Function, Biological Process, or Cellular Component."),
            c("Category", "Gene ontology Identifier."),
            c("Term", "Short definition of the category."),
            c("Over p-value", "Estimate that the set of genes provided to goseq is over-represented in the row's ontology category."),
            c("Q-value", "False discovery rate correction of the p-value."),
            c("DE genes in cat", "What genes provided to the ontology search are inside this specific ontology category?"),
            c("All genes in cat", "The full set of gene annotations included in this ontology category."),
            c("Num. de", "The number of genes in column 'F'."),
            c("Num. in cat", "The number of genes in column 'G'.")
        ))
        colnames(legend) <- c("column name", "column definition")
        xls_result <- write_xls(wb, data=legend, sheet="legend", rownames=FALSE,
                                title="Columns used in the following tables.")
    }  ## End making sure that an excel is desired.

    ## Pull out the relevant portions of the goseq data
    ## For this I am using the same (arbitrary) rules as in gather_goseq_genes()
    goseq_mf <- goseq[["mf_subset"]]
    goseq_mf <- goseq_mf[ goseq_mf[["over_represented_pvalue"]] <= pval, ]
    goseq_mf_genes <- gather_goseq_genes(goseq, ontology="MF", pval=pval)
    mf_genes <- as.data.frame(goseq_mf_genes)
    rownames(mf_genes) <- rownames(goseq_mf_genes)
    goseq_mf <- merge(goseq_mf, mf_genes, by="row.names")
    rownames(goseq_mf) <- goseq_mf[["Row.names"]]
    goseq_mf <- goseq_mf[-1]
    mf_idx <- order(goseq_mf[["qvalue"]])
    goseq_mf <- goseq_mf[mf_idx, ]

    goseq_bp <- goseq[["bp_subset"]]
    goseq_bp <- goseq_bp[ goseq_bp[["over_represented_pvalue"]] <= pval, ]
    goseq_bp_genes <- gather_goseq_genes(goseq, ontology="BP", pval=pval)
    bp_genes <- as.data.frame(goseq_bp_genes)
    rownames(bp_genes) <- rownames(goseq_bp_genes)
    goseq_bp <- merge(goseq_bp, bp_genes, by="row.names")
    rownames(goseq_bp) <- goseq_bp[["Row.names"]]
    goseq_bp <- goseq_bp[-1]
    bp_idx <- order(goseq_bp[["qvalue"]])
    goseq_bp <- goseq_bp[bp_idx, ]

    goseq_cc <- goseq[["cc_subset"]]
    goseq_cc <- goseq_cc[ goseq_cc[["over_represented_pvalue"]] <= pval, ]
    goseq_cc_genes <- gather_goseq_genes(goseq, ontology="CC", pval=pval)
    cc_genes <- as.data.frame(goseq_cc_genes)
    rownames(cc_genes) <- rownames(goseq_cc_genes)
    goseq_cc <- merge(goseq_cc, cc_genes, by="row.names")
    rownames(goseq_cc) <- goseq_cc[["Row.names"]]
    goseq_cc <- goseq_cc[-1]
    cc_idx <- order(goseq_cc[["qvalue"]])
    goseq_cc <- goseq_cc[cc_idx, ]

    goseq_mf <- goseq_mf[, c(7,1,6,2,8,10,9,4,5)]
    goseq_bp <- goseq_bp[, c(7,1,6,2,8,10,9,4,5)]
    goseq_cc <- goseq_cc[, c(7,1,6,2,8,10,9,4,5)]
    colnames(goseq_mf) <- c("Ontology","Category","Term","Over p-value", "Q-value",
                            "DE genes in cat", "All genes in cat", "Num. DE", "Num. in cat.")
    colnames(goseq_bp) <- c("Ontology","Category","Term","Over p-value", "Q-value",
                            "DE genes in cat", "All genes in cat", "Num. DE", "Num. in cat.")
    colnames(goseq_cc) <- c("Ontology","Category","Term","Over p-value", "Q-value",
                            "DE genes in cat", "All genes in cat", "Num. DE", "Num. in cat.")

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
        a_plot <- goseq[["pvalue_plots"]][["bpp_plot_over"]]
        tt <- try(print(a_plot), silent=TRUE)
        openxlsx::insertPlot(wb, sheet, width=width, height=height,
                             startCol=ncol(goseq_bp) + 2, startRow=new_row, fileType="png", units="in")
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
        a_plot <- goseq[["pvalue_plots"]][["mfp_plot_over"]]
        tt <- try(print(a_plot), silent=TRUE)
        openxlsx::insertPlot(wb, sheet, width=width, height=height,
                             startCol=ncol(goseq_mf) + 2, startRow=new_row, fileType="png", units="in")
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
        a_plot <- goseq[["pvalue_plots"]][["ccp_plot_over"]]
        tt <- try(print(a_plot), silent=TRUE)
        openxlsx::insertPlot(wb, sheet, width=width, height=height, startCol=ncol(goseq_cc) + 2, startRow=new_row, fileType="png", units="in")
    }
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
    openxlsx::setColWidths(wb, sheet=sheet, cols=6:7, widths=30)
    new_row <- new_row + nrow(goseq_cc) + 2

    res <- openxlsx::saveWorkbook(wb, excel, overwrite=TRUE)
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
#' @param ... More options, not currently used I think.
#' @export
write_gprofiler_data <- function(gprofiler_result, wb=NULL, excel="excel/gprofiler_result.xlsx",
                                 add_plots=TRUE, height=15, width=10, ...) {
    arglist <- list(...)
    table_style <- "TableStyleMedium9"
    if (!is.null(arglist[["table_style"]])) {
        table_style <- arglist[["TableStyleMedium9"]]
    }

    excel_dir <- dirname(excel)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")
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

    if (!is.null(gprofiler_result[["go"]]) & nrow(gprofiler_result[["go"]]) > 0) {
        new_row <- 1
        sheet <- "GO"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        go_data <- gprofiler_result[["go"]]
        bp_data <- go_data[go_data[["domain"]]=="BP", ]
        mf_data <- go_data[go_data[["domain"]]=="MF", ]
        cc_data <- go_data[go_data[["domain"]]=="CC", ]

        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=bp_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["bpp_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(bp_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(bp_data) + 2

        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=mf_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["mfp_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(mf_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(mf_data) + 2

        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=cc_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["ccp_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(cc_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(cc_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking if go data is null

    if (!is.null(gprofiler_result[["kegg"]]) & nrow(gprofiler_result[["kegg"]]) > 0) {
        new_row <- 1
        sheet <- "KEGG"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        kegg_data <- gprofiler_result[["kegg"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=kegg_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["kegg_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(kegg_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(kegg_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking KEGG data


    if (!is.null(gprofiler_result[["tf"]]) & nrow(gprofiler_result[["tf"]]) > 0) {
        new_row <- 1
        sheet <- "transcription_factor"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        tf_data <- gprofiler_result[["tf"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=tf_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["tf_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(tf_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(tf_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking tf data

    if (!is.null(gprofiler_result[["reactome"]]) & nrow(gprofiler_result[["reactome"]]) > 0) {
        new_row <- 1
        sheet <- "reactome"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        react_data <- gprofiler_result[["reactome"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=react_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["reactome_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(react_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(react_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking reactome data


    if (!is.null(gprofiler_result[["mi"]]) & nrow(gprofiler_result[["mi"]]) > 0) {
        new_row <- 1
        sheet <- "mirna"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        mi_data <- gprofiler_result[["mi"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=mi_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["mi_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(mi_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(mi_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking mirna data

    if (!is.null(gprofiler_result[["hp"]]) & nrow(gprofiler_result[["hp"]]) > 0) {
        new_row <- 1
        sheet <- "hp"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        hp_data <- gprofiler_result[["hp"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=hp_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["hp_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(hp_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(hp_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking corum data

    if (!is.null(gprofiler_result[["corum"]]) & nrow(gprofiler_result[["corum"]]) > 0) {
        new_row <- 1
        sheet <- "corum"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        corum_data <- gprofiler_result[["corum"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=corum_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["corum_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(corum_data) + 2, startRow=new_row, fileType="png", units="in"))
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
    return(excel_ret)
}



## EOF
