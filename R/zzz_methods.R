## Make some methods which call Biobase on expts.
## This way we don't have to do Biobase::exprs/fData/pData/notes()
## But instead R will call them for us.

## Here is a note from Hadley which is relevant here:
## Another consideration is that
## order. For example, to define the method setMethod("foo", c("bar", "baz"),
## ...) you must already have created the foo generic and the two classes. By
## default, R code is loaded in alphabetical order, but that won’t always work
## for your situation.

setGeneric("backup_expression_data", signature = signature(expt = "expt"),
           function(expt) standardGeneric("backup_expression_data"))

setGeneric("colors", signature = signature(expt = "expt"),
           function(expt) standardGeneric("colors"))

setGeneric("colors<-", signature = signature(expt = "expt"),
           function(expt, ...) standardGeneric("colors<-", ...))

setGeneric("extract_keepers", signature = c("extracted", "keepers"),
           function(extracted, keepers, ...) standardGeneric("extract_keepers"))

setGeneric("get_backup_expression_data", signature = c("expt"),
           function(expt) standardGeneric("get_backup_expression_data"))

#' Generic method to get colors from expression data.
#'
#' @param object The object from which to gather colors.
#' @param ... Additional arguments passed to function constructors.
#' @return colors!
#' @rdname methods
#' @export
setGeneric("getColors", signature = signature(expt = "expt"),
           function(expt) standardGeneric("getColors"))

#' Generic method to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setGeneric("iDA", signature = c("object"),
           function(object, ...) standardGeneric("iDA"))

setGeneric("normalizeData",
  function(expt, ...) standardGeneric("normalizeData"),
  signature = signature(expt = "expt"))

setGeneric("state",
  function(expt) standardGeneric("state"),
  signature = signature(expt = "expt"))

#setGeneric("subset_expt", signature = c("expt"),
#           function(expt, ...) standardGeneric("subset_expt"))

#' A getter for the annotation databased used to create an expt/se.
#'
#' @param object One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @importFrom Biobase annotation
#' @export
setMethod(
  "annotation", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::annotation(object[["expressionset"]])
  })

#' A setter for the annotation databased used to create an expt/se.
#'
#' @param object One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @importFrom Biobase annotation<-
#' @export
setMethod(
  "annotation<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    annotation(object[["expressionset"]]) <- value
    return(object)
  })

#' A getter to pull the assay data from an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @importFrom SummarizedExperiment assay
#' @export
setMethod(
  "assay", signature = signature(x = "expt"),
  definition = function(x, withDimnames = TRUE, ...) {
    mtrx <- Biobase::exprs(x[["expressionset"]])
    return(mtrx)
  })

#' A setter to put the assay data into an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @importFrom SummarizedExperiment assay<-
#' @export
setMethod(
  "assay<-", signature = signature(x = "expt"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::exprs(x[["expressionset"]]) <- value
    return(x)
  })

#' A getter to pull the assay data from an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @export
setMethod(
  "assay", signature = signature(x = "ExpressionSet"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::exprs(x)
  })

#' A setter to put the assay data into an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @export
setMethod(
  "assay<-", signature = signature(x = "ExpressionSet"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::exprs(x) <- value
    return(x)
  })

#' Backup the state of an expressionSet.
#'
#' @param expt An ExpressionSet.
#' @export
setMethod(
  "backup_expression_data", signature = signature(expt = "ExpressionSet"),
  definition = function(expt) {
    message("I do not know a good way to backup the expressionset data, returning.")
    return(expt)
  })

#' Backup the state of an SummarizedExperiment.
#'
#' @param expt A SummarizedExperiment.
#' @export
setMethod(
  "backup_expression_data", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt) {
    backup <- expt
    metadata(expt)[["original_se"]] <- backup
    return(expt)
  })

#' A getter to pull the sample data from an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
  "colData", signature = signature(x = "expt"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::pData(x[["expressionset"]])
  })

#' A setter to put the sample data into an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @importFrom SummarizedExperiment colData<-
#' @export
setMethod(
  "colData<-", signature = signature(x = "expt"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::pData(x[["expressionset"]]) <- value
    return(x)
  })

#' A getter to pull the sample data from an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @export
setMethod(
  "colData", signature = signature(x = "ExpressionSet"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::pData(x)
  })

#' A setter to put the sample data into an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @export
setMethod(
  "colData<-", signature = signature(x = "ExpressionSet"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::pData(x) <- value
    return(x)
  })

#' A getter to pull the colors from an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "colors", signature = signature(expt = "expt"),
  definition = function(expt) {
    expt[["colors"]]
  })

#' A getter to pull the colors from a SummarizedExperiment.
#'
#' @param expt An expt.
#' @export
setMethod(
  "colors", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt) {
    metadata(expt)[["colors"]]
  })

#' A setter to put the colors into a SummarizedExperiment.
#'
#' @param expt A SummarizedExperiment.
#' @export
setMethod(
  "colors<-", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, lst) {
    metadata(expt)[["colors"]] <- lst
    return(expt)
  })

#' A setter to put the colors into an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "colors<-", signature = signature(expt = "expt"),
  definition = function(expt, lst) {
    expt[["colors"]] <- lst
    return(expt)
  })

#' A getter to pull the expression data from an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "exprs", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::exprs(object[["expressionset"]])
  })

#' A setter to put the expression data into an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "exprs<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    testthat::expect_equal(colnames(exprs(object)),
                           colnames(value))
    if (class(value)[1] == "data.frame") {
      value <- as.matrix(value)
    }
    exprs(object[["expressionset"]]) <- value
    return(object)
  })

#' A setter to put the expression data into an expt.
#'
#' @param object ExpressionSet to modify.
#' @param value New expression data.
#' @export
setMethod(
  "exprs<-", signature = signature(object = "ExpressionSet", value = "data.frame"),
  definition = function(object, value) {
    testthat::expect_equal(colnames(exprs(object)),
                           colnames(value))
    object <- as.matrix(exprs(value))
    exprs(object) <- value
    return(object)
  })

#' A getter to pull the expression data from a SummarizedExperiment.
#'
#' @param expt A SummarizedExperiment.
#' @export
setMethod(
  "exprs", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    SummarizedExperiment::assay(object)
  })

#' A setter to put the expression data to a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @export
setMethod(
  "exprs<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, value) {
    testthat::expect_equal(colnames(exprs(object)),
                           colnames(value))
    if (class(value)[1] == "data.frame") {
      value <- as.matrix(value)
    }
    SummarizedExperiment::assay(object) <- value
    return(object)
  })

#' Pull the set of keepers from a character vector.
#'
#' Instead of a list of numerators/denominators, one might feed combine_de_tables a vector of
#' things like: 'a_vs_b'.
#' @param keepers Character vector of keepers.
#' @export
setMethod(
  "extract_keepers", signature = signature(keepers = "character"),
  definition = function(extracted, keepers, table_names,
                        all_coefficients,
                        limma, edger, ebseq, deseq, basic, noiseq,
                        adjp, annot_df,
                        include_deseq, include_edger,
                        include_ebseq, include_limma,
                        include_basic, include_noiseq,
                        excludes, padj_type,
                        fancy = FALSE, loess = FALSE,
                        lfc_cutoff = 1.0, p_cutoff = 0.05,
                        sheet_prefix = NULL, sheet_number = NULL,
                        format_sig = 4, plot_colors = plot_colors,
                        z = 1.5, alpha = 0.4, z_lines = FALSE,
                        label = 10, label_column = "hgncsymbol") {
    if (keepers[1] == "all") {
      new_keepers <- list()
      names_length <- length(table_names)
      numerators <- denominators <- c()
      for (a in seq_len(names_length)) {
        name <- table_names[a]
        splitted <- strsplit(x = name, split = "_vs_")
        denominator <- splitted[[1]][2]
        numerator <- splitted[[1]][1]
        new_keepers[[name]] <- c(numerator, denominator)
      }
    } else {
      splitted <- strsplit(x = keepers, split = "_vs_")
      numerator <- splitted[[1]][1]
      denominator <- splited[[1]][2]
      new_keepers <- list(splitted = c(numerator, denominator))
    }
    extract_keepers(extracted, new_keepers, table_names,
                    all_coefficients,
                    limma, edger, ebseq, deseq, basic, noiseq,
                    adjp, annot_df,
                    include_deseq, include_edger,
                    include_ebseq, include_limma,
                    include_basic, include_noiseq,
                    excludes, padj_type,
                    fancy = FALSE, loess = FALSE,
                    lfc_cutoff = 1.0, p_cutoff = 0.05,
                    sheet_prefix = NULL, sheet_number = NULL,
                    format_sig = 4, plot_colors = plot_colors,
                    z = 1.5, alpha = 0.4, z_lines = FALSE,
                    label = 10, label_column = "hgncsymbol")
  })

#' A getter to pull the gene annotation data from an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "fData", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::fData(object[["expressionset"]])
  })

#' A setter to put the gene annotation data into an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "fData<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    testthat::expect_equal(rownames(fData(object)),
                           rownames(value))
    fData(object[["expressionset"]]) <- value
    return(object)
  })

#' A getter to pull the gene annotation data from a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @export
setMethod(
  "fData", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    SummarizedExperiment::rowData(object)
  })

#' A setter to put the gene annotation data into a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @export
setMethod(
  "fData<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, value) {
    testthat::expect_equal(rownames(fData(object)),
                           rownames(value))
    SummarizedExperiment::rowData(object) <- value
    return(object)
  })

#' Get the backup data from a SummarizedExperiment.
#'
#' @param expt An expt containing the backup data.
#' @export
setMethod(
  "get_backup_expression_data", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt) {
    backup <- metadata(expt)[["original_se"]]
    return(backup)
  })

#' Get the backup data from an ExpressionSet.
#'
#' @param expt An ExpressionSet does not contain backup data.
#' @export
setMethod(
  "get_backup_expression_data", signature = signature(expt = "ExpressionSet"),
  definition = function(expt) {
    message("The expressionset currently does not keep a backup.")
    return(expt)
  })

#' Set method for matrix to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setMethod(
  "iDA", signature = signature(object = "matrix"),
  function(object, ...) {
    iDAoutput <- iDA::iDA_core(object, ...)
    return(iDAoutput)
  })

#setMethod(
#  "overlap_groups", signature = signature(input_mtrx = "upset", sort = "logical"),
#  definition = function(input_mtrx, sort = sort) {
#    new_mtrx <- input_mtrx == 1
#    overlap_groups(new_mtrx, sort = sort)
#  })
#
#setMethod(
#  "overlap_groups", signature = signature(input_mtrx = "list", sort = "logical"),
#  definition = function(input_mtrx, sort = sort) {
#    input_upset <- UpSetR::fromList(input_mtrx)
#    new_mtrx <- input_upset == 1
#    overlap_groups(new_mtrx, sort = sort)
#  })

setMethod(
  "normalizeData", signature = signature(expt = "expt"),
  definition = function(expt, ...) {
    normalize_expt(expt, ...)
  })

setMethod(
  "normalizeData", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, ...) {
    se <- expt
    normalize_se(se, transform = transform, norm = norm,
                 convert = convert, batch = batch, filter = filter,
                 annotations = annotations, fasta = fasta, entry_type = entry_type,
                 use_original = use_original, batch1 = batch1, batch2 = batch2,
                 batch_step = batch_step, low_to_zero = low_to_zero, thresh = thresh,
                 min_samples = min_samples, p = p, A = A, k = k, cv_min = cv_min,
                 cv_max = cv_max, na_to_zero = na_to_zero,
                 adjust_method = adjust_method, verbose = verbose, ...)
  })

#' A getter to pull the notes an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "notes", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::notes(object[["expressionset"]])
  })

#' A getter to pull the experimental metadata from an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "pData", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::pData(object[["expressionset"]])
  })

#' A setter to put the experimental metadata into an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "pData<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    testthat::expect_equal(rownames(pData(object)),
                           rownames(value))
    pData(object[["expressionset"]]) <- value
    return(object)
  })

#' A getter to pull the experimental metadata from a SummarizedExperiment.
#'
#' This is essentially synonymous with colData, except I cannot seem
#' to remember that function when I am working; so I just added
#' another signature to pData.
#'
#' @param object An expt.
#' @export
setMethod(
  "pData", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    SummarizedExperiment::colData(object)
  })

#' A setter to put the experimental metadata into a SummarizedExperiment.
#'
#' This is essentially synonymous with colData, except I cannot seem
#' to remember that function when I am working; so I just added
#' another signature to pData.
#'
#' @param object An expt.
#' @export
setMethod(
  "pData<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, value) {
    testthat::expect_equal(
      rownames(SummarizedExperiment::colData(object)),
      rownames(value))
    SummarizedExperiment::colData(object) <- value
    return(object)
  })

#' Run plot_heatmap with an expt as input.
#' @export
setMethod(
  "plot_heatmap", signature = signature(expt_data = "expt"),
  definition = function(expt_data, expt_colors = NULL, expt_design = NULL, method = "pearson",
                        expt_names = NULL, type = "correlation", batch_row = "batch",
                        plot_title = NULL, label_chars = 10, ...) {
    expt_design <- pData(expt_data)
    expt_colors <- expt_data[["colors"]]
    expt_names <- expt_data[["expt_names"]]
    expt_data <- exprs(expt_data)
    ## If plot_title is NULL, print nothing, if it is TRUE
    ## Then give some information about what happened to the data to make the plot.
    ## I tried foolishly to put this in plot_pcs(), but there is no way that receives
    ## my expt containing the normalization state of the data.
    if (isTRUE(plot_title)) {
      plot_title <- what_happened(expt_data)
    } else if (!is.null(plot_title)) {
      data_title <- what_happened(expt_data)
      plot_title <- glue("{plot_title}; {data_title}")
    } else {
      ## Leave the title blank.
    }
    plot_heatmap(expt_data, expt_colors = expt_colors, expt_design = expt_design,
      method = method, expt_names = expt_names, type = type,
      batch_row = batch_row, plot_title = plot_title,
      label_chars = label_chars, ...)
  })

#' Run plot_heatmap with a SummarizedExperiment as input.
#' @export
setMethod(
  "plot_heatmap", signature = signature(expt_data = "SummarizedExperiment"),
  definition = function(expt_data, expt_colors = NULL, expt_design = NULL, method = "pearson",
                        expt_names = NULL, type = "correlation", batch_row = "batch",
                        plot_title = NULL, label_chars = 10, ...) {
    expt_design <- pData(expt_data)
    expt_colors <- metadata(expt_data)[["colors"]]
    expt_names <- metadata(expt_data)[["expt_names"]]
    expt_data <- exprs(expt_data)
    if (isTRUE(plot_title)) {
      plot_title <- what_happened(expt_data)
    } else if (!is.null(plot_title)) {
      data_title <- what_happened(expt_data)
      plot_title <- glue("{plot_title}; {data_title}")
    } else {
      ## Leave the title blank.
    }
    plot_heatmap(expt_data, expt_colors = expt_colors, expt_design = expt_design,
      method = method, expt_names = expt_names, type = type,
      batch_row = batch_row, plot_title = plot_title,
      label_chars = label_chars, ...)
  })

#' Run plot_heatmap with a dataframe as input.
#' @export
setMethod(
  "plot_heatmap", signature = signature(expt_data = "data.frame"),
  definition = function(expt_data, expt_colors = NULL, expt_design = NULL,
                        method = "pearson", expt_names = NULL, type = "correlation",
                        batch_row = "batch", plot_title = NULL, label_chars = 10, ...) {
    expt_mtrx <- as.matrix(expt_data)
    plot_heatmap(expt_mtrx, expt_colors = expt_colors, expt_design = expt_design,
                 method = method, expt_names = expt_names, type = type,
                 batch_row = batch_row, plot_title = plot_title,
                 label_chars = label_chars, ...)
  })

#' Run plot_heatmap with an ExpressionSet as input.
#' @export
setMethod(
  "plot_heatmap", signature = signature(expt_data = "ExpressionSet"),
  definition = function(expt_data, expt_colors = NULL, expt_design = NULL, method = "pearson",
                        expt_names = NULL, type = "correlation", batch_row = "batch",
                        plot_title = NULL, label_chars = 10, ...) {
    expt_design <- pData(expt_data)
    expt_mtrx <- exprs(expt_data)
    plot_heatmap(expt_mtrx, expt_colors = expt_colors, expt_design = expt_design,
      method = method, expt_names = expt_names, type = type,
      batch_row = batch_row, plot_title = plot_title,
      label_chars = label_chars, ...)
  })

#' Send a SummarizedExperiment to plot_libsize().
#'
#' @param data SummarizedExperiment presumably created by create_se().
#' @param condition Set of conditions observed in the metadata, overriding
#'  the metadata in the SE.
#' @param colors Set of colors for the plot, overriding the SE metadata.
#' @param text Print text with the counts/sample observed at the top of the bars?
#' @param order Optionally redefine the order of the bars of the plot.
#' @param plot_title Plot title!
#' @param yscale Explicitly set the scale on the log or base10 scale.
#' @param expt_names Optionally change the names of the bars.
#' @param label_chars If the names of the bars are larger than this, abbreviate them.
#' @param ... Additonal arbitrary arguments.
#' @return Plot of library sizes and a couple tables describing the data.
#' @export
setMethod(
  "plot_libsize", signature = signature(data = "SummarizedExperiment"),
  definition = function(data, condition = NULL, colors = NULL, text = TRUE,
                        order = NULL, plot_title = NULL, yscale = NULL,
                        expt_names = NULL, label_chars = 10, ...) {
    mtrx <- as.matrix(assay(data))
    condition <- pData(data)[["condition"]]
    colors <- metadata(data)[["colors"]]
    plot_libsize(mtrx, condition = condition, colors = colors, text = text,
                 order = order, plot_title = plot_title, yscale = yscale,
                 expt_names = expt_names, label_chars = label_chars,
                 ...)
  })

#' Run plot_libsize() with a dataframe as input.
#' @export
setMethod(
  "plot_libsize", signature = signature(data = "data.frame", condition = "factor",
                                        colors = "character"),
  definition = function(data, condition, colors, text = TRUE,
                        order = NULL, plot_title = NULL, yscale = NULL,
                        expt_names = NULL, label_chars = 10, ...) {
    data <- as.matrix(data)
    plot_libsize(data, condition = condition, colors = colors,
                 text = text, order = order, plot_title = plot_title, yscale = yscale,
                 expt_names = expt_names, label_chars = label_chars, ...) # , ...)
  })

#' Run plot_libsize() with an expt as input.
#' @export
setMethod(
  "plot_libsize", signature = signature(data = "expt"),
  definition = function(data, condition = NULL, colors = NULL, text = TRUE,
                        order = NULL, plot_title = NULL, yscale = NULL,
                        expt_names = NULL, label_chars = 10, ...) {
    mtrx <- exprs(data)
    condition <- pData(data)[["condition"]]
    colors = data[["colors"]]
    plot_libsize(mtrx, condition = condition, colors = colors, text = text,
                 order = order, plot_title = plot_title, yscale = yscale,
                 expt_names = expt_names, label_chars = label_chars, ...)
  })

#' Run plot_libsize() with an ExpressionSet as input.
#' @export
setMethod(
  "plot_libsize", signature = signature(data = "ExpressionSet"),
  definition = function(data, condition = NULL, colors = NULL, text = TRUE,
                        order = NULL, plot_title = NULL, yscale = NULL,
                        expt_names = NULL, label_chars = 10, ...) {
    mtrx <- exprs(data)
    condition <- pData(data)[["condition"]]
    plot_libsize(mtrx, condition = condition, colors = colors,
                 text = text, order = order, plot_title = plot_title,
                 yscale = yscale, expt_names = expt_names, label_chars = label_chars,
                 ...)
  })

#' Feed an expt to a sankey plot.
#' @export
setMethod(
  "plot_meta_sankey", signature = signature(design = "expt"),
  definition = function(design, factors = c("condition", "batch"),
                        color_choices = NULL) {
    design <- pData(design)
    plot_meta_sankey(design, factors = factors,
                     color_choices = color_choices)
  })

#' Plot the sample heatmap of an expt.
#' @export
setMethod(
  "plot_sample_heatmap", signature = signature(data = "expt"),
  definition = function(data, colors = NULL, design = NULL, heatmap_colors = NULL,
                        expt_names = NULL, dendrogram = "column",
                        row_label = NA, plot_title = NULL, Rowv = TRUE,
                        Colv = TRUE, label_chars = 10, filter = TRUE, ...) {
    expt_design <- pData(data)
    expt_colors <- data[["colors"]]
    expt_names <- data[["expt_names"]]
    expt_data <- exprs(data)
    plot_sample_heatmap(expt_data, colors = expt_colors, design = expt_design,
      expt_names = expt_names, dendrogram = dendrogram, heatmap_colors = heatmap_colors,
      row_label = row_label, plot_title = plot_title, Rowv = Rowv,
      Colv = Colv, label_chars = label_chars, filter = filter, ...)
  })

#' Plot a sample heatmap of an ExpressionSet.
#' @export
setMethod(
  "plot_sample_heatmap", signature = (data = "ExpressionSet"),
  definition = function(data, colors = NULL, design = NULL, heatmap_colors = NULL,
                        expt_names = NULL, dendrogram = "column",
                        row_label = NA, plot_title = NULL, Rowv = TRUE,
                        Colv = TRUE, label_chars = 10, filter = TRUE, ...) {
    expt_design <- pData(data)
    expt_names <- colnames(expt_design)
    expt_data <- exprs(data)
    plot_sample_heatmap(expt_data, colors = expt_colors, design = expt_design,
      expt_names = expt_names, dendrogram = dendrogram, heatmap_colors = heatmap_colors,
      row_label = row_label, plot_title = plot_title, Rowv = Rowv,
      Colv = Colv, label_chars = label_chars, filter = filter, ...)
  })

#' Plot a sample heatmap with a SummarizedExperiment.
#' @export
setMethod(
  "plot_sample_heatmap", signature = signature(data = "SummarizedExperiment"),
  definition = function(data, colors = NULL, design = NULL, heatmap_colors = NULL,
                        expt_names = NULL, dendrogram = "column",
                        row_label = NA, plot_title = NULL, Rowv = TRUE,
                        Colv = TRUE, label_chars = 10, filter = TRUE, ...) {
    expt_design <- pData(data)
    expt_colors <- metadata(data)[["colors"]]
    expt_names <- metadata(data)[["expt_names"]]
    expt_data <- exprs(data)
    plot_sample_heatmap(expt_data, colors = expt_colors, design = expt_design,
      expt_names = expt_names, dendrogram = dendrogram, heatmap_colors = heatmap_colors,
      row_label = row_label, plot_title = plot_title, Rowv = Rowv,
      Colv = Colv, label_chars = label_chars, filter = filter, ...)
  })

#' Plot the standard median pairwise values of an expt.
#' @export
setMethod(
  "plot_sm", signature = signature(data = "expt"),
  definition = function(data, colors = NULL, method = "pearson",
                        plot_legend = FALSE, expt_names = NULL,
                        label_chars = 10, plot_title = NULL, dot_size = 5,
                        ...) {
            design <- pData(data)
            colors <- colors(data)
            conditions <- design[["condition"]]
            mtrx <- exprs(data)
            plot_sm(mtrx, design = design, colors = colors, method = method,
                    plot_legend = plot_legend, expt_names = expt_names,
                    label_chars = label_chars, plot_title = plot_title,
                    dot_size = dot_size,
                    ...)
          })

#' Plot the standard median pairwise values of a SummarizedExperiment.
#' @export
setMethod(
  "plot_sm", signature = signature(data = "SummarizedExperiment"),
  definition = function(data, colors = NULL, method = "pearson",
                        plot_legend = FALSE, expt_names = NULL,
                        label_chars = 10, plot_title = NULL, dot_size = 5,
                        ...) {
    design <- pData(data)
    colors <- colors(data)
    mtrx <- exprs(data)
    plot_sm(mtrx, design = design, colors = colors, method = method,
            plot_legend = plot_legend, expt_names = expt_names,
            label_chars = label_chars, plot_title = plot_title,
            dot_size = dot_size, ...)
  })

#' Plot the standard median pairwise values of an ExpressionSet.
#' @export
setMethod(
  "plot_sm", signature = signature(data = "ExpressionSet"),
  definition = function(data, colors = NULL, method = "pearson",
                        plot_legend = FALSE, expt_names = NULL,
                        label_chars = 10, plot_title = NULL, dot_size = 5,
                        ...) {
    design <- pData(data)
    mtrx <- exprs(data)
    plot_sm(mtrx, design = design, colors = colors, method = method,
            plot_legend = plot_legend, expt_names = expt_names,
            label_chars = label_chars, plot_title = plot_title,
            dot_size = dot_size, ...)
  })

#' Plot the standard median pairwise values of a dataframe.
#' @export
setMethod(
  "plot_sm", signature = signature(data = "data.frame"),
  definition = function(data, colors = NULL, method = "pearson",
                        plot_legend = FALSE, expt_names = NULL,
                        label_chars = 10, plot_title = NULL, dot_size = 5,
                        ...) {
    mtrx <- as.matrix(data)
    plot_sm(mtrx, colors = colors, method = method,
            plot_legend = plot_legend, expt_names = expt_names,
            label_chars = label_chars, plot_title = plot_title,
            dot_size = dot_size, ...)
  })

#' Plot the coefficient of variance values of a SummarizedExperiment.
#' @export
setMethod(
  "plot_variance_coefficients", signature = signature(data = "expt"),
  definition = function(data, design = NULL, x_axis = "condition", colors = NULL,
                        plot_title = NULL, ...) {
    design <- pData(data)
    colors <- colors(data)
    mtrx <- exprs(data)
    plot_variance_coefficients(mtrx, design = design, x_axis = x_axis,
                               colors = colors, plot_title = plot_title, ...)
  })

#' Plot the coefficient of variance values of a SummarizedExperiment.
#' @export
setMethod(
  "plot_variance_coefficients", signature = signature(data = "SummarizedExperiment"),
  definition = function(data, design = NULL, x_axis = "condition", colors = NULL,
                        plot_title = NULL, ...) {
    design <- pData(data)
    colors <- colors(data)
    mtrx <- exprs(data)
    plot_variance_coefficients(mtrx, design = design, x_axis = x_axis,
                               colors = colors, plot_title = plot_title, ...)
  })

#' Plot the coefficient of variance values of an ExpressionSet.
#' @export
setMethod(
  "plot_variance_coefficients", signature = signature(data = "ExpressionSet"),
  definition = function(data, design = NULL, x_axis = "condition", colors = NULL,
                        plot_title = NULL, ...) {
    design <- pData(data)
    mtrx <- exprs(data)
    plot_variance_coefficients(mtrx, design = design, x_axis = x_axis,
                               colors = colors, plot_title = plot_title, ...)
  })

#' A getter of the gene information from an expt, synonymous with fData().
#' @export
setMethod(
  "rowData", signature = signature(x = "expt"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::fData(x[["expressionset"]])
  })

#' A setter to put the gene information into an expt.
#' @export
setMethod(
  "rowData<-", signature = signature(x = "expt"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::fData(x[["expressionset"]]) <- value
    return(x)
  })

#' A getter of the gene information from an ExpressionSet, synonymous with fData().
#' @export
setMethod(
  "rowData", signature = signature(x = "ExpressionSet"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::fData(x)
  })

#' A getter to get the samples names from an expt.
#' @export
setMethod(
  "sampleNames", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::sampleNames(object[["expressionset"]])
  })

#' A setter to put the samples names into an expt.
#' @export
setMethod(
  "sampleNames<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    set_expt_samplenames(object, value)
  })

# #' Metadata sanitizers for an expt
# #' @export
#setMethod(
#  "sanitize_metadata", signature = signature(meta = "expt"),
#  definition = function(meta, columns = NULL, na_string = "notapplicable",
#                        lower = TRUE, punct = TRUE, factorize = "heuristic",
#                        max_levels = NULL, spaces = FALSE, numbers = NULL) {
#    old_meta <- pData(meta)
#    new_meta <- sanitize_metadata(old_meta, columns = columns, na_string = na_string,
#                                  lower = lower, punct = punct,
#                                  factorize = factorize, max_levels = max_levels,
#                                  spaces = spaces, numbers = numbers)
#    pData(meta) <- new_meta
#    return(meta)
#  })
#
# #' Metadata sanitizers for an expressionset
# #' @export
#setMethod(
#  "sanitize_metadata", signature = signature(meta = "ExpressionSet"),
#  definition = function(meta, columns = NULL, na_string = "notapplicable",
#                        lower = TRUE, punct = TRUE, spaces = FALSE,
#                        numbers = NULL) {
#    old_meta <- pData(meta)
#    new_meta <- sanitize_metadata(old_meta, columns = columns, na_string = na_string,
#                                  lower = lower, punct = punct, spaces = spaces,
#                                  numbers = numbers)
#    pData(meta) <- new_meta
#    return(meta)
#  })
#
# #' Metadata sanitizers for a Summarized Experiment.
# #' @export
#setMethod(
#  "sanitize_metadata", signature = signature(meta = "SummarizedExperiment"),
#  definition = function(meta, columns = NULL, na_string = "notapplicable",
#                        lower = TRUE, punct = TRUE, factorize = "heuristic",
#                        max_levels = NULL, spaces = FALSE, numbers = NULL) {
#    old_meta <- pData(meta)
#    new_meta <- sanitize_metadata(old_meta, columns = columns, na_string = na_string,
#                                  lower = lower, punct = punct,
#                                  factorize = factorize, max_levels = max_levels,
#                                  spaces = spaces, numbers = numbers)
#    pData(meta) <- new_meta
#    return(meta)
#  })

#' Extract the state of an expt vis a vis normalization.
#' @export
setMethod(
  "state", signature = signature(expt = "expt"),
  definition = function(expt) {
    expt[["state"]]
  })

#' Put the current state into an expt.
#' @export
setMethod(
  "state<-", signature = signature(expt = "expt"),
  definition = function(expt, value) {
    expt[["state"]] <- value
    return(expt)
  })

#' A getter to get the samples names from a SummarizedExperiment.
#' @export
setMethod(
  "sampleNames", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    BiocGenerics::colnames(object)
  })

#' A setter to put the samples names into a SummarizedExperiment.
#' @export
setMethod(
  "sampleNames<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, value) {
    BiocGenerics::colnames(object) <- value
  })

#' Get the state from a SummarizedExperiment.
#' @export
setMethod(
  "state", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt) {
    metadata(expt)[["state"]]
  })

#' Put the state into a SummarizedExperiment.
#' @export
setMethod(
  "state<-", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, value) {
    metadata(expt)[["state"]] <- value
    return(expt)
  })

#' Subset a SummarizedExperiment with some extra syntax.
#' @export
setMethod(
  "subset_expt", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, subset = NULL, ids = NULL,
                        nonzero = NULL, coverage = NULL) {
    subset_se(expt, subset = subset, ids = ids,
              nonzero = nonzero, coverage = coverage)
  })

## EOF
