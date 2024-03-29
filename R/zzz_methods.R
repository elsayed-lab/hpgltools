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

#' Set the colors for an expt.
#'
#' @param expt Expt to modify.
#' @param ... Colors!
#' @return The expt with new colors associated with each sample.
setGeneric("colors<-", signature = signature(expt = "expt"),
           function(expt, ...) standardGeneric("colors<-", ...))

## I cannot seem to define this generic correctly.
#setGeneric("extract_keepers", signature = signature(extracted = "list", keepers = "list"),
#           function(extracted, keepers, ...) standardGeneric("extract_keepers"))

setGeneric("get_backup_expression_data", signature = c("expt"),
           function(expt) standardGeneric("get_backup_expression_data"))

#' Generic method to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setGeneric("iDA", signature = c("object"),
           function(object, ...) standardGeneric("iDA"))

#' Get the state from an expt.
#'
#' @param expt One of my slightly modified ExpressionSets.
#' @return List with the methods used to modify the data (if any).
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

#' A setter for the annotation database used to create an expt/se.
#'
#' @param object One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param value New annotation slot for the expt/se.
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
#' @param i specific samples to replace the data.
#' @param withDimnames I do not know.
#' @param ... Extra args, currently unused.
#' @param value New assay values to fill in the data structure.
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
#' @param withDimnames I do not know.
#' @param ... Extra args!
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
#' @param i Subset to replace.
#' @param withDimnames I do not know, I need to look this up.
#' @param ... Extra args.
#' @param value New values for the expressionset.
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
    S4Vectors::metadata(expt)[["original_se"]] <- backup
    return(expt)
  })

#' A getter to pull the sample data from an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param withDimnames Again, haven't looked it up yet.
#' @param ... Extra args.
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
#' @param i Subset to replace.
#' @param withDimnames indeed.
#' @param ... extra args.
#' @param value New Sample data for the expt.
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
#' @param withDimnames indeed.
#' @param ... extra args.
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
#' @param i Slice to replace.
#' @param withDimnames yes
#' @param ... args for the arglist
#' @param value New values for the expressionset.
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
    S4Vectors::metadata(expt)[["colors"]]
  })

#' A setter to put the colors into a SummarizedExperiment.
#'
#' @param expt A SummarizedExperiment.
#' @param value List of new colors.
#' @export
setMethod(
  "colors<-", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, value) {
    S4Vectors::metadata(expt)[["colors"]] <- value
    return(expt)
  })

#' A setter to put the colors into an expt.
#'
#' @param expt An expt.
#' @param value List of new colors.
#' @export
setMethod(
  "colors<-", signature = signature(expt = "expt"),
  definition = function(expt, value) {
    expt[["colors"]] <- value
    return(expt)
  })

#' Count nmers given a filename instead of genome object.
#'
#' @param genome filename of the genome in question
#' @param pattern Pattern for which to search.
#' @param mismatch Number of mismatches allowed.
#' @export
setMethod(
  "count_nmer", signature = signature(genome = "character"),
  definition = function(genome, pattern = "ATG", mismatch = 0) {
    new_genome <- Rsamtools::FaFile(genome)
    count_nmer(new_genome, pattern = pattern, mismatch = mismatch)
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

#' @describeIn extract_keepers Use a character vector instead of a list.
#' @export
setMethod(
  "extract_keepers", signature = signature(extracted = "list", keepers = "character"),
  definition = function(extracted, keepers, table_names,
                        all_coefficients, apr,
                        adjp, annot_df, includes,
                        excludes, padj_type,
                        fancy = FALSE, loess = FALSE,
                        lfc_cutoff = 1.0, p_cutoff = 0.05,
                        format_sig = 4, plot_colors = plot_colors,
                        z = 1.5, alpha = 0.4, z_lines = FALSE,
                        label = 10, label_column = "hgncsymbol") {
    if (keepers[1] == "all") {
      new_keepers <- list()
      numerators <- denominators <- c()
      ## Note, I changed table_names to be sorted by method.  I can
      ## either iterate over every method, take the union of all, or
      ## arbitrarily choose a method...
      possible_names <- c()
      ## Limma and edger are the most likely to have extra contrasts,
      ## so check one of them first.
      if (!is.null(table_names[["limma"]])) {
        possible_names <- table_names[["limma"]]
      } else {
        possible_names <- table_names[[1]]
      }
      for (a in seq_along(possible_names)) {
        name <- possible_names[[a]]
        splitted <- strsplit(x = name, split = "_vs_")
        denominator <- splitted[[1]][2]
        numerator <- splitted[[1]][1]
        new_keepers[[name]] <- c(numerator, denominator)
      }
    } else {
      splitted <- strsplit(x = keepers, split = "_vs_")
      numerator <- splitted[[1]][1]
      denominator <- splitted[[1]][2]
      new_keepers <- list(splitted = c(numerator, denominator))
    }
    extract_keepers(extracted, new_keepers, table_names,
                    all_coefficients, apr,
                    adjp, annot_df, includes,
                    excludes, padj_type,
                    fancy = fancy, loess = loess,
                    lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff,
                    format_sig = format_sig, plot_colors = plot_colors,
                    z = z, alpha = alpha, z_lines = z_lines,
                    label = label, label_column = label_column)
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
    backup <- S4Vectors::metadata(expt)[["original_se"]]
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
  "normalize_expt", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, transform = "raw", norm = "raw", convert = "raw",
                        batch = "raw", filter = FALSE, annotations = NULL, fasta = NULL,
                        entry_type = "gene", use_original = FALSE, batch1 = "batch",
                        batch2 = NULL, batch_step = 4, low_to_zero = TRUE,
                        thresh = 2, min_samples = 2, p = 0.01, A = 1, k = 1,
                        cv_min = 0.01, cv_max = 1000, na_to_zero = FALSE,
                        adjust_method = "ruv", verbose = FALSE, ...) {
    normalize_se(expt, transform = transform, norm = norm, convert = convert,
                 batch = batch, filter = filter, annotations = annotations,
                 fasta = fasta, entry_type = entry_type, use_original = use_original,
                 batch1 = batch1, batch2 = batch2, batch_step = batch_step,
                 low_to_zero = low_to_zero, thresh = thresh, min_samples = min_samples,
                 p = p, A = A, k = k, cv_min = cv_min, cv_max = cv_max,
                 na_to_zero = na_to_zero, adjust_method = adjust_method,
                 verbose = verbose, ...)
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
    expt_colors <- S4Vectors::metadata(expt_data)[["colors"]]
    expt_names <- S4Vectors::metadata(expt_data)[["expt_names"]]
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
    colors <- S4Vectors::metadata(data)[["colors"]]
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

#' Make a nonzero plot given an expt.
#' @export
setMethod(
  "plot_nonzero", signature = signature(data = "expt"),
  definition = function(data, design = NULL, colors = NULL, plot_labels = "repel",
                        expt_names = NULL, max_overlaps = 5, label_chars = 10,
                        plot_legend = FALSE, plot_title = NULL, cutoff = 0.65, ...) {
    mtrx <- as.matrix(exprs(data))
    pd <- pData(data)
    condition <- pd[["condition"]]
    names <- pd[["samplenames"]]
    chosen_colors <- colors(data)
    plot_nonzero(mtrx, design = pd, colors = chosen_colors, plot_labels = plot_labels,
                 expt_names = names, max_overlaps = max_overlaps, label_chars = label_chars,
                 plot_legend = plot_legend, plot_title = plot_title, cutoff = 0.65, ...)
  })

#' Make a nonzero plot given an ExpressionSet
#' @export
setMethod(
  "plot_nonzero", signature = signature(data = "ExpressionSet"),
  definition = function(data, design = NULL, colors = NULL, plot_labels = "repel",
                        expt_names = NULL, max_overlaps = 5, label_chars = 10,
                        plot_legend = FALSE, plot_title = NULL, cutoff = 0.65, ...) {
    mtrx <- as.matrix(exprs(data))
    pd <- pData(data)
    condition <- pd[["condition"]]
    names <- pd[["samplenames"]]
    plot_nonzero(mtrx, design = pd, colors = colors, plot_labels = plot_labels,
                 expt_names = names, max_overlaps = max_overlaps,
                 label_chars = label_chars, plot_legend = plot_legend,
                 plot_title = plot_title, cutoff = 0.65, ...)
  })

#' Make a nonzero plot given a SummarizedExperiment
#' @export
setMethod(
  "plot_nonzero", signature = signature(data = "SummarizedExperiment"),
  definition = function(data, design = NULL, colors = NULL, plot_labels = "repel",
                        expt_names = NULL, max_overlaps = 5, label_chars = 10,
                        plot_legend = FALSE, plot_title = NULL, cutoff = 0.65, ...) {
    mtrx <- as.matrix(assay(data))
    pd <- SummarizedExperiment::colData(data)
    condition <- pd[["condition"]]
    names <- pd[["samplenames"]]
    colors <- S4Vectors::metadata(data)[["colors"]]
    plot_nonzero(mtrx, design = pd, colors = colors, plot_labels = plot_labels,
                 expt_names = names, max_overlaps = max_overlaps, label_chars = label_chars,
                 plot_legend = plot_legend, plot_title = plot_title, cutoff = 0.65, ...)
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
    expt_colors <- colors(data)
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
    expt_colors <- S4Vectors::metadata(data)[["colors"]]
    expt_names <- S4Vectors::metadata(data)[["expt_names"]]
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
    x <- Biobase::fData(x[["expressionset"]]) <- value
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

#' Coerce simple_topgo to accept a vector of gene IDs instead of a real dataframe of significance.
#'
#' Doing this voids the topgo warantee.
setMethod(
  "simple_topgo", signature = signature(sig_genes = "character"),
  definition = function(sig_genes, goid_map = "id2go.map", go_db = NULL,
                        pvals = NULL, limitby = "fisher", limit = 0.1,
                        signodes = 100, sigforall = TRUE, numchar = 300,
                        selector = "topDiffGenes", pval_column = "deseq_adjp",
                        overwrite = FALSE, densities = FALSE,
                        pval_plots = TRUE, excel = NULL, ...) {
    fake_df <- data.frame(row.names = sig_genes)
    fake_df[["ID"]] <- rownames(fake_df)
    fake_df[[pval_column]] <- 0.01
    warning("Faking a dataframe with significance of every gene as 0.01 because this was given a vector of gene IDs.")
    simple_topgo(fake_df, goid_map = goid_map, go_db = go_db,
                 pvals = pvals, limitby = limitby, limit = limit,
                 signodes = signodes, sigforall = sigforall,
                 numchar = numchar, selector = selector,
                 pval_column = pval_column, overwrite = overwrite,
                 densities = densities, pval_plots = pval_plots, excel = excel,
                 ...)
  })

#' Get the state from a SummarizedExperiment.
#' @export
setMethod(
  "state", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt) {
    S4Vectors::metadata(expt)[["state"]]
  })

#' Put the state into a SummarizedExperiment.
#' @export
setMethod(
  "state<-", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, value) {
    S4Vectors::metadata(expt)[["state"]] <- value
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

#' Write an xlsx file given the result of an existing xlsx write.
#' @export
setMethod(
  "write_xlsx", signature = signature(excel = "written_xlsx"),
  definition = function(data = NULL, wb = NULL, sheet = NULL, excel,
                        rownames = TRUE, start_row = 1, start_col = 1,
                        title = NULL, number_format = "0.000", data_table = TRUE,
                        freeze_first_row = TRUE, freeze_first_column = TRUE,
                        column_width = "heuristic", ...) {
    current_wb <- excel[["workbook"]]
    current_sheet <- excel[["sheet"]]
    current_row <- excel[["end_row"]]
    current_col <- excel[["end_col"]]
    current_excel <- excel[["file"]]
    if (is.null(sheet)) {
      sheet <- current_sheet
      ## You cannot have > 1 frozen first sheet, so if you are reusing the sheet,
      ## make sure freeze is off
      freeze_first_row <- FALSE
      freeze_first_column <- FALSE
      column_width <- NULL
      if (is.null(start_row)) {
        start_row <- current_row + 1
      }
      if (is.null(start_col)) {
        start_col <- 1
      }
    }
    write_xlsx(data = data, wb = current_wb, sheet = sheet, excel = current_excel,
               rownames = rownames, start_row = start_row, start_col = start_col,
               title = title, number_format = number_format, data_table = data_table,
               freeze_first_row = freeze_first_row, freeze_first_column = freeze_first_column,
               column_width = column_width, ...)
  })

#setMethod(
#  "write_gprofiler_data", signature = signature(gprofiler_result = "all_gprofiler",
#                                                contrast = "NULL"),
#  definition = function(gprofiler_result, contrast, wb = NULL, excel = "excel/profiler_result.xlsx",
#                        order_by = "recall", add_plots = TRUE, height = 15, width = 10, decreasing = FALSE, ...) {
#    outputs <- list()
#    for (contrast in names(gprofiler_result)) {
#      excel = glue("excel/gprofiler_{contrast}.xlsx")
#      input <- gprofiler_result[[contrast]]
#      outputs[[contrast]] <- write_gprofiler_data(input, excel = excel, order_by = order_by,
#                                                  add_plots = add_plots, height = height,
#                                                  width = width, decreasing = decreasing, ...)
#    }
#    return(outputs)
#  })

## EOF
