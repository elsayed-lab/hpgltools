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
setMethod("plot_libsize",
          signature = signature(data = "SummarizedExperiment"),
          definition = function(data, condition = NULL, colors = NULL, text = TRUE,
                                order = NULL, plot_title = NULL, yscale = NULL,
                                expt_names = NULL, label_chars = 10, ...) {
            mtrx <- as.matrix(assay(data))
            condition <- metadata(data)[["conditions"]]
            colors <- metadata(data)[["colors"]]
            plot_libsize(mtrx, condition = condition, colors = colors, text = text,
                         order = order, plot_title = plot_title, yscale = yscale,
                         expt_names = expt_names, label_chars = label_chars,
                         ...)
          })
setMethod("exprs", signature = "SummarizedExperiment",
          function(object) {
            SummarizedExperiment::assay(object)
          })
setMethod("exprs<-", signature = "SummarizedExperiment",
          function(object, value) {
            SummarizedExperiment::assay(object) <- value
            return(object)
          })
setMethod("fData", signature = "SummarizedExperiment",
          function(object) {
            SummarizedExperiment::rowData(object)
          })
setMethod("fData<-", signature = "SummarizedExperiment",
          function(object, value) {
            SummarizedExperiment::rowData(object) <- value
            return(object)
          })
setMethod("pData", signature = "SummarizedExperiment",
          function(object) {
            SummarizedExperiment::colData(object)
          })
setMethod("pData<-", signature = "SummarizedExperiment",
          function(object, value) {
            SummarizedExperiment::colData(object) <- value
            return(object)
          })
setMethod("sampleNames", signature = "SummarizedExperiment",
          function(object) {
            BiocGenerics::colnames(object)
          })
setMethod("sampleNames<-", signature = "SummarizedExperiment",
          function(object, value) {
            BiocGenerics::colnames(object) <- value
          })
setMethod("normalize_expt",
          signature = signature(expt = "SummarizedExperiment"),
          definition = function(expt, transform = "raw", norm = "raw", convert = "raw",
                                batch = "raw", filter = FALSE,
                                annotations = NULL, fasta = NULL, entry_type = "gene",
                                use_original = FALSE, batch1 = "batch",
                                batch2 = NULL, batch_step = 4,
                                low_to_zero = TRUE, thresh = 2, min_samples = 2,
                                p = 0.01, A = 1, k = 1, cv_min = 0.01, cv_max = 1000,
                                na_to_zero = FALSE, adjust_method = "ruv", verbose = TRUE,
                                ...) {
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
setMethod("subset_expt",
          signature = signature(expt = "SummarizedExperiment"),
          definition = function(expt, subset = NULL, ids = NULL,
                                nonzero = NULL, coverage = NULL) {
            subset_se(expt, subset = subset, ids = ids,
                      nonzero = nonzero, coverage = coverage)
          })
