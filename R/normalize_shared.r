## normalize_shared.r: Functions which bring together everything in the other
## normalize_* files.

#' Perform a default normalization of some data
#'
#' This just calls normalize expt with the most common arguments except log2
#' transformation, but that may be appended with 'transform = log2', so I don't
#' feel bad.  Indeed, it will allow you to overwrite any arguments if you wish.
#' In our work, the most common normalization is: quantile(cpm(low-filter(data))).
#'
#' @param expt An expressionset containing expt object
#' @param ... More options to pass to normalize_expt()
#' @return The normalized expt
#' @seealso [normalize_expt()]
#' @export
default_norm <- function(expt, ...) {
  arglist <- list(...)
  norm <- "quant"
  if (!is.null(arglist[["norm"]])) {
    norm <- arglist[["norm"]]
  }
  convert <- "cpm"
  if (!is.null(arglist[["convert"]])) {
    convert <- arglist[["convert"]]
  }
  filter <- TRUE
  if (!is.null(arglist[["filter"]])) {
    filter <- arglist[["filter"]]
  }

  new <- normalize_expt(expt, norm = norm, convert = convert, filter = filter, ...)
  return(new)
}

#' Normalize the data of an expt object.  Save the original data, and note what
#' was done.
#'
#' It is the responsibility of normalize_expt() to perform any arbitrary
#' normalizations desired as well as to ensure that the data integrity is
#' maintained.  In order to do this, it writes the actions performed in
#' expt$state and saves the intermediate steps of the normalization in
#' expt$intermediate_counts.  Furthermore, it should tell you every step of the
#' normalization process, from count filtering, to normalization, conversion,
#' transformation, and batch correction.
#'
#' @param expt Original expt.
#' @param transform Transformation desired, usually log2.
#' @param norm How to normalize the data? (raw, quant, sf, upperquartile, tmm, rle)
#' @param convert Conversion to perform? (raw, cpm, rpkm, cp_seq_m)
#' @param batch Batch effect removal tool to use? (limma sva fsva ruv etc)
#' @param filter Filter out low/undesired features? (cbcb, pofa, kofa, others?)
#' @param annotations Used for rpkm -- probably not needed as this is in fData now.
#' @param fasta Fasta file for cp_seq_m counting of oligos.
#' @param entry_type For getting genelengths by feature type (rpkm or cp_seq_m).
#' @param use_original Use the backup data in the expt class?
#' @param batch1 Experimental factor to extract first.
#' @param batch2 Second factor to remove (only with limma's removebatcheffect()).
#' @param batch_step From step 1-5, when should batch correction be applied?
#' @param low_to_zero When log transforming, change low numbers (< 0) to 0 to avoid NaN?
#' @param thresh Used by cbcb_lowfilter().
#' @param min_samples Also used by cbcb_lowfilter().
#' @param p Used by genefilter's pofa().
#' @param A Also used by genefilter's pofa().
#' @param k Used by genefilter's kofa().
#' @param cv_min Used by genefilter's cv().
#' @param cv_max Also used by genefilter's cv().
#' @param na_to_zero Sometimes rpkm gives some NA values for very low numbers.
#' @param verbose Print what is happening while the normalization is performed?
#'  I am not sure why, but I think they should be 0.
#' @param ... more options
#' @return Expt object with normalized data and the original data saved as
#'  'original_expressionset'
#' @seealso [convert_counts()] [normalize_counts()] [batch_counts()]
#'  [filter_counts()] [transform_counts()]
#' @examples
#' \dontrun{
#'  normed <- normalize_expt(exp, transform='log2', norm='rle', convert='cpm',
#'                           batch='raw', filter='pofa')
#'  normed_batch <- normalize_expt(exp, transform='log2', norm='rle', convert='cpm',
#'                                 batch='sva', filter='pofa')
#' }
#' @export
normalize_expt <- function(expt, ## The expt class passed to the normalizer
                           ## choose the normalization strategy
                           transform = "raw", norm = "raw", convert = "raw",
                           batch = "raw", filter = FALSE,
                           ## annotations used for rpkm/cpseqm, original may be
                           ## used to ensure double-normalization isn't
                           ## performed.
                           annotations = NULL, fasta = NULL, entry_type = "gene",
                           use_original = FALSE, batch1 = "batch", batch2 = NULL, batch_step = 4,
                           low_to_zero = TRUE, ## extra parameters for batch correction
                           thresh = 2, min_samples = 2, p = 0.01, A = 1, k = 1,
                           cv_min = 0.01, cv_max = 1000,  ## extra parameters for low-count filtering
                           na_to_zero = FALSE, adjust_method = "ruv", verbose = TRUE,
                           ...) {
  arglist <- list(...)
  expt_state <- expt[["state"]]
  new_expt <- expt
  type <- ""
  current_exprs <- expt[["expressionset"]]
  if (!is.null(arglist[["filter_low"]])) {
    ## I changed the name of this argument.
    warning("This argument has been changed to 'filter'.")
    filter <- arglist[["filter"]]
  }
  if (is.null(filter)) {
    filter <- "raw"
  }
  if (filter == FALSE) {
    filter <- "raw"
  } else if (isTRUE(filter)) {
    filter <- "cbcb"
  }
  if (is.null(convert)) {
    convert <- "raw"
  }
  if (convert == FALSE) {
    convert <- "raw"
  } else if (isTRUE(convert)) {
    convert <- "cbcbcpm"
  }
  if (is.null(norm)) {
    norm <- "raw"
  }
  if (norm == FALSE) {
    norm <- "raw"
  } else if (isTRUE(norm)) {
    norm <- "tmm"
  }
  if (is.null(transform)) {
    transform <- "raw"
  }
  if (transform == FALSE) {
    transform <- "raw"
  } else if (isTRUE(transform)) {
    transform <- "log2"
  }
  if (is.null(batch)) {
    batch <- "raw"
  }
  if (batch == FALSE) {
    batch <- "raw"
  } else if (isTRUE(batch)) {
    batch <- "sva"
  }

  if (is.null(new_expt[["original_expressionset"]])) {
    new_expt[["original_expressionset"]] <- new_expt[["expressionset"]]
  }

  mesg("This function will replace the expt$expressionset slot with:")
    operations <- what_happened(transform = transform, batch = batch, convert = convert,
                                norm = norm, filter = filter)
  mesg(operations)
  mesg("It will save copies of each step along the way
 in expt$normalized with the corresponding libsizes. Keep libsizes in mind
 when invoking limma.  The appropriate libsize is non-log(cpm(normalized)).
 This is most likely kept at:
 'new_expt$normalized$intermediate_counts$normalization$libsizes'
 A copy of this may also be found at:
 new_expt$best_libsize
")

  if (filter == "raw") {
    mesg("Filter is false, this should likely be set to something, good
 choices include cbcb, kofa, pofa (anything but FALSE).  If you want this to
 stay FALSE, keep in mind that if other normalizations are performed, then the
 resulting libsizes are likely to be strange (potentially negative!)
")
  }
  if (transform == "raw") {
    mesg("Leaving the data in its current base format, keep in mind that
 some metrics are easier to see when the data is log2 transformed, but
 EdgeR/DESeq do not accept transformed data.
")
  }
  if (convert == "raw") {
    mesg("Leaving the data unconverted.  It is often advisable to cpm/rpkm
 the data to normalize for sampling differences, keep in mind though that rpkm
 has some annoying biases, and voom() by default does a cpm (though hpgl_voom()
 will try to detect this).
")
  }
  if (norm == "raw") {
    mesg("Leaving the data unnormalized.  This is necessary for DESeq, but
 EdgeR/limma might benefit from normalization.  Good choices include quantile,
 size-factor, tmm, etc.
")
  }
  if (batch == "raw") {
    mesg("Not correcting the count-data for batch effects.  If batch is
 included in EdgerR/limma's model, then this is probably wise; but in extreme
 batch effects this is a good parameter to play with.
")
  }
  if (convert == "cpm" & transform == "tmm") {
    warning("Cpm and tmm perform similar purposes. They should not be applied to the same data.")
  }
  if (norm == "quant" & isTRUE(grepl(x = batch, pattern = "sva"))) {
    warning("Quantile normalization and sva do not always play well together.")
  }

  new_expt[["backup_expressionset"]] <- new_expt[["expressionset"]]
  data <- exprs(current_exprs)
  design <- pData(expt)
  if (is.null(annotations)) {
    annotations <- fData(current_exprs)
  }
  ## A bunch of these options should be moved into ...
  ## Having them as options to maintain is foolish
  if (isTRUE(verbose)) {
    normalized <- hpgl_norm(data, expt_state = expt_state, design = design, transform = transform,
                            norm = norm, convert = convert, batch = batch,
                            batch1 = batch1, batch2 = batch2, low_to_zero = low_to_zero,
                            filter = filter, annotations = annotations,
                            fasta = fasta, thresh = thresh, batch_step = batch_step,
                            min_samples = min_samples, p = p, A = A, k = k,
                            cv_min = cv_min, cv_max = cv_max, entry_type = entry_type,
                            adjust_method = adjust_method,
                            ...)
  } else {
    normalized <- sm(hpgl_norm(data, expt_state = expt_state, design = design, transform = transform,
                               norm = norm, convert = convert, batch = batch,
                               batch1=batch1, batch2=batch2, low_to_zero = low_to_zero,
                               filter = filter, annotations = annotations,
                               fasta = fasta, thresh = thresh, batch_step = batch_step,
                               min_samples = min_samples, p = p, A = A, k = k,
                               cv_min = cv_min, cv_max = cv_max, entry_type = entry_type,
                               adjust_method = adjust_method,
                               ...))
  }

  final_libsize <- normalized[["libsize"]]
  final_data <- as.matrix(normalized[["count_table"]])

  ## A recent update to Biobase adds a test in the function
  ## assayDataElementReplace() which no longer allows one to just
  ## replace an expressionset with a smaller version (low-filtered).
  ## Instead, one must properly subset the object first, then replace.
  ## While this is annoying, I suppose it is smart.

  ## There is another peculiarity which has popped up in my proteomics data.
  ## Some rows are all NA!  How this happened I have not yet figured out,
  ## but these rows need to go.
  row_na_idx <- !is.na(rownames(final_data))
  final_data <- final_data[row_na_idx, ]
  unfiltered_genes <- rownames(exprs(current_exprs)) %in% rownames(final_data)
  current_exprs <- current_exprs[unfiltered_genes, ]
  ## This next line was added in response to an annoying occurance when
  ## combining technical replicates of a data set for which the samples started
  ## with numbers ('2018_0315' for example).  In that instance, my
  ## concatenate_runs() does not properly check to ensure that the column names
  ## do not change. I should therefore change that, but for the moment I will
  ## add a check here.
  colnames(final_data) <- sampleNames(current_exprs)
  exprs(current_exprs) <- final_data

  ## This state slot should match the information available in
  ## new_expt$normalized$actions
  ## I am hoping this will prove a more direct place to access it and provide a
  ## chance to double-check that things match
  new_state <- list(
    "filter" = normalized[["actions"]][["filter"]],
    "normalization" = normalized[["actions"]][["normalization"]],
    "conversion" = normalized[["actions"]][["conversion"]],
    "batch" = normalized[["actions"]][["batch"]],
    "transform" = normalized[["actions"]][["transform"]])

  ## Keep in mind that low_to_zero should be ignored if transform_state is not raw.
  if (new_state[["transform"]] == "raw" & isTRUE(low_to_zero)) {
    low_idx <- exprs(current_exprs) < 0
    ## Do not forget that na does not count when looking for numbers < x.
    na_idx <- is.na(exprs(current_exprs))
    low_idx[na_idx] <- FALSE
    low_num <- sum(low_idx)
    if (low_num > 0) {
      message("Setting ", low_num, " entries to zero.")
      exprs(current_exprs)[low_idx] <- 0
    }
  }

  if (isTRUE(na_to_zero)) {
    na_idx <- is.na(exprs(current_exprs))
    exprs(current_exprs)[na_idx] <- 0
  }

  ## The original data structure contains the following slots:
  ## colors, batches, convert, conditions, design, expressionset,
  ## filtered, initial_metadata, norm, original_expressionset,
  ## original_libsize, samplenames, stages, types, transform
  ## batches, convert, transform, filtered, normalization are being
  ## replaced with 'state' containing all of them.
  ## Similarly, the multiple libsizes maintained are going to be put into
  ## libsizes = list(original, norm, best, etc)
  ## Finally, I want to remove 'stages' and 'types' those are data in design.

  ## The structure of the 'normalized' data is fairly complex and includes the following:
  ## actions_performed -- a list of the actions done by hpgl_norm()
  ## intermediate_counts -- counts after each step in the hpgl_norm() process
  ## count_table -- a dataframe count table from hpgl_norm()
  ## libsize -- a final libsize from hpgl_norm()
  new_expt[["normalized"]] <- normalized
  new_expt[["state"]] <- new_state

  ## My concept of the 'best library size' comes from Kwame's work where the
  ## libsize was kept after performing quantile normalization, but before doing
  ## a log2(cpm()) The problem with this is that, if one does a normalize() then
  ## another normalize() then the assumptions used  may get violated.
  if (!is.null(normalized[["intermediate_counts"]][["normalization"]][["libsize"]])) {
    new_expt[["best_libsize"]] <-
      normalized[["intermediate_counts"]][["normalization"]][["libsize"]]
  } else if (!is.null(normalized[["intermediate_counts"]][["filter"]][["libsize"]])) {
    new_expt[["best_libsize"]] <-
      normalized[["intermediate_counts"]][["filter"]][["libsize"]]
  } else {
    new_expt[["best_libsize"]] <- NULL
  }
  ## limma should probably use this
  new_expt[["norm_result"]] <- normalized
  new_expt[["expressionset"]] <- current_exprs
  current_notes <- glue("{new_expt[['notes']]} Normalized with {type} at {date()}.
")
  new_expt[["notes"]] <- toString(current_notes)
  return(new_expt)
}

#' Normalize a dataframe/expt, express it, and/or transform it
#'
#' There are many possible options to this function.  Refer to normalize_expt()
#' for a more complete list.
#'
#' @param data Some data as a df/expt/whatever.
#' @param ... I should put all those other options here
#' @return edgeR's DGEList expression of a count table.  This seems to
#'  me to be the easiest to deal with.
#' @seealso [edgeR] [DESeq2] [edgeR::cpm()] [filter_counts()] [batch_counts()]
#'  [convert_counts()] [transform_counts()]
#' @export
#' @examples
#' \dontrun{
#'  df_raw = hpgl_norm(expt = expt)  ## Only performs low-count filtering
#'  df_raw = hpgl_norm(df = a_df, design = a_design) ## Same, but using a df
#'  df_ql2rpkm = hpgl_norm(expt = expt, norm='quant', transform='log2',
#'                         convert='rpkm')  ## Quantile, log2, rpkm
#'  count_table = df_ql2rpkm$counts
#' }
hpgl_norm <- function(data, ...) {
  arglist <- list(...)
  batch <- arglist[["batch"]]
  filter_performed <- "raw"
  norm_performed <- "raw"
  convert_performed <- "raw"
  transform_performed <- "raw"
  batch_performed <- "raw"
  adjust_performed <- "none"
  expt_state <- list(
    "low_filter" = filter_performed,
    "normalization" = norm_performed,
    "conversion" = convert_performed,
    "batch" = batch_performed,
    "adjust" = adjust_performed,
    "transform" = transform_performed)
  data_class <- class(data)[1]
  original_counts <- NULL
  original_libsize <- NULL
  annot <- NULL
  counts <- NULL
  expt_design <- NULL
  ## I never quite realized just how nice data.tables are.  To what extent can I refactor
  ## all of my data frame usage to them?
  if (data_class == "expt") {
    original_counts <- data[["original_counts"]]
    original_libsizes <- data[["original_libsize"]]
    expt_design <- pData(data)
    annot <- fData(data)
    counts <- exprs(data)
    expt_state <- data[["state"]]
  } else if (data_class == "ExpressionSet") {
    counts <- exprs(data)
    expt_design <- pData(data)
    annot <- fData(data)
    if (!is.null(arglist[["expt_state"]])) {
      expt_state <- arglist[["expt_state"]]
    }
  } else if (data_class == "list") {
    counts <- data[["count_table"]]
    expt_design <- arglist[["design"]]
    if (is.null(data)) {
      stop("The list provided contains no count_table.")
    }
    if (!is.null(arglist[["expt_state"]])) {
      expt_state <- arglist[["expt_state"]]
    }
  } else if (data_class == "matrix" ||
             data_class == "data.frame" ||
             data_class == "data.table") {
    counts <- as.data.frame(data)  ## some functions prefer matrix, so I am
    ## keeping this explicit for the moment. In the case of data.tables, even if
    ## you set the rownames, the first column might still be rowname characters
    ## I don't yet fully understand this, so I will add an explicit test here.
    if (data_class == "data.table" & class(counts[[1]]) == "character") {
      rownames(counts) <- make.names(counts[[1]], unique = TRUE)
      counts <- counts[-1]
    }
    expt_design <- arglist[["design"]]
    if (!is.null(arglist[["expt_state"]])) {
      expt_state <- arglist[["expt_state"]]
    }
  } else {
    stop("This only understands types: expt, ExpressionSet, data.frame, and matrix.")
  }

  count_table <- as.matrix(counts)
  if (is.null(original_counts)) {
    original_counts <- counts
  }
  if (is.null(original_libsize)) {
    original_libsize <- colSums(count_table)
  }
  annotations <- NULL
  if (!is.null(arglist[["annotations"]])) {
    annotations <- arglist[["annotations"]]
  } else if (!is.null(annot)) {
    annotations <- annot
  }

  ## Make changes to this as we go.
  current_state <- expt_state

  batched_counts <- NULL
  batch_step <- 5
  if (!is.null(arglist[["batch_step"]])) {
    batch_step <- arglist[["batch_step"]]
  }
  if (!is.numeric(batch_step)) {
    batch_step <- 5
  } else if (batch_step > 5 | batch_step < 0) {
    batch_step <- 5
  }

  do_batch <- function(count_table, method = "raw", expt_design = expt_design,
                       current_state = current_state, adjust_method = adjust_method,
                       ...) {
    ##batch <- method
    ##if (!is.null(arglist[["batch"]])) {
    ##  batch <- arglist[["batch"]]
    ##}
    if (is.null(method)) {
      method <- "raw"
    }
    mesg("The batch method is: ", method, ".")

    if (method == "raw") {
      mesg("Step ", arglist[["batch_step"]], ": not doing batch correction.")
    } else {
      mesg("Step ", arglist[["batch_step"]], ": doing batch correction with ",
           arglist[["batch"]], ".")
      tmp_counts <- try(batch_counts(count_table, method = method,
                                     expt_design = expt_design,
                                     current_state = current_state,
                                     ...))
      ##tmp_counts <- batch_counts(count_table, method = method, expt_design = expt_design,
      ##                           adjust_method = adjust_method, current_state = current_state)
      if (class(tmp_counts) == "try-error") {
        warning("The batch_counts call failed.  Returning non-batch reduced data.")
        batched_counts <<- NULL
        batch_performed <- "raw"
        adjust_performed <- "none"
      } else {
        batched_counts <- tmp_counts
        batch_performed <<- batch
        adjust_performed <- adjust_method
        count_table <- batched_counts[["count_table"]]
      }
    }
    return(count_table)
  }

  if (batch_step == 1) {
    count_table <- do_batch(count_table, method = batch,
                            current_design = expt_design,
                            ...)
  }

  ## Step 1: count filtering
  filter <- FALSE
  if (!is.null(arglist[["filter"]])) {
    filter <- arglist[["filter"]]
  }
  filtered_counts <- NULL
  if (filter == FALSE | filter == "raw") {
    mesg("Step 1: not doing count filtering.")
  } else {
    if (isTRUE(filter)) {
      filter <- "cbcb"
    }
    mesg("Step 1: performing count filter with option: ", filter)
    ## All the other intermediates have a libsize slot, perhaps this should too
    filtered_counts <- filter_counts(count_table, method = filter, ...)
    ## filtered_counts <- filter_counts(count_table, method = filter)
    count_table <- filtered_counts[["count_table"]]
    filter_performed <- filter
    current_state[["filter"]] <- filter
  }

  if (batch_step == 2) {
    count_table <- do_batch(count_table, method = batch,
                            expt_design = expt_design,
                            current_state = current_state,
                            ...)
  }
  ## Step 2: Normalization
  ## This section handles the various normalization strategies
  norm <- "raw"
  if (!is.null(arglist[["norm"]])) {
    norm <- arglist[["norm"]]
  }
  normalized_counts <- NULL
  if (norm == "raw") {
    mesg("Step 2: not normalizing the data.")
  } else {
    mesg("Step 2: normalizing the data with ", norm, ".")
    if (is.null(expt_design)) {
      message("The experimental design is null.  Some normalizations will fail.")
      message("If you get an error about 'no dimensions', that is likely why.")
    }
    normalized_counts <- normalize_counts(data = count_table, method = norm, ...)
    ## normalized_counts <- normalize_counts(data = count_table, design = design, method = norm)
    count_table <- normalized_counts[["count_table"]]
    norm_performed <- norm
    current_state[["normalization"]] <- norm
  }

  ## Step 3: Convert the data to (likely) cpm
  ## The following stanza handles the three possible output types
  ## cpm and rpkm are both from edgeR
  ## They have nice ways of handling the log2 which I should consider
  if (batch_step == 3) {
    count_table <- do_batch(count_table, method = batch,
                            expt_design = expt_design,
                            current_state = current_state,
                            ...)
  }
  converted_counts <- NULL
  convert <- "raw"
  if (!is.null(arglist[["convert"]])) {
    convert <- arglist[["convert"]]
  }
  if (convert == "raw") {
    mesg("Step 3: not converting the data.")
  } else {
    mesg("Step 3: converting the data with ", convert, ".")
    converted_counts <- convert_counts(count_table, method = convert, annotations = annotations,
                                       ...)
    ## converted_counts <- convert_counts(count_table, method = convert, annotations = annotations)
    count_table <- converted_counts[["count_table"]]
    convert_performed <- convert
    current_state[["conversion"]] <- convert
  }

  ## Step 4: Transformation
  ## Finally, this considers whether to log2 the data or no
  if (batch_step == 4) {
    count_table <- do_batch(count_table, method = batch,
                            expt_design = expt_design,
                            current_state = current_state,
                            ...)
    ## count_table <- do_batch(count_table, method = batch,
    ##                         expt_design = expt_design, current_state = current_state)
  }
  transformed_counts <- NULL
  transform <- "raw"
  if (!is.null(arglist[["transform"]])) {
    transform <- arglist[["transform"]]
  }
  if (transform == "raw") {
    mesg("Step 4: not transforming the data.")
  } else {
    mesg("Step 4: transforming the data with ", transform, ".")
    transformed_counts <- transform_counts(count_table, method = transform, ...)
    ## transformed_counts <- transform_counts(count_table, transform = transform)
    count_table <- transformed_counts[["count_table"]]
    if (transform == "round") {
      transform_performed <- "raw"
      transform <- "raw"
      current_state[["rounded"]] <- TRUE
    } else {
      transform_performed <- transform
    }
    current_state[["transform"]] <- transform
  }

  if (batch_step == 5) {
    count_table <- do_batch(count_table, method = batch,
                            expt_design = expt_design,
                            current_state = current_state,
                            ...)
    ## count_table <- do_batch(count_table, arglist)
  }

  ## This list provides the list of operations performed on the data in order
  ## they were done.
  actions <- list(
    "filter" = filter_performed,
    "normalization" = norm_performed,
    "conversion" = convert_performed,
    "batch" = batch_performed,
    "adjust" = adjust_performed,
    "transform" = transform_performed)
  ## This list contains the intermediate count tables generated at each step
  ## This may be useful if there is a problem in this process.
  ## Each of them also contains the libsize at that point in the process.
  intermediate_counts <- list(
    "original" = original_counts, ## The original count table, should never
    ## change from iteration to iteration
    "input" = as.matrix(data),  ## The input provided to this function, this may
    ## diverge from original
    "filter" = filtered_counts,  ## After filtering
    "normalization" = normalized_counts,  ## and normalization
    "conversion" = converted_counts,  ## and conversion
    "batch" = batched_counts,  ## and batch correction
    "transform" = transformed_counts)  ## and finally, transformation.

  ret_list <- list(
    "actions" = actions,
    "intermediate_counts" = intermediate_counts,
    "count_table" = count_table,  ## The final count table
    "final_state" = current_state,
    "libsize" = colSums(count_table, na.rm = TRUE)  ## The final libsizes
  )
  return(ret_list)
}

normalize <- function(expt, todo = list()) {
  ## This expects a list like:
  ## list("norm" = "quant", "filter" = c(pofa, "A" = 1))
  possible_methods <- list(
    "transform" = "transform_counts",
    "norm" = "normalize_counts",
    "convert" = "convert_counts",
    "filter" = "filter_counts",
    "batch" = "batch_counts",
    "impute" = "impute_counts")
  annot <- fData(expt)
  counts <- exprs(expt)
  design <- pData(expt)
  count_table <- list("count_table" = counts,
                      "libsize" = expt[["libsize"]])
  for (i in 1:length(todo)) {
    type <- names(todo)[i]
    operation <- todo[[i]]
    method <- operation[[1]]
    args <- c()
    if (length(operation) > 1) {
       args <- operation[2:length(operation)]
    }
    if (! type %in% names(possible_methods)) {
      stop("This type of todo is not known: ", type, ".")
    }
    call <- possible_methods[[type]]
    arglist <- list()
    arglist[["count_table"]] <- count_table
    arglist[["method"]] <- method
    arglist[["design"]] <- design
    for (a in args) {
      name <- names(args)[a]
      arglist[[a]] <- args[[a]]
    }
    mesg("Invoking: ", call)
    count_table <- base::do.call(call, arglist)
  }
  return(count_table)
}

## EOF
