## Note to self, @title and @description are not needed in roxygen
## comments, the first separate #' is the title, the second the
## description, the third is the long form description.  Then add in
## the @param @return @seealso @export etc...

## My root question, to which I think I have a small idea about the answer:
## Two of the very many paths to toptable()/toptags():
##  1.  normalize data -> model(~0 + condition + batch) -> limma
###     Including batch in the model loses some power, but improves the 'truth' of the result
##  2.  normalize data -> batch correction(factor(batch)) ->  model(~0 + condition) -> limma
###     Power lost before the model, also improves the 'truth' of the result.
##  Why is #1 better than #2?
### More well understood and conservative.
##  Why have we nonetheless done #2 in a few instances?  (not only because we learned that first)

#' Perform a default normalization of some data
#'
#' This just calls normalize expt with the most common arguments except log2 transformation, but
#' that may be appended with 'transform=log2', so I don't feel bad.  Indeed, it will allow you to
#' overwrite any arguments if you wish.  In our work, the most common normalization is:
#' quantile(cpm(low-filter(data))).
#'
#' @param expt An expressionset containing expt object
#' @param ... More options to pass to normalize_expt()
#' @return The normalized expt
#' @seealso \code{\link{normalize_expt}}
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

    new <- sm(normalize_expt(expt, norm=norm, convert=convert, filter=filter, ...))
    return(new)
}

#' Normalize the data of an expt object.  Save the original data, and note what was done.
#'
#' It is the responsibility of normalize_expt() to perform any arbitrary normalizations desired as
#' well as to ensure that the data integrity is maintained.  In order to do this, it writes the
#' actions performed in expt$state and saves the intermediate steps of the normalization in
#' expt$intermediate_counts.  Furthermore, it should tell you every step of the normalization
#' process, from count filtering, to normalization, conversion, transformation, and batch
#' correction.
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
#' @param ... more options
#' @return Expt object with normalized data and the original data saved as 'original_expressionset'
#' @seealso \pkg{genefilter} \pkg{limma} \pkg{sva} \pkg{edgeR} \pkg{DESeq2}
#' @examples
#' \dontrun{
#' normed <- normalize_expt(exp, transform='log2', norm='rle', convert='cpm',
#'                          batch='raw', filter='pofa')
#' normed_batch <- normalize_expt(exp, transform='log2', norm='rle', convert='cpm',
#'                                batch='sva', filter='pofa')
#' }
#' @export
normalize_expt <- function(expt, ## The expt class passed to the normalizer
    ## choose the normalization strategy
    transform="raw", norm="raw", convert="raw", batch="raw", filter=FALSE,
    ## annotations used for rpkm/cpseqm, original may be used to ensure double-normalization isn't performed.
    annotations=NULL, fasta=NULL, entry_type="gene", use_original=FALSE,
    batch1="batch", batch2=NULL, batch_step=5, low_to_zero=FALSE, ## extra parameters for batch correction
    thresh=2, min_samples=2, p=0.01, A=1, k=1, cv_min=0.01, cv_max=1000,  ## extra parameters for low-count filtering
    ...) {
    arglist <- list(...)
    new_expt <- expt
    type <- ""
    current_exprs <- expt[["expressionset"]]
    if (!is.null(arglist[["filter_low"]])) {  ## I changed the name of this argument.
        warning("This argument has been changed to 'filter'.")
        filter <- arglist[["filter"]]
    }
    if (filter == FALSE) {
        filter <- "raw"
    } else if (isTRUE(filter)) {
        filter <- "cbcb"
    }
    if (convert == FALSE) {
        convert <- "raw"
    } else if (isTRUE(convert)) {
        convert <- "cbcbcpm"
    }
    if (norm == FALSE) {
        norm <- "raw"
    } else if (isTRUE(norm)) {
        norm <- "tmm"
    }
    if (transform == FALSE) {
        transform <- "raw"
    } else if (isTRUE(transform)) {
        transform <- "log2"
    }
    if (batch == FALSE) {
        batch <- "raw"
    } else if (isTRUE(batch)) {
        batch <- "sva"
    }

    if (is.null(new_expt[["original_expressionset"]])) {
        new_expt[["original_expressionset"]] = new_expt[["expressionset"]]
    }

<<<<<<< HEAD
=======
    message("This function will replace the expt$expressionset slot with:")
>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
    operations <- what_happened(transform=transform, batch=batch, convert=convert, norm=norm, filter=filter)
    message(operations)
    message("It backs up the current data into a slot named:
 expt$backup_expressionset. It will also save copies of each step along the way
 in expt$normalized with the corresponding libsizes. Keep the libsizes in mind
 when invoking limma.  The appropriate libsize is the non-log(cpm(normalized)).
 This is most likely kept at:
 'new_expt$normalized$intermediate_counts$normalization$libsizes'
 A copy of this may also be found at:
 new_expt$best_libsize
")


    if (filter == "raw") {
        message("Filter is false, this should likely be set to something, good
 choices include cbcb, kofa, pofa (anything but FALSE).  If you want this to
 stay FALSE, keep in mind that if other normalizations are performed, then the
 resulting libsizes are likely to be strange (potentially negative!)
")
    }
    if (transform == "raw") {
        message("Leaving the data in its current base format, keep in mind that
 some metrics are easier to see when the data is log2 transformed, but
 EdgeR/DESeq do not accept transformed data.
")
    }
    if (convert == "raw") {
        message("Leaving the data unconverted.  It is often advisable to cpm/rpkm
 the data to normalize for sampling differences, keep in mind though that rpkm
 has some annoying biases, and voom() by default does a cpm (though hpgl_voom()
 will try to detect this).
")
    }
    if (norm == "raw") {
        message("Leaving the data unnormalized.  This is necessary for DESeq, but
 EdgeR/limma might benefit from normalization.  Good choices include quantile,
 size-factor, tmm, etc.
")
    }
    if (batch == "raw") {
        message("Not correcting the count-data for batch effects.  If batch is
 included in EdgerR/limma's model, then this is probably wise; but in extreme
 batch effects this is a good parameter to play with.
")
    }
    if (convert == "cpm" & transform == "tmm") {
        warning("Cpm and tmm perform similar purposes. They should not be applied to the same data.")
    }
    if (norm == "quant" & isTRUE(grepl(x=batch, pattern="sva"))) {
        warning("Quantile normalization and sva do not always play well together.")
    }
    new_expt[["backup_expressionset"]] <- new_expt[["expressionset"]]
    current_data <- Biobase::exprs(current_exprs)
    design <- expt[["design"]]
    if (is.null(annotations)) {
        annotations <- Biobase::fData(current_exprs)
    }
    if (is.null(design)) {
        design <- Biobase::pData(current_exprs)
    }
    ## A bunch of these options should be moved into ...
    ## Having them as options to maintain is foolish
    normalized <- hpgl_norm(current_data, design=design, transform=transform,
                            norm=norm, convert=convert, batch=batch,
                            batch1=batch1, batch2=batch2, low_to_zero=low_to_zero,
                            filter=filter, annotations=annotations,
                            fasta=fasta, thresh=thresh, batch_step=batch_step,
                            min_samples=min_samples, p=p, A=A, k=k,
                            ## cv_min=cv_min, cv_max=cv_max, entry_type=entry_type)
                            cv_min=cv_min, cv_max=cv_max, entry_type=entry_type, ...)

    final_libsize <- normalized[["libsize"]]
    final_data <- as.matrix(normalized[["count_table"]])

    ## A recent update to Biobase adds a test in the function
    ## assayDataElementReplace() which no longer allows one to just
    ## replace an expressionset with a smaller version (low-filtered).
    ## Instead, one must properly subset the object first, then replace.
    ## While this is annoying, I suppose it is smart.
    unfiltered_genes <- rownames(Biobase::exprs(current_exprs)) %in% rownames(final_data)
    current_exprs <- current_exprs[unfiltered_genes, ]
    Biobase::exprs(current_exprs) <- final_data

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

    ## This state slot should match the information available in
    ## new_expt$normalized$actions
    ## I am hoping this will prove a more direct place to access it and provide a chance to double-check that things match
    new_state <- list(
        "filter" = normalized[["actions"]][["filter"]],
        "normalization" = normalized[["actions"]][["normalization"]],
        "conversion" = normalized[["actions"]][["conversion"]],
        "batch" = normalized[["actions"]][["batch"]],
        "transform" = normalized[["actions"]][["transform"]])
    new_expt[["state"]] <- new_state

    ## My concept of the 'best library size' comes from Kwame's work where the libsize was kept after performing
    ## quantile normalization, but before doing a log2(cpm())
    ## The problem with this is that, if one does a normalize() then another normalize() then the assumptions used
    ## may get violated.
    if (!is.null(normalized[["intermediate_counts"]][["normalization"]][["libsize"]])) {
        new_expt[["best_libsize"]] <- normalized[["intermediate_counts"]][["normalization"]][["libsize"]]
    } else if (!is.null(normalized[["intermediate_counts"]][["filter"]][["libsize"]])) {
        new_expt[["best_libsize"]] <- normalized[["intermediate_counts"]][["filter"]][["libsize"]]
    } else {
        new_expt[["best_libsize"]] <- NULL
    }
    ## limma should probably use this
    new_expt[["norm_result"]] <- normalized
    new_expt[["expressionset"]] <- current_exprs
    current_notes <- paste0(new_expt[["notes"]], "Normalized with ", type, " at ", date(), ".\n")
    new_expt[["notes"]] <- toString(current_notes)
    return(new_expt)
}

#' Normalize a dataframe/expt, express it, and/or transform it
#'
#' There are many possible options to this function.  Refer to normalize_expt() for a more complete list.
#'
#' @param data Some data as a df/expt/whatever.
#' @param ... I should put all those other options here
#' @return edgeR's DGEList expression of a count table.  This seems to
#' me to be the easiest to deal with.
#' @seealso \link[edgeR]{cpm} \link[edgeR]{rpkm}
#' \link{hpgl_rpkm} \link[DESeq2]{DESeqDataSetFromMatrix}
#' \link[DESeq]{estimateSizeFactors} \link[edgeR]{DGEList} \link[edgeR]{calcNormFactors}
#' @export
#' @examples
#' \dontrun{
#' df_raw = hpgl_norm(expt=expt)  ## Only performs low-count filtering
#' df_raw = hpgl_norm(df=a_df, design=a_design) ## Same, but using a df
#' df_ql2rpkm = hpgl_norm(expt=expt, norm='quant', transform='log2',
#'                        convert='rpkm')  ## Quantile, log2, rpkm
#' count_table = df_ql2rpkm$counts
#' }
hpgl_norm <- function(data, ...) {
    arglist <- list(...)
    filter_performed <- "raw"
    norm_performed <- "raw"
    convert_performed <- "raw"
    transform_performed <- "raw"
    batch_performed <- "raw"
    data_class <- class(data)[1]
    original_counts <- NULL
    original_libsize <- NULL
    annot <- NULL
    counts <- NULL
    ## I never quite realized just how nice data.tables are.  To what extent can I refactor all of my data frame usage to them?
    if (data_class == 'expt') {
        original_counts <- data[["original_counts"]]
        original_libsizes <- data[["original_libsize"]]
        design <- Biobase::pData(data[["expressionset"]])
        annot <- Biobase::fData(data[["expressionset"]])
        counts <- Biobase::exprs(data[["expressionset"]])
    } else if (data_class == 'ExpressionSet') {
        counts <- Biobase::exprs(data)
        design <- Biobase::pData(data)
        annot <- Biobase::fData(data)
    } else if (data_class == "list") {
        counts <- data[["count_table"]]
        design <- arglist[["design"]]
        if (is.null(data)) {
            stop("The list provided contains no count_table.")
        }
    } else if (data_class == "matrix" | data_class == "data.frame" | data_class == "data.table") {
        counts <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
        ## In the case of data.tables, even if you set the rownames, the first column might still be rowname characters
        ## I don't yet fully understand this, so I will add an explicit test here.
        if (data_class == "data.table" & class(counts[[1]]) == "character") {
            rownames(counts) <- make.names(counts[[1]], unique=TRUE)
            counts <- counts[-1]
        }
        design <- arglist[["design"]]
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    count_table <- as.matrix(counts)
    expt_design <- design
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

    do_batch <- function(count_table, design=design, ...) {
        batch <- "raw"
        if (!is.null(arglist[["batch"]])) {
            batch <- arglist[["batch"]]
        }
        if (batch == "raw") {
            message(paste0("Step ", arglist[["batch_step"]], ": not doing batch correction."))
        } else {
            message(paste0("Step ", arglist[["batch_step"]], ": doing batch correction with ",
                           arglist[["batch"]],"."))
            tmp_counts <- try(batch_counts(count_table, design=design, ...))
            if (class(tmp_counts) == "try-error") {
                warning("The batch_counts call failed.  Returning non-batch reduced data.")
                batched_counts <<- NULL
                batch_performed <- "raw"
            } else {
                batched_counts <- tmp_counts
                batch_performed <<- batch
                count_table <- batched_counts[["count_table"]]
            }
        }
        return(count_table)
    }

    if (batch_step == 1) {
        count_table <- do_batch(count_table, ...)
    }

    ## Step 1: count filtering
    filter <- FALSE
    if (!is.null(arglist[["filter"]])) {
        filter <- arglist[["filter"]]
    }
    filtered_counts <- NULL
    if (filter == FALSE | filter == "raw") {
        message("Step 1: not doing count filtering.")
    } else {
        if (isTRUE(filter)) {
            filter <- "cbcb"
        }
        message(paste0("Step 1: performing count filter with option: ", filter))
        ## All the other intermediates have a libsize slot, perhaps this should too
        filtered_counts <- filter_counts(count_table, ...)
        ## filtered_counts <- filter_counts(count_table, filter)
        count_table <- filtered_counts[["count_table"]]
        filter_performed <- filter
    }

    if (batch_step == 2) {
        count_table <- do_batch(count_table, ...)
    }
    ## Step 2: Normalization
    ## This section handles the various normalization strategies
    ## If nothing is chosen, then the filtering is considered sufficient
    norm <- "raw"
    if (!is.null(arglist[["norm"]])) {
        norm <- arglist[["norm"]]
    }
    normalized_counts <- NULL
    if (norm == "raw") {
        message("Step 2: not normalizing the data.")
    } else {
        message(paste0("Step 2: normalizing the data with ", arglist[["norm"]], "."))
        if (is.null(expt_design)) {
            message("The experimental design is null.  Some normalizations will therefore fail.")
            message("If you receive an error about an object with no dimensions, that is likely why.")
        }
        normalized_counts <- normalize_counts(count_table, ...)
        ## normalized_counts <- normalize_counts(count_table, norm=norm)
        count_table <- normalized_counts[["count_table"]]
        norm_performed <- norm
    }

    ## Step 3: Convert the data to (likely) cpm
    ## The following stanza handles the three possible output types
    ## cpm and rpkm are both from edgeR
    ## They have nice ways of handling the log2 which I should consider
    if (batch_step == 3) {
        count_table <- do_batch(count_table, ...)
    }
    converted_counts <- NULL
    convert <- "raw"
    if (!is.null(arglist[["convert"]])) {
        convert <- arglist[["convert"]]
    }
    if (convert == "raw") {
        message("Step 3: not converting the data.")
    } else {
        message(paste0("Step 3: converting the data with ", arglist[["convert"]], "."))
        converted_counts <- convert_counts(count_table, ...)
        ## converted_counts <- convert_counts(count_table, convert=convert)
        count_table <- converted_counts[["count_table"]]
        convert_performed <- convert
    }

    ## Step 4: Transformation
    ## Finally, this considers whether to log2 the data or no
    if (batch_step == 4) {
        count_table <- do_batch(count_table, ...)
    }
    transformed_counts <- NULL
    transform <- "raw"
    if (!is.null(arglist[["transform"]])) {
        transform <- arglist[["transform"]]
    }
    if (transform == "raw") {
        message("Step 4: not transforming the data.")
    } else {
        message(paste0("Step 4: transforming the data with ", arglist[["transform"]], "."))
        transformed_counts <- transform_counts(count_table, ...)
        ## transformed_counts <- transform_counts(count_table, transform=transform, converted=convert_performed)
        count_table <- transformed_counts[["count_table"]]
        if (transform == "round") {
            transform_performed <- "raw"
            transform <- "raw"
        } else {
            transform_performed <- transform
        }
    }

    if (batch_step == 5) {
        count_table <- do_batch(count_table, ...)
    }

    ## This list provides the list of operations performed on the data in order they were done.
    actions <- list(
        "filter" = filter_performed,
        "normalization" = norm_performed,
        "conversion" = convert_performed,
        "batch" = batch_performed,
        "transform" = transform_performed)
    ## This list contains the intermediate count tables generated at each step
    ## This may be useful if there is a problem in this process.
    ## Each of them also contains the libsize at that point in the process.
    intermediate_counts <- list(
        "original" = original_counts, ## The original count table, should never change from iteration to iteration
        "input" = as.matrix(data),  ## The input provided to this function, this may diverge from original
        "filter" = filtered_counts,  ## After filtering
        "normalization" = normalized_counts,  ## and normalization
        "conversion" = converted_counts,  ## and conversion
        "batch" = batched_counts,  ## and batch correction
        "transform" = transformed_counts)  ## and finally, transformation.

    ret_list <- list(
        "actions" = actions,
        "intermediate_counts" = intermediate_counts,
        "count_table" = count_table,  ## The final count table
        "libsize" = colSums(count_table)  ## The final libsizes
    )
    return(ret_list)
}

## EOF
