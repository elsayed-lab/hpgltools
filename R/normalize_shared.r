## Time-stamp: <Tue May  3 11:54:50 2016 Ashton Trey Belew (abelew@gmail.com)>

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

#' Normalize a dataframe/expt, express it, and/or transform it
#'
#' Sometime soon I am going to elipsis all these variables
#'
#' @param data some data
#' @param design   design dataframe must come with it
#' @param transform   defines whether to log(2|10) transform the
#' data. Defaults to raw.
#' @param norm   specify the normalization strategy.  Defaults to
#' raw.  This makes use of DESeq/EdgeR to provide: RLE, upperquartile,
#' size-factor, or tmm normalization.  I tend to like quantile, but there are
#' definitely corner-case scenarios for all strategies.
#' @param convert   defines the output type which may be raw, cpm,
#' rpkm, or cp_seq_m.  Defaults to raw.
#' @param batch   batch correction method to try out
#' @param batch1  column from design to get batch info
#' @param batch2   a second covariate to try
#' @param filter_low   choose whether to low-count filter the data.
#' @param annotations   is used for rpkm or sequence normalizations to
#' extract the lengths of sequences for normalization
#' @param entry_type   default gff entry to cull from
#' @param fasta   fasta genome for rpkm
#' @param thresh   threshold for low count filtering
#' @param min_samples   minimum samples for low count filtering
#' @param noscale   used by combatmod
#' @param p   for povera genefilter
#' @param A   for povera genefilter
#' @param k   for kovera genefilter
#' @param cv_min   for genefilter cv
#' @param cv_max   for genefilter cv
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
hpgl_norm <- function(data, design=NULL, transform="raw", norm="raw",
                      convert="raw", batch="raw", batch1="batch", batch2=NULL,
                      filter_low=FALSE, annotations=NULL, entry_type="gene",
                      fasta=NULL, thresh=2, min_samples=2,
                      noscale=TRUE, p=0.01, A=1, k=1, cv_min=0.01, cv_max=1000, ...) {
    lowfilter_performed <- FALSE
    norm_performed <- "raw"
    convert_performed <- "raw"
    transform_performed <- "raw"
    batch_performed <- "raw"
    data_class <- class(data)[1]
    original_counts <- NULL
    original_libsize <- NULL
    if (data_class == 'expt') {
        design <- data[["design"]]
        original_counts <- data[["original_counts"]]
        original_libsizes <- data[["original_libsize"]]
        data <- Biobase::exprs(data[["expressionset"]])
    } else if (data_class == 'ExpressionSet') {
        data <- Biobase::exprs(data)
    } else if (data_class == 'list') {
        data <- data$count_table
        if (is.null(data)) {
            stop("The list provided contains no count_table.")
        }
    } else if (data_class == 'matrix' | data_class == 'data.frame') {
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    count_table <- as.matrix(data)
    expt_design <- design
    if (is.null(original_counts)) {
        original_counts <- data
    }
    if (is.null(original_libsize)) {
        original_libsize <- colSums(count_table)
    }

    ## Step 1: Low-count filtering
    lowfiltered_counts <- NULL
    if (filter_low != FALSE) {
        message(paste0("Performing low-count filter with option: ", filter_low))
        ## All the other intermediates have a libsize slot, perhaps this should too
        lowfiltered_counts <- lowfilter_counts(count_table, type=filter_low, p=p, A=A, k=k, cv_min=cv_min, cv_max=cv_max, thresh=2, min_samples=2)
        count_table <- lowfiltered_counts[["count_table"]]
        lowfilter_performed <- filter_low
    }

    ## Step 2: Normalization
    ## This section handles the various normalization strategies
    ## If nothing is chosen, then the filtering is considered sufficient
    normalized_counts <- NULL
    if (norm != "raw") {
        if (is.null(expt_design)) {
            message("The experimental design is null.  Some normalizations will therefore fail.")
            message("If you receive an error about an object with no dimensions, that is likely why.")
        }
        normalized_counts <- normalize_counts(count_table, expt_design, norm=norm)
        count_table <- normalized_counts[["count_table"]]
        norm_performed <- norm
    }

    ## Step 3: Convert the data to (likely) cpm
    ## The following stanza handles the three possible output types
    ## cpm and rpkm are both from edgeR
    ## They have nice ways of handling the log2 which I should consider
    converted_counts <- NULL
    if (convert != "raw") {
        converted_counts <- convert_counts(count_table, convert=convert, annotations=annotations, fasta=fasta, entry_type=entry_type)
        count_table <- converted_counts$count_table
        convert_performed <- convert
    }

    ## Step 4: Batch correction
    batched_counts <- NULL
    if (batch != "raw") {
        ## batched_counts = batch_counts(count_table, batch=batch, batch1=batch1, batch2=batch2, design=design, ...)
        tmp_counts <- try(batch_counts(count_table, batch=batch, batch1=batch1, batch2=batch2, design=expt_design))
        if (class(tmp_counts) == 'try-error') {
            warning("The batch_counts called failed.  Returning non-batch reduced data.")
            batched_counts <- NULL
            batch_performed <- "raw"
        } else {
            batched_counts <- tmp_counts
            batch_performed <- batch
            count_table <- batched_counts[["count_table"]]
        }
    }

    ## Step 5: Transformation
    ## Finally, this considers whether to log2 the data or no
    transformed_counts <- NULL
    if (transform != "raw") {
        message(paste0("Applying: ", transform, " transformation."))
        transformed_counts <- transform_counts(count_table, transform=transform, ...)
        ##transformed_counts <- transform_counts(count_table, transform=transform, converted=convert_performed)
        count_table <- transformed_counts[["count_table"]]
        transform_performed <- transform
    }

    ## This list provides the list of operations performed on the data in order they were done.
    actions <- list(
        "lowfilter" = lowfilter_performed,
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
        "lowfilter" = lowfiltered_counts,  ## After lowfiltering
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

#'   Replace the data of an expt with normalized data.
#'
#' @param expt   The original expt
#' @param transform   The transformation desired (raw, log2, log, log10)
#' @param norm   How to normalize the data (raw, quant, sf, upperquartile, tmm, rle)
#' @param convert   Conversion to perform (raw, cpm, rpkm, cp_seq_m)
#' @param batch   Batch effect removal tool to use (limma sva fsva ruv etc)
#' @param filter_low   Filter out low sequences (cbcb, pofa, kofa, others?)
#' @param annotations  used for rpkm, a df
#' @param fasta  fasta file for cp_seq_m counting of oligos
#' @param entry_type   for getting genelengths by feature type (rpkm or cp_seq_m)
#' @param use_original   whether to use the backup data in the expt class
#' @param batch1   experimental factor to extract first
#' @param batch2   a second factor to remove (only with limma's removebatcheffect())
#' @param thresh   for cbcb_lowfilter
#' @param min_samples   for cbcb_lowfilter
#' @param p   for genefilter's pofa
#' @param A   for genefilter's pofa
#' @param k   for genefilter's kofa
#' @param cv_min   for genefilter's cv()
#' @param cv_max  for genefilter's cv()
#' @param ... more options
#' @return a new expt object with normalized data and the original data saved as 'original_expressionset'
#' @seealso \pkg{genefilter} \pkg{limma} \pkg{sva} \pkg{edgeR} \pkg{DESeq2}
#' @examples
#' \dontrun{
#' normed <- normalize_expt(exp, transform='log2', norm='rle', convert='cpm',
#'                          batch='raw', filter_low='pofa')
#' normed_batch <- normalize_expt(exp, transform='log2', norm='rle', convert='cpm',
#'                                batch='sva', filter_low='pofa')
#' }
#' @export
normalize_expt <- function(expt, ## The expt class passed to the normalizer
    ## choose the normalization strategy
    transform="raw", norm="raw", convert="raw", batch="raw", filter_low=FALSE,
    ## annotations used for rpkm/cpseqm, original may be used to ensure double-normalization isn't performed.
    annotations=NULL, fasta=NULL, entry_type="gene", use_original=FALSE,
    batch1="batch", batch2=NULL, ## extra parameters for batch correction
    thresh=2, min_samples=2, p=0.01, A=1, k=1, cv_min=0.01, cv_max=1000,  ## extra parameters for low-count filtering
    ...) {
    new_expt <- expt
    current_exprs <- expt[["expressionset"]]
    if (is.null(new_expt[["original_expressionset"]])) {
        new_expt[["original_expressionset"]] = new_expt[["expressionset"]]
    } else {
        message("This function will replace the expt$expressionset slot with:")
        type <- ""
        if (transform != "raw") {
            type <- paste0(type, transform, '(')
        }
        if (batch != "raw") {
            type <- paste0(type, 'batch-correct(')
        }
        if (convert != "raw") {
            type <- paste0(type, convert, '(')
        }
        if (norm != "raw") {
            type <- paste0(type, norm, '(')
        }
        if (filter_low != FALSE) {
            type <- paste0(type, 'low-filter(')
        }
        type <- paste0(type, 'data')
        if (transform != 'raw') {
            type <- paste0(type, ')')
        }
        if (batch != "raw") {
            type <- paste0(type, ')')
        }
        if (convert != "raw") {
            type <- paste0(type, ')')
        }
        if (norm != "raw") {
            type <- paste0(type, ')')
        }
        if (filter_low != FALSE) {
            type <- paste0(type, ')')
        }
        message(type)
        message("It backs up the current data into a slot named:
 expt$backup_expressionset. It will also save copies of each step along the way
 in expt$normalized with the corresponding libsizes. Keep the libsizes in mind
 when invoking limma.  The appropriate libsize is the non-log(cpm(normalized)).
 This is most likely kept at:
 'new_expt$normalized$intermediate_counts$normalization$libsizes'
 A copy of this may also be found at:
 new_expt$best_libsize
")
    }
    if (filter_low == FALSE) {
        message("Filter low is false, this should likely be set to something, good
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
    new_expt[["backup_expressionset"]] <- new_expt[["expressionset"]]
    current_data <- Biobase::exprs(current_exprs)
    design <- expt[["design"]]
    ## A bunch of these options should be moved into ...
    ## Having them as options to maintain is foolish
    normalized <- hpgl_norm(current_data, design=design, transform=transform,
                            norm=norm, convert=convert, batch=batch,
                            batch1=batch1, batch2=batch2,
                            filter_low=filter_low, annotations=annotations,
                            fasta=fasta, thresh=thresh,
                            min_samples=min_samples, p=p, A=A, k=k,
                            cv_min=cv_min, cv_max=cv_max, entry_type=entry_type)
    final_libsize <- normalized[["libsize"]]
    final_data <- as.matrix(normalized[["count_table"]])
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
        "lowfilter" = normalized[["actions"]][["lowfilter"]],
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
    } else if (!is.null(normalized[["intermediate_counts"]][["lowfilter"]][["libsize"]])) {
        new_expt[["best_libsize"]] <- normalized[["intermediate_counts"]][["lowfilter"]][["libsize"]]
    } else {
        new_expt[["best_libsize"]] <- NULL
    }
    ## limma should probably use this
    new_expt[["norm_result"]] <- normalized
    new_expt[["expressionset"]] <- current_exprs
    return(new_expt)
}

## EOF
