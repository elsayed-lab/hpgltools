#' Grab gene lengths from a gff file.
#'
#' This function attempts to be robust to the differences in output from importing gff2/gff3 files.
#' But it certainly isn't perfect.
#'
#' @param gff Gff file with (hopefully) IDs and widths.
#' @param type Annotation type to use (3rd column).
#' @param key Identifier in the 10th column of the gff file to use.
#' @param ... Extra arguments likely for gff2df
#' @return Data frame of gene IDs and widths.
#' @seealso \pkg{rtracklayer}
#'  \code{\link{gff2df}}
#' @examples
#' \dontrun{
#'  tt = get_genelengths('reference/fun.gff.gz')
#'  head(tt)
#' ##           ID width
#' ## 1   YAL069W   312
#' ## 2   YAL069W   315
#' ## 3   YAL069W     3
#' ## 4 YAL068W-A   252
#' ## 5 YAL068W-A   255
#' ## 6 YAL068W-A     3
#' }
#' @export
get_genelengths <- function(gff, type="gene", key="ID", ...) {
    ret <- gff2df(gff, ...)
    ret <- ret[ret[["type"]] == type, ]
    ret <- ret[, c(key, "width")]
    colnames(ret) <- c("ID", "width")
    if (dim(ret)[1] == 0) {
        stop(paste0("No genelengths were found. ",
                    "Perhaps you are using the wrong 'type' or 'key' arguments, type is: ",
                    type, ", key is: ", key))
    }
    return(ret)
}

#' Given a data frame of exon counts and annotation information, sum the exons.
#'
#' This function will merge a count table to an annotation table by the child column.
#' It will then sum all rows of exons by parent gene and sum the widths of the exons.
#' Finally it will return a list containing a df of gene lengths and summed counts.
#'
#' @param data Count tables of exons.
#' @param gff Gff filename.
#' @param annotdf Dataframe of annotations (probably from gff2df).
#' @param parent Column from the annotations with the gene names.
#' @param child Column from the annotations with the exon names.
#' @return List of 2 data frames, counts and lengths by summed exons.
#' @seealso \pkg{rtracklayer}
#'  \code{\link{gff2df}}
#' @examples
#' \dontrun{
#' summed <- sum_exons(counts, gff='reference/xenopus_laevis.gff.xz')
#' }
#' @export
sum_exons <- function(data, gff=NULL, annotdf=NULL, parent='Parent', child='row.names') {
    if (is.null(annotdf) & is.null(gff)) {
        stop("I need either a df with parents, children, and widths; or a gff filename.")
    } else if (is.null(annotdf)) {
        annotdf <- gff2df(gff)
    }

    tmp_data <- merge(data, annotdf, by=child)
    rownames(tmp_data) <- tmp_data[["Row.names"]]
    tmp_data <- tmp_data[-1]
    ## Start out by summing the gene widths
    column <- aggregate(tmp_data[, "width"], by=list(Parent=tmp_data[, parent]), FUN=sum)
    new_data <- data.frame(column[["x"]])
    rownames(new_data) <- column[["Parent"]]
    colnames(new_data) <- c("width")

    for (c in 1:length(colnames(data))) {
        column_name <- colnames(data)[[c]]
        column <- aggregate(tmp_data[, column_name], by=list(Parent=tmp_data[, parent]), FUN=sum)
        rownames(column) <- column[["Parent"]]
        new_data <- cbind(new_data, column[["x"]])
    } ## End for loop

    width_df <- data.frame(new_data[["width"]])
    rownames(width_df) <- rownames(new_data)
    colnames(width_df) <- c("width")
    new_data <- new_data[-1]
    colnames(new_data) <- colnames(data)
    rownames(new_data) <- rownames(column)
    ret <- list(
        "width" = width_df,
        "counts" = new_data)
    return(ret)
}

#' Extract annotation information from a gff file into a df
#'
#' Try to make import.gff a little more robust; I acquire (hopefully) valid gff files from various
#' sources: yeastgenome.org, microbesonline, tritrypdb, ucsc, ncbi. To my eyes, they all look like
#' reasonably good gff3 files, but some of them must be loaded with import.gff2, import.gff3, etc.
#' That is super annoying. Also, I pretty much always just do as.data.frame() when I get something
#' valid from rtracklayer, so this does that for me, I have another function which returns the
#' iranges etc.  This function wraps import.gff/import.gff3/import.gff2 calls in try() because
#' sometimes those functions fail in unpredictable ways.
#'
#' @param gff Gff filename.
#' @param type Subset the gff file for entries of a specific type.
#' @param id_col Column in a successful import containing the IDs of interest.
#' @param second_id_col Second column to check.
#' @param try Give your own function call to use for importing.
#' @return Dataframe of the annotation information found in the gff file.
#' @seealso \pkg{rtracklayer} \pkg{GenomicRanges}
#'  \code{\link[rtracklayer]{import.gff}}
#' @examples
#' \dontrun{
#'  funkytown <- gff2df('reference/gff/saccharomyces_cerevsiae.gff.xz')
#' }
#' @export
gff2df <- function(gff, type=NULL, id_col="ID", second_id_col="locus_tag", try=NULL) {
    ret <- NULL
    success <- FALSE
    attempts <- c("rtracklayer::import.gff3(gff, sequenceRegionsAsSeqinfo=TRUE)",
                  "rtracklayer::import.gff3(gff, sequenceRegionsAsSeqinfo=FALSE)",
                  "rtracklayer::import.gff2(gff, sequenceRegionsAsSeqinfo=TRUE)",
                  "rtracklayer::import.gff2(gff, sequenceRegionsAsSeqinfo=FALSE)",
                  "rtracklayer::import.gff(gff, format='gtf')",
                  "rtracklayer::import.gff(gff)")
    if (!is.null(try)) {
        attempts <- c(paste0(try, "(gff)"), attempts)
    }

    annot <- NULL
    for (att in 1:length(attempts)) {
        annotations <- NULL
        message(paste0("Trying attempt: ", attempts[[att]]))
        attempt <- attempts[[att]]
        eval_string <- paste0("annotations <- try(", attempt, ", silent=TRUE)")
        eval(parse(text=eval_string))
        if (class(annotations) == "try-error") {
            success <- FALSE
            rm(annotations)
        } else if (is.null(GenomicRanges::as.data.frame(annotations)[[id_col]]) &
                   is.null(GenomicRanges::as.data.frame(annotations)[[second_id_col]])) {
            success <- FALSE
            rm(annotations)
        } else {
            success <- TRUE
            annot <- annotations
            rm(annotations)
            message("Had a successful gff import with ", attempt)
            break
        }
    }
    ret <- NULL
    if (class(annot)[[1]] == "GRanges") {
        ret <- GenomicRanges::as.data.frame(annot)
        rm(annot)
    } else {
        stop("Unable to load gff file.")
    }
    if (!is.null(type)) {
        index <- ret[, "type"] == type
        ret <- ret[index, ]
    }
    message(paste0("Returning a df with ", ncol(ret), " columns and ", nrow(ret), " rows."))
    return(ret)
}

#' Extract annotation information from a gff file into an irange object.
#'
#' Try to make import.gff a little more robust; I acquire (hopefully) valid gff files from various
#' sources: yeastgenome.org, microbesonline, tritrypdb, ucsc, ncbi. To my eyes, they all look like
#' reasonably good gff3 files, but some of them must be loaded with import.gff2, import.gff3, etc.
#' That is super annoying. Also, I pretty much always just do as.data.frame() when I get something
#' valid from rtracklayer, so this does that for me, I have another function which returns the
#' iranges etc.  This function wraps import.gff/import.gff3/import.gff2 calls in try() because
#' sometimes those functions fail in unpredictable ways.
#'
#' This is essentially gff2df(), but returns data suitable for getSet()
#'
#' @param gff Gff filename.
#' @param type Subset to extract.
#' @return Iranges! (useful for getSeq().)
#' @seealso \pkg{rtracklayer} \link{gff2df} \pkg{Biostrings}
#'  \code{\link[rtracklayer]{import.gff}}
#' @examples
#' \dontrun{
#'  library(BSgenome.Tcruzi.clbrener.all)
#'  tc_clb_all <- BSgenome.Tcruzi.clbrener.all
#'  cds_ranges <- gff2irange('reference/gff/tcruzi_clbrener.gff.xz', type='CDS')
#'  cds_sequences <- Biostrings::getSeq(tc_clb_all, cds_ranges)
#' }
#' @export
gff2irange <- function(gff, type=NULL) {
    ret <- NULL
    annotations <- try(rtracklayer::import.gff3(gff), silent=TRUE)
    if (class(annotations) == 'try-error') {
        annotations <- try(rtracklayer::import.gff2(gff), silent=TRUE)
        if (class(annotations) == 'try-error') {
            stop("Could not extract the widths from the gff file.")
        } else {
            ret <- annotations
        }
    } else {
        ret <- annotations
    }
    ## The call to as.data.frame must be specified with the GenomicRanges namespace,
    ## otherwise one gets an error about no method to coerce an S4 class to a vector.
    if (!is.null(type)) {
        index <- ret[, "type"] == type
        ret <- ret[index, ]
    }
    return(ret)
}

#' Create a simple df from a gff which contains tooltips usable in googleVis graphs.
#'
#' The tooltip column is a handy proxy for more thorough anontations information when it would
#' otherwise be too troublesome to acquire.
#'
#' @param annotations Either a gff file or annotation data frame (which likely came from a gff file.).
#' @param desc_col Gff column from which to gather data.
#' @param type Gff type to use as the master key.
#' @param id_col Which annotation column to cross reference against?
#' @param ... Extra arguments dropped into arglist.
#' @return Df of tooltip information or name of a gff file.
#' @seealso \pkg{googleVis}
#'  \code{\link{gff2df}}
#' @examples
#' \dontrun{
#'  tooltips <- make_tooltips('reference/gff/saccharomyces_cerevisiae.gff.gz')
#' }
#' @export
make_tooltips <- function(annotations, desc_col='description', type="gene", id_col="ID", ...) {
    arglist <- list(...)
    tooltip_data <- NULL
    if (class(annotations) == "character") {
        tooltip_data <- gff2df(gff=annotations, type=type)
    } else if (class(annotations) == "data.frame") {
        tooltip_data <- annotations
    } else {
        stop("This requires either a filename or data frame.")
    }
    if (is.null(tooltip_data[[id_col]])) {
        tooltip_data[["ID"]] <- rownames(tooltip_data)
    }

    ## Attempt to use multiple columns if a c() was given
    tooltip_data[["tooltip"]] <- tooltip_data[[id_col]]
    for (col in desc_col) {
        if (is.null(tooltip_data[[col]])) {
            message(paste0("The column ", col, " is null, not using it."))
        } else {
            tooltip_data[["tooltip"]] <- paste0(tooltip_data[["tooltip"]], ": ", tooltip_data[[col]])
        }
    }
    tooltip_data[["tooltip"]] <- gsub("\\+", " ", tooltip_data[["tooltip"]])
    tooltip_data[["tooltip"]] <- gsub(": $", "", tooltip_data[["tooltip"]])
    tooltip_data[["tooltip"]] <- gsub("^: ", "", tooltip_data[["tooltip"]])
    rownames(tooltip_data) <- make.names(tooltip_data[[id_col]], unique=TRUE)
    ## Now remove extraneous columns
    tooltip_data <- tooltip_data[, c("tooltip"), drop=FALSE]
    colnames(tooltip_data) <- c("1.tooltip")
    return(tooltip_data)
}

#' Find how many times a given pattern occurs in every gene of a genome.
#'
#' There are times when knowing how many times a given string appears in a genome/CDS is helpful.
#' This function provides that information and is primarily used by cp_seq_m().
#'
#' @param fasta Genome sequence.
#' @param gff Gff of annotation information from which to acquire CDS (if not provided it will just
#'  query the entire genome).
#' @param pattern What to search for? This was used for tnseq and TA is the mariner insertion point.
#' @param type Column to use in the gff file.
#' @param key What type of entry of the gff file to key from?
#' @return Data frame of gene names and number of times the pattern appears/gene.
#' @seealso \pkg{Biostrings} \pkg{Rsamtools} \pkg{Rsamtools}
#'  \code{\link[Rsamtools]{FaFile}} \code{\link[Biostrings]{getSeq}} \code{\link[Biostrings]{PDict}}
#'  \code{\link[Biostrings]{vcountPDict}}
#' @examples
#' \dontrun{
#'  num_pattern = pattern_count_genome('mgas_5005.fasta', 'mgas_5005.gff')
#' }
#' @export
pattern_count_genome <- function(fasta, gff=NULL, pattern='TA', type='gene', key='locus_tag') {
    rawseq <- Rsamtools::FaFile(fasta)
    if (is.null(gff)) {
        entry_sequences <- rawseq
    } else {
        entries <- rtracklayer::import.gff3(gff, asRangedData=FALSE)
        ## type_entries <- subset(entries, type==type)
        type_entries <- entries[entries$type == type, ]
        names(type_entries) <- rownames(type_entries)
        entry_sequences <- Biostrings::getSeq(rawseq, type_entries)
        names(entry_sequences) <- entry_sequences[[, key]]
    }
    dict <- Biostrings::PDict(pattern, max.mismatch=0)
    result <- Biostrings::vcountPDict(dict, entry_sequences)
    num_pattern <- data.frame(
        "name" = names(entry_sequences),
        "num" = as.data.frame(t(result)))
    return(num_pattern)
}

#' Gather some simple sequence attributes.
#'
#' This extends the logic of the pattern searching in pattern_count_genome() to search on some other
#' attributes.
#'
#' @param fasta Genome encoded as a fasta file.
#' @param gff Optional gff of annotations (if not provided it will just ask the whole genome).
#' @param type Column of the gff file to use.
#' @param key What type of entry of the gff file to key from?
#' @return List of data frames containing gc/at/gt/ac contents.
#' @seealso \pkg{Biostrings} \pkg{Rsamtools}
#'  \code{\link[Rsamtools]{FaFile}} \code{\link[Biostrings]{getSeq}}
#' @examples
#' \dontrun{
#'  num_pattern = sequence_attributes('mgas_5005.fasta', 'mgas_5005.gff')
#' }
#' @export
sequence_attributes <- function(fasta, gff=NULL, type='gene', key='locus_tag') {
    rawseq <- Rsamtools::FaFile(fasta)
    if (is.null(gff)) {
        entry_sequences <- rawseq
    } else {
        entries <- rtracklayer::import.gff3(gff, asRangedData=FALSE)
        ## type_entries <- subset(entries, type==type)
        type_entries <- entries[entries[["type"]] == type, ]
        ##names(type_entries) <- rownames(type_entries)
        entry_sequences <- Biostrings::getSeq(rawseq, type_entries)
        names(entry_sequences) <- type_entries[["Name"]]
    }
    attribs <- data.frame(
        "gc" = Biostrings::letterFrequency(entry_sequences, "CG", as.prob=TRUE),
        "at" = Biostrings::letterFrequency(entry_sequences, "AT", as.prob=TRUE),
        "gt" = Biostrings::letterFrequency(entry_sequences, "GT", as.prob=TRUE),
        "ac" = Biostrings::letterFrequency(entry_sequences, "AC", as.prob=TRUE))
    rownames(attribs) <- type_entries[["locus_tag"]]
    colnames(attribs) <- c("gc","at","gt","ac")
    return(attribs)
}

## EOF
