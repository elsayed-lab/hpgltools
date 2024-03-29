## annotation_gff.r: Extract annotations from GFF files.  When mapping, we
## often/usually use the gene tagged entries from a gff file for counting along
## with one of the annotation types.  Sadly, gff files are quite a mess, so this
## tries to work around some of the more likely corner cases.  The results from
## these functions are likely to be used along with another annotation source in
## order to get a more complete view of the genome/transcriptome of interest.

#' Extract annotation information from a gff file into an irange object.
#'
#' Try to make import.gff a little more robust; I acquire (hopefully) valid gff
#' files from various sources: yeastgenome.org, microbesonline, tritrypdb, ucsc,
#' ncbi. To my eyes, they all look like reasonably good gff3 files, but some of
#' them must be loaded with import.gff2, import.gff3, etc.  That is super
#' annoying. Also, I pretty much always just do as.data.frame() when I get
#' something valid from rtracklayer, so this does that for me, I have another
#' function which returns the iranges etc.  This function wraps
#' import.gff/import.gff3/import.gff2 calls in try() because sometimes those
#' functions fail in unpredictable ways.
#'
#' This is essentially load_gff_annotations(), but returns data suitable for
#' getSet()  This is another place which should be revisited for improvements
#' via mcols().  Check snp.r. for ideas.
#'
#' @param gff Gff filename.
#' @param type Subset to extract.
#' @return Iranges! (useful for getSeq().)
#' @seealso [rtracklayer] [load_gff_annotations()]
#'  \code{\link[rtracklayer]{import.gff}}
#' @examples
#'  example_gff <- system.file("share", "gas.gff", package = "hpgldata")
#'  gas_iranges <- gff2irange(example_gff)
#'  colnames(as.data.frame(gas_iranges))
#' @export
gff2irange <- function(gff, type = NULL) {
  ret <- NULL
  annotations <- try(rtracklayer::import.gff3(gff), silent = TRUE)
  if (class(annotations) == "try-error") {
    annotations <- try(rtracklayer::import.gff2(gff), silent = TRUE)
    if (class(annotations) == "try-error") {
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

#' Extract annotation information from a gff file into a df
#'
#' Try to make import.gff a little more robust; I acquire (hopefully) valid gff
#' files from various sources: yeastgenome.org, microbesonline, tritrypdb, ucsc,
#' ncbi. To my eyes, they all look like reasonably good gff3 files, but some of
#' them must be loaded with import.gff2, import.gff3, etc. That is super
#' annoying. Also, I pretty much always just do as.data.frame() when I get
#' something valid from rtracklayer, so this does that for me, I have another
#' function which returns the iranges etc.  This function wraps
#' import.gff/import.gff3/import.gff2 calls in try() because sometimes those
#' functions fail in unpredictable ways.
#'
#' @param gff Gff filename.
#' @param type Subset the gff file for entries of a specific type.
#' @param id_col Column in a successful import containing the IDs of interest.
#' @param ret_type Return a data.frame or something else?
#' @param second_id_col Second column to check.
#' @param try Give your own function call to use for importing.
#' @param row.names Choose another column for setting the rownames of the data frame.
#' @return Dataframe of the annotation information found in the gff file.
#' @seealso [rtracklayer] [GenomicRanges]
#' @examples
#'  example_gff <- system.file("share", "gas.gff", package = "hpgldata")
#'  gas_gff_annot <- load_gff_annotations(example_gff)
#'  dim(gas_gff_annot)
#' @export
load_gff_annotations <- function(gff, type = NULL, id_col = "ID", ret_type = "data.frame",
                                 second_id_col = "locus_tag", try = NULL, row.names = NULL) {
  if (!file.exists(gff)) {
    stop("Unable to find the gff file: ", gff)
  }
  compressed <- FALSE
  newfile <- NULL
  original <- gff
  ext <- "gz"
  fun <- gzfile
  if (grepl(pattern = "\\.[x|g]z$", x = gff)) {
    compressed <- TRUE
    newfile <- gsub(pattern = "\\.[x|g]z$", replacement = "", x = gff)
    ext <- gsub(pattern = ".*\\.([x|g]z)$", replacement = "\\1", x = gff)
    if (grepl(x = ext, pattern = "^x")) {
      fun <- xzfile
    }
    uncomp <- R.utils::decompressFile(filename = gff, destname = newfile, ext = ext,
                                      FUN = fun, remove = FALSE)
    gff <- newfile
  }

  ret <- NULL
  attempts <- c("rtracklayer::import.gff3(gff, sequenceRegionsAsSeqinfo = TRUE)",
                "rtracklayer::import.gff3(gff, sequenceRegionsAsSeqinfo = FALSE)",
                "rtracklayer::import.gff2(gff, sequenceRegionsAsSeqinfo = TRUE)",
                "rtracklayer::import.gff2(gff, sequenceRegionsAsSeqinfo = FALSE)",
                "rtracklayer::import.gff(gff, format='gtf')",
                "rtracklayer::import.gff(gff)")
  if (!is.null(try)) {
    attempts <- c(glue("{try}(gff)"), attempts)
  }

  annot <- NULL
  for (att in seq_along(attempts)) {
    annotations <- NULL
    message("Trying attempt: ", attempts[[att]])
    attempt <- attempts[[att]]
    eval_string <- glue("annotations <- try({attempt}, silent = TRUE)")
    eval(parse(text = eval_string))
    if (class(annotations) == "try-error") {
      rm(annotations)
    } else if (is.null(GenomicRanges::as.data.frame(annotations)[[id_col]]) &
               is.null(GenomicRanges::as.data.frame(annotations)[[second_id_col]])) {
      rm(annotations)
    } else {
      annot <- annotations
      rm(annotations)
      message("Had a successful gff import with ", attempt)
      break
    }
  }
  ret <- NULL
  if (class(annot)[[1]] == "GRanges" && ret_type == "data.frame") {
    ret <- GenomicRanges::as.data.frame(annot)
    rm(annot)
    if (!is.null(type)) {
      index <- ret[, "type"] == type
      ret <- ret[index, ]
    }
    message("Returning a df with ", ncol(ret), " columns and ",
            nrow(ret), " rows.")
  } else if (class(annot)[[1]] == "GRanges" && ret_type == "GRanges") {
    ret <- annot
    rm(annot)
    message("Returning a GRanges with ", ncol(ret), " columns and ",
            nrow(ret), " rows.")
  } else {
    stop("Unable to load gff file.")
  }

  ## Sometimes we get some pretty weird data structures inside the gff
  ## data, the following two blocks should theoretically simplify the
  ## resulting output.
  for (col in colnames(ret)) {
    ret[[col]] <- unAsIs(ret[[col]])
  }

  for (col in colnames(ret)) {
    if (class(ret[[col]]) == "list") {
      ret[[col]] <- sapply(X = ret[[col]], FUN = unlist)
      ret[[col]] <- sapply(X = ret[[col]], FUN = toString)
    } else if (class(ret[[col]]) == "factor") {
      ret[[col]] <- as.character(ret[[col]])
    }
  }

  if (!is.null(row.names)) {
    rownames(ret) <- ret[[row.names]]
  }
  if (isTRUE(compressed)) {
    removed <- file.remove(gff)
  }

  return(ret)
}

#' Find how many times a given pattern occurs in every gene of a genome.
#'
#' There are times when knowing how many times a given string appears in a
#' genome/CDS is helpful. This function provides that information and is
#' primarily used by cp_seq_m().
#'
#' This is once again a place where mcols() usage might improve the overall
#' quality of life.
#'
#' @param fasta Genome sequence.
#' @param gff Gff of annotation information from which to acquire CDS (if not
#'  provided it will just query the entire genome).
#' @param pattern What to search for? This was used for tnseq and TA is the
#'  mariner insertion point.
#' @param type Column to use in the gff file.
#' @param key What type of entry of the gff file to key from?
#' @return Data frame of gene names and number of times the pattern appears/gene.
#' @seealso [Biostrings] [Rsamtools::FaFile()] [Biostrings::PDict()]
#' @examples
#'  pa_data <- get_paeruginosa_data()
#'  pa_fasta <- pa_data[["fasta"]]
#'  pa_gff <- pa_data[["gff"]]
#'  ta_count <- pattern_count_genome(pa_fasta, pa_gff)
#'  head(ta_count)
#' @export
pattern_count_genome <- function(fasta, gff = NULL, pattern = "TA",
                                 type = "gene", key = NULL) {
  rawseq <- Rsamtools::FaFile(fasta)
  if (is.null(key)) {
    key <- c("ID", "locus_tag")
  }
  if (is.null(gff)) {
    entry_sequences <- Biostrings::getSeq(rawseq)
  } else {
    ## entries <- rtracklayer::import.gff3(gff, asRangedData = FALSE)
    entries <- rtracklayer::import.gff(gff)
    ## keep only the ones of the type of interest (gene).
    type_entries <- entries[entries$type == type, ]
    ## Set some hopefully sensible names.
    names(type_entries) <- rownames(type_entries)
    ## Get the sequence from the genome for them.
    entry_sequences <- Biostrings::getSeq(rawseq, type_entries)
    tmp_entries <- as.data.frame(type_entries)
    ## Give them some sensible names
    for (k in key) {
      if (!is.null(tmp_entries[[k]])) {
        names(entry_sequences) <- tmp_entries[[k]]
      }
    }
  }
  dict <- Biostrings::PDict(pattern, max.mismatch = 0)
  result <- Biostrings::vcountPDict(dict, entry_sequences)
  num_pattern <- data.frame(
      "name" = names(entry_sequences),
      "num" = as.data.frame(t(result)),
      stringsAsFactors = FALSE)
  colnames(num_pattern) <- c("name", "number")
  class(num_pattern) <- "pattern_counted"
  return(num_pattern)
}

#' Given a data frame of exon counts and annotation information, sum the exons.
#'
#' This function will merge a count table to an annotation table by
#' the child column. It will then sum all rows of exons by parent gene
#' and sum the widths of the exons. Finally it will return a list
#' containing a df of gene lengths and summed counts.
#'
#' @param data Count tables of exons.
#' @param gff Gff filename.
#' @param annotdf Dataframe of annotations (probably from load_gff_annotations).
#' @param parent Column from the annotations with the gene names.
#' @param child Column from the annotations with the exon names.
#' @return List of 2 data frames, counts and lengths by summed exons.
#' @seealso [rtracklayer] [load_gff_annotations()]
#' @examples
#' \dontrun{
#'  summed <- sum_exon_widths(counts, gff = "reference/xenopus_laevis.gff.xz")
#' }
#' @author Keith Hughitt with some modifications by atb.
#' @export
sum_exon_widths <- function(data = NULL, gff = NULL, annotdf = NULL,
                            parent = "Parent", child = "row.names") {
  if (is.null(annotdf) & is.null(gff)) {
    stop("I need either a df with parents, children, and widths; or a gff filename.")
  } else if (is.null(annotdf)) {
    annotdf <- load_gff_annotations(gff)
  }

  tmp_data <- merge(data, annotdf, by = child)
  rownames(tmp_data) <- tmp_data[["Row.names"]]
  tmp_data <- tmp_data[-1]
  ## Start out by summing the gene widths
  column <- aggregate(tmp_data[, "width"],
                      by = list(Parent = tmp_data[, parent]), FUN = sum)
  new_data <- data.frame(column[["x"]], stringsAsFactors = FALSE)
  rownames(new_data) <- column[["Parent"]]
  colnames(new_data) <- c("width")

  for (c in seq_along(colnames(data))) {
    column_name <- colnames(data)[[c]]
    column <- aggregate(tmp_data[, column_name],
                        by = list(Parent = tmp_data[, parent]), FUN = sum)
    rownames(column) <- column[["Parent"]]
    new_data <- cbind(new_data, column[["x"]])
  } ## End for loop

  width_df <- data.frame(new_data[["width"]], stringsAsFactors = FALSE)
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

## EOF
