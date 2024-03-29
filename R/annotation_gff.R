## annotation_gff.r: Extract annotations from GFF files.  When mapping, we
## often/usually use the gene tagged entries from a gff file for counting along
## with one of the annotation types.  Sadly, gff files are quite a mess, so this
## tries to work around some of the more likely corner cases.  The results from
## these functions are likely to be used along with another annotation source in
## order to get a more complete view of the genome/transcriptome of interest.

#' Rewrite a gff file as a granges with full seqinfo if possible.
#'
#' @param gff Input gff file.
#' @param type Feature type to extract.
#' @param type_column Tag from the gff file to use when extracting the type.
#' @export
gff2gr <- function(gff, type = NULL, type_column = "type") {
  chr_entries <- read.delim(file = gff, header = FALSE, sep = "")
  contigs <- chr_entries[["V1"]] == "##sequence-region"
  contigs <- chr_entries[contigs, c("V2", "V3", "V4")]
  colnames(contigs) <- c("ID", "start", "end")
  contig_info <- data.frame(
    "chrom" = contigs[["ID"]],
    "length" = as.numeric(contigs[["end"]]),
    "is_circular" = NA,
    stringsAsFactors = FALSE)
  rownames(contig_info) <- make.names(contig_info[["chrom"]], unique = TRUE)

  ## Dump a granges object and save it as an rda file.
  granges_result <- rtracklayer::import.gff3(gff)
  name_order <- names(Biostrings::seqinfo(granges_result))
  contig_info <- contig_info[name_order, ]
  length_vector <- contig_info[["length"]]
  GenomeInfoDb::seqlengths(granges_result) <- length_vector
  if (!is.null(type)) {
    start <- length(granges_result)
    type_idx <- mcols(granges_result)[[type_column]] == type
    granges_result <- granges_result[type_idx]
    end <- length(granges_result)
    message("Subsetting for type ", type, " in column ", type_column,
            " reduces the grange from: ", start, " to ", end, ".")
  }
  return(granges_result)
}

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

  attempts <- list(
    "rtracklayer::import.gff3" = list(
      "con" = gff,
      "sequenceRegionsAsSeqinfo" = TRUE),
    "rtracklayer::import.gff3" = list(
      "con" = gff,
      "sequenceRegionsAsSeqinfo" = FALSE),
    "rtracklayer::import.gff2" = list(
      "con" = gff,
      "sequenceRegionsAsSeqinfo" = TRUE),
    "rtracklayer::import.gff2" = list(
      "con" = gff,
      "sequenceRegionsAsSeqinfo" = FALSE),
    "rtracklayer::import.gff" = list(
      "con" = gff,
      "format" = "gtf"),
    "rtracklayer::import.gff" = list(
      "con" = gff))

  result <- NULL
  for (a in seq_along(attempts)) {
    attempt <- names(attempts)[a]
    split_name <- strsplit(attempt, "::")[[1]]
    fun_call <- get(split_name[[2]], asNamespace(split_name[[1]]))
    args <- attempts[[a]]
    result <- try(do.call(what = fun_call, args = args, quote = TRUE), silent = TRUE)
    if ("try-error" %in% class(result)) {
      result <- NULL
      next
    } else {
      break
    }
  }
  if (is.null(result)) {
    mesg("No gff import methods were able to read this file.")
    return(NULL)
  }
  annotations <- as.data.frame(result)[[id_col]]
  if (is.null(annotations)) {
    annotations <- as.data.frame(result)[[second_id_col]]
    if (is.null(annotations)) {
      warning("Attempting to create a dataframe with ", id_col, " and ",
              second_id_col, " both failed.")
      return(NULL)
    } else {
      annot <- annotations
    }
  } else {
    annot <- annotations
  }

  ret <- NULL
  if ("GRanges" %in% class(result) && ret_type == "data.frame") {
    ret <- GenomicRanges::as.data.frame(result)
    if (!is.null(type)) {
      index <- ret[, "type"] == type
      ret <- ret[index, ]
    }
    message("Returning a df with ", ncol(ret), " columns and ",
            nrow(ret), " rows.")
  } else {
    message("Returning a GRanges with ", ncol(ret), " columns and ",
            nrow(ret), " rows.")
    return(result)
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

merge_gff_children <- function(gff, grandparent = "gene", parent = "mRNA",
                               parent_tag = "Parent", id_tag = "ID",
                               children = c("CDS", "exon", "three_prime_UTR", "five_prime_UTR")) {
  grandparents <- load_gff_annotations(gff, type = grandparent, id_col = id_tag)
  colnames(grandparents) <- paste0("toplevel_", colnames(grandparents))
  parents <- load_gff_annotations(gff, type = parent, id_col = id_tag)
  colnames(parents) <- paste0("secondlevel_", colnames(parents))
  by_gp <- paste0("toplevel_", id_tag)
  by_parent <- paste0("secondlevel_", parent_tag)
  all <- merge(grandparents, parents, by.x = by_gp,
               by.y = by_parent)
  for (child in children) {
    mesg("Merging in the data for the ", child, " type.")
    child_annot <- load_gff_annotations(gff, type = child, id_col = id_tag)
    colnames(child_annot) <- paste0(child, "_", colnames(child_annot))
    by_parent <- paste0("secondlevel_", id_tag)
    by_child <- paste0(child, "_", parent_tag)
    all <- merge(all, child_annot, by.x = by_parent, by.y = by_child)
  }
  return(all)
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
#' @param id_col Column containing the gene IDs.
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
                                 type = "gene", id_col = "ID", key = NULL) {
  rawseq <- Rsamtools::FaFile(fasta)
  if (is.null(key)) {
    key <- c("ID", "locus_tag")
  }
  if (is.null(gff)) {
    entry_sequences <- Biostrings::getSeq(rawseq)
  } else {
    ## Note, I use import.gff here instead of import.gff3 or one of
    ## the other more sensible functions because import.gff sets the
    ## seqnames to match the chromosome; I have not figured out how to
    ## do that properly using import.gff3 without some annoying shenanigans.
    entries <- rtracklayer::import.gff(gff)
    meta <- mcols(entries)
    wanted <- meta[["type"]] == type
    ## keep only the ones of the type of interest (gene).
    type_entries <- entries[wanted, ]
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
