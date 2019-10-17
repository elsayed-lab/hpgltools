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
#' @seealso \pkg{rtracklayer} \link{load_gff_annotations} \pkg{Biostrings}
#'  \code{\link[rtracklayer]{import.gff}}
#' @examples
#' \dontrun{
#'  library(BSgenome.Tcruzi.clbrener.all)
#'  tc_clb_all <- BSgenome.Tcruzi.clbrener.all
#'  cds_ranges <- gff2irange('reference/gff/tcruzi_clbrener.gff.xz', type='CDS')
#'  cds_sequences <- Biostrings::getSeq(tc_clb_all, cds_ranges)
#' }
#' @author atb
#' @export
gff2irange <- function(gff, type=NULL) {
  ret <- NULL
  annotations <- try(rtracklayer::import.gff3(gff), silent=TRUE)
  if (class(annotations) == "try-error") {
    annotations <- try(rtracklayer::import.gff2(gff), silent=TRUE)
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
#' @param ret_type  Return a data.frame or something else?
#' @param second_id_col Second column to check.
#' @param try Give your own function call to use for importing.
#' @param row.names  Choose another column for setting the rownames of the data frame.
#' @return Dataframe of the annotation information found in the gff file.
#' @seealso \pkg{rtracklayer} \pkg{GenomicRanges}
#'  \code{\link[rtracklayer]{import.gff}}
#' @examples
#' \dontrun{
#'  funkytown <- load_gff_annotations('reference/gff/saccharomyces_cerevsiae.gff.xz')
#' }
#' @author atb
#' @export
load_gff_annotations <- function(gff, type=NULL, id_col="ID", ret_type="data.frame",
                                 second_id_col="locus_tag", try=NULL, row.names=NULL) {
  if (!file.exists(gff)) {
    stop("Unable to find the gff file: ", gff)
  }
  ret <- NULL
  attempts <- c("rtracklayer::import.gff3(gff, sequenceRegionsAsSeqinfo=TRUE)",
                "rtracklayer::import.gff3(gff, sequenceRegionsAsSeqinfo=FALSE)",
                "rtracklayer::import.gff2(gff, sequenceRegionsAsSeqinfo=TRUE)",
                "rtracklayer::import.gff2(gff, sequenceRegionsAsSeqinfo=FALSE)",
                "rtracklayer::import.gff(gff, format='gtf')",
                "rtracklayer::import.gff(gff)")
  if (!is.null(try)) {
    ##  attempts <- c(paste0(try, "(gff)"), attempts)
    attempts <- c(glue("{try}(gff)"), attempts)
  }

  annot <- NULL
  for (att in 1:length(attempts)) {
    annotations <- NULL
    message("Trying attempt: ", attempts[[att]])
    attempt <- attempts[[att]]
    eval_string <- glue("annotations <- try({attempt}, silent=TRUE)")
    eval(parse(text=eval_string))
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
  if (class(annot)[[1]] == "GRanges" & ret_type == "data.frame") {
    ret <- GenomicRanges::as.data.frame(annot)
    rm(annot)
    if (!is.null(type)) {
      index <- ret[, "type"] == type
      ret <- ret[index, ]
    }
    message("Returning a df with ", ncol(ret), " columns and ", nrow(ret), " rows.")
  } else if (class(annot)[[1]] == "GRanges" & ret_type == "GRanges") {
    ret <- annot
    rm(annot)
    message("Returning a GRanges with ", ncol(ret), " columns and ", nrow(ret), " rows.")
  } else {
    stop("Unable to load gff file.")
  }

  ## Sometimes we get some pretty weird data structures inside the gff data, the
  ## following two blocks should theoretically simplify the resulting output.
  for (col in colnames(ret)) {
    ret[[col]] <- unAsIs(ret[[col]])
  }

  for (col in colnames(ret)) {
    if (class(ret[[col]]) == "list") {
      ret[[col]] <- sapply(X=ret[[col]], FUN=unlist)
      ret[[col]] <- sapply(X=ret[[col]], FUN=toString)
    } else if (class(ret[[col]]) == "factor") {
      ret[[col]] <- as.character(ret[[col]])
    }
  }

  if (!is.null(row.names)) {
    rownames(ret) <- ret[[row.names]]
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
#'   provided it will just query the entire genome).
#' @param pattern What to search for? This was used for tnseq and TA is the
#'   mariner insertion point.
#' @param type Column to use in the gff file.
#' @param key What type of entry of the gff file to key from?
#' @return Data frame of gene names and number of times the pattern appears/gene.
#' @seealso \pkg{Biostrings} \pkg{Rsamtools} \pkg{Rsamtools}
#'  \code{\link[Rsamtools]{FaFile}} \code{\link[Biostrings]{getSeq}}
#'  \code{\link[Biostrings]{PDict}} \code{\link[Biostrings]{vcountPDict}}
#' @examples
#' \dontrun{
#'  num_pattern <- pattern_count_genome('mgas_5005.fasta', 'mgas_5005.gff')
#' }
#' @author atb
#' @export
pattern_count_genome <- function(fasta, gff=NULL, pattern="TA", type="gene", key=NULL) {
  rawseq <- Rsamtools::FaFile(fasta)
  if (is.null(key)) {
    key <- c("ID", "locus_tag")
  }
  if (is.null(gff)) {
    entry_sequences <- Biostrings::getSeq(rawseq)
  } else {
    ## entries <- rtracklayer::import.gff3(gff, asRangedData=FALSE)
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
  dict <- Biostrings::PDict(pattern, max.mismatch=0)
  result <- Biostrings::vcountPDict(dict, entry_sequences)
  num_pattern <- data.frame(
    "name" = names(entry_sequences),
    "num" = as.data.frame(t(result)),
    stringsAsFactors=FALSE)
  colnames(num_pattern) <- c("name", "number")
  return(num_pattern)
}

#' Gather some simple sequence attributes.
#'
#' This extends the logic of the pattern searching in pattern_count_genome() to
#' search on some other attributes.
#'
#' @param fasta Genome encoded as a fasta file.
#' @param gff Optional gff of annotations (if not provided it will just ask the
#'   whole genome).
#' @param type Column of the gff file to use.
#' @param key What type of entry of the gff file to key from?
#' @return List of data frames containing gc/at/gt/ac contents.
#' @seealso \pkg{Biostrings} \pkg{Rsamtools}
#'  \code{\link[Rsamtools]{FaFile}} \code{\link[Biostrings]{getSeq}}
#' @examples
#' \dontrun{
#'  num_pattern = sequence_attributes('mgas_5005.fasta', 'mgas_5005.gff')
#' }
#' @author atb
#' @export
sequence_attributes <- function(fasta, gff=NULL, type="gene", key=NULL) {
  rawseq <- Rsamtools::FaFile(fasta)
  if (is.null(key)) {
    key <- c("ID", "locus_tag")
  }
  if (is.null(gff)) {
    entry_sequences <- Biostrings::getSeq(rawseq)
  } else {
    ## entries <- rtracklayer::import.gff3(gff, asRangedData=FALSE)
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
  attribs <- data.frame(
    "gc" = Biostrings::letterFrequency(entry_sequences, "CG", as.prob=TRUE),
    "at" = Biostrings::letterFrequency(entry_sequences, "AT", as.prob=TRUE),
    "gt" = Biostrings::letterFrequency(entry_sequences, "GT", as.prob=TRUE),
    "ac" = Biostrings::letterFrequency(entry_sequences, "AC", as.prob=TRUE),
    stringsAsFactors=FALSE)
  rownames(attribs) <- names(entry_sequences)
  colnames(attribs) <- c("gc", "at", "gt", "ac")
  return(attribs)
}

#' Given a data frame of exon counts and annotation information, sum the exons.
#'
#' This function will merge a count table to an annotation table by the child column.
#' It will then sum all rows of exons by parent gene and sum the widths of the exons.
#' Finally it will return a list containing a df of gene lengths and summed counts.
#'
#' @param data Count tables of exons.
#' @param gff Gff filename.
#' @param annotdf Dataframe of annotations (probably from load_gff_annotations).
#' @param parent Column from the annotations with the gene names.
#' @param child Column from the annotations with the exon names.
#' @return List of 2 data frames, counts and lengths by summed exons.
#' @seealso \pkg{rtracklayer}
#'  \code{\link{load_gff_annotations}}
#' @examples
#' \dontrun{
#' summed <- sum_exons(counts, gff='reference/xenopus_laevis.gff.xz')
#' }
#' @author Keith Hughitt with some modifications by atb.
#' @export
sum_exon_widths <- function(data=NULL, gff=NULL, annotdf=NULL,
                            parent="Parent", child="row.names") {
  if (is.null(annotdf) & is.null(gff)) {
    stop("I need either a df with parents, children, and widths; or a gff filename.")
  } else if (is.null(annotdf)) {
    annotdf <- load_gff_annotations(gff)
  }

  tmp_data <- merge(data, annotdf, by=child)
  rownames(tmp_data) <- tmp_data[["Row.names"]]
  tmp_data <- tmp_data[-1]
  ## Start out by summing the gene widths
  column <- aggregate(tmp_data[, "width"], by=list(Parent=tmp_data[, parent]), FUN=sum)
  new_data <- data.frame(column[["x"]], stringsAsFactors=FALSE)
  rownames(new_data) <- column[["Parent"]]
  colnames(new_data) <- c("width")

  for (c in 1:length(colnames(data))) {
    column_name <- colnames(data)[[c]]
    column <- aggregate(tmp_data[, column_name], by=list(Parent=tmp_data[, parent]), FUN=sum)
    rownames(column) <- column[["Parent"]]
    new_data <- cbind(new_data, column[["x"]])
  } ## End for loop

  width_df <- data.frame(new_data[["width"]], stringsAsFactors=FALSE)
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
