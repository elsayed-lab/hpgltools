#' Count n-mers in a given data set using Biostrings
#'
#' This just calls PDict() and vcountPDict() on a sequence database given a
#' pattern and number of mismatches.  This may be used by divide_seq()
#' normalization.
#'
#' @param genome  Sequence database, genome in this case.
#' @param pattern  Count off this string.
#' @param mismatch  How many mismatches are acceptable?
#' @return Set of counts by sequence.
#' @export
count_nmer <- function(genome, pattern="ATG", mismatch=0) {
  if (class(genome)[1] == "character") {
    genome <- Rsamtools::FaFile(genome)
  }
  seq_obj <- Biostrings::getSeq(genome)
  dict <- Biostrings::PDict(pattern, max.mismatch=mismatch)
  result <- as.data.frame(Biostrings::vcountPDict(dict, seq_obj))
  rownames(result) <- pattern
  colnames(result) <- names(seq_obj)
  result <- t(result)
  return(result)
}

#' Given an eupathdb species lacking UTR boundaries, extract an arbitrary region
#' before/after each gene.
#'
#' This is a very domain-specific function.
#'
#' @param species_name Species name for which to query the eupathdb.
#' @param entry EuPathDB metadatum entry.
#' @param webservice If specified, makes the query faster, I always used
#'   tritrypdb.org.
#' @param padding Number of nucleotides to gather.
#' @param ... Extra arguments for the various EuPathDB functions.
#' @return Set of padding UTR sequences/coordinates.
gather_eupath_utrs_padding <- function(species_name="Leishmania major", entry=NULL,
                                       webservice="tritrypdb", padding=200, ...) {
  if (!is.null(entry)) {
    pkg_names <- EuPathDB::get_eupath_pkgnames(entry)
  } else {
    entry <- EuPathDB::get_eupath_entry(species_name, webservice=webservice)
    pkg_names <- EuPathDB::get_eupath_pkgnames(entry)
  }
  bsgenome_name <- pkg_names[["bsgenome"]]
  orgdb_name <- pkg_names[["orgdb"]]
  if (!isTRUE(pkg_names[["bsgenome_installed"]])) {
    genome_installedp <- EuPathDB::make_eupath_bsgenome(species=species_name, entry=entry,
                                                          ...)
  }
  if (!isTRUE(pkg_names[["orgdb_installed"]])) {
    orgdb_installedp <- EuPathDB::make_eupath_orgdb(entry=entry,
                                                    ...)
  }

  ##lib_result <- sm(library(orgdb_name, character.only=TRUE))
  lib_result <- sm(requireNamespace(orgdb_name))
  att_result <- sm(try(attachNamespace(orgdb_name), silent=TRUE))
  ##lib_result <- sm(library(bsgenome_name, character.only=TRUE))
  lib_result <- sm(requireNamespace(bsgenome_name))
  att_result <- sm(try(attachNamespace(bsgenome_name), silent=TRUE))
  orgdb <- get0(orgdb_name)
  bsgenome <- get0(bsgenome_name)
  wanted_fields <- c("annot_gene_location_text",
                     "annot_gene_name",
                     "annot_gene_product",
                     "annot_gene_type")
  annot <- load_orgdb_annotations(orgdb, keytype="gid", fields=wanted_fields)
  annot <- EuPathDB::extract_gene_locations(annot[["genes"]])
  annot_complete <- complete.cases(annot)
  annot_df <- annot[annot_complete, ]
  annot_df[["chr_length"]] <- 0

  ## Gather the chromosome lengths to make sure we don't pass them.
  genome_info <- BiocGenerics::as.data.frame(BSgenome::seqinfo(bsgenome))
  chr_names <- rownames(genome_info)
  for (l in 1:length(chr_names)) {
    name <- chr_names[l]
    len <- genome_info[l, "seqlengths"]
    chr_idx <- annot_df[["chromosome"]] == name
    annot_df[chr_idx, "chr_length"] <- len
  }

  result <- gather_utrs_padding(bsgenome, annot_df, padding=padding,
                                ...)
  return(result)
}

#' Take a BSgenome and data frame of chr/start/end/strand, provide 5' and 3'
#' padded sequence.
#'
#' For some species, we do not have a fully realized set of UTR boundaries, so
#' it can be useful to query some arbitrary and consistent amount of sequence
#' before/after every CDS sequence.  This function can provide that information.
#'
#' @param bsgenome BSgenome object containing the genome of interest.
#' @param annot_df Annotation data frame containing all the entries of
#'   interest, this is generally extracted using a function in the
#'   load_something_annotations() family (load_orgdb_annotations() being the most
#'   likely).
#' @param gid Specific GID(s) to query.
#' @param name_column Give each gene a name using this column.
#' @param chr_column Column name of the chromosome names.
#' @param start_column Column name of the start information.
#' @param end_column Ibid, end column.
#' @param strand_column Ibid, strand.
#' @param type_column Subset the annotation data using this column, if not null.
#' @param gene_type Subset the annotation data using the type_column with this
#'   type.
#' @param padding Return this number of nucleotides for each gene.
#' @param ... Arguments passed to child functions (I think none currently).
#' @return Dataframe of UTR, CDS, and UTR+CDS sequences.
gather_utrs_padding <- function(bsgenome, annot_df, gid=NULL, name_column="gid",
                                chr_column="chromosome", start_column="start",
                                end_column="end", strand_column="strand",
                                type_column="annot_gene_type",
                                gene_type="protein coding", padding=120, ...) {

  if (!is.null(type_column)) {
    ## Pull the entries which are of the desired type.
    all_gene_idx <- annot_df[[type_column]] == gene_type
    annot_df <- annot_df[all_gene_idx, ]
  }
  if (!is.null(gid)) {
    annot_df <- annot_df[gid, ]
  }
  one_idx <- annot_df[[strand_column]] == 1
  if (sum(one_idx) > 0) {
    annot_df[one_idx, strand_column] <- "+"
  }
  minusone_idx <- annot_df[[strand_column]] == -1
  if (sum(minusone_idx) > 0) {
    annot_df[minusone_idx, strand_column] <- "-"
  }

  ## Here is the thing I absolutely cannot remember:
  ## start is _always_ a lower number than end.
  annot_df[["low_boundary"]] <- annot_df[[start_column]] - padding
  annot_df[["high_boundary"]] <- annot_df[[end_column]] + padding
  ## Add some logic to make sure we do not have sequences which
  ## expand beyond the chromosome boundaries.
  sub_zero_idx <- annot_df[["low_boundary"]] < 1
  if (sum(sub_zero_idx) > 0) {
    message("Found ", sum(sub_zero_idx), " genes which are less than ",
            padding, " nt. from the beginning of the chromosome.")
    annot_df[sub_zero_idx, "low_boundary"] <- 1
  }
  past_end_idx <- annot_df[["high_boundary"]] > annot_df[["chr_length"]]
  if (sum(past_end_idx) > 0) {
    message("Found ", sum(past_end_idx), " genes which are less than ",
            padding, " nt. from the end of the chromosome.")
    annot_df[past_end_idx, "high_boundary"] <- annot_df[past_end_idx, "chr_length"]
  }

  plus_idx <- annot_df[[strand_column]] == "+"
  minus_idx <- annot_df[[strand_column]] == "-"
  plus_entries <- annot_df[plus_idx, ]
  minus_entries <- annot_df[minus_idx, ]

  pluses_fivep <- GenomicRanges::GRanges(
                                   seqnames=S4Vectors::Rle(plus_entries[[chr_column]]),
                                   ranges=IRanges::IRanges(
                                                     start=plus_entries[["low_boundary"]],
                                                     end=plus_entries[[start_column]]),
                                   strand=S4Vectors::Rle(plus_entries[[strand_column]]),
                                   name=S4Vectors::Rle(plus_entries[[name_column]]))
  pluses_threep <- GenomicRanges::GRanges(
                                    seqnames=S4Vectors::Rle(plus_entries[[chr_column]]),
                                    ranges=IRanges::IRanges(
                                                      start=plus_entries[[end_column]],
                                                      end=plus_entries[["high_boundary"]]),
                                    strand=S4Vectors::Rle(plus_entries[[strand_column]]),
                                    name=S4Vectors::Rle(plus_entries[[name_column]]))
  pluses_cds <- GenomicRanges::GRanges(
                                 seqnames=S4Vectors::Rle(plus_entries[[chr_column]]),
                                 ranges=IRanges::IRanges(
                                                   start=plus_entries[[start_column]],
                                                   end=plus_entries[[end_column]]),
                                 strand=S4Vectors::Rle(plus_entries[[strand_column]]),
                                 name=S4Vectors::Rle(plus_entries[[name_column]]))
  pluses_all <- GenomicRanges::GRanges(
                                 seqnames=S4Vectors::Rle(plus_entries[[chr_column]]),
                                 ranges=IRanges::IRanges(
                                                   start=plus_entries[["low_boundary"]],
                                                   end=plus_entries[["high_boundary"]]),
                                 strand=S4Vectors::Rle(plus_entries[[strand_column]]),
                                 name=S4Vectors::Rle(plus_entries[[name_column]]))
  plus_fivep_seqstrings <- BSgenome::getSeq(bsgenome, pluses_fivep)
  names(plus_fivep_seqstrings) <- pluses_fivep$name
  plus_threep_seqstrings <- BSgenome::getSeq(bsgenome, pluses_threep)
  names(plus_threep_seqstrings) <- pluses_threep$name
  plus_cds_seqstrings <- BSgenome::getSeq(bsgenome, pluses_cds)
  names(plus_cds_seqstrings) <- pluses_cds$name
  plus_all_seqstrings <- BSgenome::getSeq(bsgenome, pluses_all)
  names(plus_all_seqstrings) <- pluses_all$name

  minuses_fivep <- GenomicRanges::GRanges(
                                    seqnames=S4Vectors::Rle(minus_entries[[chr_column]]),
                                    ranges=IRanges::IRanges(
                                                      start=minus_entries[[end_column]],
                                                      end=minus_entries[["high_boundary"]]),
                                    strand=S4Vectors::Rle(minus_entries[[strand_column]]),
                                    name=S4Vectors::Rle(minus_entries[[name_column]]))
  minuses_threep <- GenomicRanges::GRanges(
                                     seqnames=S4Vectors::Rle(minus_entries[[chr_column]]),
                                     ranges=IRanges::IRanges(
                                                       start=minus_entries[["low_boundary"]],
                                                       end=minus_entries[[start_column]]),
                                     strand=S4Vectors::Rle(minus_entries[[strand_column]]),
                                     name=S4Vectors::Rle(minus_entries[[name_column]]))
  minuses_cds <- GenomicRanges::GRanges(
                                  seqnames=S4Vectors::Rle(minus_entries[[chr_column]]),
                                  ranges=IRanges::IRanges(
                                                    start=minus_entries[[start_column]],
                                                    end=minus_entries[[end_column]]),
                                  strand=S4Vectors::Rle(minus_entries[[strand_column]]),
                                  name=S4Vectors::Rle(minus_entries[[name_column]]))
  minuses_all <- GenomicRanges::GRanges(
                                  seqnames=S4Vectors::Rle(minus_entries[[chr_column]]),
                                  ranges=IRanges::IRanges(
                                                    start=minus_entries[["low_boundary"]],
                                                    end=minus_entries[["high_boundary"]]),
                                  strand=S4Vectors::Rle(minus_entries[[strand_column]]),
                                  name=S4Vectors::Rle(minus_entries[[name_column]]))
  minus_fivep_seqstrings <- BSgenome::getSeq(bsgenome, minuses_fivep)
  names(minus_fivep_seqstrings) <- minuses_fivep$name
  minus_threep_seqstrings <- BSgenome::getSeq(bsgenome, minuses_threep)
  names(minus_threep_seqstrings) <- minuses_threep$name
  minus_cds_seqstrings <- BSgenome::getSeq(bsgenome, minuses_cds)
  names(minus_cds_seqstrings) <- minuses_cds$name
  minus_all_seqstrings <- BSgenome::getSeq(bsgenome, minuses_all)
  names(minus_all_seqstrings) <- minuses_all$name


  ## These provide data frames of the sequence lexically before/after every gene.
  plus_seqdf <- as.data.frame(plus_fivep_seqstrings)
  plus_seqdf <- cbind(plus_seqdf, as.data.frame(plus_threep_seqstrings))
  plus_seqdf <- cbind(plus_seqdf, as.data.frame(plus_cds_seqstrings))
  plus_seqdf <- cbind(plus_seqdf, as.data.frame(plus_all_seqstrings))
  colnames(plus_seqdf) <- c("fivep", "threep", "cds", "all")

  minus_seqdf <- as.data.frame(minus_fivep_seqstrings)
  minus_seqdf <- cbind(minus_seqdf, as.data.frame(minus_threep_seqstrings))
  minus_seqdf <- cbind(minus_seqdf, as.data.frame(minus_cds_seqstrings))
  minus_seqdf <- cbind(minus_seqdf, as.data.frame(minus_all_seqstrings))
  colnames(minus_seqdf) <- c("fivep", "threep", "cds", "all")

  plus_infodf <- merge(plus_entries, plus_seqdf, by="row.names")
  rownames(plus_infodf) <- plus_infodf[[1]]
  plus_infodf[[1]] <- NULL
  minus_infodf <- merge(minus_entries, minus_seqdf, by="row.names")
  rownames(minus_infodf) <- minus_infodf[[1]]
  minus_infodf[[1]] <- NULL

  allinfo_df <- rbind(plus_infodf, minus_infodf)
  return(allinfo_df)
}

#' Get UTR sequences using information provided by TxDb and fiveUTRsByTranscript
#'
#' For species like Mus musculus, load_orgdb_annotations(Mus.musculus)
#' should return a list including the requisite GRanges for the 5'/3' UTRs.
#'
#' @param bsgenome A BSGenome instance containing the encoded genome.
#' @param fivep_utr  Locations of the 5' UTRs.
#' @param threep_utr  Locations of the 3' UTRs.
#' @param start_column  What column in the annotation data contains the starts?
#' @param strand_column  What column in the annotation data contains the sequence strands?
#' @param end_column  Column in the data with the end locations.
#' @param chr_column  Column in the df with the chromosome names.
#' @param name_column  Finally, where are the gene names?
#' @param ...  Parameters passed to child functions.
#' @return UTRs!
gather_utrs_txdb <- function(bsgenome, fivep_utr=NULL, threep_utr=NULL,
                             start_column="start", end_column="end",
                             strand_column="strand", chr_column="seqnames",
                             name_column="group_name", ...) {
  fivep_df <- as.data.frame(fivep_utr)
  threep_df <- as.data.frame(threep_utr)


  ## These are named 'smaller' and 'larger' because this is extracting the ranges
  ## which are numerically smaller/larger than the CDS entry, not 5'/3' because
  ## this does not yet take into account the strand information.
  fivep_ranges <- GenomicRanges::GRanges(
                                   seqnames=S4Vectors::Rle(fivep_df[[chr_column]]),
                                   ranges=IRanges::IRanges(
                                                     start=fivep_df[[start_column]],
                                                     end=fivep_df[[end_column]] + 2),
                                   strand=S4Vectors::Rle(fivep_df[[strand_column]]),
                                   name=S4Vectors::Rle(fivep_df[[name_column]]))
  threep_ranges <- GenomicRanges::GRanges(
                                    seqnames=S4Vectors::Rle(threep_df[[chr_column]]),
                                    ranges=IRanges::IRanges(
                                                      start=threep_df[[start_column]],
                                                      end=threep_df[[end_column]] + 2),
                                    strand=S4Vectors::Rle(threep_df[[strand_column]]),
                                    name=S4Vectors::Rle(threep_df[[name_column]]))
  fivep_seqstrings <- BSgenome::getSeq(bsgenome, fivep_ranges)
  threep_seqstrings <- BSgenome::getSeq(bsgenome, threep_ranges)
  names(fivep_seqstrings) <- fivep_ranges$name
  names(threep_seqstrings) <- threep_ranges$name

  ## These provide data frames of the sequence lexically before/after every gene.
  fivep_seqdt <- data.table::as.data.table(fivep_seqstrings)
  fivep_seqdt[["name"]] <- as.character(fivep_ranges$name)
  colnames(fivep_seqdt) <- c("sequence", "name")
  threep_seqdt <- data.table::as.data.table(threep_seqstrings)
  threep_seqdt[["name"]] <- as.character(threep_ranges$name)
  colnames(threep_seqdf) <- c("sequence", "name")
  ## These provide the strand and location information for the same.
  fivep_infodt <- data.table::as.data.table(fivep_ranges)
  threep_infodt <- data.table::as.data.table(threep_ranges)
  ## Bring together the location and sequence information.
  all_fivep <- cbind(fivep_infodt, fivep_seqdt)
  all_threep <- cbind(threep_infodt, threep_seqdt)

  retlist <- list(
    "five_prime" = all_fivep,
    "three_prime" = all_threep)
  return(retlist)
}

## EOF
