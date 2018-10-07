gather_eupath_utrs_padding <- function(species_name="Leishmania major", webservice="tritrypdb", ...) {
  metadata <- sm(download_eupath_metadata(webservice=webservice))
  pkg_names <- get_eupath_pkgnames(species=species_name, metadata=metadata)
  bsgenome_name <- pkg_names[["bsgenome"]]
  orgdb_name <- pkg_names[["orgdb"]]
  if (!isTRUE(pkg_names[["bsgenome_installed"]])) {
    genome_installedp <- make_eupath_bsgenome(species=species_name, metadata=metadata, ...)
  }
  if (!isTRUE(pkg_names[["orgdb_installed"]])) {
    orgdb_installedp <- make_eupath_orgdb(species=species_name, metadata=metadata, ...)
  }

  lib_result <- sm(library(orgdb_name, character.only=TRUE))
  lib_result <- sm(library(bsgenome_name, character.only=TRUE))
  orgdb <- get0(orgdb_name)
  bsgenome <- get0(bsgenome_name)
  wanted_fields <- c("annot_gene_location_text",
                     "annot_gene_name",
                     "annot_gene_product",
                     "annot_gene_type")
  annot <- load_orgdb_annotations(orgdb, keytype="gid", fields=wanted_fields)
  annot <- extract_gene_locations(annot[["genes"]])
  annot_complete <- complete.cases(annot)
  annot_df <- annot[annot_complete, ]

  result <- gather_utrs_padding(bsgenome, annot_df, ...)
  return(result)
}

#' Take a BSgenome and data frame of chr/start/end/strand, provide 5' and 3' padded sequence.
#'
#' For some species, we do not have a fully realized set of UTR boundaries, so
#' it can be useful to query some arbitrary and consistent amount of sequence
#' before/after every CDS sequence.  This function can provide that information.
#'
#' @param bsgenome  BSgenome object containing the genome of interest.
#' @param annot_df  Annotation data frame containing all the entries of
#'   interest, this is generally extracted using a function in the
#'   load_something_annotations() family (load_orgdb_annotations() being the most
#'   likely).
#' @param name_column  Give each gene a name using this column.
#' @param chr_column  Column name of the chromosome names.
#' @param start_column  Column name of the start information.
#' @param end_column  Ibid, end column.
#' @param strand_column  Ibid, strand.
#' @param type_column  Subset the annotation data using this column, if not null.
#' @param gene_type  Subset the annotation data using the type_column with this
#'   type.
#' @param padding  Return this number of nucleotides for each gene.
#' @param ... Arguments passed to child functions (I think none currently).
#' @return List of 2 elements, the 5' and 3' regions.
gather_utrs_padding <- function(bsgenome, annot_df, name_column="gid", chr_column="chromosome",
                                start_column="start", end_column="end", strand_column="strand",
                                type_column="annot_gene_type", gene_type="protein coding",
                                padding=120, ...) {

  if (!is.null(type_column)) {
    ## Pull the entries which are of the desired type.
    all_gene_idx <- annot_df[[type_column]] == gene_type
    annot_df <- annot_df[all_gene_idx, ]
  }

  ## These are named 'smaller' and 'larger' because this is extracting the ranges
  ## which are numerically smaller/larger than the CDS entry, not 5'/3' because
  ## this does not yet take into account the strand information.
  smaller <- GenomicRanges::GRanges(
                              seqnames=S4Vectors::Rle(annot_df[[chr_column]]),
                              ranges=IRanges::IRanges(start=annot_df[[start_column]] - padding,
                                                      end=annot_df[[end_column]] + 2),
                              strand=S4Vectors::Rle(annot_df[[strand_column]]),
                              name=S4Vectors::Rle(annot_df[[name_column]]))
  larger <- GenomicRanges::GRanges(
                             seqnames=S4Vectors::Rle(annot_df[[chr_column]]),
                             ranges=IRanges::IRanges(start=annot_df[[start_column]],
                                                     end=annot_df[[end_column]] + padding),
                             strand=S4Vectors::Rle(annot_df[[strand_column]]),
                             name=S4Vectors::Rle(annot_df[[name_column]]))
  smaller_seqstrings <- BSgenome::getSeq(bsgenome, smaller)
  larger_seqstrings <- Biostrings::getSeq(bsgenome, larger)
  names(smaller_seqstrings) <- smaller$name
  names(larger_seqstrings) <- larger$name

  ## These provide data frames of the sequence lexically before/after every gene.
  smaller_seqdf <- as.data.frame(smaller_seqstrings)
  colnames(smaller_seqdf) <- "sequence"
  larger_seqdf <- as.data.frame(larger_seqstrings)
  colnames(larger_seqdf) <- "sequence"
  ## These provide the strand and location information for the same.
  smaller_infodf <- as.data.frame(smaller)
  larger_infodf <- as.data.frame(larger)
  ## Bring together the location and sequence information.
  all_smaller <- merge(smaller_infodf, smaller_seqdf, by.x="name", by.y="row.names")
  all_larger <- merge(larger_infodf, larger_seqdf, by.x="name", by.y="row.names")

  ## To pull the 5' and 3' entries for the + strand, we use smaller and larger respectively.
  fivep_plus_idx <- all_smaller[["strand"]] == "+"
  fivep_plus_entries <- all_smaller[fivep_plus_idx, ]
  threep_plus_idx <- all_larger[["strand"]] == "+"
  threep_plus_entries <- all_larger[threep_plus_idx, ]
  ## In contrast, to pull the 5' and 3' entries from the - strand, use larger and smaller.
  fivep_minus_idx <- all_larger[["strand"]] == "-"
  fivep_minus_entries <- all_larger[fivep_minus_idx, ]
  threep_minus_idx <- all_smaller[["strand"]] == "-"
  threep_minus_entries <- all_smaller[threep_minus_idx, ]

  all_fivep_entries <- rbind(fivep_plus_entries, fivep_minus_entries)
  all_threep_entries <- rbind(threep_plus_entries, threep_minus_entries)

  retlist <- list(
    "five_prime" = all_fivep_entries,
    "three_prime" = all_threep_entries)
  return(retlist)
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
                                   ranges=IRanges::IRanges(start=fivep_df[[start_column]],
                                                           end=fivep_df[[end_column]] + 2),
                                   strand=S4Vectors::Rle(fivep_df[[strand_column]]),
                                   name=S4Vectors::Rle(fivep_df[[name_column]]))
  threep_ranges <- GenomicRanges::GRanges(
                                   seqnames=S4Vectors::Rle(threep_df[[chr_column]]),
                                   ranges=IRanges::IRanges(start=threep_df[[start_column]],
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
