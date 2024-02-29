#' Simplified TM calculator
#'
#' A quick and dirty TM calculator, taken from:
#' Taken from: https://www.biostars.org/p/58437/
#'
#' @param sequence String of atgc letters (not smart enough to do RNA).
cheap_tm <- function(sequence) {
  n <- nchar(sequence)
  ng <- stringr::str_count(toupper(sequence), "G")
  nc <- stringr::str_count(toupper(sequence), "C")
  n_gc <- ng + nc
  n_at <- n - n_gc
  tm <- 0
  if (n_gc < 18) {
    method <- "2*AT + 4*GC"
    tm <- ((2 * n_at) + (4 * n_gc)) - 7
  } else {
    method <- "64.9 + (41*(GC-16.4)/n)"
    tm <- 64.9 + (41 * ((n_gc - 16.4) / n))
  }
  ret <- list(
    "method" = method,
    "tm" = tm)
  class(ret) <- "cheap_tm"
  return(ret)
}

#' Given a named vector of fun regions, make a dataframe which
#' includes putative primers and the spec strings for expected
#' variants.
#'
#' This function came out of our TMRC2 work and seeks to provide an
#' initial set of potential PCR primers which are able to distinguish
#' between different aspects of the data.  In the actual data, we were
#' looking for differences between the zymodemes 2.2 and 2.3.
#'
#' @param vector variant-based set of putative regions with variants
#'  between conditions of interest.
#' @param max_primer_length given this length as a start, whittle down
#'  to a hopefully usable primer size.
#' @param topn Choose this number of variant regions from the rather
#'  larger set of possibilities..
#' @param bin_width Separate the genome into chunks of this size when
#'  hunting for primers, this size will therefore be the approximate
#'  PCR amplicon length.
#' @param genome (BS)Genome to search.
#' @param target_temp PCR temperature to attempt to match.
choose_sequence_regions <- function(long_variant_vector, max_primer_length = 45,
                                    topn = NULL, bin_width = 600,
                                    genome = NULL, target_temp = 58,
                                    min_gc_prop = 0.25,
                                    max_nmer_run = 5) {
  recount_positions <- function(position_string, bin_width = 600) {
    positions <- strsplit(position_string, ",")[[1]]
    positions <- gsub(x = positions, pattern = " ", replacement = "")
    positions <- as.numeric(positions)
    new_positions <- toString(as.character(bin_width - positions))
    return(new_positions)
  }
  rev_variants <- function(ref_string, position_string, alt_string) {
    refs <- gsub(pattern = " ", replacement = "", x = strsplit(ref_string, ",")[[1]])
    pos <- gsub(pattern = " ", replacement = "", x = strsplit(position_string, ",")[[1]])
    alts <- gsub(pattern = " ", replacement = "", x = strsplit(alt_string, ",")[[1]])
    new <- toString(paste0(refs, pos, alts))
    return(new)
  }

  run_pattern <- paste0("A{", max_nmer_run,
                        ",}|T{", max_nmer_run,
                        ",}|G{", max_nmer_run,
                        ",}|C{", max_nmer_run, ",}")

  ## Now get the nucleotides of the first 30 nt of each window
  sequence_df <- data.frame()
  if (is.null(topn)) {
    sequence_df <- as.data.frame(long_variant_vector)
  } else {
    sequence_df <- as.data.frame(head(long_variant_vector, n = topn))
  }
  colnames(sequence_df) <- "variants"
  ## At this point we have one column which contains variant specs which start
  ## at the beginning of the region of interest.  So let us choose
  ## primers which end 1 nt before that region.
  sequence_df[["chr"]] <- gsub(x = rownames(sequence_df), pattern = "^(.*)_start_.*$",
                               replacement = "\\1")
  sequence_df[["bin_start"]] <- as.numeric(gsub(x = rownames(sequence_df), pattern = "^(.*)_start_(.*)$",
                                                replacement = "\\2"))
  ## I think it would be nice to have spec strings which count from right->left to make
  ## reading off a 3' primer easier.
  sequence_df[["fwd_references"]] <- gsub(x = sequence_df[["variants"]],
                                          pattern = "([[:alpha:]])(\\d+)([[:alpha:]])",
                                          replacement = "\\1")
  sequence_df[["rev_references"]] <- chartr("ATGC", "TACG", sequence_df[["fwd_references"]])

  sequence_df[["fwd_positions"]] <- gsub(x = sequence_df[["variants"]],
                                         pattern = "([[:alpha:]])(\\d+)([[:alpha:]])",
                                         replacement = "\\2")
  ## Getting the relative 3' position cannot be done with a simple regex/tr operation.
  ## but instead is binwidth - position
  ## which therefore requires a little apply() shenanigans
  sequence_df[["rev_positions"]] <- mapply(FUN = recount_positions,
                                           sequence_df[["fwd_positions"]])

  sequence_df[["fwd_alternates"]] <- gsub(x = sequence_df[["variants"]],
                                          pattern = "([[:alpha:]])(\\d+)([[:alpha:]])",
                                          replacement = "\\3")
  sequence_df[["rev_alternates"]] <- chartr("ATGC", "TACG", sequence_df[["fwd_alternates"]])
  ## Now lets rewrite the variants as the combinations of rev_reference/position/alternate

  sequence_df[["rev_variants"]] <- mapply(FUN = rev_variants,
                                          sequence_df[["rev_references"]],
                                          sequence_df[["rev_positions"]],
                                          sequence_df[["rev_alternates"]])

  ## I am making explicit columns for these so I can look manually
  ## and make sure I am not trying to get sequences which run off the chromosome.
  sequence_df[["fwd_superprimer_start"]] <- sequence_df[["bin_start"]] -
    max_primer_length - 1
  sequence_df[["rev_superprimer_start"]] <- sequence_df[["bin_start"]] +
    bin_width + max_primer_length + 1
  sequence_df[["fwd_superprimer_end"]] <- sequence_df[["bin_start"]] - 1
  sequence_df[["rev_superprimer_end"]] <- sequence_df[["bin_start"]] + bin_width
  ## The super_primers are the entire region of length == max_primer_length, I will
  ## iteratively remove bases from the 5' until the target is reached.
  sequence_df[["fwd_superprimer"]] <- ""
  sequence_df[["rev_superprimer"]] <- ""
  ## The fwd_primer rev_primer columns will hold the final primers.
  sequence_df[["fwd_primer"]] <- ""
  sequence_df[["rev_primer"]] <- ""
  sequence_df[["fwd_note"]] <- ""
  sequence_df[["rev_note"]] <- ""
  sequence_df[["fwd_gc_prop"]] <- 0
  sequence_df[["rev_gc_prop"]] <- 0

  ## For me, this is a little too confusing to do in an apply.
  ## I want to fill in the fwd/rev_superprimer columns with the longer
  ## sequence/revcomp regions of interest.
  ## Then I want to chip away at the beginning/end of them to get the actual fwd/rev primer.
  for (i in seq_len(nrow(sequence_df))) {
    chr_name <- sequence_df[i, "chr"]
    super_start <- sequence_df[i, "fwd_superprimer_start"]
    super_end <- sequence_df[i, "fwd_superprimer_end"]
    max_end <- length(genome[[chr_name]])
    if (super_start < 1) {
      super_start <- 1
      sequence_df[i, "fwd_note"] <- "ran over chromosome start."
    } else {
      fwd_superprimer <- Biostrings::subseq(
        genome[[chr_name]],
        start = sequence_df[i, "fwd_superprimer_start"],
        end = sequence_df[i, "fwd_superprimer_end"])
      if (grepl(x = fwd_superprimer, pattern = "N")) {
        sequence_df[i, "fwd_note"] <- "contains N"
      } else {
        fwd_superprimer <- as.character(fwd_superprimer)
        sequence_df[i, "fwd_superprimer"] <- fwd_superprimer
        fwd_primer <- try(find_subseq_target_temp(fwd_superprimer,
                                                  target = target_temp,
                                                  direction = "forward"))
        if ("try-error" %in% class(fwd_primer)) {
          sequence_df[i, "fwd_note"] <- "bad sequence for priming"
        } else {
          sequence_df[i, "fwd_primer"] <- fwd_primer
          ## Add the GC content
          fwd_gc <- Biostrings::letterFrequency(
            Biostrings::DNAString(fwd_primer), "GC", as.prob = TRUE)
          sequence_df[i, "fwd_gc_prop"] <- as.numeric(fwd_gc)
          if (fwd_gc < min_gc_prop) {
            sequence_df[i, "fwd_note"] <- "bad GC content"
          }
          over_run <- grepl(x = fwd_primer, pattern = run_pattern)
          if (isTRUE(over_run)) {
            sequence_df[i, "fwd_note"] <- paste0(sequence_df[i, "fwd_note"], " run of too many")
          }
        }
      } ## End checking for pathologically bad primers (Ns, unable to find a decent Tm)
    } ## End checking if we overran the beginning of the contig.

    ## Note, start and end here are a little confusing because I was thinking
    ## start and end in terms of 5'/3', not chromosome coordinates.
    ## The final endpoint is
    rev_end <- sequence_df[i, "rev_superprimer_start"]
    rev_start <- sequence_df[i, "rev_superprimer_end"]
    if (rev_end > max_end) {
      super_end <- max_end
      sequence_df[i, "rev_note"] <- "ran over chromosome end."
    } else {
      rev_superprimer <- Biostrings::subseq(
        genome[[chr_name]],
        start = rev_start,
        end = rev_end)
      if (grepl(x = rev_superprimer, pattern = "N")) {
        sequence_df[i, "rev_note"] <- "contains N"
      } else {
        rev_superprimer <- toupper(as.character(spgs::reverseComplement(rev_superprimer)))
        sequence_df[i, "rev_superprimer"] <- rev_superprimer
        rev_primer <- find_subseq_target_temp(rev_superprimer,
                                              target = target_temp,
                                              direction = "reverse")
        if ("try-error" %in% class(rev_primer)) {
          sequence_df[i, "rev_note"] <- "bad sequence for priming"
        } else {
          sequence_df[i, "rev_primer"] <- rev_primer
          ## Add the GC content
          rev_gc <- Biostrings::letterFrequency(
            Biostrings::DNAString(rev_primer), "GC", as.prob = TRUE)
          sequence_df[i, "rev_gc_prop"] <- as.numeric(rev_gc)
          if (rev_gc < min_gc_prop) {
            sequence_df[i, "rev_note"] <- "bad GC content"
          }
          over_run <- grepl(x = rev_primer, pattern = run_pattern)
          if (isTRUE(over_run)) {
            sequence_df[i, "fwd_note"] <- paste0(sequence_df[i, "fwd_note"], " run of too many")
          }
        }
      } ## End checking for pathologically bad primers (Ns, unable to find a decent Tm)
    } ## End checking if we over-ran the end of the contig.
  } ## End iterating over sequence_df

  ## Drop anything which is definitely useless
  ## Probably should add some logic here to choose filters.
  keepers <- sequence_df[["fwd_note"]] == ""
  filtered_sequence_df <- sequence_df[keepers, ]
  keepers <- filtered_sequence_df[["rev_note"]] == ""
  filtered_sequence_df <- filtered_sequence_df[keepers, ]

  return(filtered_sequence_df)
}

#' Find a subsequence with a target PCR temperature.
#'
#' Given a relatively large sequence, this function will iteratively
#' remove a single nucleotide and recalulate the TM until the TM falls
#' to the target temperature.
#'
#' @param sequence Starting sequence.
#' @param target Desire TM of the final sequence.
#' @param direction What strand is expected for annealing this primer?
#' @param verbose Be chatty?
find_subseq_target_temp <- function(sequence, target = 53,
                                    direction = "forward", verbose = FALSE) {
  cheapo <- cheap_tm(sequence)
  if (isTRUE(verbose)) {
    mesg("Starting with ", sequence, " and ", cheapo)
  }
  if (nchar(sequence) < 1) {
    stop("Never hit the target tm, something is wrong.")
  }
  if (cheapo[["tm"]] > target) {
    if (direction == "forward") {
      ## Then drop the first character
      subseq <- stringr::str_sub(sequence, 2, nchar(sequence))
    } else {
      ## Then drop the last character
      subseq <- stringr::str_sub(sequence, 1, nchar(sequence) - 1)
    }
    result <- find_subseq_target_temp(subseq, target = target, direction = direction)
  } else {
    mesg("Final tm is: ", cheapo[["tm"]], " from sequence: ", sequence, ".")
    return(sequence)
  }
}

#' Perform a series of tests of a putative primer.
#'
#' This function should probably replace the morass of code found in
#' snp_density_primers().  It is current used by snp_cds_primers().
#'
#' @param entry Single row of the table of potential primers.
#' @param genome bsgenome used to search against.
#' @param variant_gr GRanges of variants to xref against.
#' @param target_temp Desired Tm.
#' @param direction Either fwd or rev.
#' @param run_pattern Regex to look for bad runs of a single nt.
#' @param min_gc_prop Minimum proportion of GC content.
#' @param seq_object Used to hunt for multi-hit primers.
primer_qc <- function(entry, genome,
                      variant_gr, target_temp = 60,
                      direction = "fwd", run_pattern = "AAAA",
                      min_gc_prop = 0.3, seq_object = NULL) {
  silly <- "forward"

  ## Define the places from which to acquire the various parameters we will need.
  ## They will all be fwd_something or rev_something.
  note_name <- paste0(direction, "_note")
  super_name <- paste0(direction, "_superprimer")
  primer_name <- paste0(direction, "_primer")
  gc_name <- paste0(direction, "_gc_prop")
  primer_var_name <- paste0(direction, "_var_in_primer")
  occur_name <- paste0(direction, "_occurrences")
  ## Except the names of the start and end coordinates, they are different because of strandedness.
  start_name <- "bin_start"
  end_name <- "fwd_end"
  if (direction == "rev") {
    silly <- "reverse"
    start_name <- "rev_start"
    end_name <- "bin_end"
  }
  ## Grab the contig name
  chr_name <- entry[["bin_seqname"]]
  ## Use the contig name to get a superprimer
  super <- as.character(Biostrings::subseq(genome[[chr_name]],
                                           start = entry[[start_name]],
                                           end = entry[[end_name]]))
  if (direction == "rev") {
    ## If we are getting a reverse primer, revcomp it.
    super <- toupper(as.character(spgs::reverseComplement(super)))
  }

  ## Do not bother writing anything for sequences containing N
  if (grepl(x = super, pattern = "N")) {
    entry[[note_name]] <- "contains N"
  } else {
    entry[[super_name]] <- super
    ## find_subseq_target_temp iterates over the superprimer
    ## and removes 1 nt until it hits a Tm that is near our target.
    primer <- try(find_subseq_target_temp(super,
                                          target = target_temp,
                                          direction = silly))
    ## Given this new primer sequence, let us see if it is crappy
    if ("try-error" %in% class(primer)) {
      ## AFAIK one only gets an error if the ranges are poorly constructed.
      entry[[note_name]] <- "bad sequence for priming"
    } else if (isTRUE(grepl(x = primer, pattern = run_pattern))) {
      ## If there are runs of any nucleotide, make a note.
      entry[[note_name]] <- paste0(entry[[note_name]], " run of too many.")
    } else {
      ## Check the primer region to see if it overlaps with any variants.
      ## Just make a quick granges for the primer and
      ## see if it overlaps with any variants
      primer_length <- nchar(primer)
      primer_gr_string <- paste0(chr_name, ":", entry[[start_name]], "-", entry[[start_name]] + primer_length)
      primer_gr <- as(primer_gr_string, "GRanges")
      primer_overlap_gr <- GenomicRanges::findOverlaps(primer_gr, variant_gr)
      variant_primer_overlaps <- length(primer_overlap_gr)
      ## I figure that if this number is not too big, then the primer
      ## should be ok, as long as it doesn't turn out to be the last
      ## base of the primer; so I will not (yet) do an explicit check.
      entry[[primer_var_name]] <- variant_primer_overlaps
      entry[[primer_name]] <- primer
      ## Add the GC content
      gc <- Biostrings::letterFrequency(Biostrings::DNAString(primer),
                                        "GC", as.prob = TRUE)
      entry[[gc_name]] <- as.numeric(gc)
      ## Also add a note if the GC content is too low, I do not bother looking for a GC clamp.
      if (gc < min_gc_prop) {
        entry[[note_name]] <- paste0(entry[[note_name]], " bad GC content.")
      }
    }

    ## If this function receives a Biostrings copy of the genome, look to see if the primer matches
    ## in multiple locations and record that number.
    if (!is.null(seq_object)) {
      ## See how many times this primer is found in the genome
      ## vcount doesn't work with mismatch > 0

      dict <- NULL
      if (direction == "rev") {
        tmp_primer <- toupper(as.character(spgs::reverseComplement(primer)))
        dict <- Biostrings::PDict(tmp_primer, max.mismatch = 0)
      } else {
        dict <- Biostrings::PDict(primer, max.mismatch = 0)
      }
      num_hits <- sum(unlist(Biostrings::vcountPDict(dict, seq_obj)))
      entry[[occur_name]] <- num_hits
      if (num_hits > 1) {
        entry[[note_name]] <- paste0(entry[[note_name]], " this primer hits multiple places.")
      }
    }
  }
  return(entry)
}

#' Search a set of variants for ones which are relatively sequential.
#'
#' One potential way to screen strains is to use PCR primers which
#' should(not) anneal due to variants with respect to the genome.
#' This function seeks to find variants which are clustered
#' sufficiently close to each other that this is possible.
#'
#' @param snp_sets Result from get_snp_sets() containing the variants
#'  with respect to known conditions.
#' @param conditions Set of conditions to search against.
#' @param minimum Minimum number of variants required for a candiate.
#' @param maximum_separation How far apart from each other are these
#'  >={minimum} variants allowed to be?
#' @param one_away_file Location to write variants that are no more than 1 base apart.
#' @param two_away_file Location for those which are no more than 2 apart.
#' @param doubles_file Write out variants which are 2 in a row.
#' @param singles_file Write out the individual variants here.
sequential_variants <- function(snp_sets, conditions = NULL,
                                minimum = 3, maximum_separation = 3,
                                one_away_file = "one_away.csv",
                                two_away_file = "two_away.csv",
                                doubles_file = "doubles.csv",
                                singles_file = "singles.csv") {
  if (is.null(conditions)) {
    conditions <- 1
  }
  intersection_sets <- snp_sets[["intersections"]]
  intersection_names <- snp_sets[["set_names"]]
  chosen_intersection <- 1
  if (is.numeric(conditions)) {
    chosen_intersection <- conditions
  } else {
    intersection_idx <- intersection_names == conditions
    chosen_intersection <- names(intersection_names)[intersection_idx]
    mesg("Using intersection: ", chosen_intersection,
         ", which corresponds to: ",
         intersection_names[intersection_idx], ".")
  }

  possible_positions <- intersection_sets[[chosen_intersection]]
  position_table <- data.frame(row.names = possible_positions)
  pat <- "^chr_(.+)_pos_(.+)_ref_.*$"
  position_table[["chr"]] <- gsub(pattern = pat, replacement = "\\1",
                                  x = rownames(position_table))
  position_table[["pos"]] <- as.numeric(gsub(pattern = pat,
                                             replacement = "\\2",
                                             x = rownames(position_table)))
  position_idx <- order(position_table[, "chr"], position_table[, "pos"])
  position_table <- position_table[position_idx, ]
  position_table[["dist"]] <- 0

  last_chr <- ""
  this_chr <- position_table[1, "chr"]
  for (r in seq_len(nrow(position_table))) {
    this_chr <- position_table[r, "chr"]
    if (r == 1) {
      position_table[r, "dist"] <- position_table[r, "pos"]
      last_chr <- this_chr
      next
    }
    if (this_chr == last_chr) {
      previous_idx <- r - 1
      position_table[r, "dist"] <- position_table[r, "pos"] -
        position_table[previous_idx, "pos"]
    } else {
      position_table[r, "dist"] <- position_table[r, "pos"]
    }
    last_chr <- this_chr
  }
  message("Here is a table of the variants closest together:")
  print(head(table(position_table[["dist"]])))

  ## Working interactively here.
  doubles <- position_table[["dist"]] == 1
  doubles <- position_table[doubles, ]
  if (nrow(doubles) > 0) {
    write.csv(doubles, doubles_file)
  } else {
    message("There are no doubles to write.")
  }

  one_away <- position_table[["dist"]] == 2
  one_away <- position_table[one_away, ]
  if (nrow(one_away) > 0) {
    write.csv(one_away, one_away_file)
  } else {
    message("There are no variants with 1 nt separating them.")
  }

  two_away <- position_table[["dist"]] == 3
  two_away <- position_table[two_away, ]
  if (nrow(two_away) > 0) {
    write.csv(two_away, two_away_file)
  } else {
    message("There are no variants with 2 nt separating them.")
  }

  combined <- data.frame()
  if (nrow(doubles) > 0 && nrow(one_away) > 0) {
    combined <- rbind(doubles, one_away)
  }
  if (nrow(two_away) > 0) {
    combined <- rbind(combined, two_away)
  }

  if (nrow(combined) > 0) {
    position_idx <- order(combined[, "chr"], combined[, "pos"])
    combined <- combined[position_idx, ]
    this_chr <- ""
    last_chr <- NULL
    for (r in seq_len(nrow(combined))) {
      this_chr <- combined[r, "chr"]
      if (r == 1) {
        combined[r, "dist_pair"] <- combined[r, "pos"]
        last_chr <- this_chr
        next
      }
      if (this_chr == last_chr) {
        combined[r, "dist_pair"] <- combined[r, "pos"] - combined[r - 1, "pos"]
      } else {
        combined[r, "dist_pair"] <- combined[r, "pos"]
      }
      last_chr <- this_chr
    }

    dist_pair_maximum <- 1000
    dist_pair_minimum <- 200
    dist_pair_idx <- combined[["dist_pair"]] <= dist_pair_maximum &
      combined[["dist_pair"]] >= dist_pair_minimum
    remaining <- combined[dist_pair_idx, ]
    no_weak_idx <- grepl(pattern = "ref_(G|C)", x = rownames(remaining))
    remaining <- remaining[no_weak_idx, ]

    print(head(table(position_table[["dist"]])))
    sequentials <- position_table[["dist"]] <= maximum_separation
    message("There are ", sum(sequentials), " candidate regions.")

    ## The following can tell me how many runs of each length occurred, that is not quite what I want.
    ## Now use run length encoding to find the set of sequential sequentials!
    rle_result <- rle(sequentials)
    rle_values <- rle_result[["values"]]
    ## The following line is equivalent to just leaving values alone:
    ## true_values <- rle_result[["values"]] == TRUE
    rle_lengths <- rle_result[["lengths"]]
    true_sequentials <- rle_lengths[rle_values]
    rle_idx <- cumsum(rle_lengths)[which(rle_values)]

    position_table[["last_sequential"]] <- 0
    count <- 0
    for (r in rle_idx) {
      count <- count + 1
      position_table[r, "last_sequential"] <- true_sequentials[count]
    }
    message("The maximum sequential set is: ", max(position_table[["last_sequential"]]), ".")

    wanted_idx <- position_table[["last_sequential"]] >= minimum
    wanted <- position_table[wanted_idx, c("chr", "pos")]
  } else {
    wanted <- NULL
  }
  class(wanted) <- "sequential_variants"
  return(wanted)
}

#' Look for variants associated with CDS regions instead of high-density.
#'
#' The function snp_density_primers looks for regions with many
#' variants.  This flips the script and looks first to the set of CDS
#' regions.  It also makes heavy use of GRanges and so should prove
#' useful as a reference when looking for range examples.
#'
#' @param cds_gr GRanges of CDS features.  It does not have to be CDS,
#'  but probably should not be genes.
#' @param variant_gr GRanges of observed variants.  I get this by
#'  coercing my peculiar variant rownames into a GR.
#' @param bsgenome Genome containing all the contigs mentioned above.
#' @param amplicon_size Desired PCR amplicon for sequencing.
#' @param min_overlap Desired overlap between every genome bin and
#'  CDS.  Note I didn't say amplicon here because of the way I am
#'  making the primers.
#' @param minvar_perbin Discard bins with less than this number of
#'  variants inside them.
#' @param super_len I start out with a 'superprimer' which is assumed
#'  to be longer than needed for the target Tm.  It is this long.
#' @param target_temp Attempt to create primers with this Tm.
#' @param min_gc_prop Warn or discard primers with less than this GC content.
#' @param max_nmer_run Warn or discard primers with runs of a single
#'  base this long.
#' @param count_occurrences Count up how many times the primer is
#'  found in the genome, hopefully this is always 1.  Annoyingly, the
#'  vcountDict function does not allow mismatches.
#' @param occurrence_mismatch I cannot use this, but I want to.
snp_cds_primers <- function(cds_gr, variant_gr, bsgenome, amplicon_size = 600, min_overlap = 200,
                            minvar_perbin = 10, super_len = 40, target_temp = 60,
                            min_gc_prop = 0.3, max_nmer_run = 4, count_occurrences = TRUE,
                            occurrence_mismatch = 0) {

  ## Make a regex to search for excessive runs of a single nt.
  run_pattern <- paste0("A{", max_nmer_run,
                        ",}|T{", max_nmer_run,
                        ",}|G{", max_nmer_run,
                        ",}|C{", max_nmer_run, ",}")
  ## The amplicons will not be overlapping by the distance
  ## comprised of the superprimer - actual primer length * 2
  tile_width <- amplicon_size + (super_len * 2)
  ## because we will tile it to the amplicon + (super * 2)
  bin_gr <- GenomicRanges::tileGenome(seqlengths(cds_gr),
                                      tilewidth = tile_width,
                                      cut.last.tile.in.chrom = TRUE)
  ## Get the overlaps between the bins and CDSes
  ## These expect to have a resonable minimum overlap
  bin_cds_idx <- GenomicRanges::findOverlaps(bin_gr, cds_gr,
                                             ignore.strand = TRUE,
                                             minoverlap = min_overlap)
  bin_cds_gr <- IRanges::subsetByOverlaps(bin_gr, cds_gr,
                                          ignore.strand = TRUE,
                                          minoverlap = min_overlap)
  ## and the overlaps between the bins and variants.
  bin_var_idx <- GenomicRanges::findOverlaps(bin_gr, variant_gr, ignore.strand = TRUE)
  bin_var_gr <- IRanges::subsetByOverlaps(bin_gr, variant_gr, ignore.strand = TRUE)

  ## recast the granges to dfs so I can play with them
  ## while at it, add an index column for later.
  bin_cds_df <- as.data.frame(bin_cds_gr)
  bin_cds_df[["idx"]] <- rownames(bin_cds_df)
  bin_df <- as.data.frame(bin_gr)
  bin_df[["idx"]] <- rownames(bin_df)
  cds_df <- as.data.frame(cds_gr)
  cds_df[["idx"]] <- rownames(cds_df)
  variant_df <- as.data.frame(variant_gr)
  idx_vector <- seq_len(nrow(variant_df))
  variant_df[["idx"]] <- idx_vector

  ## Now get the overlap of overlaps, I will use this to
  ## count up the variants / bin / CDS
  both_idx <- GenomicRanges::findOverlaps(bin_cds_gr, variant_gr)
  both_gr <- IRanges::subsetByOverlaps(bin_cds_gr, variant_gr)

  ## Give column names to my various dataframes that should
  ## make it easier for me to track what is what.
  both_df <- as.data.frame(both_idx)
  colnames(both_df) <- c("bincds", "var")
  bin_cds_df <- as.data.frame(bin_cds_idx)
  colnames(bin_cds_df) <- c("bin", "cds")
  bin_var_df <- as.data.frame(bin_var_idx)
  colnames(bin_var_df) <- c("bin", "var")

  ## Group the bins+cds and count up the variants.
  variants_per_bin <- as.data.frame(both_df) %>%
    group_by(bincds) %>%
    summarise(n = n()) %>%
    arrange(desc(n))

  ## Make a big df of all these little dfs merged together.
  ## use the new column names to hopefully make it easier.
  all_merged <- merge(bin_cds_df, bin_var_df, by = "bin")
  all_merged <- merge(all_merged, both_df,
                      by.x = "var", by.y = "var",
                      all.x = TRUE)
  ## Order them by the most fun bins first, where fun is defined
  ## by bins with the most variants within them.
  all_merged <- merge(all_merged, variants_per_bin,
                      by = "bincds", all.x = TRUE) %>%
    arrange(desc(n))
  ## and ignore the bins with too few variants
  starting_rows <- nrow(all_merged)
  keepers <- all_merged[["n"]] >= minvar_perbin
  all_merged <- all_merged[keepers, ]
  ending_rows <- nrow(all_merged)
  message("Filtering to at least ", minvar_perbin, " variants per bin reduced the df from: ",
          starting_rows, " to ", ending_rows, ".")

  ## We have a whole bunch of repeated/redundant columns,
  ## but I want to keep them so that I can sanity check myself.
  colnames(bin_df) <- c("bin_seqname", "bin_start", "bin_end", "bin_width", "bin_strand", "idx")
  ## Merge in the bin data
  all_merged <- merge(all_merged, bin_df, by.x = "bin", by.y = "idx")
  ## and clean up the CDS column names before merging it in.
  colnames(cds_df)[1] <- "cds_seqname"
  colnames(cds_df)[2] <- "cds_start"
  colnames(cds_df)[3] <- "cds_end"
  colnames(cds_df)[4] <- "cds_width"
  colnames(cds_df)[5] <- "cds_strand"
  ## then merge in the CDS data
  all_merged <- merge(all_merged, cds_df, by.x = "cds", by.y = "idx")
  ## again, clean up the variant column nanes before merging.
  colnames(variant_df)[1] <- "var_seqname"
  colnames(variant_df)[2] <- "var_start"
  colnames(variant_df)[3] <- "var_end"
  variant_df[["width"]] <- NULL
  variant_df[["strand"]] <- NULL
  ## and merge the variants, note strand/width are useless here.
  all_merged <- merge(all_merged, variant_df, by.x = "var", by.y = "idx")


  ## Figure out the relative positions of the bin/variant start positions.
  ## I am not going to bother repeating this for the reverse primers with
  ## the assumption that whoever actually uses these will
  ## only read from the front.
  all_merged[["fwd_end"]] <- all_merged[["bin_start"]] + super_len
  all_merged[["rev_start"]] <- all_merged[["bin_end"]] - super_len
  all_merged[["relative"]] <- all_merged[["var_start"]] - all_merged[["fwd_end"]]

  ## Merge every row which shares a bin and format it so that the
  ## variant data is x,y,z,a,b,c for the references/alternates/positions
  final_merged <- all_merged %>%
    group_by(bincds) %>%
    distinct(relative, .keep_all = TRUE) %>%
    mutate(ref_series = paste0(reference, collapse = ","),
           alt_series = paste0(alternate, collapse = ","),
           pos_series = paste0(relative, collapse = ",")) %>%
    distinct(bincds, .keep_all = TRUE)
  final_rows <- nrow(test_merged)
  message("Collapsing variants/bin reduced the rows from: ", ending_rows, " to: ", final_rows, ".")

  ## I need to double check these number to ensure that the relative positions are correct
  ## vis a vis the beginning of the amplicon.

  ## Now I have everything in place I need to make a PCR primer?
  ## Let us try, load up the genome, and if requested a biostrings version
  ## of same...
  genome <- NULL
  if (!is.null(bsgenome)) {
    ## FIXME: Shouldn't use library()
    library(bsgenome, character.only = TRUE)
    genome <- get0(bsgenome)
  }
  seq_obj <- NULL
  if (isTRUE(count_occurrences)) {
    seq_obj <- Biostrings::getSeq(genome)
  }

  ## Fill in some empty columns which the primer_qc() function will
  ## be responsible for filling in.
  final_merged[["fwd_superprimer"]] <- ""
  final_merged[["rev_superprimer"]] <- ""
  final_merged[["fwd_primer"]] <- ""
  final_merged[["rev_primer"]] <- ""
  final_merged[["fwd_note"]] <- ""
  final_merged[["rev_note"]] <- ""
  final_merged[["fwd_gc_prop"]] <- 0
  final_merged[["rev_gc_prop"]] <- 0
  final_merged[["bin_seqname"]] <- as.character(final_merged[["bin_seqname"]])
  final_merged[["fwd_occurrences"]] <- 0
  final_merged[["rev_occurrences"]] <- 0
  final_merged[["fwd_var_in_primer"]] <- 0
  final_merged[["rev_var_in_primer"]] <- 0

  ## Iterate over every element in the merged data and go primer hunting.
  for (i in seq_len(nrow(final_merged))) {
    message("Working on ", i, " of ", nrow(final_merged))
    entry <- final_merged[i, ]
    ## note, I am doing this in two passes so I can look at the pieces
    ## before committing them to the final dataframe.
    new_entry <- primer_qc(entry, genome, variant_gr,
                           target_temp = target_temp, direction = "fwd",
                           run_pattern = run_pattern, min_gc_prop = min_gc_prop,
                           seq_object = seq_obj)
    final_entry <- primer_qc(new_entry, genome, variant_gr,
                             target_temp = target_temp, direction = "rev",
                             run_pattern = run_pattern, min_gc_prop = min_gc_prop,
                             seq_object = seq_obj)
    final_merged[i, ] <- final_entry
  }
  retlist <- list(
    "primer_table" = final_merged)
  class(retlist) <- "cds_variant_primers"
  return(retlist)
}

#' Create a density function given a variant output and some metadata
#'
#' It is hoped that this will point out regions of a genome which
#' might prove useful when designing PCR primers for a specific
#' condition in a dataset of variants.
#'
#' @param snp_count Result from count_expt_snps()
#' @param pdata_column Metadata column containing the condition of
#'  interest.
#' @param condition Chosen condition to search for variants.
#' @param cutoff Minimum number of variants in a region.
#' @param bin_width Bin size/region of genome to consider.
#' @param divide Normalize by bin width?
#' @param topn Keep only this number of candidates.
#' @param target_temp Try to get primers with this Tm.
#' @param max_primer_length Keep primers at or less than this length.
#' @param bsgenome Genome package containing the sequence of interest.
#' @param gff GFF to define regions of interest.
#' @param feature_type GFF feature type to search against.
#' @param feature_start GFF column with the starts (needed?)
#' @param feature_end GFF column with the ends (needed?)
#' @param feature_strand GFF column with strand information (needed?)
#' @param feature_chr GFF column with chromosome information.
#' @param feature_type_column GFF column with type information.
#' @param feature_id GFF tag with the ID information.
#' @param feature_name GFF tag with the names.
#' @param truncate Truncate the results to just the columns I think
#'  are useful.
#' @param min_gc_prop Minimum GC content for a suitable primer.
#' @export
snp_density_primers <- function(snp_count, pdata_column = "condition",
                                condition = NULL, cutoff = 20, bin_width = 600,
                                divide = FALSE, topn = 400,
                                target_temp = 53, max_primer_length = 50,
                                bsgenome = "BSGenome.Leishmania.panamensis.MHOMCOL81L13.v52",
                                gff = "reference/lpanamensis_col_v46.gff",
                                feature_type = "protein_coding_gene", feature_start = "start",
                                feature_end = "end", feature_strand = "strand",
                                feature_chr = "seqnames", feature_type_column = "type",
                                feature_id = "ID", feature_name = "description",
                                truncate = TRUE, xref_genes = TRUE,
                                verbose = FALSE, min_contig_length = NULL,
                                min_gc_prop = 0.25, max_nmer_run = 5) {

  ## Start out by loading the bsgenome data
  genome <- NULL
  if (!is.null(bsgenome)) {
    ## FIXME: Shouldn't use library()
    library(bsgenome, character.only = TRUE)
    genome <- get0(bsgenome)
  }

  samples_by_condition <- pData(snp_count)[[pdata_column]]
  ## Keep only those samples of the condition of interest for now,
  ## maybe make it a loop to iterate over conditions later.
  snp_table <- exprs(snp_count)
  kept_samples <- 0
  if (!is.null(condition)) {
    interest_idx <- samples_by_condition == condition
    kept_samples <- sum(interest_idx)
    if (kept_samples < 1) {
      stop("No samples were deemed interesting when looking only at ",
           condition, " in column ", pdata_column, ".")
    } else if (kept_samples == 1) {
      warning("There remains only 1 sample when subsetting on condition.")
      the_colname <- colnames(snp_table)[interest_idx]
      the_rownames <- rownames(snp_table)
      snp_vector <- snp_table[, interest_idx]
      snp_table <- as.matrix(snp_vector)
      colnames(snp_table) <- the_colname
      rownames(snp_table) <- the_rownames
    } else {
      snp_table <- snp_table[, interest_idx]
    }
  }
  na_idx <- is.na(snp_table)
  if (sum(na_idx) > 0) {
    warning("There are: ", sum(na_idx), " NA values in the data, converting them to 0.")
    snp_table[na_idx] <- 0
  }
  start_rows <- nrow(snp_table)
  if (kept_samples > 1) {
    interest_rows <- rowMeans(snp_table) >= cutoff
    snp_table <- snp_table[interest_rows, ]
  } else {
    interest_rows <- snp_table[, 1] >= cutoff
    the_colname <- colnames(snp_table)
    snp_vector <- snp_table[interest_rows, ]
    the_rownames <- names(snp_vector)
    snp_table <- as.matrix(snp_vector)
    rownames(snp_table) <- the_rownames
    colnames(snp_table) <- the_colname
  }

  mesg("Started with ", start_rows, ", candidate snps, after cutoff there are ",       sum(interest_rows), " remaining.")
  position_table <- data.frame(row.names = rownames(snp_table))
  position_table[["chromosome"]] <- gsub(x = rownames(position_table),
                                         pattern = "chr_(.*)_pos_.*",
                                         replacement = "\\1")
  position_table[["position"]] <- gsub(x = rownames(position_table),
                                       pattern = ".*_pos_(\\d+)_.*",
                                       replacement = "\\1")
  ## In the case of lpanamensis at least, chromosomes have '-' in the name, this is not allowed.
  position_table[["chromosome"]] <- gsub(pattern = "-", replacement = "_",
                                         x = position_table[["chromosome"]])
  position_table[["ref"]] <- gsub(x = rownames(position_table),
                                  pattern = ".*_pos_\\d+_ref_(\\w).*",
                                  replacement = "\\1")
  position_table[["alt"]] <- gsub(x = rownames(position_table),
                                  pattern = ".*_pos_\\d+_ref_.+_alt_(\\w)$",
                                  replacement = "\\1")
  ## Lets make a dataframe of chromosomes, their lengths, and number of variants/chromosome
  chromosomes <- levels(as.factor(position_table[["chromosome"]]))

  chromosome_df <- data.frame(row.names = chromosomes)
  chromosome_df[["length"]] <- 0
  chromosome_df[["variants"]] <- 0
  density_lst <- list()
  variant_lst <- list()
  ## This was written in an attempt to make it work if one does or does not
  ## have a BSgenome from which to get the actual chromosome lengths.

  show_progress <- TRUE
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style = 3)
  }
  end <- length(chromosomes)
  for (ch in seq_along(chromosomes)) {
    if (isTRUE(show_progress)) {
      pct_done <- ch / end
      utils::setTxtProgressBar(bar, pct_done)
    }
    chr <- chromosomes[ch]
    chr_idx <- position_table[["chromosome"]] == chr
    chr_data <- position_table[chr_idx, ]
    chr_len = NULL
    if (is.null(genome)) {
      chromosome_df[chr, "length"] <- max(chr_data[["position"]])
      chr_len <- max(chr_data[["position"]])
    } else {
      this_length <- length(genome[[chr]])
      chromosome_df[chr, "length"] <- this_length
      chr_len <- this_length
    }
    if (!is.null(min_contig_length)) {
      if (chr_len < min_contig_length) {
        next
      }
    }
    chromosome_df[chr, "variants"] <- nrow(chr_data)

    ## Stop scientific notation for this operation
    current_options <- options(scipen = 999)
    density_vector <- rep(0, ceiling(this_length / bin_width))
    variant_vector <- rep('', ceiling(this_length / bin_width))
    density_name <- 0
    for (i in seq_along(density_vector)) {
      range_min <- density_name
      range_max <- (density_name + bin_width) - 1
      ## We want to subtract the maximum primer length from the position
      ## and have each extension reaction start near/at the bin.
      chr_data[["position"]] <- as.numeric(chr_data[["position"]])
      density_idx <- chr_data[["position"]] >= range_min &
        chr_data[["position"]] <= range_max
      ## Note that if I want to completely correct in how I do this
      ## it should not be blindly bin_width, but should be the smaller of
      ## bin_width or the remainder of the chromosome.  Leaving it at bin_width
      ## will under-represent the end of each chromosome/contig, but I don't think
      ## I mind that at all.
      if (isTRUE(divide)) {
        density_vector[i] <- sum(density_idx) / bin_width
      } else {
        density_vector[i] <- sum(density_idx)
      }
      names(density_vector)[i] <- as.character(density_name)
      ## Add the observed variants to the variant_vector
      ## variant_set <- chr_data[density_idx, "mutation"]
      variant_relative_positions <- chr_data[density_idx, "position"] - range_min
      variant_set <- paste0(chr_data[density_idx, "ref"],
                            variant_relative_positions,
                            chr_data[density_idx, "alt"])
      variant_vector[i] <- toString(variant_set)
      names(variant_vector)[i] <- as.character(density_name)

      density_name <- density_name + bin_width
    } ## End iterating over every element of the density vector
    density_lst[[chr]] <- density_vector
    variant_lst[[chr]] <- variant_vector
  } ## End iterating over every chromosome
  if (isTRUE(show_progress)) {
    close(bar)
  }
  new_options <- options(current_options)

  long_density_vector <- vector()
  long_variant_vector <- vector()
  for (ch in seq_along(density_lst)) {
    chr <- names(density_lst)[ch]
    density_vec <- density_lst[[chr]]
    variant_vec <- variant_lst[[chr]]
    vecnames <- paste0(chr, "_start_", names(density_vec))
    names(density_vec) <- vecnames
    names(variant_vec) <- vecnames
    long_density_vector <- c(long_density_vector, density_vec)
    long_variant_vector <- c(long_variant_vector, variant_vec)
  }  ## End iterating over the density list
  most_idx <- order(long_density_vector, decreasing = TRUE)
  long_density_vector <- long_density_vector[most_idx]
  long_variant_vector <- long_variant_vector[most_idx]

  mesg("Extracting primer regions.")
  sequence_df <- choose_sequence_regions(long_variant_vector,
                                         max_primer_length = max_primer_length,
                                         topn = topn, bin_width = bin_width,
                                         genome = genome, target_temp = target_temp,
                                         min_gc_prop = min_gc_prop,
                                         max_nmer_run = max_nmer_run)
  if (isTRUE(xref_genes)) {
    mesg("Searching for overlapping/closest genes.")
    sequence_df <- xref_regions(sequence_df, gff, bin_width = bin_width,
                                feature_type = feature_type, feature_start = feature_start,
                                feature_end = feature_end, feature_strand = feature_strand,
                                feature_chr = feature_chr, feature_type_column = feature_type_column,
                                feature_id = feature_id, feature_name = feature_name,
                                name_type = name_type, desc_column = desc_column, desc_type = desc_type)
    ## Get rid of any regions with 'N' in them
  }

  dropme <- grepl(x = sequence_df[["fwd_superprimer"]], pattern = "N")
  mesg("Dropped ", sum(dropme), " regions with Ns in the 5' region.")
  sequence_df <- sequence_df[!dropme, ]
  dropme <- grepl(x = sequence_df[["rev_superprimer"]], pattern = "N")
  mesg("Dropped ", sum(dropme), " regions with Ns in the 3' region.")
  sequence_df <- sequence_df[!dropme, ]

  #keepers <- c("variants", "rev_variants", "fwd_superprimer", "rev_superprimer",
  #             "fwd_primer", "rev_primer", "overlap_gene_id", "overlap_gene_description",
  #             "closest_gene_before_id", "closest_gene_before_description")
  #if (isTRUE(truncate)) {
  #  sequence_df <- sequence_df[, keepers]
  #}

  ## TODO:
  ## 1. Reindex based on primer start
  ## 2. Report closest CDS regions and be able to filter on multicopies
  ## 3. Report if in CDS or not
  ## 4. Extra credit: report polyN runs in the putative amplicon
  current_options <- options(new_options)
  retlist <- list(
    "densities" = density_lst,
    "density_vector" = long_density_vector,
    "variant_vector" = long_variant_vector,
    "favorites" = sequence_df)
  class(retlist) <- "density_primers"
  return(retlist)
}

#' Write out a set of primers for testing.
#'
#' @param density_primers List containing a series of putative sequencing/PCR primers.
#' @param prefix Sequence name prefix, 'pf' meaning 'primer forward'.
#' @param column Column from the dataframe of putative primers.
#' @param fasta Output filename.
#' @return DNAStringSet of the primers with side effect of written fasta file.
#' @export
write_density_primers <- function(density_primers, prefix = "pf",
                                  column = "fwd_primer",
                                  fasta = "forward_primers.fasta") {
  df <- density_primers[["favorites"]]
  string_set <- df[[column]]
  num_strings <- numform::f_pad_zero(seq_len(nrow(df)))
  name_strings <- rownames(df)
  sequence_names <- glue("{prefix}{num_strings}_{name_strings}")
  names(string_set) <- sequence_names
  string_set <- Biostrings::DNAStringSet(x = string_set, use.names = TRUE)
  written <- Biostrings::writeXStringSet(string_set, fasta, append = FALSE,
                                         compress = FALSE, format = "fasta")
  return(string_set)
}
