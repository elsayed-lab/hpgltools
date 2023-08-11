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
choose_sequence_regions <- function(vector, max_primer_length = 45,
                                    topn = 200, bin_width = 600,
                                    genome = NULL, target_temp = 58) {
  ## Now get the nucleotides of the first 30 nt of each window
  sequence_df <- as.data.frame(head(vector, n = topn))
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
  sequence_df[["fivep_references"]] <- gsub(x = sequence_df[["variants"]],
                                            pattern = "([[:alpha:]])(\\d+)([[:alpha:]])",
                                            replacement = "\\1")
  sequence_df[["threep_references"]] <- chartr("ATGC", "TACG", sequence_df[["fivep_references"]])

  sequence_df[["fivep_positions"]] <- gsub(x = sequence_df[["variants"]],
                                           pattern = "([[:alpha:]])(\\d+)([[:alpha:]])",
                                           replacement = "\\2")
  ## Getting the relative 3' position cannot be done with a simple regex/tr operation.
  ## but instead is binwidth - position
  ## which therefore requires a little apply() shenanigans
  recount_positions <- function(position_string, bin_width = 600) {
    positions <- strsplit(position_string, ",")[[1]]
    positions <- gsub(x = positions, pattern = " ", replacement = "")
    positions <- as.numeric(positions)
    new_positions <- toString(as.character(bin_width - positions))
    return(new_positions)
  }
  sequence_df[["threep_positions"]] <- mapply(FUN = recount_positions,
                                              sequence_df[["fivep_positions"]])

  sequence_df[["fivep_alternates"]] <- gsub(x = sequence_df[["variants"]],
                                            pattern = "([[:alpha:]])(\\d+)([[:alpha:]])",
                                            replacement = "\\3")
  sequence_df[["threep_alternates"]] <- chartr("ATGC", "TACG", sequence_df[["fivep_alternates"]])
  ## Now lets rewrite the variants as the combinations of threep_reference/position/alternate
  threep_variants <- function(ref_string, position_string, alt_string) {
    refs <- gsub(pattern = " ", replacement = "", x = strsplit(ref_string, ",")[[1]])
    pos <- gsub(pattern = " ", replacement = "", x = strsplit(position_string, ",")[[1]])
    alts <- gsub(pattern = " ", replacement = "", x = strsplit(alt_string, ",")[[1]])
    new <- toString(paste0(refs, pos, alts))
    return(new)
  }
  sequence_df[["threep_variants"]] <- mapply(FUN = threep_variants,
                                             sequence_df[["threep_references"]],
                                             sequence_df[["threep_positions"]],
                                             sequence_df[["threep_alternates"]])

  ## I am making explicit columns for these so I can look manually
  ## and make sure I am not trying to get sequences which run off the chromosome.
  sequence_df[["fivep_superprimer_start"]] <- sequence_df[["bin_start"]] -
    max_primer_length - 1
  sequence_df[["threep_superprimer_start"]] <- sequence_df[["bin_start"]] +
    bin_width + max_primer_length + 1
  sequence_df[["fivep_superprimer_end"]] <- sequence_df[["bin_start"]] - 1
  sequence_df[["threep_superprimer_end"]] <- sequence_df[["bin_start"]] + bin_width
  ## The super_primers are the entire region of length == max_primer_length, I will
  ## iteratively remove bases from the 5' until the target is reached.
  sequence_df[["fivep_superprimer"]] <- ""
  sequence_df[["threep_superprimer"]] <- ""
  ## The fivep_primer threep_primer columns will hold the final primers.
  sequence_df[["fivep_primer"]] <- ""
  sequence_df[["threep_primer"]] <- ""

  ## For me, this is a little too confusing to do in an apply.
  ## I want to fill in the fivep/threep_superprimer columns with the longer
  ## sequence/revcomp regions of interest.
  ## Then I want to chip away at the beginning/end of them to get the actual fivep/threep primer.
  for (i in seq_along(nrow(sequence_df))) {
    chr <- sequence_df[i, "chr"]

    silence = TRUE
    fivep_superprimer <- Biostrings::subseq(genome[[chr]],
                                            sequence_df[i, "fivep_superprimer_start"],
                                            sequence_df[i, "fivep_superprimer_end"])
    if ("try-error" %in% class(fivep_superprimer)) {
      sequence_df[i, "fivep_superprimer"] <- "Ran over chromosome end"
    } else if (sequence_df[i, "fivep_superprimer_start"] < 0) {
      ## Added this because subseq will assume a negative value means a circular chromosome.
      sequence_df[i, "fivep_superprimer"] <- "Ran over chromosome end"
    } else {
      fivep_superprimer <- as.character(fivep_superprimer)
      sequence_df[i, "fivep_superprimer"] <- fivep_superprimer
    }
    threep_superprimer <- try(Biostrings::subseq(genome[[chr]],
                                                 sequence_df[i, "threep_superprimer_end"],
                                                 sequence_df[i, "threep_superprimer_start"]),
                              silent = silence)
    if ("try-error" %in% class(threep_superprimer)) {
      sequence_df[i, "threep_superprimer"] <- "Ran over chromosome end"
    } else {
      threep_superprimer <- toupper(as.character(spgs::reverseComplement(threep_superprimer)))
      sequence_df[i, "threep_superprimer"] <- threep_superprimer
    }

    fivep_primer <- try(find_subseq_target_temp(fivep_superprimer,
                                                target = target_temp,
                                                direction = "forward"),
                        silent = silence)
    if ("try-error" %in% class(fivep_primer)) {
      fivep_primer <- "bad sequence for priming"
    }
    sequence_df[i, "fivep_primer"] <- fivep_primer

    threep_primer <- try(find_subseq_target_temp(threep_superprimer,
                                                 target = target_temp,
                                                 direction = "reverse"),
                         silent = silence)
    if ("try-error" %in% class(threep_primer)) {
      threep_primer <- "bad sequence for priming"
    }
    sequence_df[i, "threep_primer"] <- threep_primer
  } ## End iterating over sequence_df

  ## Drop anything which is definitely useless
  dropme <- grepl(x = sequence_df[["fivep_superprimer"]], pattern = "Ran")
  sequence_df <- sequence_df[!dropme, ]
  dropme <- grepl(x = sequence_df[["threep_superprimer"]], pattern = "Ran")
  sequence_df <- sequence_df[!dropme, ]

  wanted_columns <- c("variants", "threep_variants", "fivep_superprimer", "threep_superprimer",
                      "fivep_primer", "threep_primer", "overlap_gene_description",
                      "overlap_gene_start", "overlap_gene_end", "closest_gene_before_description",
                      "closest_gene_before_start", "closest_gene_before_end",
                      "closest_gene_after_id", "closest_gene_after_description",
                      "closest_gene_after_start", "closest_gene_after_end",
                      "chr", "bin_start", "fivep_references", "fivep_positions",
                      "fivep_alternates", "threep_references", "threep_positions",
                      "threep_alternates")
  ## sequence_df <- sequence_df[, wanted_columns]
  return(sequence_df)
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
sequential_variants <- function(snp_sets, conditions = NULL,
                                minimum = 3, maximum_separation = 3,
                                one_away_file = "one_away.csv",
                                two_away_file = "two_away.csv",
                                doubles_files = "doubles.csv",
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
#' @export
snp_density_primers <- function(snp_count, pdata_column = "condition",
                                condition = "z2.3", cutoff = 20, bin_width = 600,
                                divide = FALSE, topn = 400,
                                target_temp = 53, max_primer_length = 50,
                                bsgenome = "BSGenome.Leishmania.panamensis.MHOMCOL81L13.v52",
                                gff = "reference/lpanamensis_col_v46.gff",
                                feature_type = "protein_coding_gene", feature_start = "start",
                                feature_end = "end", feature_strand = "strand",
                                feature_chr = "seqnames", feature_type_column = "type",
                                feature_id = "ID", feature_name = "description",
                                truncate = TRUE) {

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
  interest_idx <- samples_by_condition == condition
  snp_table <- exprs(snp_count)[, interest_idx]
  start_rows <- nrow(snp_table)
  interest_rows <- rowMeans(snp_table) >= cutoff
  snp_table <- snp_table[interest_rows, ]
  mesg("Started with ", start_rows, ", after cutoff there are ",
          sum(interest_rows), " rows left.")
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
  for (ch in seq_along(chromosomes)) {
    chr <- chromosomes[ch]
    mesg("Starting chromosome: ", chr, ".")
    chr_idx <- position_table[["chromosome"]] == chr
    chr_data <- position_table[chr_idx, ]
    if (is.null(genome)) {
      chromosome_df[chr, "length"] <- max(chr_data[["position"]])
    } else {
      this_length <- length(genome[[chr]])
      chromosome_df[chr, "length"] <- this_length
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
    }
    density_lst[[chr]] <- density_vector
    variant_lst[[chr]] <- variant_vector
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
  }
  most_idx <- order(long_density_vector, decreasing = TRUE)
  long_density_vector <- long_density_vector[most_idx]
  long_variant_vector <- long_variant_vector[most_idx]

  mesg("Extracting primer regions.")
  sequence_df <- choose_sequence_regions(long_variant_vector, topn = topn,
                                         max_primer_length = max_primer_length,
                                         bin_width = bin_width, genome = genome)
  mesg("Searching for overlapping/closest genes.")
  sequence_df <- xref_regions(sequence_df, gff, bin_width = bin_width,
                              feature_type = feature_type, feature_start = feature_start,
                              feature_end = feature_end, feature_strand = feature_strand,
                              feature_chr = feature_chr, feature_type_column = feature_type_column,
                              feature_id = feature_id, feature_name = feature_name)
  ## Get rid of any regions with 'N' in them

  dropme <- grepl(x = sequence_df[["fivep_superprimer"]], pattern = "N")
  mesg("Dropped ", sum(dropme), " regions with Ns in the 5' region.")
  sequence_df <- sequence_df[!dropme, ]
  dropme <- grepl(x = sequence_df[["threep_superprimer"]], pattern = "N")
  mesg("Dropped ", sum(dropme), " regions with Ns in the 3' region.")
  sequence_df <- sequence_df[!dropme, ]

  keepers <- c("variants", "threep_variants", "fivep_superprimer", "threep_superprimer",
               "fivep_primer", "threep_primer", "overlap_gene_id", "overlap_gene_description",
               "closest_gene_before_id", "closest_gene_before_description")
  if (isTRUE(truncate)) {
    sequence_df <- sequence_df[, keepers]
  }

  ## TODO:
  ## 1. Reindex based on primer start
  ## 2. Report closest CDS regions and be able to filter on multicopies
  ## 3. Report if in CDS or not
  ## 4. Extra credit: report polyN runs in the putative amplicon
  current_options <- options(new_options)
  retlist <- list(
      "density_vector" = long_density_vector,
      "variant_vector" = long_variant_vector,
      "favorites" = sequence_df)
  class(retlist) <- "density_primers"
  return(retlist)
}
