#' Pull metadata from a table (xlsx/xls/csv/whatever)
#'
#' I find that when I acquire metadata from a paper or collaborator, annoyingly
#' often there are many special characters or other shenanigans in the column
#' names.  This function performs some simple sanitizations.  In addition, if I
#' give it a filename it calls my generic 'read_metadata()' function before
#' sanitizing.
#'
#' @param metadata file or df of metadata
#' @param id_column Column in the metadat containing the sample names.
#' @param ... Arguments to pass to the child functions (read_csv etc).
#' @return Metadata dataframe hopefully cleaned up to not be obnoxious.
#' @examples
#'  \dontrun{
#'   sanitized <- extract_metadata("some_random_supplemental.xls")
#'   saniclean <- extract_metadata(some_goofy_df)
#' }
#' @export
extract_metadata <- function(metadata, id_column = "sampleid", ...) {
  arglist <- list(...)
  ## FIXME: Now that this has been yanked into its own function,
  ## Make sure it sets good, standard rownames.
  file <- NULL

  meta_dataframe <- NULL
  meta_file <- NULL
  if ("character" %in% class(metadata)) {
    ## This is a filename containing the metadata
    meta_file <- metadata
  } else if ("data.frame" %in% class(metadata)) {
    ## A data frame of metadata was passed.
    meta_dataframe <- metadata
  } else {
    stop("This requires either a file or meta data.frame.")
  }

  ## The two primary inputs for metadata are a csv/xlsx file or a dataframe, check for them here.
  if (is.null(meta_dataframe) & is.null(meta_file)) {
    stop("This requires either a csv file or dataframe of metadata describing the samples.")
  } else if (is.null(meta_file)) {
    ## punctuation is the devil
    sample_definitions <- meta_dataframe
  }  else {
    sample_definitions <- read_metadata(meta_file,
                                        ...)
    ## sample_definitions <- read_metadata(meta_file)
  }

  colnames(sample_definitions) <- gsub(pattern = "[[:punct:]]",
                                       replacement = "",
                                       x = colnames(sample_definitions))
  id_column <- tolower(id_column)
  id_column <- gsub(pattern = "[[:punct:]]",
                    replacement = "",
                    x = id_column)

  ## Get appropriate row and column names.
  current_rownames <- rownames(sample_definitions)
  bad_rownames <- as.character(1:nrow(sample_definitions))
  ## Try to ensure that we have a useful ID column by:
  ## 1. Look for data in the id_column column.
  ##  a.  If it is null, look at the rownames
  ##    i.  If they are 1...n, arbitrarily grab the first column.
  ##    ii. If not, use the rownames.
  if (is.null(sample_definitions[[id_column]])) {
    if (identical(current_rownames, bad_rownames)) {
      id_column <- colnames(sample_definitions)[1]
    } else {
      sample_definitions[[id_column]] <- rownames(sample_definitions)
    }
  }

  ## Drop empty rows in the sample sheet
  empty_samples <- which(sample_definitions[, id_column] == "" |
                         grepl(x = sample_definitions[, id_column], pattern = "^undef") |
                         is.na(sample_definitions[, id_column]) |
                         grepl(pattern = "^#", x = sample_definitions[, id_column]))
  if (length(empty_samples) > 0) {
    message("Dropped ", length(empty_samples),
            " rows from the sample metadata because they were blank.")
    sample_definitions <- sample_definitions[-empty_samples, ]
  }

  ## Drop duplicated elements.
  num_duplicated <- sum(duplicated(sample_definitions[[id_column]]))
  if (num_duplicated > 0) {
    message("There are ", num_duplicated,
            " duplicate rows in the sample ID column.")
    sample_definitions[[id_column]] <- make.names(sample_definitions[[id_column]],
                                                  unique = TRUE)
  }

  ## Now we should have consistent sample IDs, set the rownames.
  rownames(sample_definitions) <- sample_definitions[[id_column]]
  ## Check that condition and batch have been filled in.
  sample_columns <- colnames(sample_definitions)

  ## The various proteomics data I am looking at annoyingly starts with a number
  ## So make.names() prefixes it with X which is ok as far as it goes, but
  ## since it is a 's'amplename, I prefer an 's'.
  rownames(sample_definitions) <- gsub(pattern = "^X([[:digit:]])",
                                       replacement = "s\\1",
                                       x = rownames(sample_definitions))

  sample_columns_to_remove <- NULL
  for (col in 1:length(colnames(sample_definitions))) {
    sum_na <- sum(is.na(sample_definitions[[col]]))
    sum_null <- sum(is.null(sample_definitions[[col]]))
    sum_empty <- sum_na + sum_null
    if (sum_empty ==  nrow(sample_definitions)) {
      ## This column is empty.
      sample_columns_to_remove <- append(sample_columns_to_remove, col)
    }
  }
  if (length(sample_columns_to_remove) > 0) {
    sample_definitions <- sample_definitions[-sample_columns_to_remove]
  }

  ## Now check for columns named condition and batch
  found_condition <- "condition" %in% sample_columns
  if (!isTRUE(found_condition)) {
    message("Did not find the condition column in the sample sheet.")
    message("Filling it in as undefined.")
    sample_definitions[["condition"]] <- "undefined"
  } else {
    ## Make sure there are no NAs in this column.
    na_idx <- is.na(sample_definitions[["condition"]])
    sample_definitions[na_idx, "condition"] <- "undefined"
  }
  found_batch <- "batch" %in% sample_columns
  if (!isTRUE(found_batch)) {
    message("Did not find the batch column in the sample sheet.")
    message("Filling it in as undefined.")
    sample_definitions[["batch"]] <- "undefined"
  } else {
    ## Make sure there are no NAs in this column.
    na_idx <- is.na(sample_definitions[["batch"]])
    sample_definitions[na_idx, "batch"] <- "undefined"
  }

  ## Double-check that there is a usable condition column
  ## This is also an instance of simplifying subsetting, identical to
  ## sample_definitions[["condition"]] I don't think I care one way or the other which I use in
  ## this case, just so long as I am consistent -- I think because I have trouble remembering the
  ## difference between the concept of 'row' and 'column' I should probably use the [, column] or
  ## [row, ] method to reinforce my weak neurons.
  if (is.null(sample_definitions[["condition"]])) {
    ## type and stage are commonly used, and before I was consistent about always having
    ## condition, they were a proxy for it.
    sample_definitions[["condition"]] <- tolower(paste(sample_definitions[["type"]],
                                                       sample_definitions[["stage"]], sep = "_"))
  }
  ## Extract out the condition names as a factor
  condition_names <- unique(sample_definitions[["condition"]])
  if (is.null(condition_names)) {
    warning("There is no 'condition' field in the definitions, this will make many
analyses more difficult/impossible.")
  }
  ## Condition and Batch are not allowed to be numeric, so if they are just numbers,
  ## prefix them with 'c' and 'b' respectively.
  pre_condition <- unique(sample_definitions[["condition"]])
  pre_batch <- unique(sample_definitions[["batch"]])
  sample_definitions[["condition"]] <- gsub(pattern = "^(\\d+)$", replacement = "c\\1",
                                            x = sample_definitions[["condition"]])
  sample_definitions[["batch"]] <- gsub(pattern = "^(\\d+)$", replacement = "b\\1",
                                        x = sample_definitions[["batch"]])
  sample_definitions[["condition"]] <- factor(sample_definitions[["condition"]],
                                              levels = unique(sample_definitions[["condition"]]),
                                              labels = pre_condition)
  sample_definitions[["batch"]] <- factor(sample_definitions[["batch"]],
                                          levels = unique(sample_definitions[["batch"]]),
                                          labels = pre_batch)
  return(sample_definitions)
}

#' Gather metadata for trimomatic.
#'
#' Simplify the specification for extracting trimomatic data.
#'
#' @param metadata Starting metadata
#' @param file_spec Filename containing the information of interest.
#' @param columns Arbitrarily set the column names
#' @param new_metadata Filename for the new metadata.
#' @param ... extra stuff for glue.
gather_trimomatic_metadata <- function(metadata,
                                       file_spec = "preprocessing/{meta[['sampleid']]}/outputs/*-trimomatic.out",
                                       columns = NULL, new_metadata = NULL, ...) {
  specification <- list(
      "trimomatic_input" = list(
          "file" = file_spec),
      "trimomatic_output" = list(
          "file" = file_spec),
      "trimomatic_ratio" = list(
          "column" = "trimomatic_percent"))
  if (!is.null(columns)) {
    for (d in 1:length(columns)) {
      element_name <- names(columns)[d]
      element_value <- columns[[element_name]]
      specification[[element_name]][["column"]] <- element_value
    }
  }
  result <- gather_preprocessing_metadata(starting_metadata = metadata,
                                          specification = specification,
                                          new_metadata = new_metadata, ...)
  return(result)
}

#' Automagically fill in a sample sheet with the results of the
#' various preprocessing tools.
#'
#' I am hoping to fill this little function out with a bunch of useful
#' file specifications and regular expressions.  If I do a good job,
#' then it should become trivial to fill in a sample sheet with lots
#' of fun useful numbers in preparations for creating a nice table
#' S1.  I am thinking to split this up into sections for
#' trimming/mapping/etc.  But for the moment I just want to add some
#' specifications/regexes and see if it proves itself robust.  If
#' Theresa reads this, I think this is another good candidate for a
#' true OO implmentation.  E.g. make a base-class for the metadata and
#' use S4 multi-dispatch to pick up different log files.  I wrote the
#' downstream functions with this in mind already, but I am too
#' stupid/lazy to do the full implementation until I am confident that
#' these functions/ideas actually have merit.
#'
#' @param starting_metadata Already existing sample sheet.
#' @param specification List containing one element for each new
#'  column to append to the sample sheet.  Each element in turn is a
#'  list containing column names and/or input filenames (and
#' presumably other stuff as I think of it).
#' @param new_metadata Filename to which to write the new metadata
#' @param verbose Currently just used to debug the regexes.
#' @param ... This is one of the few instances where I used
#' ... intelligently.  Pass extra variables to the file specification
#' and glue will pick them up (note the {species} entries in the
#' example specifications.
#' @return For the moment it just returns the modified metadata, I
#'  suspect there is something more useful it should do.
#' @export
gather_preprocessing_metadata <- function(starting_metadata, specification = NULL,
                                          new_metadata = NULL, verbose = FALSE, ...) {
  ## I want to create a set of specifications for different tasks:
  ## tnseq/rnaseq/assembly/phage/phylogenetics/etc.
  ## For the moment, the following is specific for phage assembly.
  if (is.null(specification)) {
    specification <- list(
        ## First task performed is pretty much always trimming
        "input_r1" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/scripts/*trim_*.sh"),
        "input_r2" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/scripts/*trim_*.sh"),
        "trimomatic_input" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.out"),
        "trimomatic_output" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.out"),
        "trimomatic_ratio" = list(
            "column" = "trimomatic_percent"),
        ## Second task is likely error correction
        "racer_changed" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*racer/racer.out"),
        "host_filter_species" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/host_species.txt"),
        ## After those, things can get pretty arbitrary...
        "hisat_single_concordant" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*hisat2_*/hisat2_*.err"),
        "hisat_multi_concordant" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*hisat2_*/hisat2_*.err"),
        "hisat_single_all" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*hisat2_*/hisat2_*.err"),
        "hisat_multi_all" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*hisat2_*/hisat2_*.err"),
        "hisat_singlecon_ratio" = list(
            "column" = "hisat_single_concordant_percent"),
        "hisat_singleall_ratio" = list(
            "column" = "hisat_single_all_percent"),
        "kraken_viral_classified" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*kraken_viral*/kraken_out.txt"),
        "kraken_viral_unclassified" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*kraken_viral*/kraken_out.txt"),
        "kraken_first_viral_species" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*kraken_viral*/kraken_report.txt"),
        "kraken_first_viral_species_reads" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*kraken_viral*/kraken_report.txt"),
        "kraken_standard_classified" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*kraken_standard*/kraken_out.txt"),
        "kraken_standard_unclassified" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*kraken_standard*/kraken_out.txt"),
        "kraken_first_standard_species" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*kraken_standard*/kraken_report.txt"),
        "kraken_first_standard_species_reads" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*kraken_standard*/kraken_report.txt"),
        "salmon_mapped" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*salmon_*/salmon.err"),
        "shovill_contigs" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*shovill_*/shovill.log"),
        "shovill_length" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*shovill_*/shovill.log"),
        "shovill_estlength" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*shovill_*/shovill.log"),
        "shovill_minlength" = list(
            "file" = "preprocessing/{meta[['sampleid']]}/outputs/*shovill_*/shovill.log"),
        "unicycler_lengths" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*unicycler/*final_assembly.fasta"),
        "unicycler_relative_coverage" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*unicycler/*final_assembly.fasta"),
        "phageterm_dtr_length" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*phageterm_*/direct-terminal-repeats.fasta"),
        "prodigal_positive_strand" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*prodigal_*/predicted_cds.gff"),
        "prodigal_negative_strand" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*prodigal_*/predicted_cds.gff"),
        "interpro_signalp_hits" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
        "interpro_phobius_hits" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
        "interpro_pfam_hits" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
        "interpro_tmhmm_hits" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
        "interpro_cdd_hits" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
        "interpro_smart_hits" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
        "interpro_gene3d_hits" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
        "interpro_superfamily_hits" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
        "tRNA_hits" = list(
            file = "preprocessing/{meta[['sampleid']]}/outputs/*prokka_*/{meta[['sampleid']]}.log"),
        "notes" = list(
            file = "preprocessing/{meta[['sampleid']]}/notes.txt")        
        )
  }
  if (is.null(new_metadata)) {
    new_metadata <- gsub(x = starting_metadata, pattern = "\\.xlsx$",
                         replacement = "_modified.xlsx")
  }

  meta <- extract_metadata(starting_metadata)
  for (entry in 1:length(specification)) {
    entry_type <- names(specification[entry])
    message("Starting ", entry_type, ".")
    new_column <- entry_type
    if (!is.null(specification[[entry_type]][["column"]])) {
      new_column <- specification[[entry_type]][["column"]]
    }
    if (new_column %in% colnames(meta)) {
      warning("Column: ", new_column, " already exists, replacing it.")
    }
    input_file_spec <- specification[[entry_type]][["file"]]
    new_entries <- dispatch_metadata_extract(meta, entry_type,
                                             input_file_spec,
                                             specification,
                                             verbose = verbose,
                                             ...)
    if (is.null(new_entries)) {
      message("Not including new entries for: ", new_column, ".")
    } else {
      meta[[new_column]] <- new_entries
    }
  }
  message("Writing new metadata to: ", new_metadata)
  written <- write_xlsx(data = meta, excel = new_metadata)
  return(new_metadata)
}


#' This is basically just a switch and set of regexes for finding the
#' numbers of interest in the various log files.
#'
#' When I initially wrote this, it made sense to me to have it
#' separate from the top-level function.  I am not sure that is true
#' now, having slept on it.
#'
#' @param meta Starting metadata
#' @param entry_type String which defines the type of log entry to
#'  hunt down.  If the specification does not include a column, this
#'  will be used as the column name to write to the metadata.
#' @param input_file_spec Glue specification defining the log file for
#'  each sample to hunt down.
#' @param specification This is the reason I am thinking having this
#'  as a separate function might be stupid.  I added it to make it
#'  easier to calculate ratios of column_x/column_y; but it is a
#'  def-facto argument to either get rid of input_file_spec as an arg
#'  or to just get rid of this function.
#' @param verbose used for testing regexes.
#' @param ... passed to glue to add more variables to the file spec.
#' @return Vector of entries which will be used to populate the new
#'  column in the metadata.
dispatch_metadata_extract <- function(meta, entry_type, input_file_spec,
                                      specification, verbose = FALSE, ...) {
  switchret <- switch(
      entry_type,
      "hisat_single_concordant" = {
        search <-"^\\s+\\d+ \\(.+\\) aligned concordantly exactly 1 time" 
        replace <- "^\\s+(\\d+) \\(.+\\) aligned concordantly exactly 1 time"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "hisat_multi_concordant" = {
        search <- "^\\s+\\d+ \\(.+\\) aligned concordantly >1 times"
        replace <- "^\\s+(\\d+) \\(.+\\) aligned concordantly >1 times"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "hisat_single_all" = {
        search <- "^\\s+\\d+ \\(.+\\) aligned exactly 1 time"
        replace <- "^\\s+(\\d+) \\(.+\\) aligned exactly 1 time"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "hisat_multi_all" = {
        search <- "^\\s+\\d+ \\(.+\\) aligned concordantly >1 times"
        replace <- "^\\s+(\\d+) \\(.+\\) aligned concordantly >1 times"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "hisat_singlecon_ratio" = {
        numerator_column <- specification[["hisat_single_concordant"]][["column"]]
        denominator_column <- specification[["trimomatic_input"]][["column"]]
        entries <- dispatch_metadata_ratio(meta, numerator_column, denominator_column)
      },
      "hisat_singleall_ratio" = {
        numerator_column <- "hisat_single_all"
        if (!is.null(specification[["hisat_single_all"]][["column"]])) {
          numerator_column <- specification[["hisat_single_all"]][["column"]]
        }
        denominator_column <- "trimomatic_input"
        if (!is.null(specification[["trimomatic_input"]][["column"]])) {
          denominator_column <- specification[["trimomatic_input"]][["column"]]
        }
        entries <- dispatch_metadata_ratio(meta, numerator_column, denominator_column)
      },
      "host_filter_species" = {
        search <- "^.*$"
        replace <- "(.*)"
        entries <- dispatch_regex_search(meta, search, replace, which = "all",
                                         input_file_spec, verbose = verbose)
      },
      "input_r1" = {
        search <- "^\\s+<\\(less .+\\).*$"
        replace <- "^\\s+<\\(less (.+?)\\).*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "input_r2" = {
        search <- "^\\s+<\\(less .+\\) <\\(less .+\\).*$"
        replace <- "^\\s+<\\(less .+\\) <\\(less (.+)\\).*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "kraken_viral_classified" = {
        search <- "^\\s+\\d+ sequences classified.*$"
        replace <- "^\\s+(\\d+) sequences classified.*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "kraken_viral_unclassified" = {
        search <- "^\\s+\\d+ sequences unclassified.*$"
        replace <- "^\\s+(\\d+) sequences unclassified.*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "kraken_first_viral_species" = {
        search <- "^.*s__.*\\t\\d+$"
        replace <- "^.*s__(.*)\\t\\d+$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "kraken_first_viral_species_reads" = {
        search <- "^.*s__.*\\t\\d+$"
        replace <- "^.*s__.*\\t(\\d+)$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "kraken_standard_classified" = {
        search <- "^\\s+\\d+ sequences classified.*$"
        replace <- "^\\s+(\\d+) sequences classified.*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "kraken_standard_unclassified" = {
        search <- "^\\s+\\d+ sequences unclassified.*$"
        replace <- "^\\s+(\\d+) sequences unclassified.*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "kraken_first_standard_species" = {
        search <- "^.*s__.*\\t\\d+$"
        replace <- "^.*s__(.*)\\t\\d+$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "kraken_first_standard_species_reads" = {
        search <- "^.*s__.*\\t\\d+$"
        replace <- "^.*s__.*\\t(\\d+)$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "notes" = {
        search <- "^.*$"
        replace <- "^(.*)$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         which = "all")
      },
      "racer_changed" = {
        search <- "^Number of changed positions"
        replace <- "^Number of changed positions\\s+(\\d+)$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         which = "first",
                                         ...)
      },
      "salmon_mapped" = {
        search <- "^.* [jointLog] [info] Counted .+ total reads in the equivalence classes$"
        replace <- "^.* [jointLog] [info] Counted (.+) total reads in the equivalence classes$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "shovill_contigs" = {
        ## [shovill] It contains 1 (min=131) contigs totalling 40874 bp.
        search <- "^\\[shovill\\] It contains \\d+ .*$"
        replace <- "^\\[shovill\\] It contains (\\d+) .*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         which = "last",
                                         ...)
      },
      "shovill_length" = {
        ## [shovill] It contains 1 (min=131) contigs totalling 40874 bp.
        search <- "^\\[shovill\\] It contains \\d+ .* contigs totalling \\d+ bp\\.$"
        replace <- "^\\[shovill\\] It contains \\d+ .* contigs totalling (\\d+) bp\\.$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         which = "last",
                                         ...)
      },
      "shovill_estlength" = {
        ## [shovill] Assembly is 109877, estimated genome size was 117464 (-6.46%)
        search <- "^\\[shovill\\] Assembly is \\d+\\, estimated genome size was \\d+.*$"
        replace <- "^\\[shovill\\] Assembly is \\d+\\, estimated genome size was (\\d+).*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         which = "last",
                                         ...)
      },
      "shovill_minlength" = {
        ## [shovill] It contains 1 (min=145) contigs totalling 109877 bp.
        search <- "^\\[shovill\\] It contains \\d+ \\(min=\\d+\\).*$"
        replace <- "^\\[shovill\\] It contains \\d+ \\(min=(\\d+)\\).*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         which = "last",
                                         ...)
      },
      "trimomatic_input" = {
        search <- "^Input Read Pairs: \\d+ .*$"
        replace <- "^Input Read Pairs: (\\d+) .*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "trimomatic_output" = {
        search <- "^Input Read Pairs: \\d+ Both Surviving: \\d+ .*$"
        replace <- "^Input Read Pairs: \\d+ Both Surviving: (\\d+) .*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         ...)
      },
      "trimomatic_ratio" = {
        ## I think we can assume that the trimomatic ratio will come immediately after input/output
        numerator_column <- "trimomatic_output"
        if (!is.null(specification[["trimomatic_output"]][["column"]])) {
          numerator_column <- specification[["trimomatic_output"]][["column"]]
        }
        denominator_column <- "trimomatic_input"
        if (!is.null(specification[["trimomatic_input"]][["column"]])) {
          denominator_column <- specification[["trimomatic_input"]][["column"]]
        }
        entries <- dispatch_metadata_ratio(meta, numerator_column,
                                           denominator_column, verbose = verbose)
      },
      "unicycler_lengths" = {
        ## >1 length=40747 depth=1.00x circular=true
        search <- "^>\\d+ length=\\d+ depth.*$"
        replace <- "^>\\d+ length=(\\d+) depth.*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         which = "all",
                                         ...)
      },
      "unicycler_relative_coverage" = {
        ## >1 length=40747 depth=1.00x circular=true
        search <- "^>\\d+ length=\\d+ depth=.*x.*$"
        replace <- "^>\\d+ length=\\d+ depth=(.*)x.*$"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose,
                                         which = "all",
                                         ...)
      },
      "tRNA_hits" = {
        ## >1 length=40747 depth=1.00x circular=true
        search <- "^.*Found \\d+ tRNAs"
        replace <- "^.*Found (\\d+) tRNAs"
        entries <- dispatch_regex_search(meta, search, replace,
                                         input_file_spec, verbose = verbose)
      },
      "phageterm_dtr_length" = {
        entries <- dispatch_fasta_lengths(meta, input_file_spec, verbose = verbose)
      },
      "prodigal_positive_strand" = {
        search <- "\\t\\+\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      "prodigal_negative_strand" = {
        search <- "\\t-\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      "interpro_signalp_hits" = {
        search <- "\\tSiglapP.*\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      "interpro_phobius_hits" = {
        search <- "\\tPhobius\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      "interpro_pfam_hits" = {
        search <- "\\tPfam\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      "interpro_tmhmm_hits" = {
        search <- "\\tTMHMM\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      "interpro_cdd_hits" = {
        search <- "\\tCDD\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      "interpro_smart_hits" = {
        search <- "\\tSMART\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      "interpro_gene3d_hits" = {
        search <- "\\tGene3D\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      "interpro_superfamily_hits" = {
        search <- "\\tSUPERFAMILY\\t"
        entries <- dispatch_count_lines(meta, search, input_file_spec, verbose = verbose)
      },
      {
        stop("I do not know this spec: ", entry_type)
      })
  if (!is.null(entries)) {
    entries <- gsub(pattern = ",", replacement = "", x = entries)
  }
  return(entries)
}

dispatch_count_lines <- function(meta, search, input_file_spec, verbose = verbose) {
  filenames_with_wildcards <- glue::glue(input_file_spec)
  message("Example filename: ", filenames_with_wildcards[1], ".")
  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in 1:nrow(meta)) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    if (length(input_file) == 0) {
      warning("There is no file matching: ", filenames_with_wildcards[row],
              ".")
      next
    }
    if (is.na(input_file)) {
      warning("The input file is NA for: ", filenames_with_wildcards[row], ".")
      next
    }

    input_handle <- file(input_file, "r", blocking = FALSE)
    input_vector <- readLines(input_handle)
    last_found <- NULL
    this_found <- NULL
    all_found <- c()
    num_hits <- sum(grepl(x = input_vector, pattern = search))
    close(input_handle)
    output_entries[row] <- num_hits
  } ## End looking at every row of the metadata
  return(output_entries)
}

dispatch_fasta_lengths <- function(meta, input_file_spec, verbose = verbose) {
  filenames_with_wildcards <- glue::glue(input_file_spec)
  message("Example filename: ", filenames_with_wildcards[1], ".")
  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in 1:nrow(meta)) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    if (length(input_file) == 0) {
      warning("There is no file matching: ", filenames_with_wildcards[row],
              ".")
      next
    }
    if (is.na(input_file)) {
      warning("The input file is NA for: ", filenames_with_wildcards[row], ".")
      next
    }

    dtrs <- Biostrings::readBStringSet(input_file)
    output_entries[row] <- width(dtrs[1])
  } ## End looking at every row of the metadata
  return(output_entries)
}

#' Given two metadata columns, print a ratio.
#'
#' @param meta metadata, contains the column names!
#' @param numerator_column what it says on the tin.
#' @param denominator_column what it says on the tin.
#' @param digits Number of significant digits to keep in the output.
#' @param verbose unsed for the moment.
dispatch_metadata_ratio <- function(meta, numerator_column = NULL,
                                    denominator_column = NULL, digits = 3,
                                    verbose = FALSE) {
  column_number <- ncol(meta)
  if (is.null(numerator_column)) {
    numerator_column <- colnames(meta)[ncol(meta)]
  }
  if (is.null(denominator_column)) {
    denominator_column <- colnames(meta)[ncol(meta) - 1]
  }

  entries <- NULL
  if (is.null(meta[[numerator_column]]) | is.null(meta[[denominator_column]])) {
    message("Missing data to calculate the ratio between: ", numerator_column,
            " and ", denominator_column, ".")
  } else {
    message("The numerator column is: ", numerator_column, ".")
    message("The denominator column is: ", denominator_column, ".")
    entries <- as.numeric(meta[[numerator_column]]) / as.numeric(meta[[denominator_column]])
    if (!is.null(digits)) {
      entries <- signif(entries, digits)
    }
  }
  return(entries)
}


#' Generic dispatcher to hunt down useful information from logs.
#'
#' Given the metadata, a couple of regular expressions, and a filename
#' specification, this should be able to pull out the interesting
#' number(s) from one logfile per sample from the metadata.
#'
#' @param meta Input metadata.
#' @param search regex used to go hunting for the line of interest.
#' @param replace probably the same regex with parentheses in place
#'  for gsub().
#' @param input_file_spec filename extractor expression.
#' @param extraction the replacement portion of gsub(). I am thinking
#'  to make it possible to have this function return more interesting
#'  outputs if this changes, but for the moment I am sort of assuming
#'  \\1 will always suffice.
#' @param which Usually 'first', which means grab the first match and get out.
#' @param verbose For testing regexes.
#' @param ... Used to pass extra variables to glue for finding files.
dispatch_regex_search <- function(meta, search, replace, input_file_spec,
                                  extraction = "\\1", which = "first", verbose = FALSE,
                                  ...) {
  arglist <- list(...)
  ##if (length(arglist) > 0) {
  ##  
  ##}
  filenames_with_wildcards <- glue::glue(input_file_spec,
                                         ...)
  message("Example filename: ", filenames_with_wildcards[1], ".")
  test_file <- Sys.glob(filenames_with_wildcards[1])
  if (length(test_file) == 0) {
    message("The first filename does not exist, assuming this method was not performed.")
    return(NULL)
  }
  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in 1:nrow(meta)) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    if (is.na(input_file)) {
      ## The file did not exist.
      next
    }
    if (length(input_file) == 0) {
      warning("There is no file matching: ", filenames_with_wildcards[row],
              ".")
      next
    }

    input_handle <- file(input_file, "r", blocking = FALSE)
    input_vector <- readLines(input_handle)
    last_found <- NULL
    this_found <- NULL
    all_found <- c()
    for (i in 1:length(input_vector)) {
      if (which == "first" & found == 1) {
        output_entries[row] <- last_found
        next
      }
      input_line <- input_vector[i]
      if (grepl(x = input_line, pattern = search)) {
        if (isTRUE(verbose)) {
          message("Found the correct line: ")
          message(input_line)
        }                  
        this_found <- gsub(x = input_line,
                           pattern = replace,
                           replacement = extraction)
        found <- found + 1
      } else {
        next
      }
      last_found <- this_found
      all_found <- c(all_found, this_found)
    } ## End looking at every line of the log file specified by the input file spec for this row
    close(input_handle)
    ## Handle cases where one might want to pull only the last entry in a log, or all of them.
    if (which == "last") {
      message("TESTING: ", last_found)
      output_entries[row] <- last_found
    } else if (which == "all") {
      output_entries[row] <- toString(all_found)
    }
  } ## End looking at every row of the metadata
  return(output_entries)
}


#' Given an expressionset, sanitize pData columns of interest.
#'
#' I wrote this function after spending a couple of hours confused
#' because one cell in my metadata said 'cure ' instead of 'cure' and
#' I could not figure out why chaos reigned in my analyses.  There is
#' a sister to this somewhere else which checks that the expected
#' levels of a metadata factor are consistent; this is because in
#' another analysis we essentially had a cell which said 'cyre' and a
#' similar data explosion occurred.
#'
#' @param expt Input expressionset
#' @param columns Set of columns to check, if left NULL, all columns
#'  will be molested.
#' @param na_string Fill NA values with a string.
#' @param lower Set everything to lowercase?
#' @param punct Remove punctuation?
#' @export
sanitize_expt_metadata <- function(expt, columns = NULL, na_string = "notapplicable",
                                   lower = TRUE, punct = TRUE) {
  pd <- pData(expt)
  if (is.null(columns)) {
    columns <- colnames(pd)
  }
  for (col in 1:length(columns)) {
    todo <- columns[col]
    mesg("Sanitizing metadata column: ", todo, ".")
    if (! todo %in% colnames(pd)) {
      mesg("The column ", todo, " is missing, skipping it (also warning this).")
      warning("The column ", todo, " is missing, skipping it.")
      next
    }
    ## First get rid of trailing/leading spaces, those anger me and are crazy hard to find
    pd[[todo]] <- gsub(pattern = "^[[:space:]]", replacement = "", x = pd[[todo]])
    pd[[todo]] <- gsub(pattern = "[[:space:]]$", replacement = "", x = pd[[todo]])
    ## Set the column to lowercase, I have recently had a rash of mixed case sample sheet data.
    if (isTRUE(lower)) {
      pd[[todo]] <- tolower(pd[[todo]])
    }
    ## I think punctuation needs to go
    if (isTRUE(punct)) {
      pd[[todo]] <- gsub(pattern = "[[:punct:]]", replacement = "", x = pd[[todo]])
    }
    if (!is.null(na_string)) {
      ## Set NAs to "NotApplicable"
      na_idx <- is.na(pd[[todo]])
      pd[na_idx, todo] <- na_string
    }
  } ## End iterating over the columns of interest
  
  pData(expt[["expressionset"]]) <- pd
  expt[["design"]] <- pd
  return(expt)
}

## EOF
