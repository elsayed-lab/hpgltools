#' Create a SummarizedExperiment given some metadata
#'
#' This function was taken from create_expt() and repurposed to create SummarizedExperiments.
#'
#' @param metadata Filename or table of metadata about the samples of interest.
#' @param gene_info Annotations for the genes in the count data.
#' @param count_dataframe Optional table of counts.
#' @param sanitize_rownames Clean up unruly gene IDs?
#' @param sample_colors Specify the colors for the samples?
#' @param title Provide a title for the experiment.
#' @param notes Provide arbitrary notes.
#' @param countdir (deprecated) Directory containing count tables.
#' @param include_type Used to specify types of genes/annotations to use.
#' @param include_gff Keep a copy of the gff with the data?
#' @param file_column Metadata column containing the counts for each sample.
#' @param id_column Non-default column containing the sample IDs.
#' @param savefile Filename to which to save a rda file of the data structure.
#' @param low_files I don't remember this, I bet it is deprecated.
#' @param annotation orgDB associated with this, primarily used with gsva-like tools.
#' @param palette Color palette when auto-choosing colors for the samples.
#' @param round Round the data if/when it is not integer?
#' @param tx_gene_map When using tximport, use this to convert from transcripts to genes.
#' @param ... Extra options.
#' @export
create_se <- function(metadata = NULL, gene_info = NULL, count_dataframe = NULL,
                      sanitize_rownames = FALSE, sample_colors = NULL, title = NULL,
                      notes = NULL, countdir = NULL, include_type = "all",
                      include_gff = NULL, file_column = "file", id_column = NULL,
                      savefile = NULL, low_files = FALSE, annotation = "org.Hs.eg.db",
                      palette = "Dark2", round = FALSE, tx_gene_map = NULL,
                      ...) {
  arglist <- list(...)  ## pass stuff like sep=, header=, etc here

  if (is.null(metadata)) {
    stop("This requires some metadata at minimum.")
  }
  ## I am learning about simplifying vs. preserving subsetting
  ## This is a case of simplifying and I believe one which is good because I
  ## just want the string out from my list. Lets assume that palette is in fact
  ## an element in arglist, I really don't care that the name of the resturn is
  ## 'palette'; I already knew that by asking for it.
  if (is.null(title)) {
    title <- "This is a summarized experiment."
  }
  if (is.null(notes)) {
    notes <- glue::glue("Created on {date()}.
")
  }
  ## An expressionset needs to have a Biobase::annotation() in order for
  ## GSEABase to work with it. Reading the documentation, these are primarily
  ## used for naming the type of microarray chip used.

  ## Palette for colors when auto-chosen
  file_suffix <- ".count.gz"
  if (!is.null(arglist[["file_suffix"]])) {
    file_suffix <- arglist[["file_suffix"]]
  }
  file_prefix <- ""
  if (!is.null(arglist[["file_prefix"]])) {
    file_prefix <- arglist[["file_prefix"]]
  }
  gff_type <- "all"
  if (!is.null(arglist[["include_type"]])) {
    gff_type <- arglist[["include_type"]]
  }

  if (is.null(id_column)) {
    id_column <- "sampleid"
  } else {
    id_column <- tolower(id_column)
    id_column <- gsub(pattern = "[[:punct:]]", replacement = "", x = id_column)
  }
  file_column <- tolower(file_column)
  file_column <- gsub(pattern = "[[:punct:]]", replacement = "", x = file_column)

  ## Read in the metadata from the provided data frame, csv, or xlsx.
  message("Reading the sample metadata.")
  sample_definitions <- extract_metadata(metadata, id_column = id_column,
                                         ...)
  ## sample_definitions <- extract_metadata(metadata)
  ## Add an explicit removal of the column named 'file' file column if the option file_column is NULL.
  ## This is a just in case measure to avoid conflicts.
  if (is.null(file_column)) {
    if (!is.null(metadata[["file"]])) {
      metadata[["previous_file_column"]] <- metadata[["file"]]
      message("file_column is NULL, moving this column to 'previous_file_column'.")
      metadata[["file"]] <- NULL
    }
  }
  message("The sample definitions comprises: ", nrow(sample_definitions),
          " rows(samples) and ", ncol(sample_definitions),
          " columns(metadata fields).")
  num_samples <- nrow(sample_definitions)
  ## Create a matrix of counts with columns as samples and rows as genes
  ## This may come from either a data frame/matrix, a list of files from the metadata
  ## or it can attempt to figure out the location of the files from the sample names.
  filenames <- NULL
  all_count_tables <- NULL
  ## This set of if() statements is too complex and requires some reworking.
  if (!is.null(count_dataframe)) {
    ## Lets set the order of the count data to that of the sample definitions.
    test_col_rownames <- all.equal(sort(colnames(count_dataframe)),
                                   sort(rownames(sample_definitions)))
    if (isTRUE(test_col_rownames)) {
      count_dataframe <- count_dataframe[, rownames(sample_definitions)]
    } else {
      message("The count table column names are: ",
              toString(sort(colnames(count_dataframe))))
      message("The  meta   data  row  names are: ",
              toString(sort(rownames(sample_definitions))))
      stop("The count table column names are not the same as the sample definition row names.")
    }
    all_count_tables <- data.table::as.data.table(count_dataframe, keep.rownames = "rownames")
    ## If neither of these cases is true, start looking for the files in the
    ## processed_data/ directory
  } else if (is.null(sample_definitions[[file_column]])) {
    stop("This requires a column containing the input data.")
  }

  ## At this point sample_definitions$file should be filled in no matter what;
  ## so read the files.
  tximport_data <- NULL
  ## The count_data list should include a set of IDs and tables which are coherent
  ## Therefore, we will want to check in with it later.
  ## Notably, it has slots for: 'kept_ids' which should match 1:1 with the slot 'kept_files',
  ## 'source' which should remind us if the data came from htseq/tximport/etc.
  ## and count_table which should have one column for every kept_id/kept_file.
  count_data <- NULL
  if (is.null(all_count_tables)) {
    ## If all_count_tables does not exist, then we want to read the various files
    ## in the sample definitions to get them.
    filenames <- as.character(sample_definitions[[file_column]])
    sample_ids <- rownames(sample_definitions)
    count_data <- read_counts_expt(sample_ids, filenames, countdir = countdir,
                                   ...)
    if (count_data[["source"]] == "tximport") {
      tximport_data <- list("raw" = count_data[["tximport"]],
                            "scaled" = count_data[["tximport_scaled"]])
    }
    all_count_tables <- count_data[["count_table"]]
  } else {
    ## if all_count_tables _did_ exist, then we already had the count tables and so
    ## count_data should have them and all ids as 'kept'.
    count_data <- list(
      "source" = "dataframe",
      "raw" = all_count_tables,
      "kept_ids" = rownames(sample_definitions))
    ## Remember that R does not like rownames to start with a number, and if they do
    ## I already changed the count table rownames to begin with 's'.
    count_data[["kept_ids"]] <- gsub(pattern = "^([[:digit:]])",
                                     replacement = "s\\1",
                                     x = count_data[["kept_ids"]])
  }
  ## Here we will prune the metadata for any files/ids which were dropped
  ## when reading in the count tables.
  kept_definitions_idx <- rownames(sample_definitions) %in% count_data[["kept_ids"]]
  if (sum(kept_definitions_idx) < length(kept_definitions_idx)) {
      warning("Some samples were removed when cross referencing the samples against the count data.")
  }
  sample_definitions <- sample_definitions[kept_definitions_idx, ]
  ## While we are removing stuff...
  ## I have had a couple data sets with incomplete counts, get rid of those rows
  ## before moving on.
  all_count_tables <- all_count_tables[complete.cases(all_count_tables), ]

  numeric_columns <- colnames(all_count_tables) != "rownames"
  for (col in colnames(all_count_tables)[numeric_columns]) {
    ## Ensure there are no stupid entries like target_id est_counts
    all_count_tables[[col]] <- as.numeric(all_count_tables[[col]])
  }
  ## Features like exon:alicethegene-1 are annoying and entirely too common in TriTrypDB data
  if (isTRUE(sanitize_rownames)) {
      all_count_tables[["rownames"]] <- gsub(pattern = "^exon:", replacement = "",
                                             x = all_count_tables[["rownames"]])
      all_count_tables[["rownames"]] <- make.names(gsub(pattern = ":\\d+", replacement = "",
                                                        x = all_count_tables[["rownames"]]),
                                                   unique = TRUE)
  }

  ## There is an important caveat here!!
  ## data.table::as.data.table(stuff, keep.rownames='column') will change the
  ## rownames to remove punctuation including ':'!  Which means that if I have a
  ## rowname that looks like 'LmjF.01.0010:mRNA', it will get changed to
  ## 'LmjF.01.0010.mRNA' Which will of course kill any downstream analyses which
  ## depend on consistent rownames between the count table and any tximport
  ## data, since the tximport data  will still have the ':'...
  ## I am not certain what the best solution is, I am thinking perhaps to recast
  ## the tximport data as a set of data tables so that whatever silly stuff it
  ## does it will at least do consistently.  That of course will have unintended
  ## consequences for other tools which use tximport data (DESeq2), but if my
  ## rownames are consistent, then other problems will be easier to handle via
  ## recasting the data to a matrix or df or whatever the downstream tool
  ## requires.
  ## In contrast, I can take a simpler but less transparent route and change the
  ## rownames of all the tximport imported data to match that returned by
  ## as.data.table().
  if (!is.null(tximport_data[["raw"]])) {
    rownames(tximport_data[["raw"]][["abundance"]]) <- gsub(
      pattern = ":", replacement = "\\.",
      x = rownames(tximport_data[["raw"]][["abundance"]]))
    rownames(tximport_data[["raw"]][["counts"]]) <- gsub(
      pattern = ":", replacement = "\\.",
      x = rownames(tximport_data[["raw"]][["counts"]]))
    rownames(tximport_data[["raw"]][["length"]]) <- gsub(
      pattern = ":", replacement = "\\.",
      x = rownames(tximport_data[["raw"]][["length"]]))
  }
  if (!is.null(tximport_data[["scaled"]])) {
    rownames(tximport_data[["scaled"]][["abundance"]]) <- gsub(
      pattern = ":", replacement = "\\.",
      x = rownames(tximport_data[["scaled"]][["abundance"]]))
    rownames(tximport_data[["scaled"]][["counts"]]) <- gsub(
      pattern = ":", replacement = "\\.",
      x = rownames(tximport_data[["scaled"]][["counts"]]))
    rownames(tximport_data[["scaled"]][["length"]]) <- gsub(
      pattern = ":", replacement = "\\.",
      x = rownames(tximport_data[["scaled"]][["length"]]))
  }

  ## Try a couple different ways of getting gene-level annotations into the expressionset.
  annotation <- NULL
  if (is.null(gene_info)) {
    ## Including, if all else fails, just grabbing the gene names from the count tables.
    if (is.null(include_gff)) {
      gene_info <- data.table::as.data.table(all_count_tables[["rownames"]],
                                             keep.rownames = "rownames")
      names(gene_info) <- "rownames"
    } else {
      ## Or reading a gff file.
      message("create_expt(): Reading annotation gff, this is slow.")
      annotation <- load_gff_annotations(gff = include_gff, type = gff_type)
      gene_info <- data.table::as.data.table(annotation, keep.rownames = "rownames")
    }
  } else if (class(gene_info)[[1]] == "list" & !is.null(gene_info[["genes"]])) {
    ## In this case, it is using the output of reading a OrgDB instance
    gene_info <- data.table::as.data.table(gene_info[["genes"]], keep.rownames = "rownames")
  } else if (class(gene_info)[[1]] == "data.table" || class(gene_info)[[1]] == "tbl_df") {
    ## Try to make the data table usage consistent by rownames.
    ## Sometimes we take these from data which did "keep.rownames='some_column'"
    ## Sometimes we take these from data which set rownames(dt)
    ## And sometimes the rownames were never set.
    ## Therefore I will use rownames(dt) as the master, dt$rownames as secondary, and
    ## as a fallback take the first column in the data.
    if (is.null(rownames(gene_info)) & is.null(gene_info[["rownames"]])) {
      gene_info[["rownames"]] <- make.names(rownames[[1]], unique = TRUE)
      message("Both rownames() and $rownames were null.")
    }
  } else {
    gene_info <- data.table::as.data.table(gene_info, keep.rownames = "rownames")
  }

  ## It turns out that loading the annotation information from orgdb/etc may not set the
  ## row names. Perhaps I should do that there, but I will add a check here, too.
  found_sum <- sum(gene_info[["rownames"]] %in% all_count_tables[["rownames"]])
  if (found_sum == 0) {
    if (!is.null(gene_info[["geneid"]])) {
      gene_info[["rownames"]] <- gene_info[["geneid"]]
      found_sum <- sum(gene_info[["rownames"]] %in% all_count_tables[["rownames"]])
    }
  }
  if (found_sum == 0) {
    warning("Even after changing the rownames in gene info, they do not match the count table.")
    message("Even after changing the rownames in gene info, they do not match the count table.")
    message("Here are the first few rownames from the count tables:")
    message(toString(head(all_count_tables[["rownames"]])))
    message("Here are the first few rownames from the gene information table:")
    message(toString(head(gene_info[["rownames"]])))
  } else {
    message("Matched ", found_sum, " annotations and counts.")
  }

  ## Take a moment to remove columns which are blank
  columns_to_remove <- NULL
  for (col in 1:length(colnames(gene_info))) {
    sum_na <- sum(is.na(gene_info[[col]]))
    sum_null <- sum(is.null(gene_info[[col]]))
    sum_empty <- sum_na + sum_null
    if (sum_empty ==  nrow(gene_info)) {
      ## This column is empty.
      columns_to_remove <- append(columns_to_remove, col)
    }
    ## While we are looping through the columns,
    ## Make certain that no columns in gene_info are lists or factors.
    ## FIXME: 202104: I am no longer sure why I changed factors.
    if (class(gene_info[[col]]) == "factor" |
        class(gene_info[[col]]) == "AsIs" |
        class(gene_info[[col]]) == "list") {
      gene_info[[col]] <- as.character(gene_info[[col]])
    }
  }
  if (length(columns_to_remove) > 0) {
    gene_info <- gene_info[-columns_to_remove]
  }

  ## There should no longer be blank columns in the annotation data.
  ## Maybe I will copy/move this to my annotation collection toys?
  ## This temporary id number will be used to ensure that the order of features in everything
  ## will remain consistent, as we will call order() using it later.
  all_count_tables[["temporary_id_number"]] <- 1:nrow(all_count_tables)
  message("Bringing together the count matrix and gene information.")
  ## The method here is to create a data.table of the counts and annotation data,
  ## merge them, then split them apart.

  ## Made a small change to check for new tximport rownames in the gene information.
  ## This should automagically check and fix rownames when they would otherwise
  ## not match after using tximport.
  if (!is.null(tx_gene_map)) {
    matched_rows <- sum(rownames(gene_info) %in% tx_gene_map[[2]])
    if (matched_rows < 1) {
      message("The mapped IDs are not the rownames of your gene information, changing them now.")
      if (names(tx_gene_map)[2] %in% colnames(gene_info)) {
        new_name <- names(tx_gene_map)[2]
        rownames(gene_info) <- make.names(tx_gene_map[[new_name]], unique = TRUE)
      } else {
        warning("Cannot find an appropriate column in gene_info, refusing to use the tx_map.")
      }
    }
  }

  counts_and_annotations <- merge(all_count_tables, gene_info, by = "rownames", all.x = TRUE)
  ## In some cases, the above merge will result in columns being set to NA
  ## We should set all the NA fields to something I think.
  na_entries <- is.na(counts_and_annotations)
  if (sum(na_entries) > 0) {
    message("Some annotations were lost in merging, setting them to 'undefined'.")
  }
  counts_and_annotations[na_entries] <- "undefined"
  ## Set an incrementing id number to make absolutely paranoidly certain the
  ## order stays constant.
  counts_and_annotations <- counts_and_annotations[
    order(counts_and_annotations[["temporary_id_number"]]), ]
  ## Pull out the annotation data and convert to data frame.
  kept_columns <- colnames(counts_and_annotations) %in% colnames(gene_info)
  final_annotations <- counts_and_annotations[, kept_columns, with = FALSE]
  final_annotations <- as.data.frame(final_annotations, stringsAsFactors = FALSE)
  rownames(final_annotations) <- final_annotations[["rownames"]]
  final_kept <- colnames(final_annotations) != "rownames"
  final_annotations <- final_annotations[, final_kept]

  ## There are some shenanigans, Maddy is getting an error on countsdt...
  final_counts <- counts_and_annotations
  kept_columns <- colnames(counts_and_annotations) %in% colnames(all_count_tables) &
    colnames(counts_and_annotations) != "temporary_id_number"
  final_counts <- final_counts[, kept_columns, with = FALSE]
  final_counts <- as.data.frame(final_counts)
  rownames(final_counts) <- final_counts[["rownames"]]
  final_kept <- colnames(final_counts) != "rownames"
  final_counts <- final_counts[, final_kept]
  final_counts <- as.matrix(final_counts)

  ## I found a non-bug but utterly obnoxious behaivor in R
  ## Imagine a dataframe with 2 entries: TcCLB.511511.3 and TcCLB.511511.30
  ## Then imagine that TcCLB.511511.3 gets removed because it is low abundance.
  ## Then imagine what happens if I go to query 511511.3...
  ## Here are some copy/pasted lines illustrating it:
  ## > find_fiveeleven["TcCLB.511511.3", ]
  ##                 logFC AveExpr    t   P.Value adj.P.Val     B    qvalue
  ## TcCLB.511511.30  5.93   6.315 69.6 1.222e-25 7.911e-25 48.86 9.153e-27
  ## Here is the line in the dataframe documentation explaining this nonsense:
  ## https://stat.ethz.ch/R-manual/R-devel/library/base/html/Extract.data.frame.html
  ## Both [ and [[ extraction methods partially match row names. By default neither partially
  ## match column names, but [[ will if exact = FALSE (and with a warning if exact = NA). If you
  ## want to exact matching on row names use match, as in the examples
  ## How about you go eff yourself?  If you then look carefully at the match help, you will see
  ## that this is a feature, not a bug and that if you want truly exact matches, then the string
  ## must not end with a numeric value... oooo....kkk....
  ## Therefore, the following line replaces a terminal numeric rowname with the number and .
  ## > test_df = data.frame(a = c(1,1,1), b = c(2,2,2))
  ## > rownames(test_df) = c("TcCLB.511511.3","TcCLB.511511.30","bob")
  ## > test_df
  ##                 a b
  ## TcCLB.511511.3  1 2
  ## TcCLB.511511.30 1 2
  ## bob             1 2
  ## > rownames(test_df) <- gsub(pattern = "(\\d)$", replacement = "\\1\\.", x = rownames(test_df))
  ## > test_df
  ##                  a b
  ## TcCLB.511511.3.  1 2
  ## TcCLB.511511.30. 1 2
  ## bob              1 2
  ## This is so stupid I think I am going to go and cry.
  ##rownames(final_annotations) <- gsub(pattern = "(\\d)$", replacement = "\\1\\.",
  ##                                    x = rownames(final_annotations), perl = TRUE)
  ##rownames(final_counts) <- gsub(pattern = "(\\d)$", replacement = "\\1\\.",
  ##                                    x = rownames(final_counts), perl = TRUE)

  ##final_counts <- final_counts[, -1, drop = FALSE]

  ## If the user requests input of non-int counts, fix that here.
  if (isTRUE(round)) {
    final_counts <- round(final_counts)
    less_than <- final_counts < 0
    final_counts[less_than] <- 0
  }

  ## I moved the color choices to this area pretty late in the process to make sure that
  ## there was time to remove unused samples.
  ## Make sure we have a viable set of colors for plots
  chosen_colors <- generate_expt_colors(sample_definitions, sample_colors = sample_colors,
                                        chosen_palette = palette)

  ## Fill in incomplete tables.
  if (is.null(sample_definitions[["condition"]])) {
    sample_definitions[["condition"]] <- "unknown"
  }
  if (is.null(sample_definitions[["batch"]])) {
    sample_definitions[["batch"]] <- "unknown"
  }
  if (is.null(sample_definitions[["file"]])) {
    sample_definitions[["file"]] <- "null"
  }

  requireNamespace("SummarizedExperiment")
  ## SummarizedExperiments vs. ExpressionSets:
  ## assays() vs. exprs()
  ## rowData()/rowRanges() vs. fData()
  ## colData() vs. pData()
  ## Samples metadata access via $ accessor
  ## Experimental metadata (e.g. publication, lab, sra, whatever) via metadata()
  ## I may need to do some reorganizing to avoid confusion between my single experimental metadata
  ## and the metadata() provided by SummarizedExperiment
  ## Note that metadata() is just a list, so anything may be dumped here.
  se <- SummarizedExperiment(assays = final_counts,
                             rowData = final_annotations,
                             colData = metadata)
  metadata(se)[["notes"]] <- notes
  metadata(se)[["title"]] <- title
  metadata(se)[["annotation"]] <- annotation
  metadata(se)[["gff_file"]] <- include_gff
  ## the 'state' slot in the expt is used to keep track of how the data is modified over time.
  starting_state <- list(
    "filter" = "raw",
    "normalization" = "raw",
    "conversion" = "raw",
    "batch" = "raw",
    "transform" = "raw")
  metadata(se)[["state"]] <- starting_state
  se_conditions <- sample_definitions[["condition"]]
  names(se_conditions) <- rownames(sample_definitions)
  metadata(se)[["conditions"]] <- se_conditions
  se_batches <- sample_definitions[["batch"]]
  names(se_batches) <- rownames(sample_definitions)
  metadata(se)[["batches"]] <- se_batches
  se_libsizes <- colSums(final_counts)
  names(se_libsizes) <- rownames(sample_definitions)
  metadata(se)[["libsize"]] <- se_libsizes

  if (sum(se_libsizes == 0) > 0) {
    zero_idx <- se_libsizes == 0
    zero_samples <- names(se_libsizes)[zero_idx]
    warning("The following samples have no counts! ", zero_samples)
  }

  ## Save the chosen colors
  names(chosen_colors) <- rownames(sample_definitions)
  metadata(se)[["colors"]] <- chosen_colors
  metadata(se)[["tximport"]] <- tximport_data

  ## Save an rdata file of the expressionset.
  if (is.null(savefile)) {
    if ("character" %in% class(metadata)) {
      name <- paste0(gsub(x = basename(metadata), pattern = "^(.*)\\..*", replacement = "\\1"), ".rda")
    } else {
      message("Saving the expressionset to 'expt.rda'.")
      savefile <- "expt.rda"
    }
  }
  save_result <- try(save(expt, file = savefile), silent = TRUE)
  if (class(save_result) == "try-error") {
    warning("Saving the expt object failed, perhaps you do not have permissions?")
  }
  message("The final expressionset has ", nrow(assays(se)),
          " rows and ", ncol(colData(se)), " columns.")
  return(se)
}

#' Analagous function to make_pombe_expt()
#'
#' @export
make_pombe_se <- function(annotation = TRUE) {
  fission <- new.env()
  tt <- sm(please_install("fission"))
  tt <- sm(requireNamespace("fission"))
  tt <- sm(try(attachNamespace("fission"), silent = TRUE))
  tt <- data(fission, envir = fission)
  ## some minor shenanigans to get around the oddities of loading from data()
  fission <- fission[["fission"]]
  meta <- as.data.frame(fission@colData)
  meta[["condition"]] <- glue::glue("{meta[['strain']]}.{meta[['minute']]}")
  meta[["batch"]] <- meta[["replicate"]]
  meta[["sample.id"]] <- rownames(meta)
  meta <- meta[, c("sample.id", "id", "strain", "minute",
                   "replicate", "condition", "batch")]

  fission_data <- fission@assays$data[["counts"]]

  annotations <- NULL
  if (isTRUE(annotation)) {
    ## Neat, it works, and even figures out that the default mart is incorrect by itself.
    pombe_annotations <- load_biomart_annotations(
        host = "fungi.ensembl.org", trymart = "fungi_mart",
        trydataset = "spombe_eg_gene",
        gene_requests = c("pombase_transcript", "ensembl_gene_id", "ensembl_transcript_id",
                        "hgnc_symbol", "description", "gene_biotype"),
        species = "spombe", overwrite = TRUE)
    pombe_mart <- pombe_annotations[["mart"]]
    annotations <- pombe_annotations[["annotation"]]
    rownames(annotations) <- make.names(gsub(pattern = "\\.\\d+$",
                                             replacement = "",
                                             x = rownames(annotations)), unique = TRUE)
  }
  pombe_se <- sm(create_se(metadata = meta,
                           count_dataframe = fission_data,
                           gene_info = annotations))
  detach("package:fission")
  return(pombe_se)
}


setMethod("exprs", signature = "SummarizedExperiment",
          function(object) {
            SummarizedExperiment::assay(object)
          })
setMethod("exprs<-", signature = "SummarizedExperiment",
          function(object, value) {
            SummarizedExperiment::assay(object) <- value
            return(object)
          })
setMethod("fData", signature = "SummarizedExperiment",
          function(object) {
            SummarizedExperiment::rowData(object)
          })
setMethod("fData<-", signature = "SummarizedExperiment",
          function(object, value) {
            SummarizedExperiment::rowData(object) <- value
            return(object)
          })
setMethod("pData", signature = "SummarizedExperiment",
          function(object) {
            SummarizedExperiment::colData(object)
          })
setMethod("pData<-", signature = "SummarizedExperiment",
          function(object, value) {
            SummarizedExperiment::rowData(object) <- value
            return(object)
          })
