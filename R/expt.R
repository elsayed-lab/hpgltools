## expt.r: Create and play with expressionsets.  This actually is an
## implementation of a non-S4 superclass of the expressionSet.  I did this
## because I was having annoying problems creating expressionSets.  Thus, >90%
## of the logic here is intended to simplify and standardize that process.

#' Take two expressionsets and smoosh them together.
#'
#' Because of the extra sugar I added to expressionSets, the combine() function
#' needs a little help when combining expts.  Notably, the information from
#' tximport needs some help.
#'
#' @param expt1 First expt object.
#' @param expt2 Second expt object.
#' @param condition Column with which to reset the conditions.
#' @param all_x Keep all of the first expt's annotations/counts if there are mismatches?
#' @param all_y Keep all the second expt's annotations/counts if there are mismatches?
#' @param batch Column with which to reset the batches.
#' @param merge_meta Merge the metadata when they mismatch?  This should perhaps default to TRUE.
#' @return Larger expt.
#' @seealso [set_expt_batches()] [set_expt_conditions()] [set_expt_colors()]
#'  [set_expt_genenames()] [set_expt_samplenames()] [subset_expt()] [create_expt()]
#' @examples
#'  \dontrun{
#'   ## I am trying to get rid of all my dontrun sections, but I don't have two
#'   ## expressionsets to combine.
#'   expt1 <- create_expt(first_meta)
#'   expt2 <- create_expt(second_meta)
#'   combined <- combine_expts(expt1, expt2, merge_meta = TRUE)
#' }
#' @export
combine_expts <- function(expt1, expt2, condition = "condition", all_x = TRUE, all_y = TRUE,
                          batch = "batch", merge_meta = TRUE) {
  exp1 <- expt1[["expressionset"]]
  exp2 <- expt2[["expressionset"]]
  fData(exp2) <- fData(exp1)

  if (isTRUE(merge_meta)) {
    design1 <- pData(exp1)
    d1_rows <- 1:nrow(design1)
    design2 <- pData(exp2)
    both <- as.data.frame(data.table::rbindlist(list(design1, design2), fill = TRUE))
    ## na_idx <- is.na(both)
    ## both[na_idx] <- ""
    d2_rows <- (nrow(design1) + 1):nrow(both)
    new_design1 <- both[d1_rows, ]
    rownames(new_design1) <- rownames(design1)
    new_design2 <- both[d2_rows, ]
    rownames(new_design2) <- rownames(design2)
    pData(exp1) <- new_design1
    pData(exp2) <- new_design2
  }

  new <- combine(exp1, exp2)
  expt1[["expressionset"]] <- new
  expt1[["design"]] <- pData(new)
  expt1[["conditions"]] <- pData(expt1)[["condition"]]
  names(expt1[["conditions"]]) <- rownames(pData(expt1))
  expt1[["batches"]] <- pData(expt1)[["batch"]]
  names(expt1[["batches"]]) <- rownames(pData(expt1))
  expt1[["colors"]] <- c(expt1[["colors"]], expt2[["colors"]])
  expt1 <- set_expt_conditions(expt1, fact = condition)

  if (!is.null(expt1[["tximport"]])) {
    raw1 <- expt1[["tximport"]][["raw"]]
    raw2 <- expt2[["tximport"]][["raw"]]
    merged <- merge(raw1[["abundance"]], raw2[["abundance"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged <- merged[, -1]
    expt1[["tximport"]][["raw"]][["abundance"]] <- merged
    merged <- merge(raw1[["counts"]], raw2[["counts"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged <- merged[, -1]
    expt1[["tximport"]][["raw"]][["counts"]] <- merged
    merged <- merge(raw1[["length"]], raw2[["length"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged <- merged[, -1]
    expt1[["tximport"]][["raw"]][["length"]] <- merged
    scaled1 <- expt1[["tximport"]][["scaled"]]
    scaled2 <- expt2[["tximport"]][["scaled"]]
    merged <- merge(scaled1[["abundance"]], scaled2[["abundance"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged <- merged[, -1]
    expt1[["tximport"]][["scaled"]][["abundance"]] <- merged
    merged <- merge(scaled1[["counts"]], scaled2[["counts"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged <- merged[, -1]
    expt1[["tximport"]][["scaled"]][["counts"]] <- merged
    merged <- merge(scaled1[["length"]], scaled2[["length"]], by = "row.names")
    rownames(merged) <- merged[["Row.names"]]
    merged <- merged[, -1]
    expt1[["tximport"]][["scaled"]][["length"]] <- merged
  }

  ## It appears combining expressionsets can lead to NAs?
  new_exprs <- exprs(expt1)
  na_idx <- is.na(new_exprs)
  new_exprs[na_idx] <- 0
  tmp <- expt1[["expressionset"]]
  exprs(tmp) <- as.matrix(new_exprs)
  expt1[["expressionset"]] <- tmp

  if (isFALSE(all_x) | isFALSE(all_y)) {
    found_first <- rownames(exprs(expt1)) %in% rownames(exprs(exp1))
    shared_first_ids <- rownames(exprs(expt1))[found_first]
    expt1 <- exclude_genes_expt(expt1, ids = shared_first_ids, method = "keep")

    found_second <- rownames(exprs(expt1)) %in% rownames(exprs(exp2))
    shared_second_ids <- rownames(exprs(expt1))[found_second]
    expt1 <- exclude_genes_expt(expt1, ids = shared_second_ids, method = "keep")
  }

  return(expt1)
}

#' Sum the reads/gene for multiple sequencing runs of a single condition/batch.
#'
#' On occasion we have multiple technical replicates of a sequencing run.  This
#' can use a column in the experimental design to identify those replicates and
#' sum the counts into a single column in the count tables.
#'
#' Untested as of 2016-12-01, but used in a couple of projects where sequencing
#' runs got repeated.
#'
#' @param expt Experiment class containing the requisite metadata and count tables.
#' @param column Column of the design matrix used to specify which samples are replicates.
#' @return Expt with the concatenated counts, new design matrix, batches, conditions, etc.
#' @seealso [Biobase] [exprs()] [fData()] [pData()] [create_expt()]
#' @examples
#' \dontrun{
#'  compressed <- concatenate_runs(expt)
#' }
#' @export
concatenate_runs <- function(expt, column = "replicate") {
  design <- pData(expt)
  message("The original expressionset has ", nrow(design), " samples.")
  replicates <- levels(as.factor(design[[column]]))
  final_expt <- expt
  final_data <- NULL
  final_design <- NULL
  column_names <- list()
  colors <- list()
  conditions <- list()
  batches <- list()
  samplenames <- list()
  for (rep in replicates) {
    ## expression <- paste0(column, "=='", rep, "'")
    expression <- glue::glue("{column} == '{rep}'")
    tmp_expt <- sm(subset_expt(expt, expression))
    tmp_data <- rowSums(exprs(tmp_expt))
    tmp_design <- pData(tmp_expt)[1, ]
    final_data <- cbind(final_data, tmp_data)
    final_design <- rbind(final_design, tmp_design)
    column_names[[rep]] <- as.character(tmp_design[, "sampleid"])
    colors[[rep]] <- as.character(tmp_expt[["colors"]][1])
    batches[[rep]] <- as.character(pData(tmp_expt)[["batch"]][1])
    conditions[[rep]] <- as.character(pData(tmp_expt)[["condition"]][1])
    samplenames[[rep]] <- paste(conditions[[rep]], batches[[rep]], sep = "-")
    colnames(final_data) <- column_names
  }
  final_expt[["design"]] <- final_design
  metadata <- new("AnnotatedDataFrame", final_design)
  sampleNames(metadata) <- colnames(final_data)
  feature_data <- new("AnnotatedDataFrame", fData(expt))
  featureNames(feature_data) <- rownames(final_data)
  experiment <- new("ExpressionSet", exprs = final_data,
                    phenoData = metadata, featureData = feature_data)
  final_expt[["expressionset"]] <- experiment
  final_expt[["samples"]] <- final_design
  final_expt[["colors"]] <- as.character(colors)
  final_expt[["batches"]] <- as.character(batches)
  final_expt[["conditions"]] <- as.character(conditions)
  final_expt[["samplenames"]] <- as.character(samplenames)
  message("The final expressionset has ", nrow(pData(final_expt)), " samples.")
  return(final_expt)
}

#' Wrap bioconductor's expressionset to include some extra information.
#'
#' The primary innovation of this function is that it will check the metadata
#' for columns containing filenames for the count tables, thus hopefully making
#' the collation and care of metadata/counts easier.  For example, I have some
#' data which has been mapped against multiple species.  I can use this function
#' and just change the file_column argument to pick up each species' tables.
#'
#' @param metadata Comma separated file (or excel) describing the samples with
#'  information like condition, batch, count_filename, etc.
#' @param gene_info Annotation information describing the rows of the data set,
#'  this often comes from a call to import.gff() or biomart or organismdbi.
#' @param count_dataframe If one does not wish to read the count tables from the
#'  filesystem, they may instead be fed as a data frame here.
#' @param sanitize_rownames Clean up weirdly written gene IDs?
#' @param sample_colors List of colors by condition, if not provided it will
#'  generate its own colors using colorBrewer.
#' @param title Provide a title for the expt?
#' @param notes Additional notes?
#' @param countdir Directory containing count tables.
#' @param include_type I have usually assumed that all gff annotations should be
#'  used, but that is not always true, this allows one to limit to a specific
#'  annotation type.
#' @param include_gff Gff file to help in sorting which features to keep.
#' @param file_column Column to use in a gene information dataframe for
#' @param id_column Column which contains the sample IDs.
#' @param savefile Rdata filename prefix for saving the data of the resulting
#'  expt.
#' @param low_files Explicitly lowercase the filenames when searching the
#'  filesystem?
#' @param handle_na How does one wish to deal with NA values in the data?
#' @param researcher Used to make the creation of gene sets easier, set the researcher tag.
#' @param study_name Ibid, but set the study tag.
#' @param annotation_name Ibid, but set the orgdb (or other annotation) instance.
#' @param ... More parameters are fun!
#' @return experiment an expressionset
#' @seealso [Biobase] [cdm_expt_rda] [example_gff] [sb_annot] [sb_data] [extract_metadata()]
#'  [set_expt_conditions()] [set_expt_batches()] [set_expt_samplenames()] [subset_expt()]
#'  [set_expt_colors()] [set_expt_genenames()] [tximport] [load_annotations()]
#' @examples
#'  cdm_expt_rda <- system.file("share", "cdm_expt.rda", package = "hpgltools")
#'  load(file = cdm_expt_rda)
#'  head(cdm_counts)
#'  head(cdm_metadata)
#'  ## The gff file has differently labeled locus tags than the count tables, also
#'  ## the naming standard changed since this experiment was performed, therefore I
#'  ## downloaded a new gff file.
#'  example_gff <- system.file("share", "gas.gff", package = "hpgltools")
#'  gas_gff_annot <- load_gff_annotations(example_gff)
#'  rownames(gas_gff_annot) <- make.names(gsub(pattern = "(Spy)_", replacement = "\\1",
#'                                             x = gas_gff_annot[["locus_tag"]]), unique = TRUE)
#'  mgas_expt <- create_expt(metadata = cdm_metadata, gene_info = gas_gff_annot,
#'                           count_dataframe = cdm_counts)
#'  head(pData(mgas_expt))
#'  ## An example using count tables referenced in the metadata.
#'  sb_annot <- system.file("share", "sb", "trinotate_head.csv.xz", package = "hpgltools")
#'  sb_annot <- load_trinotate_annotations(trinotate = sb_annot)
#'  sb_annot <- as.data.frame(sb_annot)
#'  rownames(sb_annot) <- make.names(sb_annot[["transcript_id"]], unique = TRUE)
#'  sb_annot[["rownames"]] <- NULL
#'  sb_data <- system.file("share", "sb", "preprocessing.tar.xz", package = "hpgltools")
#'  untarred <- utils::untar(tarfile = sb_data)
#'  sb_expt <- create_expt(metadata = "preprocessing/kept_samples.xlsx",
#'                         gene_info = sb_annot)
#'  dim(exprs(sb_expt))
#'  dim(fData(sb_expt))
#'  pData(sb_expt)
#'  ## There are lots of other ways to use this, for example:
#'  \dontrun{
#'   new_experiment <- create_expt(metadata = "some_csv_file.csv", gene_info = gene_df)
#'   ## Remember that this depends on an existing data structure of gene annotations.
#'   meta <- extract_metadata("some_supplementary_materials_xls_file_I_downloaded.xls")
#'   another_expt <- create_expt(metadata = meta, gene_info = annotations, count_dataframe = df_I_downloaded)
#'  }
#' @import Biobase
#' @export
create_expt <- function(metadata = NULL, gene_info = NULL, count_dataframe = NULL,
                        sanitize_rownames = FALSE, sample_colors = NULL, title = NULL,
                        notes = NULL, countdir = NULL, include_type = "all",
                        include_gff = NULL, file_column = "file", id_column = NULL,
                        savefile = NULL, low_files = FALSE, handle_na = "drop",
                        researcher = "elsayed", study_name = NULL,
                        annotation_name = "org.Hs.eg.db", tx_gene_map = NULL, ...) {
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
    title <- "This is an expt class."
  }
  if (is.null(notes)) {
    notes <- glue::glue("Created on {date()}.
")
  }

  ## Palette for colors when auto-chosen
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }
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

  round <- FALSE
  if (!is.null(arglist[["round"]])) {
    round <- arglist[["round"]]
  }

  ## Read in the metadata from the provided data frame, csv, or xlsx.
  message("Reading the sample metadata.")
  sample_definitions <- extract_metadata(metadata, id_column = id_column,
                                         ...)
  ## sample_definitions <- extract_metadata(metadata)
  ## Add an explicit removal of the file column if the option file_column is NULL.
  ## This is a just in case measure to avoid conflicts.
  if (is.null(file_column)) {
    if (!is.null(metadata[["file"]])) {
      message("file_column is NULL, removing the column named 'file'.")
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
      message("The metadata row names are: ",
              toString(sort(rownames(sample_definitions))))
      stop("The count table column names are not the same as the sample definition row names.")
    }
    all_count_tables <- data.table::as.data.table(count_dataframe, keep.rownames = "rownames")
    ## If neither of these cases is true, start looking for the files in the
    ## processed_data/ directory
  } else if (is.null(sample_definitions[[file_column]])) {
    stop("This requires either a count dataframe/matrix or column containing the filenames.")
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
                                   tx_gene_map = tx_gene_map, ...)
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

  if (handle_na == "drop") {
    all_count_tables <- all_count_tables[complete.cases(all_count_tables), ]
  } else {
    na_idx <- is.na(all_count_tables)
    all_count_tables[na_idx] <- 0
  }

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
  for (col in seq_along(colnames(gene_info))) {
    sum_na <- sum(is.na(gene_info[[col]]))
    sum_null <- sum(is.null(gene_info[[col]]))
    sum_empty <- sum_na + sum_null
    if (sum_empty ==  nrow(gene_info)) {
      ## This column is empty.
      columns_to_remove <- append(columns_to_remove, col)
    }
    ## While we are looping through the columns,
    ## Make certain that no columns in gene_info are lists or factors.
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
  if (!is.null(arglist[["tx_gene_map"]])) {
    tx_gene_map <- arglist[["tx_gene_map"]]
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
                                        chosen_palette = chosen_palette)

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

  ## Finally, create the ExpressionSet using the counts, annotations, and metadata.
  requireNamespace("Biobase")  ## AnnotatedDataFrame is from Biobase
  metadata <- methods::new("AnnotatedDataFrame", sample_definitions)
  new_samplenames <- try(Biobase::sampleNames(metadata) <- colnames(final_counts))
  if (class(new_samplenames) == "try-error") {
    message("Something is wrong with the sample names for the experimental metadata.")
    message("They must be set to the column names of the count table, which are: ",
            toString(colnames(final_counts)))
    message("If the above column names have .x/.y in them, they may be duplicated,",
            " check your sample sheet.")
    message("The sample definitions are:")
    print(knitr::kable(sample_definitions))
    stop()
  }

  feature_data <- methods::new("AnnotatedDataFrame", final_annotations)
  Biobase::featureNames(feature_data) <- rownames(final_counts)

  ## An expressionset needs to have a Biobase::annotation() in order for
  ## GSEABase to work with it. Reading the documentation, these are primarily
  ## used for naming the type of microarray chip used.
  experiment <- methods::new("ExpressionSet",
                             exprs = final_counts,
                             phenoData = metadata,
                             annotation = annotation_name,
                             featureData = feature_data)

  Biobase::notes(experiment) <- toString(notes)
  ## These entries in new_expt are intended to maintain a record of
  ## the transformation status of the data, thus if we now call
  ## normalize_expt() it should change these.
  ## Therefore, if we call a function like DESeq() which requires
  ## non-log2 counts, we can check these values and convert accordingly

  ## Now that the expressionset has been created, pack it into an expt object so that I
  ## can keep backups etc.
  expt <- sm(subset_expt(experiment)) ## I think this is spurious now.

  if (is.null(study_name)) {
    study_name <- basename(getwd())
  }
  expt[["study"]] <- study_name
  expt[["researcher"]] <- researcher
  expt[["title"]] <- title
  expt[["notes"]] <- toString(notes)
  expt[["design"]] <- sample_definitions
  expt[["annotation"]] <- annotation
  expt[["gff_file"]] <- include_gff
  ## the 'state' slot in the expt is used to keep track of how the data is modified over time.
  starting_state <- list(
      "filter" = "raw",
      "normalization" = "raw",
      "conversion" = "raw",
      "batch" = "raw",
      "transform" = "raw")
  expt[["state"]] <- starting_state
  ## Just in case there are condition names which are not used.
  ## Ditto here, this should not be needed.
  ## expt[["conditions"]] <- droplevels(as.factor(sample_definitions[, "condition"]))
  expt[["conditions"]] <- sample_definitions[["condition"]]
  ## This might be redundant, but it ensures that no all-numeric conditions exist.
  ## expt[["conditions"]] <- gsub(pattern = "^(\\d+)$", replacement = "c\\1", x = expt[["conditions"]])
  names(expt[["conditions"]]) <- rownames(sample_definitions)
  ## Ditto for batches
  expt[["batches"]] <- sample_definitions[["batch"]]
  ## expt[["batches"]] <- droplevels(as.factor(sample_definitions[, "batch"]))
  ## expt[["batches"]] <- gsub(pattern = "^(\\d+)$", replacement = "b\\1", x = expt[["batches"]])
  names(expt[["batches"]]) <- rownames(sample_definitions)
  ## Keep a backup of the metadata in case we do semantic filtering or somesuch.
  ## Keep a backup of the library sizes for limma.
  expt[["libsize"]] <- colSums(Biobase::exprs(experiment))
  if (sum(expt[["libsize"]] == 0) > 0) {
    zero_idx <- expt[["libsize"]] == 0
    zero_samples <- names(expt[["libsize"]])[zero_idx]
    warning("The following samples have no counts! ", zero_samples)
  }
  names(expt[["libsize"]]) <- rownames(sample_definitions)
  ## Save the chosen colors
  expt[["colors"]] <- chosen_colors
  names(expt[["colors"]]) <- rownames(sample_definitions)
  expt[["tximport"]] <- tximport_data
  ## Save an rdata file of the expressionset.
  if (is.null(savefile)) {
    if ("character" %in% class(metadata)) {
      name <- paste0(gsub(x = basename(metadata), pattern = "^(.*)\\..*", replacement = "\\1"), ".rda")
    } else {
      message("Saving the expressionset to 'expt.rda'.")
      savefile <- "expt.rda"
    }
  }
  if (!dir.exists(dirname(savefile))) {
    created <- dir.create(dirname(savefile), recursive = TRUE)
  }
  save_result <- try(save(expt, file = savefile), silent = TRUE)
  if (class(save_result) == "try-error") {
    warning("Saving the expt object failed, perhaps you do not have permissions?")
  }
  message("The final expressionset has ", nrow(exprs(expt)),
          " features and ", ncol(exprs(expt)), " samples.")
  return(expt)
}

#' Exclude some genes given a pattern match
#'
#' Because I am too lazy to remember that expressionsets use matrix subsets for
#' gene and sample.  Also those methods lead to shenanigans when I want to know
#' what happened to the data over the course of the subset.
#'
#' @param expt Expressionset containing expt object.
#' @param column fData column to use for subsetting.
#' @param method Either remove explicit rows, or keep them.
#' @param ids Specific IDs to exclude.
#' @param warning_cutoff Print the sample IDs for anything which has less than this percent left.
#' @param meta_column Save the amount of data lost to this metadata column when not null.
#' @param patterns Character list of patterns to remove/keep
#' @param ... Extra arguments are passed to arglist, currently unused.
#' @return A smaller expt
#' @seealso [create_expt()] [Biobase]
#' @examples
#'  \dontrun{
#'   all_expt <- create_expt(metadata)
#'   ## This assumes a column in the metadata named 'txtype' containing the
#'   ## information telling us what type of transcript each gene is.
#'   no_ribosomes <- exclude_genes_expt(all_expt, column = "txtype",
#'                                      patterns = c("snRNA", "tRNA", "rRNA"))
#'   i_hate_these_genes <- exclude_genes_expt(all_expt, ids = c("gene1", "gene2"))
#'   only_ribosomes <- exclude_genes_expt(all_expt, method = "keep")
#' }
#' @export
exclude_genes_expt <- function(expt, column = "txtype", method = "remove", ids = NULL,
                               warning_cutoff = 90, meta_column = NULL,
                               patterns = c("snRNA", "tRNA", "rRNA"), ...) {
  arglist <- list(...)
  ex <- expt[["expressionset"]]
  annotations <- Biobase::fData(ex)
  expression <- Biobase::exprs(ex)
  if (is.null(ids) & is.null(annotations[[column]])) {
    message("The ", column, " column is null, doing nothing.")
    return(expt)
  }

  pattern_string <- ""
  for (pat in patterns) {
    pattern_string <- glue::glue("{pattern_string}{pat}|")
  }
  silly_string <- gsub(pattern = "\\|$", replacement = "", x = pattern_string)
  idx <- rep(x = TRUE, times = nrow(annotations))
  if (is.null(ids)) {
    idx <- grepl(pattern = silly_string, x = annotations[[column]], perl = TRUE)
  } else if (is.logical(ids)) {
    idx <- ids
  } else {
    idx <- rownames(expression) %in% ids
  }
  kept <- NULL
  removed <- NULL
  kept_sums <- NULL
  removed_sums <- NULL
  if (method == "remove") {
    kept <- ex[!idx, ]
    removed <- ex[idx, ]
    if (!is.null(expt[["tximport"]])) {
      df <- expt[["tximport"]][["raw"]][["abundance"]][idx, ]
      expt[["tximport"]][["raw"]][["abundance"]] <- df
      df <- expt[["tximport"]][["raw"]][["counts"]][idx, ]
      expt[["tximport"]][["raw"]][["counts"]] <- df
      df <- expt[["tximport"]][["raw"]][["length"]][idx, ]
      expt[["tximport"]][["raw"]][["length"]] <- df

      df <- expt[["tximport"]][["scaled"]][["abundance"]][idx, ]
      expt[["tximport"]][["scaled"]][["abundance"]] <- df
      df <- expt[["tximport"]][["scaled"]][["counts"]][idx, ]
      expt[["tximport"]][["scaled"]][["counts"]] <- df
      df <- expt[["tximport"]][["scaled"]][["length"]][idx, ]
      expt[["tximport"]][["scaled"]][["length"]] <- df
    }
  } else {
    kept <- ex[idx, ]
    removed <- ex[!idx, ]
    if (!is.null(expt[["tximport"]])) {
      df <- expt[["tximport"]][["raw"]][["abundance"]][idx, ]
      expt[["tximport"]][["raw"]][["abundance"]] <- df
      df <- expt[["tximport"]][["raw"]][["counts"]][idx, ]
      expt[["tximport"]][["raw"]][["counts"]] <- df
      df <- expt[["tximport"]][["raw"]][["length"]][idx, ]
      expt[["tximport"]][["raw"]][["length"]] <- df

      df <- expt[["tximport"]][["scaled"]][["abundance"]][idx, ]
      expt[["tximport"]][["scaled"]][["abundance"]] <- df
      df <- expt[["tximport"]][["scaled"]][["counts"]][idx, ]
      expt[["tximport"]][["scaled"]][["counts"]] <- df
      df <- expt[["tximport"]][["scaled"]][["length"]][idx, ]
      expt[["tximport"]][["scaled"]][["length"]] <- df
    }
  }

  message("remove_genes_expt(), before removal, there were ",
          nrow(Biobase::fData(ex)), " genes, now there are ",
          nrow(Biobase::fData(kept)), ".")
  all_tables <- Biobase::exprs(ex)
  all_sums <- colSums(all_tables)
  kept_tables <- Biobase::exprs(kept)
  kept_sums <- colSums(kept_tables)
  removed_tables <- Biobase::exprs(removed)
  removed_sums <- colSums(removed_tables)
  pct_kept <- (kept_sums / all_sums) * 100.0
  pct_na <- is.na(pct_kept)
  pct_kept[pct_na] <- 0
  pct_removed <- (removed_sums / all_sums) * 100.0
  pct_na <- is.na(pct_removed)
  pct_removed[pct_na] <- 0
  summary_table <- rbind(kept_sums, removed_sums, all_sums,
                         pct_kept, pct_removed)
  rownames(summary_table) <- c("kept_sums", "removed_sums", "all_sums",
                               "pct_kept", "pct_removed")
  mesg("Percent of the counts kept after filtering: ", toString(sprintf(fmt = "%.3f", pct_kept)))

  expt[["expressionset"]] <- kept
  expt[["summary_table"]] <- summary_table
  if (!is.null(meta_column)) {
    pData(expt[["expressionset"]])[meta_column] <- summary_table["pct_removed", ]
  }

  warning_idx <- summary_table["pct_kept", ] < warning_cutoff
  if (sum(warning_idx) > 0) {
    message("There are ", sum(warning_idx), " samples which kept less than ",
            warning_cutoff, " percent counts.")
    print(summary_table["pct_kept", warning_idx])
  }
  return(expt)
}

#' Do features_greater_than() inverted!
#'
#' @param  ... Arguments passed to features_greather_than()
#' @return The set of features less than whatever you would have done with
#'   features_greater_than().
#' @seealso [features_greater_than()]
#' @export
features_less_than <- function(...) {
  features_greater_than(..., inverse = TRUE)
}

#' Count the number of features(genes) greater than x in a data set.
#'
#' Sometimes I am asked how many genes have >= x counts.  Well, here you go.
#'
#' Untested as of 2016-12-01 but used with Lucia.  I think it would be interesting to iterate
#' this function from small to large cutoffs and plot how the number of kept genes decreases.
#'
#' @param data Dataframe/exprs/matrix/whatever of counts.
#' @param cutoff Minimum number of counts.
#' @param hard Greater-than is hard, greater-than-equals is not.
#' @param inverse when inverted, this provides features less than the cutoff.
#' @return A list of two elements, the first comprised of the number of genes
#'   greater than the cutoff, the second with the identities of said genes.
#' @seealso [Biobase]
#' @examples
#' \dontrun{
#'  features <- features_greater_than(expt)
#'  fewer <- features_greater_than(expt, cutoff = 100)
#' }
#' @export
features_greater_than <- function(data, cutoff = 1, hard = TRUE, inverse = FALSE) {
  if ("expt" %in% class(data) | "ExpressionSet" %in% class(data) |
      "SummarizedExperiment" %in% class(data)) {
    data <- as.data.frame(exprs(data))
  } else {
    data <- as.data.frame(data)
  }
  number_table <- numeric(length = ncol(data))
  names(number_table) <- colnames(data)
  feature_tables <- list()
  for (col in seq_along(colnames(data))) {
    column_name <- colnames(data)[col]
    column_data <- data[[column_name]]
    num_features <- NULL
    if (isTRUE(hard)) {
      if (isTRUE(inverse)) {
        feature_idx <- column_data < cutoff
      } else {
        feature_idx <- column_data > cutoff
      }
    } else {
      if (isTRUE(inverse)) {
        feature_idx <- column_data <= cutoff
      } else {
        feature_idx <- column_data >= cutoff
      }
    }
    num_features <- sum(feature_idx)
    number_table[[column_name]] <- num_features
    passed_filter <- rownames(data)[feature_idx]
    feature_tables[[column_name]] <- passed_filter
  }
  result <- list(
      "number" = number_table,
      "features" = feature_tables)
  return(result)
}

#' I want an easy way to answer the question: what features are in only condition x?
#'
#' The answer to this lies in a combination of subset_expt() and
#' features_greater_than().
#'
#' @param expt An experiment to query.
#' @param cutoff What is the minimum number of counts required to define
#'  'included.'
#' @param factor What metadata factor to query?
#' @param chosen Either choose a subset or all conditions to query.
#' @return A set of features.
#' @seealso [subset_expt()]
#' @examples
#'  \dontrun{
#'   unique_genes
#' }
#' @export
features_in_single_condition <- function(expt, cutoff = 2, factor = "condition", chosen = NULL) {
  condition_set <- levels(as.factor(pData(expt)[[factor]]))
  solo_this_list <- list()
  solo_other_list <- list()
  shared_list <- list()
  neither_list <- list()
  for (cond in condition_set) {
    if (!is.null(chosen)) {
      if (! cond %in% chosen) {
        next
      }
    }
    extract_string <- glue::glue("{factor} == '{cond}'")
    single_expt <- sm(subset_expt(expt = expt, subset = extract_string))
    extract_string <- glue::glue("{factor} != '{cond}'")
    others_expt <- sm(subset_expt(expt = expt, subset = extract_string))
    single_data <- exprs(single_expt)
    others_data <- exprs(others_expt)

    single_positive_idx <- rowSums(single_data) >= cutoff
    single_negative_idx <- rowSums(single_data) < cutoff
    single_positive <- single_data[single_positive_idx, ]
    single_negative <- single_data[single_negative_idx, ]

    others_positive_idx <- rowSums(others_data) >= cutoff
    others_negative_idx <- rowSums(others_data) < cutoff
    others_positive <- others_data[others_positive_idx, ]
    others_negative <- others_data[others_negative_idx, ]

    in_both_sets_idx <- rownames(single_positive) %in% rownames(others_positive)
    in_both_sets <- rownames(single_positive)[in_both_sets_idx]

    only_this_set_idx <- rownames(single_positive) %in% rownames(others_negative)
    only_this_set <- rownames(single_positive)[only_this_set_idx]

    only_other_set_idx <- rownames(others_positive) %in% rownames(single_negative)
    only_other_set <- rownames(others_positive)[only_other_set_idx]

    neither_set_this_idx <- rownames(single_negative) %in% rownames(others_negative)
    neither_set_this <- rownames(single_negative)[neither_set_this_idx]
    neither_set_that_idx <- rownames(others_negative) %in% rownames(single_negative)
    neither_set_that <- rownames(others_negative)[neither_set_that_idx]

    shared_list[[cond]] <- in_both_sets
    solo_this_list[[cond]] <- only_this_set
    solo_other_list[[cond]] <- only_other_set
    neither_list[[cond]] <- neither_set_this
  }
  retlist <- list(
      "shared" = shared_list,
      "solo_this" = solo_this_list,
      "solo_other" = solo_other_list,
      "neither" = neither_list
  )
  return(retlist)
}

#' Set up default colors for a data structure containing usable metadata
#'
#' In theory this function should be useful in any context when one has a blob
#' of metadata and wants to have a set of colors.  Since my taste is utterly
#' terrible, I rely entirely upon RColorBrewer, but also allow one to choose
#' his/her own colors.
#'
#' @param sample_definitions Metadata, presumably containing a 'condition'
#'  column.
#' @param cond_column Which column in the sample data provides the set of
#'  'conditions' used to define the colors?
#' @param by Name the factor of colors according to this column.
#' @param ... Other arguments like a color palette, etc.
#' @return  Colors!
#' @seealso [create_expt()]
generate_expt_colors <- function(sample_definitions, cond_column = "condition",
                                 by = "sampleid", ...) {
  arglist <- list(...)
  ## First figure out how many conditions we have
  colnames(sample_definitions) <- tolower(colnames(sample_definitions))
  ## If there is no condition to start, then it may be NA
  chosen_colors <- as.character(sample_definitions[[cond_column]])
  na_idx <- is.na(chosen_colors)
  chosen_colors[na_idx] <- "undefined"
  num_conditions <- length(levels(as.factor(chosen_colors)))

  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }
  sample_colors <- NULL
  if (!is.null(arglist[["sample_colors"]])) {
    sample_colors <- arglist[["sample_colors"]]
  }
  ## And also the number of samples
  num_samples <- nrow(sample_definitions)
  if (!is.null(sample_colors) & length(sample_colors) == num_samples) {
    ## Thus if we have a numer of colors == the number of samples, set each sample
    ## with its own color
    chosen_colors <- sample_colors
  } else if (!is.null(sample_colors) & length(sample_colors) == num_conditions) {
    ## If instead there are colors == number of conditions, set them appropriately.
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else if (is.null(sample_colors)) {
    ## If nothing is provided, let RColorBrewer do it.
    sample_colors <- sm(
        grDevices::colorRampPalette(
                       RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else {
    ## If none of the above are true, then warn the user and let RColorBrewer do it.
    warning("The number of colors provided does not match the number of conditions nor samples.")
    warning("Unsure of what to do, so choosing colors with RColorBrewer.")
    sample_colors <- sm(
        grDevices::colorRampPalette(
                       RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  }
  ## Set the color names
  names(chosen_colors) <- sample_definitions[[by]]
  return(chosen_colors)
}

#' Small hack of limma's exampleData() to allow for arbitrary data set
#' sizes.
#'
#' exampleData has a set number of genes/samples it creates. This
#' relaxes that restriction.
#'
#' @param ngenes How many genes in the fictional data set?
#' @param columns How many samples in this data set?
#' @return Matrix of pretend counts.
#' @seealso [limma] [DESeq2] [stats]
#' @examples
#' \dontrun{
#'  pretend = make_exampledata()
#' }
#' @export
make_exampledata <- function(ngenes = 1000, columns = 5) {
  q0 <- stats::rexp(ngenes, rate = 1/250)
  is_DE <- stats::runif(ngenes) < 0.3
  lfc <- stats::rnorm(ngenes, sd = 2)
  q0A <- ifelse(is_DE, q0 * 2^ (lfc / 2), q0)
  q0B <- ifelse(is_DE, q0 * 2^ (-lfc / 2), q0)
  ##    true_sf <- c(1, 1.3, 0.7, 0.9, 1.6)
  true_sf <- abs(stats::rnorm(columns, mean = 1, sd = 0.4))
  cond_types <- ceiling(sqrt(columns))
  conds <- sample(LETTERS[1:cond_types], columns, replace = TRUE)
  m <- t(sapply(seq_len(ngenes),
                function(i) sapply(1:columns,
                                   function(j) rnbinom(1,
                                                       mu = true_sf[j] * ifelse(conds[j] == "A",
                                                                                q0A[i], q0B[i]),
                                                       size = 1/0.2))))
  rownames(m) <- paste("gene", seq_len(ngenes), ifelse(is_DE, "T", "F"), sep = "_")
  colnames(m) <- paste0("sample", 1:ncol(m))
  design <- as.data.frame(conds)
  rownames(design) <- colnames(m)
  colnames(design) <- "condition"
  for (i in seq_len(ncol(m))) {
    m[, i] <- as.integer(m[, i])
  }
  design[["condition"]] <- as.factor(design[["condition"]])
  example <- sm(DESeq2::DESeqDataSetFromMatrix(countData = m,
                                               colData = design, design = ~condition))
  return(example)
}

#' Create a data frame of the medians of rows by a given factor in the data.
#'
#' This assumes of course that (like expressionsets) there are separate columns
#' for each replicate of the conditions.  This will just iterate through the
#' levels of a factor describing the columns, extract them, calculate the
#' median, and add that as a new column in a separate data frame.
#'
#' Used in write_expt() as well as a few random collaborations.
#'
#' @param data Data frame, presumably of counts.
#' @param fact Factor describing the columns in the data.
#' @param fun Optionally choose mean or another function.
#' @return Data frame of the medians.
#' @seealso [Biobase] [matrixStats]
#' @examples
#' \dontrun{
#'  compressed = median_by_factor(data, experiment$condition)
#' }
#' @export
median_by_factor <- function(data, fact = "condition", fun = "median") {
  if (length(fact) == 1) {
    design <- pData(data)
    fact <- design[[fact]]
    na_idx <- is.na(fact)
    if (sum(na_idx) > 0) {
      fact[na_idx] <- "unknown"
    }
    names(fact) <- rownames(design)
  }
  if ("expt" %in% class(data) | "ExpressionSet" %in% class(data) |
      "SummarizedExperiment" %in% class(data)) {
    data <- exprs(data)
  }

  medians <- data.frame("ID" = rownames(data), stringsAsFactors = FALSE)
  cvs <- data.frame("ID" = rownames(data), stringsAsFactors = FALSE)
  data <- as.matrix(data)
  rownames(medians) <- rownames(data)
  rownames(cvs) <- rownames(data)
  fact <- as.factor(fact)
  used_columns <- c()
  group_indexes <- list()
  for (type in levels(fact)) {
    ## columns <- grep(pattern = type, x = fact)
    columns <- as.character(fact) == type
    group_indexes[[type]] <- columns
    med <- NULL
    if (sum(columns) < 1) {
      warning("The level ", type, " of the factor has no columns.")
      next
    }
    used_columns <- c(used_columns, type)
    if (sum(columns) == 1) {
      message("The factor ", type, " has only 1 row.")
      med <- as.data.frame(data[, columns], stringsAsFactors = FALSE)
      cv <- as.data.frame(data[, columns], stringsAsFactors = FALSE)
    } else {
      if (fun == "median") {
        message("The factor ", type, " has ", sum(columns), " rows.")
        med <- matrixStats::rowMedians(data[, columns], na.rm = TRUE)
        cv <- matrixStats::rowMads(data[, columns], na.rm = TRUE)
        cv <- cv / med
        ## I am not really sure if this is appropriate.
        nan_idx <- is.nan(cv)
        cv[nan_idx] <- 0
      } else if (fun == "mean") {
        message("The factor ", type, " has ", sum(columns), " rows.")
        med <- BiocGenerics::rowMeans(data[, columns], na.rm = TRUE)
        cv <- matrixStats::rowSds(data[, columns], na.rm = TRUE)
        cv <- cv / med
        nan_idx <- is.nan(cv)
        cv[nan_idx] <- 0
      } else {
        stop("I do not understand that function.")
      }
    }
    medians <- cbind(medians, med)
    cvs <- cbind(cvs, cv)
  }
  medians <- medians[, -1, drop = FALSE]
  cvs <- cvs[, -1, drop = FALSE]
  ## Sometimes not all levels of the original experimental design are used.
  ## Thus lets make sure to use only those which appeared.
  colnames(medians) <- used_columns
  colnames(cvs) <- used_columns
  retlist <- list(
      "method" = fun,
      "medians" = medians,
      "cvs" = cvs,
      "indexes" = group_indexes)
  return(retlist)
}

#' Create a Schizosaccharomyces cerevisiae expt.
#'
#' This just saves some annoying typing if one wishes to make a standard
#' expressionset superclass out of the publicly available fission data set.
#'
#' @param annotation Add annotation data?
#' @return Expressionset/expt of fission.
#' @seealso [fission] [create_expt()]
#' @export
make_pombe_expt <- function(annotation = TRUE) {
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
    pombe_annotations <- try(load_biomart_annotations(
        host = "fungi.ensembl.org", trymart = "fungi_mart",
        trydataset = "spombe_eg_gene",
        gene_requests = c("pombase_transcript", "ensembl_gene_id", "ensembl_transcript_id",
                          "hgnc_symbol", "description", "gene_biotype"),
        species = "spombe", overwrite = TRUE))
    if ("try-error" %in% class(pombe_annotations)) {
      warning("There was an error downloading the pombe annotations, this will still return.")
    } else {
      pombe_mart <- pombe_annotations[["mart"]]
      annotations <- pombe_annotations[["annotation"]]
      rownames(annotations) <- make.names(gsub(pattern = "\\.\\d+$",
                                               replacement = "",
                                               x = rownames(annotations)), unique = TRUE)
    }
  }
  pombe_expt <- sm(create_expt(metadata = meta,
                               count_dataframe = fission_data,
                               gene_info = annotations))
  detach("package:fission")
  return(pombe_expt)
}

#' Read a bunch of count tables and create a usable data frame from them.
#'
#' It is worth noting that this function has some logic intended for the elsayed
#' lab's data storage structure. It shouldn't interfere with other usages, but
#' it attempts to take into account different ways the data might be stored.
#'
#' Used primarily in create_expt()
#' This is responsible for reading count tables given a list of filenames.  It
#' tries to take into account upper/lowercase filenames and uses data.table to
#' speed things along.
#'
#' @param ids List of experimental ids.
#' @param files List of files to read.
#' @param header Whether or not the count tables include a header row.
#' @param include_summary_rows Whether HTSeq summary rows should be included.
#' @param all.x When merging (as opposed to join), choose the x data column.
#' @param all.y When merging (as opposed to join), choose the y data column.
#' @param merge_type Choose one, merge or join.
#' @param suffix Optional suffix to add to the filenames when reading them.
#' @param countdir Optional count directory to read from.
#' @param tx_gene_map Dataframe which provides a mapping between transcript IDs and gene IDs.
#' @param ... More options for happy time!
#' @return Data frame of count tables.
#' @seealso [data.table] [create_expt()] [tximport]
#' @examples
#' \dontrun{
#'  count_tables <- hpgl_read_files(as.character(sample_ids), as.character(count_filenames))
#' }
#' @export
read_counts_expt <- function(ids, files, header = FALSE, include_summary_rows = FALSE,
                             all.x = TRUE, all.y = FALSE, merge_type = "merge",
                             suffix = NULL, countdir = NULL, tx_gene_map = NULL, ...) {
  ## load first sample
  arglist <- list(...)
  retlist <- list()
  ## Add an optional directory if I don't feel like specifying in the sample sheet.
  if (!is.null(countdir)) {
    files <- file.path(countdir, files)
  }
  skippers <- (files == "" | files == "undef" | is.null(files) | is.na(files))
  files <- files[!skippers]
  ids <- ids[!skippers]
  skippers <- (ids == "" | ids == "undef" | is.null(ids) | is.na(ids))
  files <- files[!skippers]
  ids <- ids [!skippers]
  retlist[["kept_ids"]] <- ids
  retlist[["kept_files"]] <- files
  ## lower_filenames <- files
  dirs <- dirname(files)

  count_table <- NULL
  for (f in seq_along(files)) {
    ## Get rid of any lurking spaces
    files[f] <- gsub(pattern = " ", replacement = "", x = files[f])
    ## Check if these are relative or absolute paths and standardize them to absolute.
    if (!grepl(pattern = "^\\/", x = files[f])) {
      files[f] <- file.path(getwd(), files[f])
    }
  }

  ## When I start using sailfish/salmon/etc, I will need to add more conditions to this
  ## and probably change it to a switch to more prettily take them into account.

  ## Likely important options:
  ## txOut: When true, do not back-convert to gene-level.
  ## Otherwise, tximport requires a tx2gene data frame with 2 columns, TXNAME and GENEID.
  ## These columns should be definition be available in my fData annotation.
  ## Therefore, I will set the flags tx2gene and txOut accordingly.
  mesg("Reading count tables.")
  txout <- TRUE
  if (!is.null(tx_gene_map)) {
    mesg("Using the transcript to gene mapping.")
    message("In some cases, (notably salmon) the format of the IDs used by this can be tricky.
It is likely to require the transcript ID followed by a '.' and the ensembl column:
'transcript_version', which is explicitly different than the gene version column.
If this is not correctly performed, very few genes will be observed")
    txout <- FALSE
  }
  if (grepl(pattern = "\\.tsv|\\.h5", x = files[1])) {
    mesg("Reading kallisto data with tximport.")
    ## This hits if we are using the kallisto outputs.
    names(files) <- ids
    if (!all(file.exists(files))) {
      warning(files)
    }
    import <- NULL
    import_scaled <- NULL
    if (is.null(tx_gene_map)) {
      import <- sm(tximport::tximport(files = files, type = "kallisto", txOut = txout))
      import_scaled <- sm(tximport::tximport(
                                        files = files, type = "kallisto",
                                        txOut = txout, countsFromAbundance = "lengthScaledTPM"))
    } else {
      import <- sm(tximport::tximport(
                                 files = files, type = "kallisto", tx2gene = tx_gene_map, txOut = txout))
      import_scaled <- sm(tximport::tximport(
                                        files = files, type = "kallisto", tx2gene = tx_gene_map,
                                        txOut = txout, countsFromAbundance = "lengthScaledTPM"))
    }
    retlist[["count_table"]] <- data.table::as.data.table(
                                                import[["counts"]], keep.rownames = "rownames")
    retlist[["count_table"]] <- data.table::setkey(retlist[["count_table"]], rownames)
    retlist[["tximport"]] <- import
    retlist[["tximport_scaled"]] <- import_scaled
    retlist[["source"]] <- "tximport"
  } else if (grepl(pattern = "\\.genes\\.results", x = files[1])) {
    mesg("Reading rsem data with tximport.")
    names(files) <- ids
    import <- NULL
    import_scaled <- NULL
    if (is.null(tx_gene_map)) {
      import <- tximport::tximport(files = files, type = "rsem", txOut = txout)
      import_scaled <- tximport::tximport(files = files, type = "rsem",
                                          txOut = txout, countsFromAbundance = "lengthScaledTPM")
    } else {
      import <- tximport::tximport(files = files, type = "rsem", tx2gene = tx_gene_map, txOut = txout)
      import_scaled <- tximport::tximport(files = files, type = "rsem", tx2gene = tx_gene_map,
                                          txOut = txout, countsFromAbundance = "lengthScaledTPM")
    }
    retlist[["count_table"]] <- data.table::as.data.table(import[["counts"]],
                                                          keep.rownames = "rownames")
    retlist[["tximport"]] <- import
    retlist[["tximport_scaled"]] <- import_scaled
    retlist[["source"]] <- "tximport"
  } else if (grepl(pattern = "\\.sf", x = files[1])) {
    mesg("Reading salmon data with tximport.")
    ## This hits if we are using the salmon outputs.
    names(files) <- ids
    missing_idx <-!file.exists(files)
    if (sum(missing_idx) > 0) {
      missing_files <- files[missing_idx]
      warning("Not all the files exist: ", toString(missing_files))
    }
    import <- NULL
    import_scaled <- NULL
    if (is.null(tx_gene_map)) {
      import <- sm(tximport::tximport(files = files, type = "salmon", txOut = txout))
      import_scaled <- sm(tximport::tximport(
                                        files = files, type = "salmon",
                                        txOut = txout, countsFromAbundance = "lengthScaledTPM"))
    } else {
      import <- sm(tximport::tximport(
                                 files = files, type = "salmon", tx2gene = tx_gene_map, txOut = txout))
      import_scaled <- sm(tximport::tximport(files = files, type = "salmon", tx2gene = tx_gene_map,
                                             txOut = txout, countsFromAbundance = "lengthScaledTPM"))
    }
    retlist[["count_table"]] <- data.table::as.data.table(
                                                import[["counts"]], keep.rownames = "rownames")
    retlist[["count_table"]] <- data.table::setkey(retlist[["count_table"]], rownames)
    retlist[["tximport"]] <- import
    retlist[["tximport_scaled"]] <- import_scaled
    retlist[["source"]] <- "tximport"
  } else {
    mesg("Reading count files with read.table().")
    ## Use this codepath when we are working with htseq
    count_table <- read.table(files[1], header = header, stringsAsFactors = FALSE)
    colnames(count_table) <- c("rownames", ids[1])
    count_table[, 2] <- as.numeric(count_table[, 2])
    na_idx <- is.na(count_table[[2]])
    ## This is a bit more circuituous than I would like.
    ## I want to make sure that na_rownames does not evaluate to something
    ## annoying like 'logical(0)' which it will if there are no nas in the count
    ## table and you just ask for is.na(count_table).
    na_rownames <- ""
    if (sum(na_idx) > 0) {
      na_rownames <- count_table[na_idx, "rownames"]
    }
    keepers_idx <- count_table[["rownames"]] != na_rownames
    count_table <- count_table[keepers_idx, ]
    count_table <- data.table::as.data.table(count_table)
    count_table <- data.table::setkey(count_table, rownames)
    first_rownames <- sort(count_table[["rownames"]])
    if (class(count_table)[1] == "try-error") {
      stop("There was an error reading: ", files[1])
    }
    mesg(files[1], " contains ", length(rownames(count_table)), " rows.")
    ## iterate over and append remaining samples
    for (num in seq(from = 2, to = length(files))) {
      table <- files[num]
      if (file.exists(tolower(table))) {
        table <- tolower(table)
      }
      tmp_count <- try(read.table(table, header = header))
      ## Drop the rows with NAs before coercing to numeric.
      keepers_idx <- tmp_count[[1]] != na_rownames
      tmp_count <- tmp_count[keepers_idx, ]
      current_rownames <- sort(tmp_count[[1]])
      mismatched_rownames <- sum(first_rownames != current_rownames)
      if (mismatched_rownames > 0) {
        warning("The file: ", table, " has mismatched rownames.")
      }
      tmp_count[, 2] <- as.numeric(tmp_count[, 2])
      if (class(tmp_count)[1] == "try-error") {
        stop("There was an error reading: ", table)
      }
      colnames(tmp_count) <- c("rownames", ids[num])
      tmp_count <- data.table::as.data.table(tmp_count)
      pre_merge <- nrow(tmp_count)
      if (merge_type == "merge") {
        count_table <- merge(count_table, tmp_count, by = "rownames", all.x = all.x, all.y = all.y)
      } else {
        count_table <- plyr::join(count_table, tmp_count)
      }
      ## rownames(count_table) <- count_table[, "Row.names"]
      ## count_table <- count_table[, -1, drop = FALSE]
      ## post_merge <- length(rownames(count_table))
      post_merge <- nrow(count_table)
      mesg(table, " contains ", pre_merge, " rows and merges to ", post_merge, " rows.")
    } ## End for loop

    ## remove summary fields added by HTSeq
    if (!isTRUE(include_summary_rows)) {
      ## Depending on what happens when the data is read in, these rows may get prefixed with 'X'
      ## In theory, only 1 of these two cases should ever be true.
      htseq_meta_rows <- c("__no_feature", "__ambiguous", "__too_low_aQual",
                           "__not_aligned", "__alignment_not_unique",
                           "X__no_feature", "X__ambiguous", "X__too_low_aQual",
                           "X__not_aligned", "X__alignment_not_unique")
      kept_rows <- !count_table[["rownames"]] %in% htseq_meta_rows
      count_table <- count_table[kept_rows, ]
      retlist[["count_table"]] <- count_table
      retlist[["source"]] <- "htseq"
    } ## End the difference between tximport and reading tables.
    retlist[["count_table"]] <- data.table::setkey(retlist[["count_table"]], rownames)
  }
  mesg("Finished reading count data.")
  return(retlist)
}

#' Given a table of meta data, read it in for use by create_expt().
#'
#' Reads an experimental design in a few different formats in preparation for
#' creating an expt.
#'
#' @param file Csv/xls file to read.
#' @param ... Arguments for arglist, used by sep, header and similar
#'  read_csv/read.table parameters.
#' @return Df of metadata.
#' @seealso [openxlsx] [readODS]
#' @export
read_metadata <- function(file, ...) {
  arglist <- list(...)
  if (is.null(arglist[["sep"]])) {
    arglist[["sep"]] <- ","
  }
  if (is.null(arglist[["header"]])) {
    arglist[["header"]] <- TRUE
  }
  extension <- tools::file_ext(file)
  if (extension == "csv") {
    definitions <- read.csv(file = file, comment.char = "#",
                            sep = arglist[["sep"]], header = arglist[["header"]])
  } else if (extension == "tsv") {
    definitions <- try(readr::read_tsv(file, ...))
  } else if (extension == "xlsx") {
    ## xls = loadWorkbook(file, create = FALSE)
    ## tmp_definitions = readWorksheet(xls, 1)
    definitions <- try(openxlsx::read.xlsx(xlsxFile = file, sheet = 1))
    if (class(definitions)[1] == "try-error") {
      stop("Unable to read the metadata file: ", file)
    }
  } else if (extension == "xls") {
    ## This is not correct, but it is a start
    definitions <- readxl::read_excel(path = file, sheet = 1)
  } else if (extension == "ods") {
    sheet <- 1
    if (!is.null(arglist[["sheet"]])) {
      sheet <- arglist[["sheet"]]
    }
    definitions <- readODS::read_ods(path = file, sheet = sheet)
  } else {
    definitions <- read.table(file = file, sep = arglist[["sep"]], header = arglist[["header"]])
  }
  colnames(definitions) <- tolower(gsub(pattern = "[[:punct:]]",
                                        replacement = "",
                                        x = colnames(definitions)))
  ## I recently received a sample sheet with a blank sample ID column name...
  empty_idx <- colnames(definitions) == ""
  colnames(definitions)[empty_idx] <- "empty"
  colnames(definitions) <- make.names(colnames(definitions), unique = TRUE)
  return(definitions)
}

#' Remove/keep specifically named genes from an expt.
#'
#' I find subsetting weirdly confusing.  Hopefully this function will allow one
#' to include/exclude specific genes/families based on string comparisons.
#'
#' @param input Expt to filter.
#' @param invert The default is to remove the genes with the semantic strings.
#'  Keep them when inverted.
#' @param topn Take the topn most abundant genes rather than a text based heuristic.
#' @param semantic Character list of strings to search for in the annotation
#'  data.
#' @param semantic_column Column in the annotations to search.
#' @return A presumably smaller expt.
#' @seealso [Biobase]
#' @export
semantic_expt_filter <- function(input, invert = FALSE, topn = NULL,
                                 semantic = c("mucin", "sialidase", "RHS", "MASP", "DGF", "GP63"),
                                 semantic_column = "description") {
  mtrx <- exprs(input)
  annots <- fData(input)
  if (isTRUE(invert)) {
    new_annots <- data.frame()
    new_mtrx <- data.frame()
  } else {
    new_annots <- annots
    new_mtrx <- mtrx
  }
  start_rows <- nrow(mtrx)
  numbers_removed <- 0
  if (is.null(topn)) {
    for (string in semantic) {
      pre_remove_size <- nrow(new_annots)
      idx <- NULL
      if (isTRUE(invert)) {
        ## Keep the rows which match the ~7 strings above.
        ## For these, we will re-grep the full table each time and just add the matches.
        type <- "Kept"
        if (semantic_column == "rownames") {
          idx <- grepl(pattern = string, x = rownames(annots))
        } else {
          idx <- grepl(pattern = string, x = annots[, semantic_column])
        }
        message("Hit ", sum(idx), " genes for term ", string, ".")
        ## Then, after grepping, just append the matched rows to the new annotations and matrix.
        tmp_annots <- annots[idx, ]
        tmp_mtrx <- mtrx[idx, ]
        new_annots <- rbind(new_annots, tmp_annots)
        new_mtrx <- rbind(new_mtrx, tmp_mtrx)
      } else {
        type <- "Removed"
        ## In the case of removals, I need to only grep what is left after each iteration.
        if (semantic_column == "rownames") {
          idx <- grepl(pattern = string, x = rownames(new_annots))
        } else {
          idx <- grepl(pattern = string, x = new_annots[, semantic_column])
        }
        mesg("Hit ", sum(idx), " genes for term ", string, ".")
        idx <- ! idx
        ## So, we take the index of stuff to keep, and just subset on that index.
        new_annots <- new_annots[idx, ]
        new_mtrx <- new_mtrx[idx, ]
      }
    } ## End for loop
    end_rows <- nrow(new_mtrx)
    lost_rows <- start_rows - end_rows
    message("semantic_expt_filter(): Removed ", lost_rows, " genes.")
  } else {
    ## Instead of a string based sematic filter, take the topn most abundant
    medians <- rowMedians(mtrx)
    new_order <- order(medians, decreasing = TRUE)
    reordered <- mtrx[new_order, ]
    subset <- rownames(head(reordered, n = topn))
    new_annots <- annots[subset, ]
  }

  expressionset <- input[["expressionset"]]
  keepers <- rownames(new_annots)
  new_expressionset <- expressionset[keepers, ]
  new_libsizes <- colSums(exprs(new_expressionset))
  input[["expressionset"]] <- new_expressionset
  input[["libsize"]] <- new_libsizes
  return(input)
}

#' Change the batches of an expt.
#'
#' When exploring differential analyses, it might be useful to play with the
#' conditions/batches of the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify.
#' @param fact Batches to replace using this factor.
#' @param ids Specific samples to change.
#' @param ... Extra options are like spinach.
#' @return The original expt with some new metadata.
#' @seealso [create_expt()] [set_expt_conditions()] [Biobase]
#' @examples
#' \dontrun{
#'  expt = set_expt_batches(big_expt, factor = c(some,stuff,here))
#' }
#' @export
set_expt_batches <- function(expt, fact, ids = NULL, ...) {
  arglist <- list(...)
  original_batches <- pData(expt)[["batch"]]
  original_length <- length(original_batches)
  if (length(fact) == 1) {
    ## Assume it is a column in the design
    if (fact %in% colnames(pData(expt))) {
      fact <- pData(expt)[[fact]]
    } else {
      stop("The provided factor is not in the design matrix.")
    }
  }

  if (length(fact) != original_length) {
    stop("The new factor of batches is not the same length as the original.")
  }
  expt[["batches"]] <- fact
  pData(expt[["expressionset"]])[["batch"]] <- fact
  expt[["design"]][["batch"]] <- fact
  return(expt)
}


#' Change the colors of an expt
#'
#' When exploring differential analyses, it might be useful to play with the
#' conditions/batches of the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify
#' @param colors colors to replace
#' @param chosen_palette I usually use Dark2 as the RColorBrewer palette.
#' @param change_by Assuming a list is passed, cross reference by condition or sample?
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_conditions()] [set_expt_batches()] [RColorBrewer]
#' @examples
#' \dontrun{
#' unique(esmer_expt$design$conditions)
#' chosen_colors <- list(
#'    "cl14_epi" = "#FF8D59",
#'    "clbr_epi" = "#962F00",
#'    "cl14_tryp" = "#D06D7F",
#'    "clbr_tryp" = "#A4011F",
#'    "cl14_late" = "#6BD35E",
#'    "clbr_late" = "#1E7712",
#'    "cl14_mid" = "#7280FF",
#'    "clbr_mid" = "#000D7E")
#' esmer_expt <- set_expt_colors(expt = esmer_expt, colors = chosen_colors)
#' }
#' @export
set_expt_colors <- function(expt, colors = TRUE, chosen_palette = "Dark2", change_by = "condition") {
  condition_factor <- as.factor(pData(expt)[["condition"]])

  ## Since I already have logic for named characters, just convert a list to one...
  if ("list" %in% class(colors)) {
    new_colors <- as.character(colors)
    names(new_colors) <- names(colors)
    colors <- new_colors
  }

  num_conditions <- length(levels(condition_factor))
  num_samples <- nrow(expt[["design"]])
  sample_ids <- expt[["design"]][["sampleid"]]
  ## chosen_colors <- expt[["conditions"]]
  chosen_colors <- condition_factor
  chosen_names <- names(chosen_colors)
  sample_colors <- NULL
  if (is.null(colors) | isTRUE(colors)) {
    sample_colors <- sm(
        grDevices::colorRampPalette(
                       RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else if (class(colors) == "factor") {
    if (change_by == "condition") {
      mesg("The new colors are a factor, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      colors_allocated <- names(colors) %in% levels(pData(expt)[["condition"]])
      if (sum(colors_allocated) < length(colors)) {
        missing_colors <- colors[!colors_allocated]
        stop("Colors for the following categories are not being used: ",
                names(missing_colors), ".")
      }
      possible_conditions <- levels(pData(expt)[["condition"]])
      conditions_allocated <- possible_conditions %in% names(colors)
      if (sum(conditions_allocated) < length(possible_conditions)) {
        missing_conditions <- possible_conditions[!conditions_allocated]
        missing_samples <- c()
        for (cond in missing_conditions) {
          missing_by_condition <- pData(expt)[["condition"]] == cond
          missing_samples_by_cond <- rownames(pData(expt))[missing_by_condition]
          missing_samples <- c(missing_samples, missing_samples_by_cond)
        }
        warning("Some conditions do not have a color: ", missing_conditions, ".")
        warning("These samples are: ", missing_samples, ".")
      }
      mapping <- colors
      chosen_colors <- mapping[as.character(chosen_colors)]
      names(chosen_colors) <- chosen_names
    } else if (change_by == "sample") {
      mesg("The new colors are a factor, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- expt[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (class(colors) == "character") {
    if (is.null(names(colors))) {
      names(colors) <- levels(as.factor(expt[["conditions"]]))
    }
    if (change_by == "condition") {
      mesg("The new colors are a character, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- colors
      pd_factor <- as.factor(pData(expt)[["condition"]])
      possible_conditions <- levels(pd_factor)
      colors_allocated <- names(colors) %in% possible_conditions
      if (sum(colors_allocated) < length(colors)) {
        missing_colors <- colors[!colors_allocated]
        warning("Colors for the following categories are not being used: ",
                names(missing_colors), ".")
      }
      conditions_allocated <- possible_conditions %in% names(colors)
      if (sum(conditions_allocated) < length(possible_conditions)) {
        missing_conditions <- possible_conditions[!conditions_allocated]
        missing_samples <- c()
        for (cond in missing_conditions) {
          missing_by_condition <- pData(expt)[["condition"]] == cond
          missing_samples_by_cond <- rownames(pData(expt))[missing_by_condition]
          missing_samples <- c(missing_samples, missing_samples_by_cond)
        }
        warning("Some conditions do not have a color: ", missing_conditions, ".")
        warning("These samples are: ", missing_samples, ".")
      }
      chosen_colors <- mapping[as.character(chosen_colors)]
      names(chosen_colors) <- chosen_names
    } else if (change_by == "sample") {
      mesg("The new colors are a character, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- expt[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (class(colors) == "list") {
    if (change_by == "condition") {
      mesg("The new colors are a list, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- as.character(colors)
      names(mapping) <- names(colors)
      chosen_colors <- mapping[chosen_colors]
    } else if (change_by == "sample") {
      mesg("The new colors are a list, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- expt[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
        ## Set the condition for the changed samples to something unique.
        original_condition <- expt[["design"]][sampleid, "condition"]
        changed_condition <- glue::glue("{original_condition}{snum}")
        expt[["design"]][sampleid, "condition"] <- changed_condition
        tmp_pdata <- pData(expt)
        old_levels <- levels(tmp_pdata[["condition"]])
        new_levels <- c(old_levels, changed_condition)
        levels(tmp_pdata[["condition"]]) <- new_levels
        tmp_pdata[sampleid, "condition"] <- changed_condition
        pData(expt[["expressionset"]]) <- tmp_pdata
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (is.null(colors)) {
    mesg("Setting colors according to a color ramp.")
    colors <- sm(
        grDevices::colorRampPalette(
                       RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    ## Check that all conditions are named in the color list:
    mapping <- setNames(colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else {
    warning("Number of colors provided does not match the number of conditions nor samples.")
    warning("Unsure of what to do, so choosing colors with RColorBrewer.")
    sample_colors <- suppressWarnings(
        grDevices::colorRampPalette(
                       RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  }

  ## Catchall in case I forgot to set the names before now.
  names(chosen_colors) <- chosen_names
  expt[["colors"]] <- chosen_colors
  return(expt)
}

#' Change the condition of an expt
#'
#' When exploring differential analyses, it might be useful to play with the
#' conditions/batches of the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify
#' @param fact Conditions to replace
#' @param ids Specific sample IDs to change.
#' @param null_cell How to fill elements of the design which are null?
#' @param ... Extra arguments are given to arglist.
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_batches()] [create_expt()]
#' @examples
#' \dontrun{
#'  expt = set_expt_conditions(big_expt, factor = c(some,stuff,here))
#' }
#' @export
set_expt_conditions <- function(expt, fact = NULL, ids = NULL, null_cell = "null", ...) {
  arglist <- list(...)
  original_conditions <- pData(expt)[["condition"]]
  original_length <- length(original_conditions)
  original_num_conditions <- length(levels(as.factor(original_conditions)))
  new_expt <- expt  ## Explicitly copying expt to new_expt
  ## because when I run this as a function call() it seems to be not properly setting
  ## the conditions and I do not know why.
  if (!is.null(ids)) {
    ## Change specific id(s) to given condition(s).
    mesg("Setting condition for ids ", toString(ids), " to ", fact, ".")
    old_pdata <- pData(expt)
    old_cond <- as.character(old_pdata[["condition"]])
    names(old_cond) <- rownames(old_pdata)
    new_cond <- old_cond
    new_cond[ids] <- fact
    new_pdata <- old_pdata
    new_pdata[["condition"]] <- as.factor(new_cond)
    pData(new_expt[["expressionset"]]) <- new_pdata
    new_conditions <- as.character(new_expt[["conditions"]])
    names(new_conditions) <- names(new_expt[["conditions"]])
    new_conditions[ids] <- fact
    new_expt[["conditions"]] <- as.factor(new_conditions)
    new_expt[["design"]][["condition"]] <- new_cond
  } else if (length(fact) == 1) {
    ## Assume it is a column in the design
    if (fact %in% colnames(pData(expt))) {
      new_fact <- pData(expt)[[fact]]
      null_ids <- is.na(new_fact) | is.null(new_fact)
      ## Only do this if there are some null entries.
      if (sum(null_ids) > 0) {
        new_fact[null_ids] <- null_cell
      }
      new_expt[["conditions"]] <- new_fact
      pData(new_expt[["expressionset"]])[["condition"]] <- new_fact
      new_expt[["design"]][["condition"]] <- new_fact
    } else {
      stop("The provided factor is not in the design matrix.")
    }
  } else if (length(fact) != original_length) {
    stop("The new factor of conditions is not the same length as the original.")
  } else {
    new_expt[["conditions"]] <- fact
    pData(new_expt[["expressionset"]])[["condition"]] <- fact
    new_expt[["design"]][["condition"]] <- fact
  }

  new_expt <- set_expt_colors(new_expt)
  return(new_expt)
}

#' Change the factors (condition and batch) of an expt
#'
#' When exploring differential analyses, it might be useful to play with the
#' conditions/batches of the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify
#' @param condition New condition factor
#' @param batch New batch factor
#' @param ids Specific sample IDs to change.
#' @param table When set to 'metadata', use pData, otherwise fData.
#' @param class Set the data to this class by default.
#' @param columns Change these columns.
#' @param ... Arguments passed along (likely colors)
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_conditions()] [set_expt_batches()]
#' @examples
#' \dontrun{
#'  expt = set_expt_factors(big_expt, condition = "column", batch = "another_column")
#' }
#' @export
set_expt_factors <- function(expt, condition = NULL, batch = NULL, ids = NULL,
                             table = "metadata", class = "factor", columns = NULL, ...) {
  arglist <- list(...)
  if (!is.null(condition)) {
    expt <- set_expt_conditions(expt, fact = condition, ...)
  }
  if (!is.null(batch)) {
    expt <- set_expt_batches(expt, fact = batch, ...)
  }
  fd <- fData(expt)
  pd <- pData(expt)
  if (!is.null(columns)) {
    if (is.null(class)) {
      stop("If columns is set, then this assumes you want to set those columns to a given class.")
    }
    meta_columns <- colnames(pd)
    for (col in columns) {
      if (! col %in% meta_columns) {
        warning("The column: ", col, " is not in the metadata, skipping it.")
        next
      }
      mesg("Setting ", col, " to type ", class, ".")
      if (class == "factor") {
        if (table == "metadata") {
          pd[[col]] <- as.factor(pd[[col]])
        } else {
          fd[[col]] <- as.factor(fd[[col]])
        }
      } else if (class == "character") {
        if (table == "metadata") {
          pd[[col]] <- as.character(pd[[col]])
        } else {
          fd[[col]] <- as.character(fd[[col]])
        }
      } else if (class == "numeric") {
        if (table == "metadata") {
          pd[[col]] <- as.numeric(pd[[col]])
        } else {
          fd[[col]] <- as.numeric(fd[[col]])
        }
      } else {
        stop("I do not know this class.")
      }
    }
  }
  pData(expt[["expressionset"]]) <- pd
  fData(expt[["expressionset"]]) <- fd
  return(expt)
}

#' Change the gene names of an expt.
#'
#' I want to change all the gene names of a big expressionset to the
#' ortholog groups.  But I want to also continue using my expts.
#' Ergo this little function.
#'
#' @param expt Expt to modify
#' @param ids Specific sample IDs to change.
#' @param ... Extra arguments are given to arglist.
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_conditions()] [create_expt()]
#' @examples
#' \dontrun{
#'  expt = set_expt_conditions(big_expt, factor = c(some,stuff,here))
#' }
#' @export
set_expt_genenames <- function(expt, ids = NULL, ...) {
  arglist <- list(...)
  expr <- expt[["expressionset"]]

  ## Make sure the order of the IDs stays consistent.
  current_ids <- rownames(exprs(expr))
  if (class(ids) == "data.frame") {
    ## our_column contains the IDs from my species.
    our_column <- NULL
    ## their_column contains the IDs from the species we want to map against.
    their_column <- NULL
    ## Grab the first ID in the first column.
    ## We will explicitly assume there are 2 columns in the data frame.
    test_first <- sum(ids[[1]] %in% current_ids)
    test_second <- sum(ids[[2]] %in% current_ids)
    if (test_first > 0) {
      mesg("Found: ", test_first, " ids in common using the first column of the IDs.")
      our_column <- colnames(ids)[1]
      their_column <- colnames(ids)[2]
    } else if (test_second > 0) {
      mesg("Found: ", test_second, " ids in common using the second column of the IDs.")
      our_column <- colnames(ids)[2]
      their_column <- colnames(ids)[1]
    } else {
      stop("Unable to match the IDs.")
    }
    mesg("Our column is: ", our_column, ", their column is: ", their_column, ".")
    ## Now the job is to ensure that the ordering is maintained.
    ## We need therefore to merge the rownames of the current IDs into the
    ## data frame of the ID mapping between our species.
    ## In addition, we must keep all of the IDs from our species (all.x = TRUE).
    exprs_id_df <- as.data.frame(current_ids)
    reordered <- merge(exprs_id_df, ids, by.x = "current_ids", by.y = our_column, all.x = TRUE)
    ## This merge should give a NA in the other column if there is no gene in the
    ## other species for a gene in our species.  In addition, it should give us
    ## duplicate IDs where a single gene from our species maps against multiple
    ## genes from the other species.
    ## Take care of the first scenario:
    na_ids <- is.na(reordered[[their_column]])
    reordered[na_ids, ] <- reordered[na_ids, "current_ids"]
    ## Then get rid of the duplicates.
    dup_ids <- duplicated(reordered[["current_ids"]])
    reordered <- reordered[!dup_ids, ]
    ## Now we should have a data frame with the same number of rows as our expressionset.
    ## The second column of this should contain the IDs from the other species when possible
    ## and a copy of the ID from our species when it is not.

    ## One final caveat: some of our new IDs may be duplicated in this (multigene families present
    ## in the other species),
    ## I will make.names() them and report how many there are, this might not be the correct way
    ## to handle this scenario!
    dup_ids <- sum(duplicated(reordered[[their_column]]))
    mesg("There are ", dup_ids, " duplicated IDs in the ", colnames(reordered)[2], " column.")
    ids <- make.names(reordered[[their_column]], unique = TRUE)
  }

  rownames(expr) <- ids
  expt[["expressionset"]] <- expr

  if (!is.null(expt[["tximport"]])) {
    rownames(expt[["tximport"]][["raw"]][["abundance"]]) <- ids
    rownames(expt[["tximport"]][["raw"]][["counts"]]) <- ids
    rownames(expt[["tximport"]][["raw"]][["length"]]) <- ids
    rownames(expt[["tximport"]][["scaled"]][["abundance"]]) <- ids
    rownames(expt[["tximport"]][["scaled"]][["counts"]]) <- ids
    rownames(expt[["tximport"]][["scaled"]][["length"]]) <- ids
  }
  return(expt)
}

#' Change the sample names of an expt.
#'
#' Sometimes one does not like the hpgl identifiers, so provide a way to change
#' them on-the-fly.
#'
#' @param expt Expt to modify
#' @param newnames New names, currently only a character vector.
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_conditions()] [set_expt_batches()]
#' @examples
#' \dontrun{
#'  expt = set_expt_samplenames(expt, c("a","b","c","d","e","f"))
#' }
#' @export
set_expt_samplenames <- function(expt, newnames) {
  if (length(newnames) == 1) {
    ## assume it is a factor in the metadata and act accordingly.
    newer_names <- make.names(pData(expt)[[newnames]], unique=TRUE)
    result <- set_expt_samplenames(expt, newer_names)
    return(result)
  }

  new_expt <- expt
  ## oldnames <- rownames(new_expt[["design"]])
  oldnames <- sampleNames(new_expt)
  newnames <- make.unique(newnames)
  newnote <- glue::glue("Sample names changed from: {toString(oldnames)} \\
                   to: {toString(newnames)} at: {date()}
")
  ## Things to modify include: batches, conditions
  names(new_expt[["batches"]]) <- newnames
  names(new_expt[["colors"]]) <- newnames
  names(new_expt[["conditions"]]) <- newnames
  newdesign <- new_expt[["design"]]
  newdesign[["oldnames"]] <- rownames(newdesign)
  rownames(newdesign) <- newnames
  newdesign[["sampleid"]] <- newnames
  new_expt[["design"]] <- newdesign
  new_expressionset <- new_expt[["expressionset"]]
  Biobase::sampleNames(new_expressionset) <- newnames
  pData(new_expressionset)[["sampleid"]] <- newnames
  new_expt[["expressionset"]] <- new_expressionset
  names(new_expt[["libsize"]]) <- newnames
  new_expt[["samplenames"]] <- newnames
  return(new_expt)
}

#' Extract a subset of samples following some rule(s) from an
#' experiment class.
#'
#' Sometimes an experiment has too many parts to work with conveniently, this
#' operation allows one to break it into smaller pieces.
#'
#' @param expt Expt chosen to extract a subset of data.
#' @param subset Valid R expression which defines a subset of the design to keep.
#' @param nonzero Look for a minimal number of nonzero genes.
#' @param ids List of sample IDs to extract.
#' @param coverage Request a minimum coverage/sample rather than text-based subset.
#' @return metadata Expt class which contains the smaller set of data.
#' @seealso [Biobase] [pData()] [exprs()] [fData()]
#' @examples
#' \dontrun{
#'  smaller_expt <- expt_subset(big_expt, "condition=='control'")
#'  all_expt <- expt_subset(expressionset, "")  ## extracts everything
#' }
#' @export
subset_expt <- function(expt, subset = NULL, ids = NULL,
                        nonzero = NULL, coverage = NULL) {
  starting_expressionset <- NULL
  starting_metadata <- NULL
  starting_samples <- sampleNames(expt)
  if ("ExpressionSet" %in% class(expt)) {
    starting_expressionset <- expt
    starting_metadata <- pData(starting_expressionset)
  } else if ("expt" %in% class(expt)) {
    starting_expressionset <- expt[["expressionset"]]
    starting_metadata <- pData(expt)
  } else {
    stop("expt is neither an expt nor ExpressionSet")
  }

  if (!is.null(ids)) {
    string <- ""
    for (id in ids) {
      string <- glue::glue("{string}|sampleid=='{id}'")
    }
    ## Remove the leading |
    subset <- substring(string, 2)
  }

  note_appended <- NULL
  subset_design <- NULL
  if (is.null(coverage) & is.null(nonzero)) {
    if (is.null(subset)) {
      subset_design <- starting_metadata
    } else {
      mesg("Using a subset expression.")
      r_expression <- paste("subset(starting_metadata,", subset, ")")
      subset_design <- eval(parse(text = r_expression))
      note_appended <- glue::glue("Subsetted with {subset} on {date()}.
")
    }
    if (nrow(subset_design) == 0) {
      stop("When the subset was taken, the resulting design has 0 members.")
    }
    subset_design <- as.data.frame(subset_design, stringsAsFactors = FALSE)
  } else if (is.null(nonzero)) {
    ## If coverage is defined, then use it to subset based on the minimal desired coverage
    ## Perhaps in a minute I will make this work for strings like '1z' to get the lowest
    ## standard deviation or somesuch...
    mesg("Subsetting given a minimal number of counts/sample.")
    coverages <- colSums(exprs(expt))
    subset_idx <- coverages >= as.numeric(coverage) ## In case I quote it on accident.
    subset_design <- starting_metadata[subset_idx, ]
    subset_design <- as.data.frame(subset_design, stringsAsFactors = FALSE)
    message("The samples removed (and read coverage) when filtering samples with less than ",
            coverage, " reads are: ")
    print(colSums(exprs(expt))[!subset_idx])
  } else if (is.null(coverage)) {
    ## Remove samples with less than this number of non-zero genes.
    nonzero_idx <- exprs(expt) != 0
    remove_idx <- colSums(nonzero_idx) < nonzero
    subset_design <- starting_metadata[!remove_idx, ]
    subset_design <- as.data.frame(subset_design, stringsAsFactors = FALSE)
    message("The samples (and read coverage) removed when filtering ",
            nonzero, " non-zero genes are: ")
    print(colSums(exprs(expt))[remove_idx])
  } else {
    stop("Unable to determine what is being subset.")
  }

  ## This is to get around stupidity with respect to needing all factors to be
  ## in a DESeqDataSet
  starting_ids <- rownames(starting_metadata)
  subset_ids <- rownames(subset_design)
  subset_positions <- starting_ids %in% subset_ids
  starting_colors <- expt[["colors"]]
  subset_colors <- starting_colors[subset_positions, drop = TRUE]
  starting_conditions <- expt[["conditions"]]
  subset_conditions <- starting_conditions[subset_positions, drop = TRUE]
  starting_batches <- expt[["batches"]]
  subset_batches <- starting_batches[subset_positions, drop = TRUE]
  current_libsize <- expt[["libsize"]]
  subset_current_libsize <- current_libsize[subset_positions, drop = TRUE]
  subset_expressionset <- starting_expressionset[, subset_positions]

  notes <- expt[["notes"]]
  if (!is.null(note_appended)) {
    notes <- glue::glue("{notes}{note_appended}")
  }

  for (col in seq_len(ncol(subset_design))) {
    if (class(subset_design[[col]])[1] == "factor") {
      subset_design[[col]] <- droplevels(subset_design[[col]])
    }
  }
  pData(subset_expressionset) <- subset_design

  new_expt <- list(
      "title" = expt[["title"]],
      "notes" = toString(notes),
      "initial_metadata" = subset_design,
      "expressionset" = subset_expressionset,
      "design" = subset_design,
      "conditions" = subset_conditions,
      "batches" = subset_batches,
      "samplenames" = subset_ids,
      "colors" = subset_colors,
      "state" = expt[["state"]],
      "libsize" = subset_current_libsize)
  class(new_expt) <- "expt"
  final_samples <- sampleNames(new_expt)
  message("subset_expt(): There were ", length(starting_samples), ", now there are ",
          length(final_samples), " samples.")
  return(new_expt)
}

#' Try a very literal subtraction
#'
#' @param expt Input expressionset.
#' @param new_meta dataframe containing the new metadata.
#' @param sample_column Column in the sample sheet to use to acquire the sample IDs given the
#'  subtractions.
#' @param convert_state Expected state of the input data vis a vis conversion (rpkm/cpm).
#' @param transform_state Expected state of the input data vis a vis transformation (log/linear).
#' @param handle_negative Set negative subtracted values to zero?
#' @param savefile Save the new expt data to this file.
#' @param ... Parameters to pass to normalize_expt()
#' @return New expt
#' @export
subtract_expt <- function(expt, new_meta, sample_column = "sample",
                          convert_state = "cpm", transform_state = "raw",
                          handle_negative = "zero", savefile = "subtracted.rda",
                          ...) {
  arglist <- list(...)
  if (expt[["state"]][["conversion"]] != convert_state) {
    expt <- normalize_expt(expt, convert = convert_state)
  }
  if (expt[["state"]][["transform"]] != transform_state) {
    expt <- normalize_expt(expt, transform = transform_state)
  }

  meta <- pData(expt)
  mtrx <- as.data.frame(exprs(expt))
  samples <- colnames(mtrx)
  new_exprs <- data.frame(row.names = rownames(mtrx))
  new_pdata <- data.frame()
  sub_names <- rownames(new_meta)
  for (s in seq_len(nrow(new_meta))) {
    sub_name <- sub_names[s]
    numerator <- new_meta[sub_name, "numerator"]
    denominator <- new_meta[sub_name, "denominator"]
    s1_idx <- meta[[sample_column]] == numerator
    s1 <- samples[s1_idx]
    s2_idx <- meta[[sample_column]] == denominator
    s2 <- samples[s2_idx]
    if (sum(s1_idx) != 1) {
      stop("Do not have 1 sample for subtraction: ", sub_name, ", ", numerator, "." )
    }
    if (sum(s2_idx) != 1) {
      stop("Do not have 1 sample for subtraction: ", sub_name, ", ", denominator, ".")
    }
    new_exprs[[sub_name]] <- mtrx[[s1]] - mtrx[[s2]]
    new_pdatum <- meta[s1, ]
    new_pdatum[["condition"]] <- new_meta[s, "condition"]
    new_pdatum[["batch"]] <- new_meta[s, "batch"]
    new_pdatum[["numerator"]] <- numerator
    new_pdatum[["denominator"]] <- denominator
    new_pdata <- rbind(new_pdata, new_pdatum)
  }
  rownames(new_pdata) <- sub_names
  colnames(new_exprs) <- sub_names

  negative_idx <- new_exprs < 0
  negative_pct <- (sum(negative_idx) / (nrow(new_exprs) * ncol(new_exprs))) * 100.0
  message("There are ", sum(negative_idx), " elements which are less than 0, (",
          signif(negative_pct, 3), "%)")
  if (is.null(handle_negative)) {
    message("Leaving negative values alone.")
  } else if (handle_negative[1] == "zero") {
    message("Setting negative values to zero.")
    new_exprs[negative_idx] <- 0
  } else if (handle_negative[1] == "na") {
    message("Setting negative values to NA.")
    new_exprs[negative_idx] <- NA
  } else {
    message("I do not understand this option, leaving negative values alone.")
  }

  new_pdata[["subtracted_samplenames"]] <- sub_names
  ## Now add the number of negative values observed.
  new_pdata[["negative_values"]] <- colSums(negative_idx)
  new_expt <- sm(create_expt(metadata = as.data.frame(new_pdata),
                             gene_info = fData(expt),
                             count_dataframe = as.matrix(new_exprs),
                             savefile = savefile,
                             id_column = "subtracted_samplenames"))
  new_expt[["na_values"]] <- as.data.frame(negative_idx)
  return(new_expt)
}

#' I want an easy way to sum counts in eupathdb-derived data sets.
#' These have a few things which should make this relatively easy.
#' Notably: The gene IDs look like: "exon_ID-1 exon_ID-2 exon_ID-3"
#' Therefore we should be able to quickly merge these.
#'
#' @param counts Matrix/df/dt of count data.
#' @return The same data type but with the exons summed.
sum_eupath_exon_counts <- function(counts) {
  rownames(counts) <- gsub(pattern = "^exon_", replacement = "", x = rownames(counts))
  rownames(counts) <- gsub(pattern = "\\-1$", replacement = "", x = rownames(counts))
  multi_exon_idx <- grep(pattern = "\\-\\d+$", x = rownames(counts))
  for (idx in multi_exon_idx) {
    gene <- gsub(pattern = "\\-\\d+$", replacement = "", x = rownames(counts)[idx])
    exon_counts <- counts[idx, ]
    gene_counts <- counts[gene, ]
    total_counts <- gene_counts + exon_counts
    counts[gene, ] <- total_counts
  }
  counts <- counts[-multi_exon_idx, ]
  counts <- as.matrix(counts)
  return(counts)
}

#' Print a string describing what happened to this data.
#'
#' Sometimes it is nice to have a string like: log2(cpm(data)) describing what
#' happened to the data.
#'
#' @param expt The expressionset.
#' @param transform How was it transformed?
#' @param convert How was it converted?
#' @param norm How was it normalized?
#' @param filter How was it filtered?
#' @param batch How was it batch-corrected?
#' @return An expression describing what has been done to this data.
#' @seealso [create_expt()] [normalize_expt()]
#' @export
what_happened <- function(expt = NULL, transform = "raw", convert = "raw",
                          norm = "raw", filter = "raw", batch = "raw") {
  if (! "expt" %in% class(expt)) {
    return("")
  }
  if (!is.null(expt)) {
    transform <- expt[["state"]][["transform"]]
    if (is.null(transform)) {
      transform <- "raw"
    }
    batch <- expt[["state"]][["batch"]]
    if (is.null(batch)) {
      batch <- "raw"
    }
    convert <- expt[["state"]][["conversion"]]
    if (is.null(convert)) {
      convert <- "raw"
    }
    norm <- expt[["state"]][["normalization"]]
    if (is.null(norm)) {
      norm <- "raw"
    }
    filter <- expt[["state"]][["filter"]]
    if (is.null(filter)) {
      filter <- "raw"
    }
  }
  ## Short circuit if nothing was done.
  if (transform == "raw" & batch == "raw" &
      convert == "raw" & norm == "raw" &
      filter == "raw") {
    what <- "raw(data)"
    return(what)
  }

  what <- ""
  if (transform != "raw") {
    what <- glue::glue("{what}{transform}(")
  }
  if (batch != "raw") {
    if (isTRUE(batch)) {
      what <- glue::glue("{what}batch-correct(")
    } else {
      what <- glue::glue("{what}{batch}(")
    }
  }
  if (convert != "raw") {
    what <- glue::glue("{what}{convert}(")
  }
  if (norm != "raw") {
    what <- glue::glue("{what}{norm}(")
  }
  if (filter != "raw") {
    what <- glue::glue("{what}{filter}(")
  }
  what <- glue::glue("{what}data")
  if (transform != "raw") {
    what <- glue::glue("{what})")
  }
  if (batch != "raw") {
    what <- glue::glue("{what})")
  }
  if (convert != "raw") {
    what <- glue::glue("{what})")
  }
  if (norm != "raw") {
    what <- glue::glue("{what})")
  }
  if (filter != "raw") {
    what <- glue::glue("{what})")
  }

  return(what)
}

#' Make pretty xlsx files of count data.
#'
#' Some folks love excel for looking at this data.  ok.
#'
#' Tested in test_03graph_metrics.R
#' This performs the following:  Writes the raw data, graphs the raw data,
#' normalizes the data, writes it, graphs it, and does a median-by-condition and
#' prints that.  I replaced the openxlsx function which writes images into xlsx
#' files with one which does not require an opening of a pre-existing plotter.
#' Instead it (optionally)opens a pdf device, prints the plot to it, opens a png
#' device, prints to that, and inserts the resulting png file.  Thus it
#' sacrifices some flexibility for a hopefully more consistent behaivor.  In
#' addition, one may use the pdfs as a set of images importable into illustrator
#' or whatever.
#'
#' @param expt An expressionset to print.
#' @param excel Filename to write.
#' @param norm Normalization to perform.
#' @param violin Include violin plots?
#' @param sample_heat Include sample heatmaps?
#' @param convert Conversion to perform.
#' @param transform Transformation used.
#' @param batch Batch correction applied.
#' @param filter Filtering method used.
#' @param med_or_mean When printing mean by condition, one may want median.
#' @param color_na Color cells which were NA before imputation this color.
#' @param merge_order Used to decide whether to put the counts or annotations first when
#'  printing count tables.
#' @param ... Parameters passed down to methods called here (graph_metrics, etc).
#' @return A big honking excel file and a list including the dataframes and images created.
#' @seealso [openxlsx] [Biobase] [normalize_expt()] [graph_metrics()]
#' @examples
#' \dontrun{
#'  excel_sucks <- write_expt(expt)
#' }
#' @export
write_expt <- function(expt, excel = "excel/pretty_counts.xlsx", norm = "quant",
                       violin = TRUE, sample_heat = TRUE, convert = "cpm", transform = "log2",
                       batch = "svaseq", filter = TRUE, med_or_mean = "mean",
                       color_na = "#DD0000", merge_order = "counts_first", ...) {
  arglist <- list(...)
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  new_row <- 1
  new_col <- 1

  ## Set up a vector of images to clean up when finished.
  image_files <- c()

  ## I think the plot dimensions should increase as the number of samples increase
  plot_dim <- 6
  num_samples <- ncol(exprs(expt))
  if (num_samples > 12) {
    plot_dim <- ceiling(num_samples / 4)
  }
  plot_cols <- floor(plot_dim * 1.5)
  plot_rows <- ceiling(plot_dim * 5.0)

  ## Write an introduction to this foolishness.
  message("Writing the first sheet, containing a legend and some summary data.")
  sheet <- "legend"
  norm_state <- glue::glue("{transform}({convert}({norm}({batch}({filter}(counts)))))")
  legend <- data.frame(
      "sheet" = c("1.", "2.", "3.", "4.", "5.", "6."),
      "sheet_definition" = c(
          "This sheet, including the experimental design.",
          "The raw counts and annotation data on worksheet 'raw_data'.",
          "Some graphs describing the distribution of raw data in worksheet 'raw_plots'.",
          glue::glue("The counts normalized with: {norm_state}"),
          "Some graphs describing the distribution of the normalized data on 'norm_plots'.",
          "The median normalized counts by condition factor on 'median_data'."),
      stringsAsFactors = FALSE)
  colnames(legend) <- c("Worksheets", "Contents")
  xls_result <- write_xlsx(data = legend, wb = wb, sheet = sheet, rownames = FALSE,
                           title = "Columns used in the following tables.")
  rows_down <- nrow(legend)
  new_row <- new_row + rows_down + 3
  annot <- as.data.frame(pData(expt), strinsAsFactors = FALSE)
  xls_result <- write_xlsx(data = annot, wb = wb, start_row = new_row, rownames = FALSE,
                           sheet = sheet, start_col = 1, title = "Experimental Design.")

  ## Get the library sizes and other raw plots before moving on...
  metrics <- sm(graph_metrics(expt, qq = TRUE,
                              ...))
  new_row <- new_row + nrow(pData(expt)) + 3
  libsizes <- as.data.frame(metrics[["libsizes"]])[, c("id", "sum", "condition")]
  xls_result <- write_xlsx(data = libsizes, wb = wb, start_row = new_row,
                           rownames = FALSE, sheet = sheet, start_col = 1,
                           title = "Library sizes.")

  new_row <- new_row + nrow(libsizes) + 3
  libsize_summary <- as.data.frame(metrics[["libsize_summary"]])
  xls_result <- write_xlsx(data = libsize_summary, wb = wb, start_row = new_row,
                           rownames = FALSE, sheet = sheet, start_col = 1,
                           title = "Library size summary.")

  ## Write the raw read data and gene annotations
  mesg("Writing the raw reads.")
  sheet <- "raw_reads"
  new_row <- 1
  new_col <- 1
  reads <- exprs(expt)
  info <- fData(expt)

  if (!is.null(info[["Row.names"]])) {
    ridx <- colnames(info) == "Row.names"
    message("Hey, you merged the annotation data and did not reset the column names!")
    colnames(info)[ridx] <- "old_row_names"
  }
  read_info <- merge(info, reads, by = "row.names")
  xls_result <- write_xlsx(data = read_info, wb = wb, sheet = sheet, rownames = FALSE,
                           start_row = new_row, start_col = new_col, title = "Raw Reads.")

  ## Potentially useful for proteomics data and subtracted data.
  if (!is.null(color_na)) {
    na_style <- openxlsx::createStyle(fontColour = color_na)
    nas <- expt[["na_values"]]
    for (col in colnames(nas)) {
      row_definition <- which(nas[[col]]) + 2
      col_idx <- colnames(read_info) == col
      col_definition <- which(col_idx)
      colored <- openxlsx::addStyle(wb, sheet = sheet, na_style,
                                    rows = row_definition, cols = col_definition)
    }
  }

  ## Write some graphs for the raw data
  mesg("Graphing the raw reads.")
  sheet <- "raw_graphs"
  newsheet <- try(openxlsx::addWorksheet(wb, sheetName = sheet))
  if (class(newsheet) == "try-error") {
    warning("Failed to add the sheet: ", sheet)
  }

  ## Start with library sizes.
  openxlsx::writeData(wb, sheet = sheet, x = "Legend.", startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw library sizes.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Non-zero genes.",
                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  new_col <- 1
  legend_plot <- metrics[["legend"]]
  try_result <- xlsx_insert_png(legend_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "01_legend", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  libsize_plot <- metrics[["libsize"]]
  try_result <- xlsx_insert_png(libsize_plot, wb = wb, sheet = sheet,
                              width = plot_dim, height = plot_dim,
                              start_col = new_col, start_row = new_row,
                              plotname = "02_libsize", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }

  ## Same row, non-zero plot
  new_col <- new_col + plot_cols + 1
  nonzero_plot <- metrics[["nonzero"]]
  try_result <- xlsx_insert_png(nonzero_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "03_nonzero", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1

  ## Visualize distributions
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw data density plot.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw Boxplot.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  density_plot <- metrics[["density"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(density_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "04_density", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  boxplot_plot <- metrics[["boxplot"]]
  try_result <- xlsx_insert_png(boxplot_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "05_boxplot", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  topn_plot <- metrics[["topnplot"]]
  try_result <- xlsx_insert_png(topn_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "06_topnplot", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  cv_plot <- metrics[["cvplot"]]
  try_result <- xlsx_insert_png(cv_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "07_cvplot", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  ## Move down next set of rows, heatmaps
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw correlation heatmap.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw distance heatmap.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  new_row <- new_row + 1
  corheat_plot <- metrics[["corheat"]]
  try_result <- xlsx_insert_png(corheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "08_corheat", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  disheat_plot <- metrics[["disheat"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(disheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "09_disheat", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  if (isTRUE(sample_heat)) {
    tmp_expt <- sm(normalize_expt(expt, transform = "log2", filter = TRUE))
    sampleheat_plot <- plot_sample_heatmap(tmp_expt)
    new_col <- new_col + plot_cols + 1
    try_result <- xlsx_insert_png(sampleheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "09a_sampleheat", savedir = excel_basename)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
  }
  new_col <- 1

  ## SM plots
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw standard median correlation.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw standard distance correlation.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  new_row <- new_row + 1
  smc_plot <- metrics[["smc"]]
  try_result <- xlsx_insert_png(smc_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "10_smc", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  smd_plot <- metrics[["smd"]]
  try_result <- xlsx_insert_png(smd_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "11_smd", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  ## PCA, PCA(l2cpm) and qq_log
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw PCA.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "top 40 PC1 loadings.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "PCA(log2(cpm())).",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw TSNE.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "TSNE(log2(cpm())).",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw QQ, log scale.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  new_row <- new_row + 1
  pca_plot <- metrics[["pc_plot"]]
  pca_topn <- metrics[["pc_loadplot"]]
  pca_table <- metrics[["pc_table"]]
  tsne_plot <- metrics[["tsne_plot"]]
  tsne_table <- metrics[["tsne_table"]]
  try_result <- xlsx_insert_png(pca_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "12_pcaplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(pca_topn, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "13_pctopn", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  tmp_data <- sm(normalize_expt(expt, transform = "log2", convert = "cpm", filter = TRUE,
                                ...))
  rpca <- plot_pca(tmp_data,
                   ...)
  rtsne <- plot_tsne(tmp_data,
                     ...)
  rspca_plot <- rpca[["plot"]]
  rtsne_plot <- rtsne[["plot"]]
  rpca_table <- rpca[["residual_df"]]
  rtsne_table <- rtsne[["residual_df"]]
  rm(tmp_data)
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(rspca_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "14_norm_pcaplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(tsne_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "15_tsneplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(rtsne_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "16_rtsneplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  qq_plot <- metrics[["qqlog"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(qq_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "17_qqlog", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  violin_plot <- NULL
  pct_plot <- NULL
  ## Violin plots
  if (isTRUE(violin)) {
    filt <- sm(normalize_expt(expt, filter = "simple"))
    varpart_raw <- suppressWarnings(try(simple_varpart(
        filt, predictor = NULL, factors = c("condition", "batch"))))
    if (class(varpart_raw) != "try-error") {
      violin_plot <- varpart_raw[["partition_plot"]]
      new_row <- new_row + plot_rows + 2
      new_col <- 1
      try_result <- xlsx_insert_png(violin_plot, wb = wb, sheet = sheet, width = plot_dim,
                                  height = plot_dim, start_col = new_col, start_row = new_row,
                                  plotname = "18_violin", savedir = excel_basename)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      new_col <- new_col + plot_cols + 1

      pct_plot <- varpart_raw[["percent_plot"]]
      try_result <- xlsx_insert_png(pct_plot, wb = wb, sheet = sheet, width = plot_dim,
                                  height = plot_dim, start_col = new_col, start_row = new_row,
                                  plotname = "19_pctvar", savedir = excel_basename)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
    }
  }

  ## PCA table
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Raw PCA res.",
                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  xls_result <- write_xlsx(data = metrics[["pc_summary"]], wb = wb, rownames = FALSE,
                           sheet = sheet, start_col = new_col, start_row = new_row)
  new_col <- xls_result[["end_col"]] + 6
  new_row <- new_row - 1
  openxlsx::writeData(wb, sheet, "Raw PCA table.",
                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  xls_result <- write_xlsx(data = metrics[["pc_table"]], wb = wb, rownames = FALSE,
                           sheet = sheet, start_row = new_row, start_col = new_col)

  ## Move on to the next sheet, normalized data
  mesg("Writing the normalized reads.")
  sheet <- "norm_data"
  new_col <- 1
  new_row <- 1
  norm_data <- sm(normalize_expt(expt = expt, transform = transform,
                                 convert = convert, batch = batch, filter = filter,
                                 ...))
  norm_reads <- exprs(norm_data)
  info <- fData(norm_data)
  read_info <- merge(norm_reads, info, by = "row.names")
  title <- what_happened(norm_data)
  xls_result <- write_xlsx(wb = wb, data = read_info, rownames = FALSE,
                           start_row = new_row, start_col = new_col, sheet = sheet, title = title)

  ## Potentially useful for proteomics data and subtracted data.
  if (!is.null(color_na)) {
    na_style <- openxlsx::createStyle(fontColour = color_na)
    nas <- expt[["na_values"]]
    for (col in colnames(nas)) {
      row_definition <- which(nas[[col]]) + 2
      col_idx <- colnames(read_info) == col
      col_definition <- which(col_idx)
      colored <- openxlsx::addStyle(wb, sheet = sheet, na_style,
                                    rows = row_definition, cols = col_definition)
    }
  }


  ## Graphs of the normalized data
  mesg("Graphing the normalized reads.")
  sheet <- "norm_graphs"
  newsheet <- try(openxlsx::addWorksheet(wb, sheetName = sheet))
  norm_metrics <- sm(graph_metrics(norm_data, qq = TRUE,
                                   ...))
  openxlsx::writeData(wb, sheet = sheet, x = "Raw PCA res.",
                      startRow = new_row, startCol = new_col)
  ## Start with library sizes.
  openxlsx::writeData(wb, sheet, "Legend.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet, "Normalized library sizes.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Non-zero genes.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  new_row <- new_row + 1
  new_plot <- norm_metrics[["legend"]]
  try_result <- xlsx_insert_png(new_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  nlibsize_plot <- norm_metrics[["libsize"]]
  try_result <- xlsx_insert_png(nlibsize_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "20_nlibsize", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  ## Same row, non-zero plot
  new_col <- new_col + plot_cols + 1
  nnzero_plot <- norm_metrics[["nonzero"]]
  try_result <- xlsx_insert_png(nnzero_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "21_nnzero", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1

  ## Visualize distributions
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized data density plot.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized Boxplot.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  ndensity_plot <- norm_metrics[["density"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(ndensity_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "22_ndensity", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  nboxplot_plot <- norm_metrics[["boxplot"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(nboxplot_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "23_nboxplot", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  ntopn_plot <- norm_metrics[["topnplot"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(ntopn_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "24_nboxplot", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  ## Move down next set of rows, heatmaps
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized correlation heatmap.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized distance heatmap.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  ncorheat_plot <- norm_metrics[["corheat"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(ncorheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "25_ncorheat", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  ndisheat_plot <- norm_metrics[["disheat"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(ndisheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "26_ndisheat", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  if (isTRUE(sample_heat)) {
    sampleheat_plot <- plot_sample_heatmap(norm_data)
    new_col <- new_col + plot_cols + 1
    try_result <- xlsx_insert_png(sampleheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "26a_sampleheat", savedir = excel_basename)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
  }
  new_col <- 1

  ## SM plots
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized standard median correlation.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized standard distance correlation.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  nsmc_plot <- norm_metrics[["smc"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(nsmc_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "27_nsmc", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  nsmd_plot <- norm_metrics[["smd"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(nsmd_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "28_nsmd", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  ## PCA and qq_log
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized PCA.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Top 40 PC1 loadings.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized TSNE.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized QQ, log scale.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  npca_plot <- norm_metrics[["pc_plot"]]
  npc_topnplot <- norm_metrics[["pc_loadplot"]]
  ntsne_plot <- norm_metrics[["tsne_plot"]]
  npca_table <- norm_metrics[["pc_table"]]
  ntsne_table <- norm_metrics[["tsne_table"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(npca_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "29_npcaplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(npc_topnplot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "30_npcloadplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(ntsne_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "31_ntsneplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  nqq_plot <- norm_metrics[["qqlog"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(nqq_plot, wb = wb, sheet = sheet, width = plot_dim,
                              height = plot_dim, start_col = new_col, start_row = new_row,
                              plotname = "32_nqqplot", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }

  new_col <- 1

  ## Violin plots
  nvarpart_plot <- NULL
  npct_plot <- NULL
  if (isTRUE(violin)) {
    varpart_norm <- suppressWarnings(try(simple_varpart(norm_data, predictor = NULL,
                                                        factors = c("condition", "batch"))))
    if (class(varpart_norm) != "try-error") {
      nvarpart_plot <- varpart_norm[["partition_plot"]]
      new_row <- new_row + plot_rows + 2
      new_col <- 1
      try_result <- xlsx_insert_png(nvarpart_plot, wb = wb, sheet = sheet, width = plot_dim,
                                  height = plot_dim, start_col = new_col, start_row = new_row,
                                  plotname = "33_nviolin", savedir = excel_basename)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      new_col <- new_col + plot_cols + 1
      npct_plot <- varpart_norm[["percent_plot"]]
      try_result <- xlsx_insert_png(npct_plot, wb = wb, sheet = sheet, width = plot_dim,
                                  height = plot_dim, start_col = new_col, start_row = new_row,
                                  plotname = "34_npctplot", savedir = excel_basename)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
    }
  }

  ## PCA table
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized PCA res.",
                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  xls_result <- write_xlsx(data = norm_metrics[["pc_summary"]], wb = wb, rownames = FALSE,
                           sheet = sheet, start_col = new_col, start_row = new_row)
  new_col <- xls_result[["end_col"]] + 6
  new_row <- new_row - 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized PCA table.",
                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  xls_result <- write_xlsx(data = norm_metrics[["pc_table"]], wb = wb, sheet = sheet,
                           rownames = FALSE, start_col = new_col, start_row = new_row)


  ## Give a median-by-factor accounting of the data
  mesg("Writing the median reads by factor.")
  sheet <- "median_data"
  new_col <- 1
  new_row <- 1
  median_data <- sm(median_by_factor(exprs(norm_data), fun = med_or_mean,
                                     fact = norm_data[["conditions"]]))
  med <- median_data[["medians"]]
  colnames(med) <- paste0(med_or_mean, "_", colnames(med))
  cv <- median_data[["cvs"]]
  colnames(cv) <- paste0("cv_", colnames(cv))
  median_data <- merge(med, cv, by = "row.names")
  rownames(median_data) <- median_data[["Row.names"]]
  median_data[["Row.names"]] <- NULL
  median_data_merged <- data.frame()
  if (merge_order == "annot_first") {
    median_data_merged <- merge(info, median_data, by.x = "row.names", by.y = "row.names")
  } else {
    median_data_merged <- merge(median_data, info, by.x = "row.names", by.y = "row.names")
  }
  rownames(median_data_merged) <- median_data_merged[["Row.names"]]
  median_data_merged[["Row.names"]] <- NULL
  xls_result <- write_xlsx(wb, data = median_data_merged, start_row = new_row, start_col = new_col,
                           rownames = FALSE, sheet = sheet, title = "Median Reads by factor.")

  ## Save the result
  save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  retlist <- list(
      "save" = save_result,
      "legend" = legend,
      "annotations" = info,
      "raw_reads" = reads,
      "design" = annot,
      "legend_plot" = legend_plot,
      "raw_libsize" = libsize_plot,
      "raw_nonzero" = nonzero_plot,
      "raw_density" = density_plot,
      "raw_cv" = cv_plot,
      "raw_boxplot" = boxplot_plot,
      "raw_corheat" = corheat_plot,
      "raw_disheat" = disheat_plot,
      "raw_smc" = smc_plot,
      "raw_smd" = smd_plot,
      "raw_pctopn" = pca_topn,
      "raw_pca" = pca_plot,
      "raw_pca_table" = pca_table,
      "raw_tsne" = tsne_plot,
      "raw_tsne_table" = tsne_table,
      "raw_scaled_pca" = rspca_plot,
      "raw_scaled_pca_table" = rpca_table,
      "raw_scaled_tsne" = rtsne_plot,
      "raw_scaled_tsne_table" = rtsne_table,
      "raw_qq" = qq_plot,
      "raw_violin" = violin_plot,
      "raw_percent" = pct_plot,
      "norm_reads" = norm_reads,
      "norm_libsize" = nlibsize_plot,
      "norm_nonzero" = nnzero_plot,
      "norm_density" = ndensity_plot,
      "norm_boxplot" = nboxplot_plot,
      "norm_corheat" = ncorheat_plot,
      "norm_disheat" = ndisheat_plot,
      "norm_smc" = nsmc_plot,
      "norm_smd" = nsmd_plot,
      "norm_pca" = npca_plot,
      "norm_pctopn" = npc_topnplot,
      "norm_pca_table" = npca_table,
      "norm_tsne" = ntsne_plot,
      "norm_tsne_table" = ntsne_table,
      "norm_qq" = nqq_plot,
      "norm_violin" = nvarpart_plot,
      "norm_pct" = npct_plot,
      "medians" = median_data
  )
  for (img in image_files) {
    removed <- try(suppressWarnings(file.remove(img)), silent=TRUE)
  }
  return(retlist)
}

## Make some methods which call Biobase on expts.
## This way we don't have to do Biobase::exprs/fData/pData/notes()
## But instead R will call them for us.

## Here is a note from Hadley which is relevant here:
## Another consideration is that
## order. For example, to define the method setMethod("foo", c("bar", "baz"),
## ...) you must already have created the foo generic and the two classes. By
## default, R code is loaded in alphabetical order, but that won’t always work
## for your situation.

#' An expt is an ExpressionSet superclass with a shorter name.
#'
#' It is also a simple list so that one may summarize it more simply,
#' provides colors and some slots to make one's life easier.
#' It is created via the function create_expt() which perhaps should be changed.
#'
#' Another important caveat: expressionSets and their methods are all S4; but I
#' did not want to write S4 methods, so I made my expt a S3 class.  As a result,
#' in order to make use of exprs, notes, pData, fData, and friends, I made use
#' of setMethod() to set up calls for the expressionSet portion of the expt
#' objects.
#'
#' @param ... Parameters for create_expt()
#' @slot title Title for the expressionSet.
#' @slot notes Notes for the expressionSet (redundant with S4 notes()).
#' @slot design Copy of the experimental metadata (redundant with pData()).
#' @slot annotation Gene annotations (redundant with fData()).
#' @slot gff_file filename of a gff file which feeds this data.
#' @slot state What is the state of the data vis a vis normalization,
#'  conversion, etc.
#' @slot conditions Usually the condition column from pData.
#' @slot batches Usually the batch column from pData.
#' @slot libsize Library sizes of the data in its current state.
#' @slot colors Chosen colors for plotting the data.
#' @slot tximport Data provided by tximport() to create the exprs() data.
#' @export expt
expt <- function(...) {
  create_expt(...)
}
expt_set <- setOldClass("expt")

## > showMethods("assay")
## Function: assay (package SummarizedExperiment)
## x="SummarizedExperiment", i="character"
## x="SummarizedExperiment", i="missing"
## x="SummarizedExperiment", i="numeric"

#' A series of setMethods for expts, ExpressionSets, and SummarizedExperiments.
#' @importFrom SummarizedExperiment assay assay<- colData colData<- rowData rowData<-

setMethod("annotation", signature = "expt",
          function(object) {
            Biobase::annotation(object[["expressionset"]])
          })
setMethod("annotation<-", signature = "expt",
          function(object, value) {
            fData(object[["expressionset"]]) <- value
            return(object)
          })
setMethod("assay", signature = "expt",
          function(x, withDimnames = TRUE, ...) {
            mtrx <- Biobase::exprs(x[["expressionset"]])
            return(mtrx)
          })
setMethod("assay", "ExpressionSet",
          function(x, withDimnames = TRUE, ...) {
            Biobase::exprs(x)
          })
setMethod("assay<-", signature = "expt",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::exprs(x[["expressionset"]]) <- value
            return(x)
          })
setMethod("assay<-", signature = "ExpressionSet",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::exprs(x) <- value
            return(x)
          })
setMethod("colData", "expt",
          function(x, withDimnames = TRUE, ...) {
            Biobase::pData(x[["expressionset"]])
          })
setMethod("colData", "ExpressionSet",
          function(x, withDimnames = TRUE, ...) {
            Biobase::pData(x)
          })
setMethod("colData<-", signature = "expt",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::pData(x[["expressionset"]]) <- value
            return(x)
          })
setMethod("colData<-", signature = "ExpressionSet",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::pData(x) <- value
            return(x)
          })
setMethod("exprs", signature = "expt",
          function(object) {
            Biobase::exprs(object[["expressionset"]])
          })
setMethod("exprs<-", signature = "expt",
          function(object, value) {
            exprs(object[["expressionset"]]) <- value
            return(object)
          })
setMethod("fData", signature = "expt",
          function(object) {
            Biobase::fData(object[["expressionset"]])
          })
setMethod("fData<-", signature = "expt",
          function(object, value) {
            fData(object[["expressionset"]]) <- value
            return(object)
          })
setMethod("notes", signature = "expt",
          function(object) {
            Biobase::notes(object[["expressionset"]])
          })
setMethod("pData", signature = "expt",
          function(object) {
            Biobase::pData(object[["expressionset"]])
          })
setMethod("pData<-", signature = "expt",
          function(object, value) {
            pData(object[["expressionset"]]) <- value
            return(object)
          })
setMethod("rowData", "expt",
          function(x, withDimnames = TRUE, ...) {
            Biobase::fData(x[["expressionset"]])
          })
setMethod("rowData", "ExpressionSet",
          function(x, withDimnames = TRUE, ...) {
            Biobase::fData(x)
          })
setMethod("rowData<-", signature = "expt",
          function(x, i, withDimnames = TRUE, ..., value) {
            Biobase::fData(x[["expressionset"]]) <- value
            return(x)
          })
##setMethod("rowData<-", signature = "Expressionset",
##          function(x, i, withDimnames = TRUE, ..., value) {
##            Biobase::fData(x) <- value
##            return(x)
##          })
setMethod("sampleNames", signature = "expt",
          function(object) {
            Biobase::sampleNames(object[["expressionset"]])
          })
setMethod("sampleNames<-", signature = "expt",
          function(object, value) {
            set_expt_samplenames(object, value)
          })

## EOF
