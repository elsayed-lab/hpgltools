#' Sum the reads/gene for multiple sequencing runs of a single condition/batch.
#'
#' On occasion we have multiple technical replicates of a sequencing run.  This can use a column in
#' the experimental design to identify those replicates and sum the counts into a single column in
#' the count tables.
#'
#' Untested as of 2016-12-01, but used in a couple of projects where sequencing runs got repeated.
#'
#' @param expt Experiment class containing the requisite metadata and count tables.
#' @param column Column of the design matrix used to specify which samples are replicates.
#' @return Expt with the concatenated counts, new design matrix, batches, conditions, etc.
#' @seealso \pkg{Biobase}
#'  \code{\link[Biobase]{exprs}} \code{\link[Biobase]{fData}} \code{\link[Biobase]{pData}}
#' @examples
#' \dontrun{
#'  compressed <- concatenate_runs(expt)
#' }
#' @export
concatenate_runs <- function(expt, column="replicate") {
  design <- expt[["design"]]
  replicates <- levels(as.factor(design[, column]))
  final_expt <- expt
  final_data <- NULL
  final_design <- NULL
  column_names <- list()
  colors <- list()
  conditions <- list()
  batches <- list()
  samplenames <- list()
  for (rep in replicates) {
    expression <- paste0(column, "=='", rep, "'")
    tmp_expt <- subset_expt(expt, expression)
    tmp_data <- rowSums(exprs(tmp_expt))
    tmp_design <- tmp_expt[["design"]][1, ]
    final_data <- cbind(final_data, tmp_data)
    final_design <- rbind(final_design, tmp_design)
    column_names[[rep]] <- as.character(tmp_design[, "sampleid"])
    colors[[rep]] <- as.character(tmp_expt[["colors"]][1])
    batches[[rep]] <- as.character(tmp_expt[["batches"]][1])
    conditions[[rep]] <- as.character(tmp_expt[["conditions"]][1])
    samplenames[[rep]] <- paste(conditions[[rep]], batches[[rep]], sep="-")
    colnames(final_data) <- column_names
  }
  final_expt[["design"]] <- final_design
  metadata <- new("AnnotatedDataFrame", final_design)
  sampleNames(metadata) <- colnames(final_data)
  feature_data <- new("AnnotatedDataFrame", fData(expt))
  featureNames(feature_data) <- rownames(final_data)
  experiment <- new("ExpressionSet", exprs=final_data,
                    phenoData=metadata, featureData=feature_data)
  final_expt[["expressionset"]] <- experiment
  final_expt[["original_expressionset"]] <- experiment
  final_expt[["samples"]] <- final_design
  final_expt[["colors"]] <- as.character(colors)
  final_expt[["batches"]] <- as.character(batches)
  final_expt[["conditions"]] <- as.character(conditions)
  final_expt[["samplenames"]] <- as.character(samplenames)
  return(final_expt)
}

#' Wrap bioconductor's expressionset to include some other extraneous
#' information.
#'
#' It is worth noting that this function has a lot of logic used to
#' find the count tables in the local filesystem.  This logic has been
#' superceded by simply adding a field to the .csv file called
#' 'file'.  create_expt() will then just read that filename, it may be
#' a full pathname or local to the cwd of the project.
#'
#' @param metadata Comma separated file (or excel) describing the samples with information like
#'     condition, batch, count_filename, etc.
#' @param gene_info Annotation information describing the rows of the data set, this often comes
#'     from a call to import.gff() or biomart or organismdbi.
#' @param count_dataframe If one does not wish to read the count tables from the filesystem, they
#'     may instead be fed as a data frame here.
#' @param sample_colors List of colors by condition, if not provided it will generate its own colors
#'     using colorBrewer.
#' @param title Provide a title for the expt?
#' @param notes Additional notes?
#' @param include_type I have usually assumed that all gff annotations should be used, but that is
#'     not always true, this allows one to limit to a specific annotation type.
#' @param include_gff Gff file to help in sorting which features to keep.
#' @param savefile Rdata filename prefix for saving the data of the resulting expt.
#' @param low_files Explicitly lowercase the filenames when searching the filesystem?
#' @param ... More parameters are fun!
#' @return  experiment an expressionset
#' @seealso \pkg{Biobase}
#'  \code{\link[Biobase]{pData}} \code{\link[Biobase]{fData}} \code{\link[Biobase]{exprs}}
#'  \code{\link{read_counts_expt}} \code{\link[hash]{as.list.hash}}
#' @examples
#' \dontrun{
#'  new_experiment <- create_expt("some_csv_file.csv", color_hash)
#'  ## Remember that this depends on an existing data structure of gene annotations.
#' }
#' @export
create_expt <- function(metadata=NULL, gene_info=NULL, count_dataframe=NULL,
                        sample_colors=NULL, title=NULL, notes=NULL,
                        include_type="all", include_gff=NULL,
                        savefile="expt", low_files=FALSE, ...) {
  arglist <- list(...)  ## pass stuff like sep=, header=, etc here

  if (is.null(metadata)) {
    stop("This requires some metadata at minimum.")
  }
  ## Palette for colors when auto-chosen
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }
  ## I am learning about simplifying vs. preserving subsetting
  ## This is a case of simplifying and I believe one which is good because I just want
  ## the string out from my list. Lets assume that palette is in fact an element in arglist,
  ## I really don't care that the name of the resturn is 'palette'; I already
  ## knew that by asking for it.
  if (is.null(title)) {
    title <- paste0("This is an expt class.")
  }
  if (is.null(notes)) {
    notes <- paste0("Created on ", date(), ".\n")
  }
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
  file_column <- "file"
  if (!is.null(arglist[["file_column"]])) {
    file_column <- arglist[["file_column"]]  ## Make it possible to have multiple count
    file_column <- tolower(file_column)
    file_column <- gsub(pattern="[[:punct:]]", replacement="", x=file_column)
    ## tables / sample in one sheet.
  }
  round <- FALSE
  if (!is.null(arglist[["round"]])) {
    round <- arglist[["round"]]
  }

  ## Read in the metadata from the provided data frame, csv, or xlsx.
  message("Reading the sample metadata.")
  sample_definitions <- data.frame()
  file <- NULL
  meta_dataframe <- NULL
  if (class(metadata) == "character") {
    ## This is a filename containing the metadata
    file <- metadata
  } else if (class(metadata) == "data.frame") {
    ## A data frame of metadata was passed.
    meta_dataframe <- metadata
  } else {
    stop("This requires either a file or meta data.frame.")
  }
  ## The two primary inputs for metadata are a csv/xlsx file or a dataframe, check for them here.
  if (is.null(meta_dataframe) & is.null(file)) {
    stop("This requires either a csv file or dataframe of metadata describing the samples.")
  } else if (is.null(file)) {
    ## punctuation is the devil
    sample_definitions <- meta_dataframe
    colnames(sample_definitions) <- tolower(colnames(sample_definitions))
    colnames(sample_definitions) <- gsub(pattern="[[:punct:]]", replacement="",
                                         x=colnames(sample_definitions))
  }  else {
    sample_definitions <- read_metadata(file, ...)
    ## sample_definitions <- read_metadata(file)
  }

  colnames(sample_definitions) <- tolower(colnames(sample_definitions))
  colnames(sample_definitions) <- gsub(pattern="[[:punct:]]", replacement="",
                                       x=colnames(sample_definitions))
  ## In case I am a doofus and repeated some column names.
  colnames(sample_definitions) <- make.names(colnames(sample_definitions), unique=TRUE)
  ## Check that condition and batch have been filled in.
  sample_columns <- colnames(sample_definitions)
  sample_column <- NULL
  ## The sample ID column should have the word 'sample' in it, otherwise this will fail.
  found_sample <- grepl(pattern="sample", x=sample_columns)
  if (sum(found_sample) == 0) {
    sample_definitions[["sampleid"]] <- make.names(rownames(sample_definitions), unique=TRUE)
  } else {
    ## Take the first column with the word 'sample' in it as the sampleid column
    sample_column <- sample_columns[found_sample][[1]]
  }
  rownames(sample_definitions) <- make.names(sample_definitions[[sample_column]], unique=TRUE)
  empty_samples <- which(sample_definitions[, sample_column] == "" |
                         is.na(sample_definitions[, sample_column]) |
                         grepl(pattern="^#", x=sample_definitions[, sample_column]))
                           if (length(empty_samples) > 0) {
                             sample_definitions <- sample_definitions[-empty_samples, ]
                           }
  ## The folks at Biobase, have changed new("AnnotatedDataFrame") so that if a column has
  ## the same information as the rownames(), then one gets the following error:
  ## new("AnnotatedDataFrame", ttt), which is utter BS.
  ## Ergo, if I have my sample IDs as the rownames _and_ as a column titled 'sampleid', it
  ## will now fail to create the ExpressionSet.
  ## I can therefore either remove the sampleid column or change it in some way.
  ## sample_definitions[[sample_column]] <- paste0("s", sample_definitions[[sample_column]])
  ## This error seems to have stopped happening, so I commented out the previous line.
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
  }
  found_batch <- "batch" %in% sample_columns
  if (!isTRUE(found_batch)) {
    message("Did not find the batch column in the sample sheet.")
    message("Filling it in as undefined.")
    sample_definitions[["batch"]] <- "undefined"
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
                                                       sample_definitions[["stage"]], sep="_"))
  }
  ## Extract out the condition names as a factor
  condition_names <- unique(sample_definitions[["condition"]])
  if (is.null(condition_names)) {
    warning("There is no 'condition' field in the definitions, this will make many
analyses more difficult/impossible.")
  }
  ## Condition and Batch are not allowed to be numeric, so if they are just numbers,
  ## prefix them with 'c' and 'b' respectively.
  sample_definitions[["condition"]] <- gsub(pattern="^(\\d+)$", replacement="c\\1",
                                            x=sample_definitions[["condition"]])
  sample_definitions[["batch"]] <- gsub(pattern="^(\\d+)$", replacement="b\\1",
                                        x=sample_definitions[["batch"]])

  ## Explicitly skip those samples which are "", null, or "undef" for the filename.
  if (is.null(count_dataframe)) {
    skippers <- (sample_definitions[[file_column]] == "" |
                 is.null(sample_definitions[[file_column]]) |
                 is.na(sample_definitions[[file_column]]) |
                 sample_definitions[[file_column]] == "undef")
    if (length(skippers) > 0) {
      ## If there is nothing to skip, do not try.
      sample_definitions <- sample_definitions[!skippers, ]
    }
  }
  num_samples <- nrow(sample_definitions)
  ## Create a matrix of counts with columns as samples and rows as genes
  ## This may come from either a data frame/matrix, a list of files from the metadata
  ## or it can attempt to figure out the location of the files from the sample names.
  filenames <- NULL
  all_count_tables <- NULL
  if (!is.null(count_dataframe)) {
    all_count_tables <- data.table::as.data.table(count_dataframe, keep.rownames="rownames")
    testthat::expect_equal(colnames(all_count_tables), c("rownames", rownames(sample_definitions)))
    ## If neither of these cases is true, start looking for the files in the
    ## processed_data/ directory
  } else if (is.null(sample_definitions[[file_column]])) {
    success <- 0
    ## There are two main organization schemes I have used in the past, the following
    ## checks for both in case I forgot to put a file column in the metadata.
    ## Look for files organized by sample
    test_filenames <- file.path("preprocessing", "count_tables",
                                as.character(sample_definitions[[sample_column]]),
                                paste0(file_prefix,
                                       as.character(sample_definitions[[sample_column]]),
                                       file_suffix))
    num_found <- sum(file.exists(test_filenames))
    if (num_found == num_samples) {
      success <- success + 1
      sample_definitions[[file_column]] <- test_filenames
    } else {
      lower_test_filenames <- tolower(test_filenames)
      num_found <- sum(file.exists(lower_test_filenames))
      if (num_found == num_samples) {
        success <- success + 1
        sample_definitions[[file_column]] <- lower_test_filenames
      }
    }
    if (success == 0) {
      ## Did not find samples by id, try them by type
      test_filenames <- file.path("preprocessing", "count_tables",
                                  tolower(as.character(sample_definitions[["type"]])),
                                  tolower(as.character(sample_definitions[["stage"]])),
                                  paste0(sample_definitions[[sample_column]], file_suffix))
      num_found <- sum(file.exists(test_filenames))
      if (num_found == num_samples) {
        success <- success + 1
        sample_definitions[[file_column]] <- test_filenames
      } else {
        test_filenames <- tolower(test_filenames)
        num_found <- sum(file.exists(test_filenames))
        if (num_found == num_samples) {
          success <- success + 1
          sample_definitions[[file_column]] <- test_filenames
        }
      }
    } ## tried by type
    if (success == 0) {
      stop("I could not find your count tables by sample nor type, uppercase nor lowercase.")
    }
  }

  ## At this point sample_definitions$file should be filled in no matter what;
  ## so read the files.
  tximport_data <- NULL
  if (is.null(all_count_tables)) {
    filenames <- as.character(sample_definitions[[file_column]])
    sample_ids <- as.character(sample_definitions[[sample_column]])
    all_count_tables <- read_counts_expt(sample_ids, filenames, ...)
    ## all_count_tables <- read_counts_expt(sample_ids, filenames, arglist)
    if (all_count_tables[["source"]] == "tximport") {
      tximport_data <- list("raw" = all_count_tables[["tximport"]],
                            "scaled" = all_count_tables[["tximport_scaled"]])
    }
    all_count_tables <- all_count_tables[["count_table"]]
  }

  ## I have had a couple data sets with incomplete counts, get rid of those rows before moving on.
  all_count_tables <- all_count_tables[complete.cases(all_count_tables), ]
  numeric_columns <- colnames(all_count_tables) != "rownames"
  for (col in colnames(all_count_tables)[numeric_columns]) {
    ## Ensure there are no stupid entries like target_id est_counts
    all_count_tables[[col]] <- as.numeric(all_count_tables[[col]])
  }
  ## Features like exon:alicethegene-1 are annoying and entirely too common in TriTrypDB data
  all_count_tables[["rownames"]] <- gsub(pattern="^exon:", replacement="",
                                         x=all_count_tables[["rownames"]])
  all_count_tables[["rownames"]] <- make.names(gsub(pattern=":\\d+", replacement="",
                                                    x=all_count_tables[["rownames"]]),
                                               unique=TRUE)

  ## Try a couple different ways of getting gene-level annotations into the expressionset.
  annotation <- NULL
  if (is.null(gene_info)) {
    ## Including, if all else fails, just grabbing the gene names from the count tables.
    if (is.null(include_gff)) {
      gene_info <- data.table::as.data.table(all_count_tables[["rownames"]], keep.rownames="rownames")
      names(gene_info) <- "rownames"
    } else {
      ## Or reading a gff file.
      message("create_expt(): Reading annotation gff, this is slow.")
      annotation <- load_gff_annotations(gff=include_gff, type=gff_type)
      gene_info <- data.table::as.data.table(annotation, keep.rownames="rownames")
    }
  } else if (class(gene_info)[[1]] == "list" & !is.null(gene_info[["genes"]])) {
    ## In this case, it is using the output of reading a OrgDB instance
    gene_info <- data.table::as.data.table(gene_info[["genes"]], keep.rownames="rownames")
  } else if (class(gene_info)[[1]] == "data.table") {
    ## Try to make the data table usage consistent by rownames.
    ## Sometimes we take these from data which did "keep.rownames='some_column'"
    ## Sometimes we take these from data which set rownames(dt)
    ## And sometimes the rownames were never set.
    ## Therefore I will use rownames(dt) as the master, dt$rownames as secondary, and
    ## as a fallback take the first column in the data.
    if (is.null(rownames(gene_info)) & is.null(gene_info[["rownames"]])) {
      gene_info[["rownames"]] <- make.names(rownames[[1]], unique=TRUE)
    } else if (!is.null(rownames(gene_info))) {
      gene_info[["rownames"]] <- rownames(gene_info)
    }
  } else {
    gene_info <- data.table::as.data.table(gene_info, keep.rownames="rownames")
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
      message(paste0("Matched ", found_sum, " annotations and counts."))
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
  counts_and_annotations <- merge(all_count_tables, gene_info, by="rownames", all.x=TRUE)
  ## In some cases, the above merge will result in columns being set to NA
  ## We should set all the NA fields to something I think.
  na_entries <- is.na(counts_and_annotations)
  if (sum(na_entries) > 0) {
    message("Some annotations were lost in merging, setting them to 'undefined'.")
  }
  counts_and_annotations[na_entries] <- "undefined"
  ## Set an incrementing id number to make absolutely paranoidly certain the order stays constant.
  counts_and_annotations <- counts_and_annotations[order(counts_and_annotations[["temporary_id_number"]]), ]
  ## Pull out the annotation data and convert to data frame.
  kept_columns <- colnames(counts_and_annotations) %in% colnames(gene_info)
  final_annotations <- counts_and_annotations[, kept_columns, with=FALSE]
  final_annotations <- as.data.frame(final_annotations, stringsAsFactors=FALSE)
  rownames(final_annotations) <- final_annotations[["rownames"]]
  final_kept <- colnames(final_annotations) != "rownames"
  final_annotations <- final_annotations[, final_kept]

  ## There are some shenanigans, Maddy is getting an error on countsdt...
  final_counts <- counts_and_annotations
  kept_columns <- colnames(counts_and_annotations) %in% colnames(all_count_tables) &
    colnames(counts_and_annotations) != "temporary_id_number"
  final_counts <- final_counts[, kept_columns, with=FALSE]
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
  ## > test_df = data.frame(a=c(1,1,1), b=c(2,2,2))
  ## > rownames(test_df) = c("TcCLB.511511.3","TcCLB.511511.30","bob")
  ## > test_df
  ##                 a b
  ## TcCLB.511511.3  1 2
  ## TcCLB.511511.30 1 2
  ## bob             1 2
  ## > rownames(test_df) <- gsub(pattern="(\\d)$", replacement="\\1\\.", x=rownames(test_df))
  ## > test_df
  ##                  a b
  ## TcCLB.511511.3.  1 2
  ## TcCLB.511511.30. 1 2
  ## bob              1 2
  ## This is so stupid I think I am going to go and cry.
  ##rownames(final_annotations) <- gsub(pattern="(\\d)$", replacement="\\1\\.",
  ##                                    x=rownames(final_annotations), perl=TRUE)
  ##rownames(final_counts) <- gsub(pattern="(\\d)$", replacement="\\1\\.",
  ##                                    x=rownames(final_counts), perl=TRUE)

  ##final_counts <- final_counts[, -1, drop=FALSE]

  ## If the user requests input of non-int counts, fix that here.
  if (isTRUE(round)) {
    final_counts <- round(final_counts)
    less_than <- final_counts < 0
    final_counts[less_than] <- 0
  }

  ## I moved the color choices to this area pretty late in the process to make sure that
  ## there was time to remove unused samples.
  ## Make sure we have a viable set of colors for plots
  ## First figure out how many conditions we have
  chosen_colors <- as.character(sample_definitions[["condition"]])
  num_conditions <- length(condition_names)
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
    sample_colors <- suppressWarnings(grDevices::colorRampPalette(
                                                   RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else {
    ## If none of the above are true, then warn the user and let RColorBrewer do it.
    warning("The number of colors provided does not match either the number of conditions nor samples.")
    warning("Unsure of what to do, so choosing colors with RColorBrewer.")
    sample_colors <- suppressWarnings(grDevices::colorRampPalette(
                                                   RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  }
  ## Set the color names
  names(chosen_colors) <- sample_definitions[[sample_column]]

  ## Perhaps I do not understand something about R's syntactic sugar
  ## Given a data frame with columns bob, jane, alice -- but not foo
  ## I can do df[["bob"]]) or df[, "bob"] to get the column bob
  ## however df[["foo"]] gives me null while df[, "foo"] gives an error.
  if (is.null(sample_definitions[["condition"]])) {
    sample_definitions[["condition"]] <- "unknown"
  }
  if (is.null(sample_definitions[["batch"]])) {
    sample_definitions[["batch"]] <- "unknown"
  }
  if (is.null(sample_definitions[["file"]])) {
    sample_definitions[["file"]] <- "null"
  }

  ## Adding these so that deseq does not complain about characters when
  ## calling DESeqDataSetFromMatrix()
  sample_definitions[["condition"]] <- as.factor(sample_definitions[["condition"]])
  sample_definitions[["batch"]] <- as.factor(sample_definitions[["batch"]])

  ## Finally, create the ExpressionSet using the counts, annotations, and metadata.
  requireNamespace("Biobase")  ## AnnotatedDataFrame is from Biobase
  metadata <- methods::new("AnnotatedDataFrame", sample_definitions)
  Biobase::sampleNames(metadata) <- colnames(final_counts)

  feature_data <- methods::new("AnnotatedDataFrame", final_annotations)
  Biobase::featureNames(feature_data) <- rownames(final_counts)
  experiment <- methods::new("ExpressionSet",
                             exprs=final_counts,
                             phenoData=metadata,
                             featureData=feature_data)
  Biobase::notes(experiment) <- toString(notes)

  ## These entries in new_expt are intended to maintain a record of
  ## the transformation status of the data, thus if we now call
  ## normalize_expt() it should change these.
  ## Therefore, if we call a function like DESeq() which requires
  ## non-log2 counts, we can check these values and convert accordingly

  ## Now that the expressionset has been created, pack it into an expt object so that I
  ## can keep backups etc.
  expt <- subset_expt(experiment) ## I think this is spurious now.
  expt[["original_expressionset"]] <- experiment
  expt[["original_metadata"]] <- Biobase::pData(experiment)

  ## I only leared fairly recently that there is quite a bit of redundancy between my expt
  ## and ExpressionSets. I do not mind this. Yet.
  expt[["title"]] <- title
  expt[["notes"]] <- toString(notes)
  expt[["design"]] <- sample_definitions
  expt[["annotation"]] <- annotation
  expt[["gff_file"]] <- include_gff
  ## the 'state' slot in the expt is used to keep track of how the data is modified over time.
  starting_state <- list(
    "lowfilter" = "raw",
    "normalization" = "raw",
    "conversion" = "raw",
    "batch" = "raw",
    "transform" = "raw")
  expt[["state"]] <- starting_state
  ## Just in case there are condition names which are not used.
  expt[["conditions"]] <- droplevels(as.factor(sample_definitions[, "condition"]))
  ## This might be redundant, but it ensures that no all-numeric conditions exist.
  expt[["conditions"]] <- gsub(pattern="^(\\d+)$", replacement="c\\1", x=expt[["conditions"]])
  names(expt[["conditions"]]) <- rownames(sample_definitions)
  ## Ditto for batches
  expt[["batches"]] <- droplevels(as.factor(sample_definitions[, "batch"]))
  expt[["batches"]] <- gsub(pattern="^(\\d+)$", replacement="b\\1", x=expt[["batches"]])
  names(expt[["batches"]]) <- rownames(sample_definitions)
  ## Keep a backup of the metadata in case we do semantic filtering or somesuch.
  expt[["original_metadata"]] <- metadata
  ## Keep a backup of the library sizes for limma.
  expt[["original_libsize"]] <- colSums(Biobase::exprs(experiment))
  names(expt[["original_libsize"]]) <- rownames(sample_definitions)
  expt[["libsize"]] <- expt[["original_libsize"]]
  names(expt[["libsize"]]) <- rownames(sample_definitions)
  ## Save the chosen colors
  expt[["colors"]] <- chosen_colors
  names(expt[["colors"]]) <- rownames(sample_definitions)
  expt[["tximport"]] <- tximport_data
  ## Save an rdata file of the expressionset.
  if (!is.null(savefile)) {
    save_result <- try(save(list = c("expt"), file=paste(savefile, ".Rdata", sep="")), silent=TRUE)
  }
  if (class(save_result) == "try-error") {
    warning("Saving the expt object failed, perhaps you do not have permissions?")
  }
  return(expt)
}

#' Exclude some genes given a pattern match
#'
#' Because I am too lazy to remember that expressionsets use matrix subsets for [gene,sample]
#'
#' @param expt  Expressionset containing expt object.
#' @param column  fData column to use for subsetting.
#' @param method  Either remove explicit rows, or keep them.
#' @param ids  Specific IDs to exclude.
#' @param patterns  Character list of patterns to remove/keep
#' @param ...  Extra arguments are passed to arglist, currently unused.
#' @return  A smaller expt
#' @seealso \code{\link{create_expt}}
#' @export
exclude_genes_expt <- function(expt, column="txtype", method="remove", ids=NULL,
                               patterns=c("snRNA","tRNA","rRNA"), ...) {
  arglist <- list(...)
  ex <- expt[["expressionset"]]
  annotations <- Biobase::fData(ex)
  if (is.null(ids) & is.null(annotations[[column]])) {
    message(paste0("The ", column, " column is null, doing nothing."))
    return(expt)
  }
  pattern_string <- ""
  for (pat in patterns) {
    pattern_string <- paste0(pattern_string, pat, "|")
  }
  silly_string <- gsub(pattern="\\|$", replacement="", x=pattern_string)
  idx <- rep(x=TRUE, times=nrow(annotations))
  if (is.null(ids)) {
    idx <- grepl(pattern=silly_string, x=annotations[[column]], perl=TRUE)
  } else if (is.logical(ids)) {
    idx <- ids
  } else {
    idx <- rownames(annotations) %in% ids
  }
  kept <- NULL
  removed <- NULL
  kept_sums <- NULL
  removed_sums <- NULL
  if (method == "remove") {
    kept <- ex[!idx, ]
    removed <- ex[idx, ]
  } else {
    kept <- ex[idx, ]
    removed <- ex[!idx, ]
  }

  message(paste0("Before removal, there were ", nrow(Biobase::fData(ex)), " entries."))
  message(paste0("Now there are ", nrow(Biobase::fData(kept)), " entries."))
  all_tables <- Biobase::exprs(ex)
  all_sums <- colSums(all_tables)
  kept_tables <- Biobase::exprs(kept)
  kept_sums <- colSums(kept_tables)
  removed_tables <- Biobase::exprs(removed)
  removed_sums <- colSums(removed_tables)
  pct_kept <- (kept_sums / all_sums) * 100.0
  pct_removed <- (removed_sums / all_sums) * 100.0
  summary_table <- rbind(kept_sums, removed_sums, all_sums,
                         pct_kept, pct_removed)
  rownames(summary_table) <- c("kept_sums", "removed_sums", "all_sums",
                               "pct_kept", "pct_removed")
  message(paste0("Percent kept: ", toString(sprintf(fmt="%.3f", pct_kept))))
  message(paste0("Percent removed: ", toString(sprintf(fmt="%.3f", pct_removed))))
  expt[["expressionset"]] <- kept
  expt[["summary_table"]] <- summary_table
  return(expt)
}

#' Count the number of features(genes) greater than x in a data set.
#'
#' Sometimes I am asked how many genes have >= x counts.  Well, here you go.
#'
#' Untested as of 2016-12-01 but used with Lucia.  I think it would be interesting to iterate
#' this function from small to large cutoffs and plot how the number of kept genes decreases.
#'
#' @param data  A dataframe/exprs/matrix/whatever of counts.
#' @param cutoff  Minimum number of counts.
#' @param hard  Greater-than is hard, greater-than-equals is not.
#' @return  Number of genes.
#' @seealso \pkg{Biobase}
#' @examples
#' \dontrun{
#'  features <- features_greater_than(expt)
#' }
#' @export
features_greater_than <- function(data, cutoff=1, hard=TRUE) {
  if (class(data) == "expt" | class(data) == "ExpressionSet") {
    data <- as.data.frame(exprs(data))
  } else {
    data <- as.data.frame(data)
  }
  result <- numeric(length=ncol(data))
  names(result) <- colnames(data)
  for (col in 1:length(colnames(data))) {
    column_name <- colnames(data)[[col]]
    column_data <- data[[column_name]]
    num_features <- NULL
    if (isTRUE(hard)) {
      num_features <- sum(column_data > cutoff)
    } else {
      num_features <- sum(column_data >= cutoff)
    }
    result[[column_name]] <- num_features
  }
  return(result)
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
#' @seealso \pkg{limma} \pkg{stats} \pkg{DESeq}
#' @examples
#' \dontrun{
#'  pretend = make_exampledata()
#' }
#' @export
make_exampledata <- function (ngenes=1000, columns=5) {
  q0 <- stats::rexp(ngenes, rate = 1/250)
  is_DE <- stats::runif(ngenes) < 0.3
  lfc <- stats::rnorm(ngenes, sd = 2)
  q0A <- ifelse(is_DE, q0 * 2^ (lfc / 2), q0)
  q0B <- ifelse(is_DE, q0 * 2^ (-lfc / 2), q0)
  ##    true_sf <- c(1, 1.3, 0.7, 0.9, 1.6)
  true_sf <- abs(stats::rnorm(columns, mean=1, sd=0.4))
  cond_types <- ceiling(sqrt(columns))
  ##    conds <- c("A", "A", "B", "B", "B")
  ##x <- sample( LETTERS[1:4], 10000, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05) )
  conds <- sample(LETTERS[1:cond_types], columns, replace=TRUE)
  m <- t(sapply(seq_len(ngenes),
                function(i) sapply(1:columns,
                                   function(j) rnbinom(1,
                                                       mu = true_sf[j] * ifelse(conds[j] == "A",
                                                                                q0A[i], q0B[i]),
                                                       size = 1/0.2))))
  rownames(m) <- paste("gene", seq_len(ngenes), ifelse(is_DE, "T", "F"), sep = "_")
  example <- DESeq::newCountDataSet(m, conds)
  return(example)
}

#' Create a data frame of the medians of rows by a given factor in the data.
#'
#' This assumes of course that (like expressionsets) there are separate columns for each replicate
#' of the conditions.  This will just iterate through the levels of a factor describing the columns,
#' extract them, calculate the median, and add that as a new column in a separate data frame.
#'
#' Used in write_expt() as well as a few random collaborations.
#'
#' @param data Data frame, presumably of counts.
#' @param fact Factor describing the columns in the data.
#' @return Data frame of the medians.
#' @seealso \pkg{Biobase} \pkg{matrixStats}
#' @examples
#' \dontrun{
#'  compressed = median_by_factor(data, experiment$condition)
#' }
#' @export
median_by_factor <- function(data, fact="condition") {
  if (length(fact) == 1) {
    design <- pData(data)
    fact <- design[[fact]]
    names(fact) <- rownames(design)
  }
  if (class(data) == "expt" | class(data) == "ExpressionSet") {
    data <- exprs(data)
  }

  medians <- data.frame("ID"=rownames(data), stringsAsFactors=FALSE)
  data <- as.matrix(data)
  rownames(medians) <- rownames(data)
  fact <- as.factor(fact)
  for (type in levels(fact)) {
    columns <- grep(pattern=type, x=fact)
    med <- NULL
    if (length(columns) < 1) {
      warning("This level of the factor has no columns.")
      next
    } else if (length(columns) == 1) {
      message(paste0("The factor ", type, " has only 1 row."))
      med <- as.data.frame(data[, columns], stringsAsFactors=FALSE)
    } else {
      message(paste0("The factor ", type, " has ", length(columns), " rows."))
      med <- matrixStats::rowMedians(data[, columns])
    }
    medians <- cbind(medians, med)
  }
  medians <- medians[, -1, drop=FALSE]
  colnames(medians) <- levels(fact)
  return(medians)
}

#' Create a Schizosaccharomyces cerevisiae expt.
#'
#' @param annotation  Add annotation data?
#' @export
make_pombe_expt <- function(annotation=TRUE) {
  tt <- sm(require.auto("fission"))
  tt <- sm(requireNamespace("fission"))
  tt <- sm(try(attachNamespace("fission"), silent=TRUE))
  tt <- data(fission)
  meta <- as.data.frame(fission@colData)
  meta[["condition"]] <- paste0(meta[["strain"]], ".", meta[["minute"]])
  meta[["batch"]] <- meta[["replicate"]]
  meta[["sample.id"]] <- rownames(meta)
  fission_data <- fission@assays$data[["counts"]]

  annotations <- NULL
  if (isTRUE(annotation)) {
    ## Neat, it works, and even figures out that the default mart is incorrect by itself.
    pombe_annotations <- sm(load_biomart_annotations(
      host="fungi.ensembl.org",
      trymart="fungal_mart",
      trydataset="spombe_eg_gene",
      gene_requests=c("pombase_transcript", "ensembl_gene_id", "ensembl_transcript_id",
                      "hgnc_symbol", "description", "gene_biotype"),
      species="spombe", overwrite=TRUE))
    pombe_mart <- pombe_annotations[["mart"]]
    annotations <- pombe_annotations[["annotation"]]
    rownames(annotations) <- make.names(gsub(pattern="\\.\\d+$",
                                             replacement="",
                                             x=rownames(annotations)), unique=TRUE)
  }
  pombe_expt <- sm(create_expt(metadata=meta,
                               count_dataframe=fission_data,
                               gene_info=annotations))
  detach("package:fission")
  return(pombe_expt)
}

#' Read a bunch of count tables and create a usable data frame from them.
#'
#' It is worth noting that this function has some logic intended for the elsayed lab's data
#' storage structure. It shouldn't interfere with other usages, but it attempts to take into
#' account different ways the data might be stored.
#'
#' Used primarily in create_expt()
#' This is responsible for reading count tables given a list of filenames.  It tries to take into
#' account upper/lowercase filenames and uses data.table to speed things along.
#'
#' @param ids List of experimental ids.
#' @param files List of files to read.
#' @param header Whether or not the count tables include a header row.
#' @param include_summary_rows Whether HTSeq summary rows should be included.
#' @param suffix Optional suffix to add to the filenames when reading them.
#' @param ... More options for happy time!
#' @return Data frame of count tables.
#' @seealso \pkg{data.table}
#'  \code{\link{create_expt}}
#' @examples
#' \dontrun{
#'  count_tables <- hpgl_read_files(as.character(sample_ids), as.character(count_filenames))
#' }
#' @export
read_counts_expt <- function(ids, files, header=FALSE, include_summary_rows=FALSE,
                             suffix=NULL, ...) {
  ## load first sample
  arglist <- list(...) 
  skippers <- (files == "" | files == "undef" | is.null(files))
  files <- files[!skippers]
  lower_filenames <- files
  dirs <- dirname(lower_filenames)
  low_files <- tolower(basename(files))
  if (!is.null(suffix)) {
    low_hpgl <- gsub(suffix, "", basename(files))
    low_hpgl <- tolower(low_hpgl)
    low_hpgl <- paste(low_hpgl, suffix, sep="")
  } else {
    low_hpgl <- gsub("HPGL", "hpgl", basename(files))
  }
  lower_filenames <- paste(dirs, low_files, sep="/")
  lowhpgl_filenames <- paste(dirs, low_hpgl, sep="/")
  if (file.exists(tolower(files[1]))) {
    files[1] <- tolower(files[1])
  } else if (file.exists(lowhpgl_filenames[1])) {
    files[1] <- lowhpgl_filenames[1]
  } else if (file.exists(lower_filenames[1])) {
    files[1] <- lower_filenames[1]
  }

  ## Add an optional directory if I don't feel like specifying in the sample sheet.
  if (!is.null(arglist[["countdir"]])) {
    files <- file.path(arglist[["countdir"]], files)
  }
  count_table <- NULL
  for (f in 1:length(files)) {
    files[f] <- gsub(pattern=" ", replacement="", x=files[f])
    if (!grepl(pattern="^\\/", x=files[f])) {
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
  message("Reading count tables.")
  retlist <- list()
  txout <- TRUE
  tx_gene_map <- NULL
  if (!is.null(arglist[["tx_gene_map"]])) {
    message("Using the transcript<->gene mapping.")
    txout <- FALSE
    tx_gene_map <- arglist[["tx_gene_map"]]
  }
  if (grepl(pattern="\\.tsv|\\.h5", x=files[1])) {
    ## This hits if we are using the kallisto outputs.
    names(files) <- ids
    if (!all(file.exists(files))) {
      warning(files)
    }
    import <- NULL
    import_scaled <- NULL
    if (is.null(tx_gene_map)) {
      import <- sm(tximport::tximport(files=files, type="kallisto", txOut=txout))
      import_scaled <- sm(tximport::tximport(files=files, type="kallisto",
                                             txOut=txout, countsFromAbundance="lengthScaledTPM"))
    } else {
      import <- sm(tximport::tximport(files=files, type="kallisto", tx2gene=tx_gene_map, txOut=txout))
      import_scaled <- sm(tximport::tximport(files=files, type="kallisto", tx2gene=tx_gene_map,
                                             txOut=txout, countsFromAbundance="lengthScaledTPM"))
    }
    retlist[["count_table"]] <- data.table::as.data.table(import[["counts"]], keep.rownames="rownames")
    tt <- setkey(retlist[["count_table"]], rownames)
    retlist[["tximport"]] <- import
    retlist[["tximport_scaled"]] <- import_scaled
    retlist[["source"]] <- "tximport"
  } else if (grepl(pattern="\\.genes\\.results", x=files[1])) {
    names(files) <- ids
    import <- NULL
    import_scaled <- NULL
    if (is.null(tx_gene_map)) {
      import <- tximport::tximport(files=files, type="rsem", txOut=txout)
      import_scaled <- tximport::tximport(files=files, type="rsem",
                                          txOut=txout, countsFromAbundance="lengthScaledTPM")
    } else {
      import <- tximport::tximport(files=files, type="rsem", tx2gene=tx_gene_map, txOut=txout)
      import_scaled <- tximport::tximport(files=files, type="rsem", tx2gene=tx_gene_map,
                                          txOut=txout, countsFromAbundance="lengthScaledTPM")
    }
    retlist[["count_table"]] <- data.table::as.data.table(import[["counts"]], keep.rownames="rownames")
    retlist[["tximport"]] <- import
    retlist[["tximport_scaled"]] <- import_scaled
    retlist[["source"]] <- "tximport"
  } else {
    ## This is used when 'normal' htseq-based counts were generated.
    count_table <- try(read.table(files[1], header=header))
    colnames(count_table) <- c("rownames", ids[1])
    count_table <- data.table::as.data.table(count_table)
    tt <- data.table::setkey(count_table, rownames)
    if (class(count_table)[1] == "try-error") {
      stop(paste0("There was an error reading: ", files[1]))
    }
    message(paste0(files[1], " contains ", length(rownames(count_table)), " rows."))
        ## Following lines not needed for data.table
    ## rownames(count_table) <- make.names(count_table[, "ID"], unique=TRUE)
    ## count_table <- count_table[, -1, drop=FALSE]
    ## iterate over and append remaining samples
    for (table in 2:length(files)) {
      if (file.exists(tolower(files[table]))) {
        files[table] <- tolower(files[table])
      } else if (file.exists(lowhpgl_filenames[table])) {
        files[table] <- lowhpgl_filenames[table]
      } else if (file.exists(lower_filenames[table])) {
        files[table] <- lower_filenames[table]
      }
      tmp_count <- try(read.table(files[table], header=header))
      if (class(tmp_count)[1] == "try-error") {
        stop(paste0("There was an error reading: ", files[table]))
      }
      colnames(tmp_count) <- c("rownames", ids[table])
      tmp_count <- data.table::as.data.table(tmp_count)
      pre_merge <- nrow(tmp_count)
      count_table <- merge(count_table, tmp_count, by="rownames", all.x=TRUE)
      ## rownames(count_table) <- count_table[, "Row.names"]
      ## count_table <- count_table[, -1, drop=FALSE]
      ## post_merge <- length(rownames(count_table))
      post_merge <- nrow(count_table)
      message(paste0(files[table], " contains ", pre_merge,
                     " rows and merges to ", post_merge, " rows."))
    }
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
      tt <- setkey(retlist[["count_table"]], rownames)
  }
  return(retlist)
}

#' Given a table of meta data, read it in for use by create_expt().
#'
#' Reads an experimental design in a few different formats in preparation for creating an expt.
#'
#' @param file Csv/xls file to read.
#' @param ... Arguments for arglist, used by sep, header and similar read.csv/read.table parameters.
#' @return Df of metadata.
#' @seealso \pkg{tools} \pkg{openxlsx} \pkg{XLConnect}
read_metadata <- function(file, ...) {
  arglist <- list(...)
  if (is.null(arglist[["sep"]])) {
    arglist[["sep"]] <- ","
  }
  if (is.null(arglist[["header"]])) {
    arglist[["header"]] <- TRUE
  }

  if (tools::file_ext(file) == "csv") {
    definitions <- read.csv(file=file, comment.char="#",
                            sep=arglist[["sep"]], header=arglist[["header"]])
  } else if (tools::file_ext(file) == "xlsx") {
    ## xls = loadWorkbook(file, create=FALSE)
    ## tmp_definitions = readWorksheet(xls, 1)
    definitions <- try(openxlsx::read.xlsx(xlsxFile=file, sheet=1))
    if (class(definitions) == "try-error") {
      stop(paste0("Unable to read the metadata file: ", file))
    }
  } else if (tools::file_ext(file) == "xls") {
    ## This is not correct, but it is a start
    definitions <- readxl::read_xls(path=file, sheet=1)
  } else {
    definitions <- read.table(file=file, sep=arglist[["sep"]], header=arglist[["header"]])
  }
  return(definitions)
}

#' Change the batches of an expt.
#'
#' When exploring differential analyses, it might be useful to play with the conditions/batches of
#' the experiment.  Use this to make that easier.
#'
#' @param expt  Expt to modify.
#' @param fact  Batches to replace using this factor.
#' @param ids  Specific samples to change.
#' @param ...  Extra options are like spinach.
#' @return  The original expt with some new metadata.
#' @seealso \code{\link{create_expt}} \code{\link{set_expt_conditions}}
#' @examples
#' \dontrun{
#'  expt = set_expt_batches(big_expt, factor=c(some,stuff,here))
#' }
#' @export
set_expt_batches <- function(expt, fact, ids=NULL, ...) {
  arglist <- list(...)
  original_batches <- expt[["batches"]]
  original_length <- length(original_batches)
  if (length(fact) == 1) {
    ## Assume it is a column in the design
    if (fact %in% colnames(expt[["design"]])) {
      fact <- expt[["design"]][[fact]]
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
#' When exploring differential analyses, it might be useful to play with the conditions/batches of
#' the experiment.  Use this to make that easier.
#'
#' @param expt  Expt to modify
#' @param colors  colors to replace
#' @param chosen_palette  I usually use Dark2 as the RColorBrewer palette.
#' @param change_by  Assuming a list is passed, cross reference by condition or sample?
#' @return expt Send back the expt with some new metadata
#' @seealso \code{\link{set_expt_conditions}} \code{\link{set_expt_batches}}
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
#' esmer_expt <- set_expt_colors(expt=esmer_expt, colors=chosen_colors)
#' }
#' @export
set_expt_colors <- function(expt, colors=TRUE, chosen_palette="Dark2", change_by="condition") {
  num_conditions <- length(levels(as.factor(expt[["conditions"]])))
  num_samples <- nrow(expt[["design"]])
  sample_ids <- expt[["design"]][["sampleid"]]
  chosen_colors <- expt[["conditions"]]
  sample_colors <- NULL
  if (is.null(colors) | isTRUE(colors)) {
    sample_colors <- sm(grDevices::colorRampPalette(RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else if (class(colors) == "character") {
    if (is.null(names(colors))) {
      names(colors) <- levels(as.factor(expt[["conditions"]]))
    }
    if (change_by == "condition") {
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- colors
      chosen_colors <- mapping[chosen_colors]
    } else if (change_by == "sample") {
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- expt[["colors"]]
      for (snum in 1:length(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (class(colors) == "list") {
    if (change_by == "condition") {
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- as.character(colors)
      names(mapping) <- names(colors)
      chosen_colors <- mapping[chosen_colors]
    } else if (change_by == "sample") {
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- expt[["colors"]]
      for (snum in 1:length(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
        ## Set the condition for the changed samples to something unique.
        original_condition <- expt[["design"]][sampleid, "condition"]
        changed_condition <- paste0(original_condition, snum)
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
    colors <- sm(grDevices::colorRampPalette(RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    ## Check that all conditions are named in the color list:
    mapping <- setNames(colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else {
    warning("The number of colors provided does not match either the number of conditions nor samples.")
    warning("Unsure of what to do, so choosing colors with RColorBrewer.")
    sample_colors <- suppressWarnings(grDevices::colorRampPalette(
                                                   RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  }
  names(chosen_colors) <- sample_ids

  expt[["colors"]] <- chosen_colors
  return(expt)
}

#' Change the condition of an expt
#'
#' When exploring differential analyses, it might be useful to play with the conditions/batches of
#' the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify
#' @param fact Conditions to replace
#' @param ids Specific sample IDs to change.
#' @param ...  Extra arguments are given to arglist.
#' @return expt Send back the expt with some new metadata
#' @seealso \code{\link{set_expt_batches}} \code{\link{create_expt}}
#' @examples
#' \dontrun{
#'  expt = set_expt_conditions(big_expt, factor=c(some,stuff,here))
#' }
#' @export
set_expt_conditions <- function(expt, fact=NULL, ids=NULL, ...) {
  arglist <- list(...)
  original_conditions <- expt[["conditions"]]
  original_length <- length(original_conditions)
  new_expt <- expt  ## Explicitly copying expt to new_expt
  ## because when I run this as a function call() it seems to be not properly setting
  ## the conditions and I do not know why.
  if (!is.null(ids)) {
    ## Change specific id(s) to given condition(s).
    old_pdata <- pData(expt)
    old_cond <- as.character(old_pdata[["condition"]])
    names(old_cond) <- rownames(old_pdata)
    new_cond <- old_cond
    new_cond[ids] <- fact
    new_pdata <- old_pdata
    new_pdata[["condition"]] <- as.factor(new_cond)
    pData(expt[["expressionset"]]) <- new_pdata
    new_expt[["conditions"]][ids] <- fact
    new_expt[["design"]][["condition"]] <- new_cond
  } else if (length(fact) == 1) {
    ## Assume it is a column in the design
    if (fact %in% colnames(expt[["design"]])) {
      new_fact <- expt[["design"]][[fact]]
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

  tmp_expt <- set_expt_colors(new_expt)
  rm(new_expt)
  return(tmp_expt)
}

#' Change the factors (condition and batch) of an expt
#'
#' When exploring differential analyses, it might be useful to play with the conditions/batches of
#' the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify
#' @param condition New condition factor
#' @param batch New batch factor
#' @param ids Specific sample IDs to change.
#' @param ... Arguments passed along (likely colors)
#' @return expt Send back the expt with some new metadata
#' @seealso \code{\link{set_expt_conditions}} \code{\link{set_expt_batches}}
#' @examples
#' \dontrun{
#'  expt = set_expt_factors(big_expt, condition="column", batch="another_column")
#' }
#' @export
set_expt_factors <- function(expt, condition=NULL, batch=NULL, ids=NULL, ...) {
  arglist <- list(...)
  if (!is.null(condition)) {
    expt <- set_expt_conditions(expt, fact=condition, ...)
  }
  if (!is.null(batch)) {
    expt <- set_expt_batches(expt, fact=batch, ...)
  }
  return(expt)
}

#' Change the sample names of an expt.
#'
#' Sometimes one does not like the hpgl identifiers, so provide a way to change them on-the-fly.
#'
#' @param expt Expt to modify
#' @param newnames New names, currently only a character vector.
#' @return expt Send back the expt with some new metadata
#' @seealso \code{\link{set_expt_conditions}} \code{\link{set_expt_batches}}
#' @examples
#' \dontrun{
#'  expt = set_expt_samplenames(expt, c("a","b","c","d","e","f"))
#' }
#' @export
set_expt_samplenames <- function(expt, newnames) {
  new_expt <- expt
  oldnames <- rownames(new_expt[["design"]])
  newnames <- make.unique(newnames)
  newnote <- paste0("Sample names changed from: ", toString(oldnames),
                    " to: ", toString(newnames), " at: ", date(), ".\n")
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
  new_expt[["expressionset"]] <- new_expressionset
  names(new_expt[["libsize"]]) <- newnames
  new_expt[["samplenames"]] <- newnames
  return(new_expt)
}

#' Extract a subset of samples following some rule(s) from an
#' experiment class.
#'
#' Sometimes an experiment has too many parts to work with conveniently, this operation allows one
#' to break it into smaller pieces.
#'
#' @param expt Expt chosen to extract a subset of data.
#' @param subset Valid R expression which defines a subset of the design to keep.
#' @return metadata Expt class which contains the smaller set of data.
#' @seealso \pkg{Biobase}
#'  \code{\link[Biobase]{pData}} \code{\link[Biobase]{exprs}} \code{\link[Biobase]{fData}}
#' @examples
#' \dontrun{
#'  smaller_expt = expt_subset(big_expt, "condition=='control'")
#'  all_expt = expt_subset(expressionset, "")  ## extracts everything
#' }
#' @export
subset_expt <- function(expt, subset=NULL) {
  starting_expressionset <- NULL
  starting_metadata <- NULL
  if (class(expt)[[1]] == "ExpressionSet") {
    starting_expressionset <- expt
    starting_metadata <- pData(starting_expressionset)
  } else if (class(expt)[[1]] == "expt") {
    starting_expressionset <- expt[["expressionset"]]
    starting_metadata <- pData(expt)
  } else {
    stop("expt is neither an expt nor ExpressionSet")
  }

  note_appended <- NULL
  if (is.null(subset)) {
    subset_design <- starting_metadata
  } else {
    r_expression <- paste("subset(starting_metadata,", subset, ")")
    subset_design <- eval(parse(text=r_expression))
    ## design = data.frame(sample=samples$sample, condition=samples$condition, batch=samples$batch)
    note_appended <- paste0("Subsetted with ", subset, " on ", date(), ".\n")
  }
  if (nrow(subset_design) == 0) {
    stop("When the subset was taken, the resulting design has 0 members, check your expression.")
  }
  subset_design <- as.data.frame(subset_design, stringsAsFactors=FALSE)
  ## This is to get around stupidity with respect to needing all factors to be in a DESeqDataSet
  starting_ids <- rownames(starting_metadata)
  subset_ids <- rownames(subset_design)
  subset_positions <- starting_ids %in% subset_ids
  starting_colors <- expt[["colors"]]
  subset_colors <- starting_colors[subset_positions]
  starting_conditions <- expt[["conditions"]]
  subset_conditions <- starting_conditions[subset_positions, drop=TRUE]
  starting_batches <- expt[["batches"]]
  subset_batches <- starting_batches[subset_positions, drop=TRUE]
  original_libsize <- expt[["original_libsize"]]
  subset_original_libsize <- original_libsize[subset_positions, drop=TRUE]
  current_libsize <- expt[["libsize"]]
  subset_current_libsize <- current_libsize[subset_positions, drop=TRUE]
  subset_expressionset <- starting_expressionset[, subset_positions]
  original_expressionset <- expt[["original_expressionset"]]
  subset_original_expressionset <- original_expressionset[, subset_positions]

  notes <- expt[["notes"]]
  if (!is.null(note_appended)) {
    notes <- paste0(notes, note_appended)
  }

  for (col in 1:ncol(subset_design)) {
    if (class(subset_design[[col]]) == "factor") {
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
    "original_expressionset" = subset_original_expressionset,
    "original_libsize" = subset_original_libsize,
    "libsize" = subset_current_libsize)
  class(new_expt) <- "expt"
  return(new_expt)
}

#' Print a string describing what happened to this data.
#'
#' Sometimes it is nice to have a string like: log2(cpm(data)) describing what happened to the data.
#'
#' @param expt  The expressionset.
#' @param transform  How was it transformed?
#' @param convert  How was it converted?
#' @param norm  How was it normalized?
#' @param filter  How was it filtered?
#' @param batch  How was it batch-corrected?
#' @return An expression describing what has been done to this data.
#' @seealso \code{\link{create_expt}}
#' @export
what_happened <- function(expt=NULL, transform="raw", convert="raw",
                          norm="raw", filter="raw", batch="raw") {
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
    what <- paste0(what, transform, '(')
  }
  if (batch != "raw") {
    if (isTRUE(batch)) {
      what <- paste0(what, 'batch-correct(')
    } else {
      what <- paste0(what, batch, '(')
    }
  }
  if (convert != "raw") {
    what <- paste0(what, convert, '(')
  }
  if (norm != "raw") {
    what <- paste0(what, norm, '(')
  }
  if (filter != "raw") {
    what <- paste0(what, filter, '(')
  }
  what <- paste0(what, 'data')
  if (transform != 'raw') {
    what <- paste0(what, ')')
  }
  if (batch != "raw") {
    what <- paste0(what, ')')
  }
  if (convert != "raw") {
    what <- paste0(what, ')')
  }
  if (norm != "raw") {
    what <- paste0(what, ')')
  }
  if (filter != "raw") {
    what <- paste0(what, ')')
  }

  return(what)
}

#' Make pretty xlsx files of count data.
#'
#' Some folks love excel for looking at this data.  ok.
#'
#' Tested in test_03graph_metrics.R
#' This performs the following:  Writes the raw data, graphs the raw data, normalizes the data,
#' writes it, graphs it, and does a median-by-condition and prints that.  I replaced the openxlsx
#' function which writes images into xlsx files with one which does not require an opening of a
#' pre-existing plotter.  Instead it (optionally)opens a pdf device, prints the plot to it, opens a
#' png device, prints to that, and inserts the resulting png file.  Thus it sacrifices some
#' flexibility for a hopefully more consistent behaivor.  In addition, one may use the pdfs as
#' a set of images importable into illustrator or whatever.
#'
#' @param expt  An expressionset to print.
#' @param excel  Filename to write.
#' @param norm  Normalization to perform.
#' @param violin  Include violin plots?
#' @param convert  Conversion to perform.
#' @param transform  Transformation used.
#' @param batch  Batch correction applied.
#' @param filter  Filtering method used.
#' @return  A big honking excel file and a list including the dataframes and images created.
#' @seealso \pkg{openxlsx} \pkg{Biobase}
#'  \code{\link{normalize_expt}} \code{\link{graph_metrics}}
#' @examples
#' \dontrun{
#'  excel_sucks <- write_expt(expt)
#' }
#' @export
write_expt <- function(expt, excel="excel/pretty_counts.xlsx", norm="quant", violin=FALSE,
                       convert="cpm", transform="log2", batch="sva", filter="cbcb") {
  wb <- openxlsx::createWorkbook(creator="hpgltools")
  plot_dim <- 6
  plot_cols <- floor(plot_dim * 1.5)
  plot_rows <- ceiling(plot_dim * 5.0)
  new_row <- 1
  new_col <- 1
  excel_basename <- gsub(pattern="\\.xlsx", replacement="", x=excel)

  ## Write an introduction to this foolishness.
  message("Writing the legend.")
  sheet <- "legend"
  norm_state <- paste0(transform, "(", convert, "(", norm, "(", batch, "(", filter, "(counts)))))")
  legend <- data.frame(
    "sheet" = c("1.", "2.", "3.", "4.", "5.", "6."),
    "sheet_definition" = c("This sheet, including the experimental design.",
                           "The raw counts and annotation data on worksheet 'raw_data'.",
                           "Some graphs describing the distribution of raw data in worksheet 'raw_plots'.",
                           paste0("The counts normalized with: ", norm_state),
                           "Some graphs describing the distribution of the normalized data on 'norm_plots'.",
                           "The median normalized counts by condition factor on 'median_data'."),
    stringsAsFactors=FALSE)
  colnames(legend) <- c("Worksheets", "Contents")
  xls_result <- write_xls(wb, data=legend, sheet=sheet, rownames=FALSE,
                          title="Columns used in the following tables.")
  rows_down <- nrow(legend)
  new_row <- new_row + rows_down + 3
  annot <- as.data.frame(pData(expt), stringsAsFactors=FALSE)
  xls_result <- write_xls(data=annot, wb=wb, start_row=new_row, rownames=FALSE,
                          sheet=sheet, start_col=1, title="Experimental Design.")

  ## Write the raw read data and gene annotations
  message("Writing the raw reads.")
  sheet <- "raw_reads"
  new_row <- 1
  new_col <- 1
  reads <- exprs(expt)
  info <- fData(expt)
  read_info <- merge(info, reads, by="row.names")
  xls_result <- write_xls(data=read_info, wb=wb, sheet=sheet, rownames=FALSE,
                          start_row=new_row, start_col=new_col, title="Raw Reads.")

  ## Write some graphs for the raw data
  message("Graphing the raw reads.")
  sheet <- "raw_graphs"
  newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet))
  if (class(newsheet) == "try-error") {
    warning(paste0("Failed to add the sheet: ", sheet))
  }
  metrics <- sm(graph_metrics(expt, qq=TRUE))
  ## Start with library sizes.
  openxlsx::writeData(wb, sheet=sheet, x="Legend.", startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw library sizes.", startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Non-zero genes.", startRow=new_row, startCol=new_col)
  new_row <- new_row + 1
  new_col <- 1
  legend_plot <- metrics[["legend"]][["plot"]]
  try_result <- xlsx_plot_png(legend_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="01_legend", savedir=excel_basename, fancy_type="svg")
  new_col <- new_col + plot_cols + 1
  libsize_plot <- metrics[["libsize"]]
  try_result <- xlsx_plot_png(libsize_plot, wb=wb, sheet=sheet, width=plot_dim, height=plot_dim,
                              start_col=new_col, start_row=new_row,
                              plotname="02_libsize", savedir=excel_basename)
  ## Same row, non-zero plot
  new_col <- new_col + plot_cols + 1
  nonzero_plot <- metrics[["nonzero"]]
  try_result <- xlsx_plot_png(nonzero_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="03_nonzero", savedir=excel_basename)
  new_col <- new_col + plot_cols + 1

  ## Visualize distributions
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw data density plot.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw Boxplot.",
                      startRow=new_row, startCol=new_col)
  new_col <- 1
  density_plot <- metrics[["density"]]
  new_row <- new_row + 1
  try_result <- xlsx_plot_png(density_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="04_density", savedir=excel_basename)
  new_col <- new_col + plot_cols + 1
  boxplot_plot <- metrics[["boxplot"]]
  try_result <- xlsx_plot_png(boxplot_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="05_boxplot", savedir=excel_basename)
  new_col <- 1

  ## Move down next set of rows, heatmaps
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw correlation heatmap.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw distance heatmap.",
                      startRow=new_row, startCol=new_col)
  new_col <- 1
  new_row <- new_row + 1
  corheat_plot <- metrics[["corheat"]]
  try_result <- xlsx_plot_png(corheat_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="06_corheat", savedir=excel_basename)
  disheat_plot <- metrics[["disheat"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(disheat_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="07_disheat", savedir=excel_basename)
  new_col <- 1

  ## SM plots
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw standard median correlation.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw standard distance correlation.",
                      startRow=new_row, startCol=new_col)
  new_col <- 1
  new_row <- new_row + 1
  smc_plot <- metrics[["smc"]]
  try_result <- xlsx_plot_png(smc_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="08_smc", savedir=excel_basename, fancy_type="svg")
  new_col <- new_col + plot_cols + 1
  smd_plot <- metrics[["smd"]]
  try_result <- xlsx_plot_png(smd_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="09_smd", savedir=excel_basename, fancy_type="svg")
  new_col <- 1

  ## PCA, PCA(l2cpm) and qq_log
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw PCA.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="PCA(log2(cpm())).",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw TSNE.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="TSNE(log2(cpm())).",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw QQ, log scale.",
                      startRow=new_row, startCol=new_col)
  new_col <- 1
  new_row <- new_row + 1
  pca_plot <- metrics[["pcaplot"]]
  pca_table <- metrics[["pcatable"]]
  tsne_plot <- metrics[["tsneplot"]]
  tsne_table <- metrics[["tsnetable"]]
  try_result <- xlsx_plot_png(pca_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="10_pcaplot", savedir=excel_basename, fancy_type="svg")
  tmp_data <- sm(normalize_expt(expt, transform="log2", convert="cpm"))
  rpca <- plot_pca(tmp_data)
  rtsne <- plot_tsne(tmp_data)
  rspca_plot <- rpca[["plot"]]
  rtsne_plot <- rtsne[["plot"]]
  rpca_table <- rpca[["table"]]
  rtsne_table <- rtsne[["table"]]
  rm(tmp_data)
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(rspca_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="11_norm_pcaplot", savedir=excel_basename, fancy_type="svg")
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(tsne_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="12_tsneplot", savedir=excel_basename, fancy_type="svg")
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(rtsne_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="13_rtsneplot", savedir=excel_basename, fancy_type="svg")
  qq_plot <- metrics[["qqlog"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(qq_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="14_qqlog", savedir=excel_basename)
  new_col <- 1

  violin_plot <- NULL
  pct_plot <- NULL
  ## Violin plots
  if (isTRUE(violin)) {
    varpart_raw <- try(varpart(expt, predictor=NULL, factors=c("condition", "batch")))
    if (class(varpart_raw) != "try-error") {
      violin_plot <- varpart_raw[["partition_plot"]]
      new_row <- new_row + plot_rows + 2
      new_col <- 1
      try_result <- xlsx_plot_png(violin_plot, wb=wb, sheet=sheet, width=plot_dim,
                                  height=plot_dim, start_col=new_col, start_row=new_row,
                                  plotname="15_violin", savedir=excel_basename)
      new_col <- new_col + plot_cols + 1

      pct_plot <- varpart_raw[["percent_plot"]]
      try_result <- xlsx_plot_png(pct_plot, wb=wb, sheet=sheet, width=plot_dim,
                                  height=plot_dim, start_col=new_col, start_row=new_row,
                                  plotname="16_pctvar", savedir=excel_basename)
    }
  }

  ## PCA table
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Raw PCA res.",
                      startRow=new_row, startCol=new_col)
  new_row <- new_row + 1
  xls_result <- write_xls(data=metrics[["pcares"]], wb=wb, rownames=FALSE,
                          sheet=sheet, start_col=new_col, start_row=new_row)
  new_col <- xls_result[["end_col"]] + 6
  new_row <- new_row - 1
  openxlsx::writeData(wb, sheet, "Raw PCA table.",
                      startRow=new_row, startCol=new_col)
  new_row <- new_row + 1
  xls_result <- write_xls(data=metrics[["pcatable"]], wb=wb, rownames=FALSE,
                          sheet=sheet, start_row=new_row, start_col=new_col)

  ## Move on to the next sheet, normalized data
  message("Writing the normalized reads.")
  sheet <- "norm_data"
  new_col <- 1
  new_row <- 1
  norm_data <- sm(normalize_expt(expt=expt, transform=transform, norm=norm,
                                 convert=convert, batch=batch, filter=filter))
  norm_reads <- exprs(norm_data)
  info <- fData(norm_data)
  read_info <- merge(norm_reads, info, by="row.names")
  title <- what_happened(norm_data)
  xls_result <- write_xls(wb=wb, data=read_info, rownames=FALSE,
                          start_row=new_row, start_col=new_col, sheet=sheet, title=title)

  ## Graphs of the normalized data
  message("Graphing the normalized reads.")
  sheet <- "norm_graphs"
  newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet))
  norm_metrics <- sm(graph_metrics(norm_data, qq=TRUE))
  ## Start with library sizes.
  openxlsx::writeData(wb, sheet, "Legend.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet, "Normalized library sizes.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Non-zero genes.",
                      startRow=new_row, startCol=new_col)
  new_col <- 1
  new_row <- new_row + 1
  new_plot <- norm_metrics[["legend"]][["plot"]]
  try_result <- xlsx_plot_png(new_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row)
  new_col <- new_col + plot_cols + 1
  nlibsize_plot <- norm_metrics[["libsize"]]
  try_result <- xlsx_plot_png(nlibsize_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="17_nlibsize", savedir=excel_basename)
  ## Same row, non-zero plot
  new_col <- new_col + plot_cols + 1
  nnzero_plot <- norm_metrics[["nonzero"]]
  try_result <- xlsx_plot_png(nnzero_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="18_nnzero", savedir=excel_basename)
  new_col <- new_col + plot_cols + 1

  ## Visualize distributions
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized data density plot.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized Boxplot.",
                      startRow=new_row, startCol=new_col)
  new_col <- 1
  ndensity_plot <- norm_metrics[["density"]]
  new_row <- new_row + 1
  try_result <- xlsx_plot_png(ndensity_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="19_ndensity", savedir=excel_basename)
  nboxplot_plot <- norm_metrics[["boxplot"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(nboxplot_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="20_nboxplot", savedir=excel_basename)
  new_col <- 1

  ## Move down next set of rows, heatmaps
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized correlation heatmap.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized distance heatmap.",
                      startRow=new_row, startCol=new_col)
  new_col <- 1
  ncorheat_plot <- norm_metrics[["corheat"]]
  new_row <- new_row + 1
  try_result <- xlsx_plot_png(ncorheat_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="21_ncorheat", savedir=excel_basename)
  ndisheat_plot <- norm_metrics[["disheat"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(ndisheat_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="22_ndisheat", savedir=excel_basename)
  new_col <- 1

  ## SM plots
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized standard median correlation.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized standard distance correlation.",
                      startRow=new_row, startCol=new_col)
  new_col <- 1
  nsmc_plot <- norm_metrics[["smc"]]
  new_row <- new_row + 1
  try_result <- xlsx_plot_png(nsmc_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="23_nsmc", savedir=excel_basename, fancy_type="svg")
  nsmd_plot <- norm_metrics[["smd"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(nsmd_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="24_nsmd", savedir=excel_basename, fancy_type="svg")
  new_col <- 1

  ## PCA and qq_log
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized PCA.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized TSNE.",
                      startRow=new_row, startCol=new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized QQ, log scale.",
                      startRow=new_row, startCol=new_col)
  new_col <- 1
  npca_plot <- norm_metrics[["pcaplot"]]
  ntsne_plot <- norm_metrics[["tsneplot"]]
  npca_table <- norm_metrics[["pcatable"]]
  ntsne_table <- norm_metrics[["tsnetable"]]
  new_row <- new_row + 1
  try_result <- xlsx_plot_png(npca_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="25_npcaplot", savedir=excel_basename, fancy_type="svg")
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(ntsne_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="26_ntsneplot", savedir=excel_basename, fancy_type="svg")
  nqq_plot <- norm_metrics[["qqlog"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_plot_png(nqq_plot, wb=wb, sheet=sheet, width=plot_dim,
                              height=plot_dim, start_col=new_col, start_row=new_row,
                              plotname="27_nqqplot", savedir=excel_basename)

  new_col <- 1

  ## Violin plots
  nvarpart_plot <- NULL
  npct_plot <- NULL
  if (isTRUE(violin)) {
    varpart_norm <- try(varpart(norm_data, predictor=NULL, factors=c("condition", "batch")))
    if (class(varpart_norm) != "try-error") {
      nvarpart_plot <- varpart_norm[["partition_plot"]]
      new_row <- new_row + plot_rows + 2
      new_col <- 1
      try_result <- xlsx_plot_png(nvarpart_plot, wb=wb, sheet=sheet, width=plot_dim,
                                  height=plot_dim, start_col=new_col, start_row=new_row,
                                  plotname="28_nviolin", savedir=excel_basename)
      new_col <- new_col + plot_cols + 1
      npct_plot <- varpart_norm[["percent_plot"]]
      try_result <- xlsx_plot_png(npct_plot, wb=wb, sheet=sheet, width=plot_dim,
                                  height=plot_dim, start_col=new_col, start_row=new_row,
                                  plotname="29_npctplot", savedir=excel_basename)
    }
  }

  ## PCA table
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized PCA res.", startRow=new_row, startCol=new_col)
  new_row <- new_row + 1
  xls_result <- write_xls(data=norm_metrics[["pcares"]], wb=wb, rownames=FALSE,
                          sheet=sheet, start_col=new_col, start_row=new_row)
  new_col <- xls_result[["end_col"]] + 6
  new_row <- new_row - 1
  openxlsx::writeData(wb, sheet=sheet, x="Normalized PCA table.", startRow=new_row, startCol=new_col)
  new_row <- new_row + 1
  xls_result <- write_xls(data=norm_metrics[["pcatable"]], wb=wb, sheet=sheet,
                          rownames=FALSE, start_col=new_col, start_row=new_row)


  ## Give a median-by-factor accounting of the data
  message("Writing the median reads by factor.")
  sheet <- "median_data"
  new_col <- 1
  new_row <- 1
  median_data <- median_by_factor(exprs(norm_data),
                                  norm_data[["conditions"]])
  median_data_merged <- merge(median_data, info, by="row.names")
  xls_result <- write_xls(wb, data=median_data_merged, start_row=new_row, start_col=new_col,
                          rownames=FALSE, sheet=sheet, title="Median Reads by factor.")

  ## Save the result
  save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite=TRUE))
  retlist <- list(
    "save" = save_result,
    "legend" = legend,
    "annotations" = annot,
    "raw_reads" = reads,
    "design" = info,
    "legend" = legend_plot,
    "raw_libsize" = libsize_plot,
    "raw_nonzero" = nonzero_plot,
    "raw_density" = density_plot,
    "raw_boxplot" = boxplot_plot,
    "raw_corheat" = corheat_plot,
    "raw_disheat" = disheat_plot,
    "raw_smc" = smc_plot,
    "raw_smd" = smd_plot,
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
    "norm_pca_table" = npca_table,
    "norm_tsne" = ntsne_plot,
    "norm_tsne_table" = ntsne_table,
    "norm_qq" = nqq_plot,
    "norm_violin" = nvarpart_plot,
    "norm_pct" = npct_plot,
    "medians" = median_data
  )
  return(retlist)
}

## Make some methods which call Biobase on expts.
## This way we don't have to do Biobase::exprs/fData/pData/notes()
## But instead R will call them for us.

## Here is a note from Hadley which is relevant here:
## Another consideration is that S4 code often needs to run in a certain
## order. For example, to define the method setMethod("foo", c("bar", "baz"),
## ...) you must already have created the foo generic and the two classes. By
## default, R code is loaded in alphabetical order, but that won’t always work
## for your situation.

#' Extend Biobase::exprs to handle expt ojects.
#'
#' @name exprs
#' @aliases exprs, exprs-methods
#' @param object  The expt object from which to extract the expressionset.
#' @importFrom Biobase exprs
#' @docType methods
#' @rdname exprs-methods
#' @export exprs
## setOldClass("expt")
setMethod("exprs", signature="expt",
          function(object) {
            Biobase::exprs(object[["expressionset"]])
          })

#' Extend Biobase::fData to handle expt objects.
#'
#' @name fData
#' @aliases fData, fData-methods
#' @param object  An expt from which to extract the expressionset.
#' @importFrom Biobase fData
#' @docType methods
#' @rdname fData-methods
#' @export fData
setMethod("fData", signature="expt",
          function(object) {
            Biobase::fData(object[["expressionset"]])
          })

#' Extend Biobase::pData to handle expt objects.
#'
#' @name pData
#' @aliases pData, pData-methods
#' @param object  The expt object from which to extract the expressionset.
#' @importFrom Biobase pData
#' @docType methods
#' @rdname pData-methods
#' @export pData
setMethod("pData", signature="expt",
          function(object) {
            Biobase::pData(object[["expressionset"]])
          })

#' Extend Biobase::sampleNames to handle expt objects.
#'
#' @name sampleNames
#' @aliases sampleNames, sampleNames-methods
#' @param object  The expt object from which to extract the expressionset.
#' @importFrom Biobase sampleNames
#' @docType methods
#' @rdname sampleNames-methods
#' @export sampleNames
setMethod("sampleNames", signature="expt",
          function(object) {
            Biobase::sampleNames(object[["expressionset"]])
          })

#' Extend Biobase::sampleNames<- to handle expt objects.
#'
#' @name sampleNames<-
#' @aliases sampleNames<-, sampleNames<--methods
#' @param object  The expt object from which to extract the expressionset.
#' @importFrom Biobase sampleNames<-
#' @docType methods
#' @rdname sampleNamesto-methods
#' @export sampleNames<-
setMethod("sampleNames<-", signature="expt",
          function(object) {
            set_expt_samplenames(object)
          })

#' Extend Biobase::notes to handle expt objects.
#'
#' @name notes
#' @aliases notes, notes-methods
#' @param object  The expt object from which to extract the expressionset.
#' @importFrom Biobase notes
#' @docType methods
#' @rdname notes-methods
#' @export notes
setMethod("notes", signature="expt",
          function(object) {
            Biobase::notes(object[["expressionset"]])
          })
## EOF
