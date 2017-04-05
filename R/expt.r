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
#'  \code{\link{expt_read_counts}} \code{\link[hash]{as.list.hash}}
#' @examples
#' \dontrun{
#'  new_experiment = create_expt("some_csv_file.csv", color_hash)
#'  ## Remember that this depends on an existing data structure of gene annotations.
#' }
#' @export
create_expt <- function(metadata, gene_info=NULL, count_dataframe=NULL,
                        sample_colors=NULL, title=NULL, notes=NULL,
                        include_type="all", include_gff=NULL,
                        savefile="expt", low_files=FALSE, ...) {
    arglist <- list(...)  ## pass stuff like sep=, header=, etc here
    ## Palette for colors when auto-chosen
    chosen_palette <- "Dark2"
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
    round <- FALSE
    if (!is.null(arglist[["round"]])) {
        round <- arglist[["round"]]
    }

    ## Read in the metadata from the provided data frame, csv, or xlsx.
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
    }

    colnames(sample_definitions) <- tolower(colnames(sample_definitions))
    colnames(sample_definitions) <- gsub(pattern="[[:punct:]]", replacement="",
                                         x=colnames(sample_definitions))
    ## Check that condition and batch have been filled in.
    sample_columns <- colnames(sample_definitions)
    sample_column <- NULL
    ## The sample ID column should have the word 'sample' in it, otherwise this will fail.
    found_sample <- grepl(pattern="sample", x=sample_columns)
    if (sum(found_sample) == 0) {
        message("Did not find the sample column in the sample sheet.")
        message("Was it perhaps saved as a .xls?")
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
        message("Was it perhaps saved as a .xls?")
    }
    found_batch <- "batch" %in% sample_columns
    if (!isTRUE(found_batch)) {
        message("Did not find the batch column in the sample sheet.")
        message("Was it perhaps saved as a .xls?")
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

    ## Create a matrix of counts with columns as samples and rows as genes
    ## This may come from either a data frame/matrix, a list of files from the metadata
    ## or it can attempt to figure out the location of the files from the sample names.
    filenames <- NULL
    all_count_tables <- NULL
    if (!is.null(count_dataframe)) {
        all_count_tables <- count_dataframe
        testthat::expect_equal(colnames(all_count_tables), rownames(sample_definitions))
        ## If neither of these cases is true, start looking for the files in the
        ## processed_data/ directory
    } else if (is.null(sample_definitions[[file_column]])) {
        success <- 0
        ## There are two main organization schemes I have used in the past, the following
        ## checks for both in case I forgot to put a file column in the metadata.
        ## Look for files organized by sample
        test_filenames <- paste0("processed_data/count_tables/",
                                 as.character(sample_definitions[[sample_column]]), "/",
                                 file_prefix,
                                 as.character(sample_definitions[[sample_column]]),
                                 file_suffix)
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
            test_filenames <- paste0("processed_data/count_tables/",
                                     tolower(as.character(sample_definitions[["type"]])), "/",
                                     tolower(as.character(sample_definitions[["stage"]])), "/",
                                     sample_definitions[[sample_column]], file_suffix)
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
    if (is.null(all_count_tables)) {
        filenames <- as.character(sample_definitions[[file_column]])
        sample_ids <- as.character(sample_definitions[[sample_column]])
        all_count_tables <- expt_read_counts(sample_ids, filenames, ...)
    }

    ## Recast the data as a data.frame and make sure everything is numeric
    all_count_tables <- as.data.frame(all_count_tables)
    for (col in colnames(all_count_tables)) {
        ## Ensure there are no stupid entries like target_id est_counts
        all_count_tables[[col]] <- as.numeric(all_count_tables[[col]])
    }
    ## I have had a couple data sets with incomplete counts, get rid of those rows before moving on.
    all_count_tables <- all_count_tables[complete.cases(all_count_tables), ]
    ## Features like exon:alicethegene-1 are annoying and entirely too common in TriTrypDB data
    rownames(all_count_tables) <- gsub(pattern="^exon:", replacement="",
                                       x=rownames(all_count_tables))
    rownames(all_count_tables) <- make.names(gsub(pattern=":\\d+", replacement="",
                                                  x=rownames(all_count_tables)), unique=TRUE)

    ## Try a couple different ways of getting gene-level annotations into the expressionset.
    annotation <- NULL
    tooltip_data <- NULL
    if (is.null(gene_info)) {
        ## Including, if all else fails, just grabbing the gene names from the count tables.
        if (is.null(include_gff)) {
            gene_info <- as.data.frame(rownames(all_count_tables))
            rownames(gene_info) <- rownames(all_count_tables)
            colnames(gene_info) <- "name"
        } else {
            ## Or reading a gff file.
            message("create_expt(): Reading annotation gff, this is slow.")
            annotation <- gff2df(gff=include_gff, type=gff_type)
            tooltip_data <- make_tooltips(annotations=annotation, type=gff_type, ...)
            gene_info <- annotation
        }
    } else if (class(gene_info) == "list" & !is.null(gene_info[["genes"]])) {
        ## In this case, it is using the output of reading a OrgDB instance
        gene_info <- as.data.frame(gene_info[["genes"]])
    }

    ## It turns out that loading the annotation information from orgdb/etc may not set the
    ## row names. Perhaps I should do that there, but I will add a check here, too.
    if (sum(rownames(gene_info) %in% rownames(all_count_tables)) == 0) {
        if (!is.null(gene_info[["geneid"]])) {
            rownames(gene_info) <- gene_info[["geneid"]]
        }
        if (sum(rownames(gene_info) %in% rownames(all_count_tables)) == 0) {
            warning("Even after changing the rownames in gene info, they do not match the count table.")
        }
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
    tmp_countsdt <- data.table::as.data.table(all_count_tables, keep.rownames="rownames")
    ##tmp_countsdt[["rownames"]] <- rownames(all_count_tables)
    ## This temporary id number will be used to ensure that the order of features in everything
    ## will remain consistent, as we will call order() using it later.
    tmp_countsdt[["temporary_id_number"]] <- 1:nrow(tmp_countsdt)
    gene_infodt <- data.table::as.data.table(gene_info, keep.rownames="rownames")
    ##gene_infodt[["rownames"]] <- rownames(gene_info)

    message("Bringing together the count matrix and gene information.")
    ## The method here is to create a data.table of the counts and annotation data,
    ## merge them, then split them apart.
    counts_and_annotations <- merge(tmp_countsdt, gene_infodt, by="rownames", all.x=TRUE)
    counts_and_annotations <- counts_and_annotations[order(counts_and_annotations[["temporary_id_number"]]), ]
    counts_and_annotations <- as.data.frame(counts_and_annotations)
    final_annotations <- counts_and_annotations[, colnames(counts_and_annotations) %in% colnames(gene_infodt) ]
    final_annotations <- final_annotations[, -1, drop=FALSE]
    ##colnames(final_annotations) <- colnames(gene_info)
    ##rownames(final_annotations) <- counts_and_annotations[["rownames"]]
    final_countsdt <- counts_and_annotations[, colnames(counts_and_annotations) %in% colnames(all_count_tables) ]
    final_counts <- as.data.frame(final_countsdt)
    rownames(final_counts) <- counts_and_annotations[["rownames"]]

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
    rm(counts_and_annotations)
    rm(tmp_countsdt)
    rm(gene_infodt)
    rm(final_countsdt)

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
                               exprs=as.matrix(final_counts),
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
    expt <- expt_subset(experiment) ## I think this is spurious now.
    expt[["original_expressionset"]] <- experiment
    expt[["original_metadata"]] <- Biobase::pData(experiment)

    ## I only leared fairly recently that there is quite a bit of redundancy between my expt
    ## and ExpressionSets. I do not mind this. Yet.
    expt[["title"]] <- title
    expt[["notes"]] <- toString(notes)
    expt[["design"]] <- sample_definitions
    expt[["annotation"]] <- annotation
    expt[["gff_file"]] <- include_gff
    expt[["tooltip"]] <- tooltip_data
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
    ## Save an rdata file of the expressionset.
    if (!is.null(savefile)) {
        save_result <- try(save(list = c("expt"), file=paste(savefile, ".Rdata", sep="")))
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
#' @param patterns  Character list of patterns to remove/keep
#' @param ...  Extra arguments are passed to arglist, currently unused.
#' @return  A smaller expt
#' @seealso \code{\link{create_expt}}
#' @export
expt_exclude_genes <- function(expt, column="txtype", method="remove",
                               patterns=c("snRNA","tRNA","rRNA"), ...) {
    arglist <- list(...)
    ex <- expt[["expressionset"]]
    annotations <- Biobase::fData(ex)
    pattern_string <- ""
    for (pat in patterns) {
        pattern_string <- paste0(pattern_string, pat, "|")
    }
    silly_string <- gsub(pattern="\\|$", replacement="", x=pattern_string)
    idx <- grepl(pattern=silly_string, x=annotations[[column]])
    ex2 <- NULL
    if (method == "remove") {
        ex2 <- ex[!idx, ]
    } else {
        ex2 <- ex[idx, ]
    }
    message(paste0("Before removal, there were ", nrow(Biobase::fData(ex)), " entries."))
    message(paste0("Now there are ", nrow(Biobase::fData(ex2)), " entries."))
    expt[["expressionset"]] <- ex2
    return(expt)
}

#' An alias to expt_subset, because it is stupid to have something start with verbs
#' and others start with nouns.
#'
#' This just calls expt_subset.
#'
#' @param ...  All arguments are passed to expt_subset.
#' @export
subset_expt <- function(...) {
    expt_subset(...)
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
expt_subset <- function(expt, subset=NULL) {
    starting_expressionset <- NULL
    starting_metadata <- NULL
    if (class(expt)[[1]] == "ExpressionSet") {
        starting_expressionset <- expt
        starting_metadata <- Biobase::pData(starting_expressionset)
    } else if (class(expt)[[1]] == "expt") {
        starting_expressionset <- expt[["expressionset"]]
        starting_metadata <- Biobase::pData(expt[["expressionset"]])
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
    subset_design <- as.data.frame(subset_design)
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
    Biobase::pData(subset_expressionset) <- subset_design

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
        definitions <- XLConnect::read.xls(xlsFile=file, sheet=1)
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
#' @seealso \code{\link{create_expt}} \code{\link{set_expt_condition}}
#' @examples
#' \dontrun{
#'  expt = set_expt_batch(big_expt, factor=c(some,stuff,here))
#' }
#' @export
set_expt_batch <- function(expt, fact, ids=NULL, ...) {
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
    Biobase::pData(expt[["expressionset"]])[["batch"]] <- fact
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
#' @seealso \code{\link{set_expt_condition}} \code{\link{set_expt_batch}}
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
#' @seealso \code{\link{set_expt_batch}} \code{\link{create_expt}}
#' @examples
#' \dontrun{
#'  expt = set_expt_condition(big_expt, factor=c(some,stuff,here))
#' }
#' @export
set_expt_condition <- function(expt, fact=NULL, ids=NULL, ...) {
    arglist <- list(...)
    original_conditions <- expt[["conditions"]]
    original_length <- length(original_conditions)
    new_expt <- expt  ## Explicitly copying expt to new_expt
    ## because when I run this as a function call() it seems to be not properly setting
    ## the conditions and I do not know why.
    if (!is.null(ids)) {
        ## Change specific id(s) to given condition(s).
        old_pdata <- Biobase::pData(expt[["expressionset"]])
        old_cond <- as.character(old_pdata[["condition"]])
        names(old_cond) <- rownames(old_pdata)
        new_cond <- old_cond
        new_cond[ids] <- fact
        new_pdata <- old_pdata
        new_pdata[["condition"]] <- as.factor(new_cond)
        Biobase::pData(expt[["expressionset"]]) <- new_pdata
        new_expt[["conditions"]][ids] <- fact
        new_expt[["design"]][["condition"]] <- new_cond
    } else if (length(fact) == 1) {
        ## Assume it is a column in the design
        if (fact %in% colnames(expt[["design"]])) {
            new_fact <- expt[["design"]][[fact]]
            new_expt[["conditions"]] <- new_fact
            Biobase::pData(new_expt[["expressionset"]])[["condition"]] <- new_fact
            new_expt[["design"]][["condition"]] <- new_fact
        } else {
            stop("The provided factor is not in the design matrix.")
        }
    } else if (length(fact) != original_length) {
            stop("The new factor of conditions is not the same length as the original.")
    } else {
        new_expt[["conditions"]] <- fact
        Biobase::pData(new_expt[["expressionset"]])[["condition"]] <- fact
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
#' @seealso \code{\link{set_expt_condition}} \code{\link{set_expt_batch}}
#' @examples
#' \dontrun{
#'  expt = set_expt_factors(big_expt, condition="column", batch="another_column")
#' }
#' @export
set_expt_factors <- function(expt, condition=NULL, batch=NULL, ids=NULL, ...) {
    arglist <- list(...)
    if (!is.null(condition)) {
        expt <- set_expt_condition(expt, fact=condition, ...)
    }
    if (!is.null(batch)) {
        expt <- set_expt_batch(expt, fact=batch, ...)
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
#' @seealso \code{\link{set_expt_condition}} \code{\link{set_expt_batch}}
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

## EOF
