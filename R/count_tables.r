## Time-stamp: <Fri Mar 25 13:38:44 2016 Ashton Trey Belew (abelew@gmail.com)>


#' Given a table of meta data, read it in for use by create_expt()
#'
#' @param file a csv/xls file to read
#' @param header does the table have a header (usually for csv)
#' @param sep separator for csv files
#' @return a df of metadata
read_metadata <- function(file, header=FALSE, sep=",") {
    if (tools::file_ext(file) == 'csv') {
        definitions <- read.csv(file=file, comment.char="#", sep=sep)
    } else if (tools::file_ext(file) == 'xlsx') {
        ## xls = loadWorkbook(file, create=FALSE)
        ## tmp_definitions = readWorksheet(xls, 1)
        definitions <- openxlsx::read.xlsx(xlsxFile=file, sheet=1)
    } else if (tools::file_ext(file) == 'xls') {
        ## This is not correct, but it is a start
        definitions <- XLConnect::read.xls(xlsFile=file, sheet=1)
    } else {
        definitions <- read.table(file=file)
    }

    colnames(definitions) <- tolower(colnames(definitions))
    rownames(definitions) <- make.names(definitions[["sample.id"]], unique=TRUE)
    ## "no visible binding for global variable 'sample.id'"  ## hmm sample.id is a column from the csv file.
    ## tmp_definitions <- subset(tmp_definitions, sample.id != "")
    empty_samples <- which(definitions$sample.id == "" | is.na(definitions$sample.id) | grepl(pattern="^#", x=definitions$sample.id))
    if (length(empty_samples) > 0) {
        definitions <- definitions[-empty_samples, ]
    }
    return(definitions)
}


#' Wrap bioconductor's expressionset to include some other extraneous
#' information.  This simply calls create_experiment and then does
#' expt_subset for everything
#'
#' this is relevant because the ceph object storage by default lowercases filenames.
#'
#' It is worth noting that this function has a lot of logic used to
#' find the count tables in the local filesystem.  This logic has been
#' superceded by simply adding a field to the .csv file called
#' 'file'.  create_expt() will then just read that filename, it may be
#' a full pathname or local to the cwd of the project.
#'
#' Also, the logic of this and create_experiment are a bit of a mess and should be redone!
#'
#' @param file   a comma separated file describing the samples with
#' information like condition,batch,count_filename,etc
#' @param sample_colors   a list of colors by condition, if not provided
#' it will generate its own colors using colorBrewer
#' @param gene_info   annotation information describing the rows of the data set, usually
#' this comes from a call to import.gff()
#' @param by_type   when looking for count tables, are they organized by type?
#' @param by_sample   or by sample?  I do all mine by sample, but others do by type...
#' @param include_type   I have usually assumed that all gff annotations should be used,
#' but that is not always true, this allows one to limit.
#' @param include_gff   A gff file to help in sorting which features to keep
#' @param count_dataframe   If one does not wish to read the count tables from processed_data/
#' they may instead be fed here
#' @param meta_dataframe   an optional dataframe containing the metadata rather than a file
#' @param savefile   an Rdata filename prefix for saving the data of the resulting expt.
#' @param low_files   whether or not to explicitly lowercase the filenames when searching in processed_data/
#' @param ... more parameters are fun
#' @return  experiment an expressionset
#' @seealso \pkg{Biobase} \link[Biobase]{pData} \link[Biobase]{fData} \link[Biobase]{exprs}
#' \link{hpgl_read_files} \link[hash]{as.list.hash}
#' @examples
#' \dontrun{
#' new_experiment = create_experiment("some_csv_file.csv", color_hash)
#' ## Remember that this depends on an existing data structure of gene annotations.
#' }
#' @export
create_expt <- function(file=NULL, sample_colors=NULL, gene_info=NULL, by_type=FALSE,
                        by_sample=FALSE, include_type="all", include_gff=NULL, count_dataframe=NULL,
                        meta_dataframe=NULL, savefile="expt", low_files=FALSE, ...) {
    arglist <- list(...)  ## pass stuff like sep=, header=, etc here

    ## Palette for colors when auto-chosen
    chosen_palette <- "Dark2"
    if (!is.null(arglist[["palette"]])) {
        chosen_palette <- arglist[["palette"]]
    }
    file_suffix <- ".count.gz"
    if (!is.null(arglist[["suffix"]])) {
        file_suffix <- arglist[["suffix"]]
    }
    gff_type <- "all"
    if (!is.null(arglist[["include_type"]])) {
        gff_type <- arglist[["include_type"]]
    }

    ## Read in the metadata from the provided data frame, csv, or xlsx.
    sample_definitions <- data.frame()
    if (is.null(meta_dataframe) & is.null(file)) {
        stop("This requires either a csv file or dataframe of metadata describing the samples.")
    } else if (is.null(file)) {
        sample_definitions <- meta_dataframe
    }  else {
        sample_definitions <- read_metadata(file, header=header, sep=sep)
    }

    ## Double-check that there is a usable condition column
    if (is.null(sample_definitions[["condition"]])) {
        sample_definitions[, "condition"] <- tolower(paste(sample_definitions[, "type"], sample_definitions[, "stage"], sep="_"))
    }
    condition_names <- unique(sample_definitions[, "condition"])
    if (is.null(condition_names)) {
        warning("There is no 'condition' field in the definitions, this will make many analyses more difficult/impossible.")
    }

    ## Make sure we have a viable set of colors for plots
    chosen_colors <- as.character(sample_definitions[, "condition"])
    num_conditions <- length(condition_names)
    num_samples <- nrow(sample_definitions)
    if (!is.null(sample_colors) & length(sample_colors) == num_samples) {
        chosen_colors <- sample_colors
    } else if (!is.null(sample_colors) & length(sample_colors) == num_conditions) {
        mapping <- setNames(sample_colors, unique(chosen_colors))
        chosen_colors <- mapping[chosen_colors]
        names(chosen_colors) <- sample_definitions[, "sample.id"]
    } else if (is.null(sample_colors)) {
        sample_colors <- suppressWarnings(grDevices::colorRampPalette(
            RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
        mapping <- setNames(sample_colors, unique(chosen_colors))
        chosen_colors <- mapping[chosen_colors]
        names(chosen_colors) <- sample_definitions[, "sample.id"]
    } else {
        warning("The number of colors provided does not match either the number of conditions nor samples.")
        warning("Unsure of what to do, so choosing colors with RColorBrewer.")
        sample_colors <- suppressWarnings(grDevices::colorRampPalette(
            RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
        mapping <- setNames(sample_colors, unique(chosen_colors))
        chosen_colors <- mapping[chosen_colors]
        names(chosen_colors) <- sample_definitions[, "sample.id"]
    }

    ## Create a matrix of counts with columns as samples and rows as genes
    ## This may come from either a data frame/matrix, a list of files from the metadata
    ## or it can attempt to figure out the location of the files from the sample names.
    filenames <- NULL
    found_counts <- NULL
    all_count_tables <- NULL
    if (!is.null(count_dataframe)) {
        all_count_tables <- count_dataframe
        colnames(all_count_tables) <- rownames(sample_definitions)
        ## If neither of these cases is true, start looking for the files in the processed_data/ directory
    } else if (is.null(sample_definitions[, "file"])) {
        success <- 0
        ## Look for files organized by sample
        test_filenames <- paste0("processed_data/count_tables/", as.character(sample_definitions[, "sample.id"]), "/",
                                 as.character(sample_definitions[, "sample.id"]), file_suffix)
        num_found <- sum(file.exists(test_filenames))
        if (num_found == num_samples) {
            success <- success + 1
            sample_definitions[, "file"] <- test_filenames
        } else {
            test_filenames <- tolower(test_filenames)
            num_found <- sum(file.exists(test_filenames))
            if (num_found == num_samples) {
                success <- success + 1
                sample_definitions[, "file"] <- test_filenames
            }
        }
        if (success == 0) {
            ## Did not find samples by id, try them by type
            test_filenames <- paste0("processed_data/count_tables/", tolower(as.character(sample_definitions[, "type"])), "/",
                                     tolower(as.character(sample_definitions[, "stage"])), "/",
                                     sample_definitions[, "sample.id"], file_suffix)
            num_found <- sum(file.exists(test_filenames))
            if (num_found == num_samples) {
                success <- success + 1
                sample_definitions[, "file"] <- test_filenames
            } else {
                test_filenames <- tolower(test_filenames)
                num_found <- sum(file.exists(test_filenames))
                if (num_found == num_samples) {
                    success <- success + 1
                    sample_definitions[, "file"] <- test_filenames
                }
            }
        } ## tried by type
        if (success == 0) {
            stop("I could not find your count tables organised either by sample nor by type, uppercase nor lowercase.")
        }
    }

    ## At this point sample_definitions$file should be filled in no matter what
    if (is.null(all_count_tables)) {
        filenames <- as.character(sample_definitions[, "file"])
        sample_ids <- as.character(sample_definitions[, "sample.id"])
        all_count_tables <- hpgl_read_files(sample_ids, filenames, ...)
    }

    all_count_matrix <- as.matrix(all_count_tables)
    rownames(all_count_matrix) <- gsub("^exon:", "", rownames(all_count_matrix))
    rownames(all_count_matrix) <- make.names(gsub(":\\d+", "", rownames(all_count_matrix)), unique=TRUE)
    if (is.null(gene_info)) {
        gene_info <- data.frame(all_count_matrix)
    } else {
        if (is.null(gene_info[, "ID"])) {
            gene_info[, "ID"] <- rownames(gene_info)
        }
        gene_info <- gene_info[gene_info[, "ID"] %in% rownames(all_count_matrix),]
        all_count_matrix <- all_count_matrix[rownames(all_count_matrix) %in% gene_info[, "ID"], ]
    }
    ## Make sure that all columns have been filled in for every gene.
    complete_index <- complete.cases(all_count_matrix)
    all_count_matrix <- all_count_matrix[complete_index, ]

    annotation <- NULL
    tooltip_data <- NULL
    if (!is.null(include_gff)) {
        message("create_experiment(): Reading annotation gff, this is slow.")
        annotation <- gff2df(gff=include_gff, type=gff_type)
        tooltip_data <- make_tooltips(annotations=annotation, type=gff_type, ...)
    }

    ## Perhaps I do not understand something about R's syntactic sugar
    ## Given a data frame with columns bob, jane, alice -- but not foo
    ## I can do df[["bob"]]) or df[, "bob"] to get the column bob
    ## however df[["foo"]] gives me null while df[, "foo"] gives an error.
    if (is.null(sample_definitions[["stage"]])) {
        sample_definitions[, "stage"] <- "unknown"
    }
    if (is.null(sample_definitions[["type"]])) {
        sample_definitions[, "type"] <- "unknown"
    }
    if (is.null(sample_definitions[["condition"]])) {
        sample_definitions[, "condition"] <- "unknown"
    }
    if (is.null(sample_definitions[["batch"]])) {
        sample_definitions[, "batch"] <- "unknown"
    }
    if (is.null(sample_definitions[["intercounts"]])) {
        sample_definitions[, "intercounts"] <- "unknown"
    }

    meta_frame <- data.frame(
        sample=as.character(sample_definitions[, "sample.id"]),
        stage=as.character(sample_definitions[, "stage"]),
        type=as.character(sample_definitions[, "type"]),
        condition=as.character(sample_definitions[, "condition"]),
        batch=as.character(sample_definitions[, "batch"]),
        colors=chosen_colors,
        counts=sample_definitions[, "file"],
        intercounts=sample_definitions[, "intercounts"])
    requireNamespace("Biobase")
    metadata <- methods::new("AnnotatedDataFrame", meta_frame)  ## AnnotatedDataFrame is from Biobase
    Biobase::sampleNames(metadata) <- colnames(all_count_matrix)
    feature_data <- methods::new("AnnotatedDataFrame", gene_info)
    Biobase::featureNames(feature_data) <- rownames(all_count_matrix)
    experiment <- methods::new("ExpressionSet", exprs=all_count_matrix,
                               phenoData=metadata, featureData=feature_data)

    ## These entries in new_expt are intended to maintain a record of
    ## the transformation status of the data, thus if we now call
    ## normalize_expt() it should change these.
    ## Therefore, if we call a function like DESeq() which requires
    ## non-log2 counts, we can check these values and convert accordingly
    expt <- expt_subset(experiment)
    expt[["design"]] <- sample_definitions
    expt[["annotation"]] <- annotation
    expt[["gff_file"]] <- include_gff
    expt[["tooltip"]] <- tooltip_data
    expt[["filtered"]] <- FALSE
    expt[["transform"]] <- "raw"
    expt[["norm"]] <- "raw"
    expt[["convert"]] <- "raw"
    expt[["original_libsize"]] <- colSums(Biobase::exprs(experiment))
    if (!is.null(savefile)) {
        save(list = c("expt"), file=paste(savefile, ".Rdata", sep=""))
    }
    return(expt)
}

#'  Extract a subset of samples following some rule(s) from an
#' experiment class
#'
#' @param expt  an expt which is a home-grown class containing an
#' expressionSet, design, colors, etc.
#' @param subset  a valid R expression which defines a subset of the
#' design to keep.
#' @return  metadata an expt class which contains the smaller set of
#' data
#' @seealso \pkg{Biobase} \link[Biobase]{pData}
#' \link[Biobase]{exprs} \link[Biobase]{fData}
#' @examples
#' \dontrun{
#'  smaller_expt = expt_subset(big_expt, "condition=='control'")
#'  all_expt = expt_subset(expressionset, "")  ## extracts everything
#' }
#' @export
expt_subset <- function(expt, subset=NULL) {
    if (class(expt) == "ExpressionSet") {
        expressionset <- expt
    } else if (class(expt) == "expt") {
        expressionset <- expt[["expressionset"]]
    } else {
        stop("expt is neither an expt nor ExpressionSet")
    }
    if (is.null(expt[["design"]])) {
        ## warning("There is no expt$definitions, using the expressionset.")
        initial_metadata <- Biobase::pData(expressionset)
    } else {
        initial_metadata <- expt[["design"]]
    }
    
    if (is.null(subset)) {
        subset_design <- initial_metadata
    } else {
        r_expression <- paste("subset(initial_metadata,", subset, ")")
        subset_design <- eval(parse(text=r_expression))
        ## design = data.frame(sample=samples$sample, condition=samples$condition, batch=samples$batch)
    }
    subset_design <- as.data.frame(subset_design)
    ## This is to get around stupidity with respect to needing all factors to be in a DESeqDataSet
    conditions <- as.factor(as.character(subset_design[["condition"]]))
    batches <- as.factor(as.character(subset_design[["batch"]]))
    subset_design[["condition"]] <- conditions
    subset_design[["batch"]] <- batches
    if (is.null(subset_design[["sample.id"]])) {
        samplenames <- as.character(subset_design[["sample"]])
    } else {
        samplenames <- as.character(subset_design[["sample.id"]])
    }
    colors <- as.character(subset_design[["colors"]])
    names <- paste(conditions, batches, sep="-")
    ##subset_definitions <- expt$definitions[rownames(expt$definitions) %in% samplenames, ]
    subset_libsize <- expt$original_libsize[names(expt$original_libsize) %in% samplenames]
    expressionset <- expressionset[, Biobase::sampleNames(expressionset) %in% samplenames]
    columns <- data.frame(sample=colnames(Biobase::exprs(expressionset)))
    rownames(columns) <- colnames(Biobase::exprs(expressionset))
    metadata <- list(initial_metadata=initial_metadata,
                     original_expressionset=expressionset,
                     expressionset=expressionset,
                     design=subset_design,
                     stages=initial_metadata$stage,
                     types=initial_metadata$type,
                     conditions=conditions,
                     batches=batches,
                     samplenames=samplenames,
                     colors=colors,
                     names=names,
                     filtered=expt$filtered,
                     transform=expt$transform,
                     norm=expt$norm,
                     convert=expt$convert,
                     original_libsize=subset_libsize,
                     columns=columns)
    class(metadata) <- "expt"
    return(metadata)
}

#'  Read a bunch of count tables and create a usable data frame from
#' them.
#' It is worth noting that this function has some logic intended for the elsayed lab's data storage structure.
#' It shouldn't interfere with other usages, but it attempts to take into account different ways the data might be stored.
#'
#' @param ids  a list of experimental ids
#' @param files  a list of files to read
#' @param header   whether or not the count tables include a header row.
#' @param include_summary_rows   whether HTSeq summary rows should be included.
#' @param suffix  an optional suffix to add to the filenames when reading them.
#' @param ... more options for happy time
#' @return  count_table a data frame of count tables
#' @seealso \link{create_experiment}
#' @examples
#' \dontrun{
#'  count_tables = hpgl_read_files(as.character(sample_ids), as.character(count_filenames))
#' }
#' @export
hpgl_read_files <- function(ids, files, header=FALSE, include_summary_rows=FALSE, suffix=NULL, ...) {
    ## load first sample
    lower_filenames <- files
    dirs <- dirname(lower_filenames)
    low_files <- tolower(basename(files))
    if (!is.null(suffix)) {
        low_hpgl <- gsub(suffix, "", basename(files))
        low_hpgl <- tolower(low_hpgl)
        low_hpgl <- paste(low_hpgl, suffix, sep="")
    } else {
        low_hpgl <- gsub("HPGL","hpgl", basename(files))
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
    ##count_table = read.table(files[1], header=header, ...)
    count_table <- try(read.table(files[1], header=header))
    if (class(count_table)[1] == 'try-error') {
        stop(paste0("There was an error reading: ", files[1]))
    }
    message(paste0(files[1], " contains ", length(rownames(count_table)), " rows."))
    colnames(count_table) <- c("ID", ids[1])
    rownames(count_table) <- make.names(count_table$ID, unique=TRUE)
    count_table <- count_table[-1]
    ## iterate over and append remaining samples
    for (table in 2:length(files)) {
        if (file.exists(tolower(files[table]))) {
            files[table] <- tolower(files[table])
        } else if (file.exists(lowhpgl_filenames[table])) {
            files[table] <- lowhpgl_filenames[table]
        } else if (file.exists(lower_filenames[table])) {
            files[table] <- lower_filenames[table]
        }
        tmp_count = try(read.table(files[table], header=header))
        if (class(tmp_count)[1] == 'try-error') {
            stop(paste0("There was an error reading: ", files[table]))
        }
        colnames(tmp_count) <- c("ID", ids[table])
        ##tmp_count <- tmp_count[, c("ID", ids[table])]
        rownames(tmp_count) <- make.names(tmp_count$ID, unique=TRUE)
        tmp_count <- tmp_count[-1]
        pre_merge <- length(rownames(tmp_count))
        count_table <- merge(count_table, tmp_count, by.x="row.names", by.y="row.names", all.x=TRUE)
        rownames(count_table) <- count_table$Row.names
        count_table <- count_table[-1]
        post_merge <- length(rownames(count_table))
        message(paste0(files[table], " contains ", pre_merge, " rows and merges to ", post_merge, " rows."))
    }

    rm(tmp_count)
    ## set row and columns ids
    ## rownames(count_table) <- make.names(count_table$ID, unique=TRUE)
    ## count_table <- count_table[-1]
    colnames(count_table) <- ids

    ## remove summary fields added by HTSeq
    if (!include_summary_rows) {
        htseq_meta_rows <- c('__no_feature', '__ambiguous', '__too_low_aQual',
                            '__not_aligned', '__alignment_not_unique')
        count_table <- count_table[!rownames(count_table) %in% htseq_meta_rows,]
    }
    return(count_table)
}

#' \code{concatenate_runs()}  Sum the reads/gene for multiple sequencing runs of a single condition/batch
#'
#' @param expt  an experiment class containing the requisite metadata and count tables
#' @param column  a column of the design matrix used to specify which samples are replicates
#' @return the input expt with the new design matrix, batches, conditions, colors, and count tables.
#' @seealso
#' \pkg{Biobase}
#' @examples
#' \dontrun{
#'  compressed = concatenate_runs(expt)
#' }
#' @export
concatenate_runs <- function(expt, column='replicate') {
    design <- expt$definitions
    replicates <- levels(as.factor(design[,column]))
    final_expt <- expt
    final_data <- NULL
    final_design <- NULL
    final_definitions <- NULL
    column_names <- list()
    colors <- list()
    conditions <- list()
    batches <- list()
    names <- list()
    for (rep in replicates) {
        expression <- paste0(column, "=='", rep, "'")
        tmp_expt <- expt_subset(expt, expression)
        tmp_data <- rowSums(Biobase::exprs(tmp_expt$expressionset))
        tmp_design <- tmp_expt$design[1, ]
        tmp_definitions <- tmp_expt$definitions[1, ]
        final_data <- cbind(final_data, tmp_data)
        final_design <- rbind(final_design, tmp_design)
        final_definitions <- rbind(final_definitions, tmp_definitions)
        column_names[[rep]] <- as.character(tmp_design$sample.id)
        colors[[rep]] <- as.character(tmp_design$color)
        batches[[rep]] <- as.character(tmp_design$batch)
        conditions[[rep]] <- as.character(tmp_design$condition)
        names[[rep]] <- paste(conditions[[rep]], batches[[rep]], sep='-')
        colnames(final_data) <- column_names
    }
    final_expt$design <- final_design
    final_expt$definitions <- final_definitions
    metadata <- new("AnnotatedDataFrame", final_design)
    Biobase::sampleNames(metadata) <- colnames(final_data)
    feature_data <- new("AnnotatedDataFrame", Biobase::fData(expt$expressionset))
    Biobase::featureNames(feature_data) <- rownames(final_data)
    experiment <- new("ExpressionSet", exprs=final_data,
                      phenoData=metadata, featureData=feature_data)
    final_expt$expressionset <- experiment
    final_expt$original_expressionset <- experiment
    final_expt$samples <- final_design
    final_expt$colors <- as.character(colors)
    final_expt$batches <- as.character(batches)
    final_expt$conditions <- as.character(conditions)
    final_expt$names <- as.character(names)
    return(final_expt)
}

#' Create a data frame of the medians of rows by a given factor in the data
#'
#' This assumes of course that (like expressionsets) there are separate columns for each replicate
#' of the conditions.  This will just iterate through the levels of a factor describing the columns,
#' extract them, calculate the median, and add that as a new column in a separate data frame.
#'
#' @param data  a data frame, presumably of counts.
#' @param fact  a factor describing the columns in the data.
#' @return a data frame of the medians
#' @examples
#' \dontrun{
#'  compressed = hpgltools:::median_by_factor(data, experiment$condition)
#' }
#' @export
median_by_factor <- function(data, fact) {
    medians <- data.frame("ID"=rownames(data))
    rownames(medians) = rownames(data)
    for (type in levels(fact)) {
        columns <- grep(pattern=type, fact)
        med <- matrixStats::rowMedians(data[, columns])
        medians <- cbind(medians, med)
    }
    medians <- medians[-1]
    colnames(medians) <- levels(fact)
    return(medians)
}

## EOF
