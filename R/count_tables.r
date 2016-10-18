#' Read a bunch of count tables and create a usable data frame from them.
#'
#' It is worth noting that this function has some logic intended for the elsayed lab's data storage structure.
#' It shouldn't interfere with other usages, but it attempts to take into account different ways the data might be stored.
#'
#' @param ids List of experimental ids.
#' @param files List of files to read.
#' @param header Whether or not the count tables include a header row.
#' @param include_summary_rows Whether HTSeq summary rows should be included.
#' @param suffix Optional suffix to add to the filenames when reading them.
#' @param ... More options for happy time!
#' @return Data frame of count tables.
#' @seealso \link{create_expt}
#' @examples
#' \dontrun{
#'  count_tables = hpgl_read_files(as.character(sample_ids), as.character(count_filenames))
#' }
#' @export
expt_read_counts <- function(ids, files, header=FALSE, include_summary_rows=FALSE, suffix=NULL, ...) {
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
    count_dt <- data.table::as.data.table(count_table)
    if (class(count_table)[1] == "try-error") {
        stop(paste0("There was an error reading: ", files[1]))
    }
    message(paste0(files[1], " contains ", length(rownames(count_dt)), " rows."))
    colnames(count_dt) <- c("rownames", ids[1])
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
        tmp_count = try(read.table(files[table], header=header))
        if (class(tmp_count)[1] == "try-error") {
            stop(paste0("There was an error reading: ", files[table]))
        }
        colnames(tmp_count) <- c("rownames", ids[table])
        tmp_count <- data.table::as.data.table(tmp_count)
        ##tmp_count <- tmp_count[, c("ID", ids[table])]
        ##rownames(tmp_count) <- make.names(tmp_count[, "ID"], unique=TRUE)
        ##tmp_count <- tmp_count[, -1, drop=FALSE]
        pre_merge <- length(rownames(tmp_count))
        count_dt <- merge(count_dt, tmp_count, by="rownames", all.x=TRUE)
        ## rownames(count_table) <- count_table[, "Row.names"]
        ## count_table <- count_table[, -1, drop=FALSE]
        ## post_merge <- length(rownames(count_table))
        post_merge <- nrow(count_dt)
        message(paste0(files[table], " contains ", pre_merge, " rows and merges to ", post_merge, " rows."))
    }
    count_table <- as.data.frame(count_dt)
    rownames(count_table) <- count_table[["rownames"]]
    count_table <- count_table[, -1, drop=FALSE]
    rm(count_dt)
    rm(tmp_count)
    ## set row and columns ids
    ## rownames(count_table) <- make.names(count_table$ID, unique=TRUE)
    ## count_table <- count_table[, -1, drop=FALSE]
    ## colnames(count_table) <- ids

    ## remove summary fields added by HTSeq
    if (!include_summary_rows) {
        ## Depending on what happens when the data is read in, these rows may get prefixed with 'X'
        ## In theory, only 1 of these two cases should ever be true.
        htseq_meta_rows <- c("__no_feature", "__ambiguous", "__too_low_aQual",
                             "__not_aligned", "__alignment_not_unique",
                             "X__no_feature", "X__ambiguous", "X__too_low_aQual",
                             "X__not_aligned", "X__alignment_not_unique")

        count_table <- count_table[!rownames(count_table) %in% htseq_meta_rows, ]
    }
    return(count_table)
}

#' Sum the reads/gene for multiple sequencing runs of a single condition/batch.
#'
#' On occasion we have multiple technical replicates of a sequencing run.  This can use a column in
#' the experimental design to identify those replicates and sum the counts into a single column in
#' the count tables.
#'
#' @param expt Experiment class containing the requisite metadata and count tables.
#' @param column Column of the design matrix used to specify which samples are replicates.
#' @return Expt with the concatenated counts, new design matrix, batches, conditions, etc.
#' @seealso
#' \pkg{Biobase}
#' @examples
#' \dontrun{
#'  compressed = concatenate_runs(expt)
#' }
#' @export
concatenate_runs <- function(expt, column='replicate') {
    design <- expt[["design"]]
    replicates <- levels(as.factor(design[, column]))
    final_expt <- expt
    final_data <- NULL
    final_design <- NULL
    final_definitions <- NULL
    column_names <- list()
    colors <- list()
    conditions <- list()
    batches <- list()
    samplenames <- list()
    for (rep in replicates) {
        expression <- paste0(column, "=='", rep, "'")
        tmp_expt <- expt_subset(expt, expression)
        tmp_data <- rowSums(Biobase::exprs(tmp_expt[["expressionset"]]))
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
    Biobase::sampleNames(metadata) <- colnames(final_data)
    feature_data <- new("AnnotatedDataFrame", Biobase::fData(expt[["expressionset"]]))
    Biobase::featureNames(feature_data) <- rownames(final_data)
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

#' Create a data frame of the medians of rows by a given factor in the data.
#'
#' This assumes of course that (like expressionsets) there are separate columns for each replicate
#' of the conditions.  This will just iterate through the levels of a factor describing the columns,
#' extract them, calculate the median, and add that as a new column in a separate data frame.
#'
#' @param data Data frame, presumably of counts.
#' @param fact Factor describing the columns in the data.
#' @return Data frame of the medians.
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
    medians <- medians[, -1, drop=FALSE]
    colnames(medians) <- levels(fact)
    return(medians)
}

#' Count the number of features(genes) greater than x in a data set.
#'
#' Sometimes I am asked how many genes have >= x counts.  Well, here you go.
#'
#' @param data  A dataframe/exprs/matrix/whatever of counts.
#' @param cutoff  Minimum number of counts.
#' @param hard  Greater-than is hard, greater-than-equals is not.
#' @return  Number of genes.
#' @export
features_greater_than <- function(data, cutoff=1, hard=TRUE) {
    if (class(data) == "expt") {
        data <- as.data.frame(Biobase::exprs(data[["expressionset"]]))
    } else if (class(data) == "ExpressionSet") {
        data <- as.data.frame(Biobase::exprs(data))
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

write_pretty_counts <- function(expt, excel="excel/pretty_counts.xlsx", norm="quant",
                                convert="cpm", transform="log2", model_batch="sva", filter=TRUE) {
    wb <- openxlsx::createWorkbook(creator="hpgltools")
    plot_dim <- 6
    plot_cols <- floor(plot_dim * 1.5)
    plot_rows <- ceil(plot_dim * 0.8)
    new_row <- 1
    new_col <- 1


    ## Write an introduction to this foolishness.
    message("Writing the legend.")
    sheet <- "legend"
    legend <- data.frame(rbind(
        c("This excel workbook the following worksheets:", ""),
        c("1.", "This sheet, including the experimental design."),
        c("2.", "The raw counts and annotation data on worksheet 'raw_data'."),
        c("3.", "Some graphs describing the distribution of raw data in worksheet 'raw_plots'."),
        c("4.", paste0("The counts normalized with: ", transform, "(", convert, "(", norm, "(", batch, "(", filter, "(counts))))) on 'norm_data'.")),
        c("5.", "Some graphs describing the distribution of the normalized data on 'norm_plots'."),
        c("6.", "The median normalized counts by condition factor on 'median_data'.")))
    colnames(legend) <- c("sheet", "sheet definition")
    xls_result <- write_xls(wb, data=legend, sheet=sheet, rownames=FALSE,
                            title="Columns used in the following tables.")
    rows_down <- nrow(legend)
    new_row <- new_row + rows_down + 1
    annot <- Biobase::fData(expt[["expressionset"]])
    xls_result <- write_xls(wb, data=annot, start_row=new_row, sheet=sheet, start_col=1, title="Experimental Design.")


    ## Write the raw read data and gene annotations
    message("Writing the raw reads.")
    sheet <- "raw_reads"
    new_row <- 1
    reads <- Biobase::exprs(expt[["expressionset"]])
    info <- Biobase::pData(expt[["expressionset"]])
    read_info <- merge(reads, info, by="row.names")
    xls_result <- write_xls(wb, data=read_info, start_row=new_row, start_col=new_col, sheet=sheet, title="Raw Reads.")


    ## Write some graphs for the raw data
    message("Writing the raw reads.")
    sheet <- "raw_graphs"
    metrics <- graph_metrics(expt[["expressionset"]], qq=TRUE)
    ## Start with library sizes.
    openxlsx::writeData(wb, sheet=sheet, x="Legend.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw library sizes.", startRow=new_row, startCol=new_col + plot_cols + 1)
    openxlsx::writeData(wb, sheet=sheet, x="Non-zero genes.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_row <- new_row + 1
    new_plot <- metrics[["legend"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- new_col + plot_cols + 1
    new_plot <- metrics[["libsize"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    ## Same row, non-zero plot
    new_col <- new_col + plot_cols + 1
    new_plot <- metrics[["non_zero"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- new_col + plot_cols + 1

    ## Visualize distributions
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw data density plot.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw Boxplot.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_plot <- metrics[["density"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_plot <- metrics[["boxplot"]]
    new_col <- new_col + plot_cols + 1
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- 1

    ## Move down next set of rows, heatmaps
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw correlation heatmap.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw distance heatmap.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_plot <- metrics[["corheat"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_plot <- metrics[["disheat"]]
    new_col <- new_col + plot_cols + 1
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- 1

    ## SM plots
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw standard median correlation.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw standard distance correlation.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_plot <- metrics[["smc"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_plot <- metrics[["smd"]]
    new_col <- new_col + plot_cols + 1
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- 1

    ## PCA and qq_log
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw PCA.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw QQ, log scale.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_plot <- metrics[["pcaplot"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_plot <- metrics[["qqlog"]]
    new_col <- new_col + plot_cols + 1
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- 1

    ## PCA table
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw PCA res.", startRow=new_row, startCol=new_col)
    xls_result <- write_xls(data=metrics[["pcares"]], wb=wb, sheet="raw_graphs")
    new_col <- xls_result[["end_col"]] + 2
    openxlsx::writeData(wb, sheet=sheet, x="Raw PCA table.", startRow=new_row, startCol=new_col)
    xls_result <- write_xls(data=metrics[["pcatable"]], wb=wb, sheet="raw_graphs")


    ## Move on to the next sheet, normalized data
    message("Writing the normalized reads.")
    sheet <- "norm_data"
    new_col <- 1
    new_row <- 1
    norm <- normalize_expt(expt=expt, transform=transform, norm=norm, convert=convert, batch=batch, filter=filter, ...)
    reads <- Biobase::exprs(norm[["expressionset"]])
    info <- Biobase::pData(norm[["expressionset"]])
    read_info <- merge(reads, info, by="row.names")
    xls_result <- write_xls(wb, data=read_info, start_row=new_row, start_col=new_col, sheet=sheet, title="Raw Reads.")


    ## Graphs of the normalized data
    message("Writing the normalized graphs.")
    sheet <- "norm_graphs"
    metrics <- graph_metrics(norm[["expressionset"]], qq=TRUE)
    ## Start with library sizes.
    openxlsx::writeData(wb, sheet=sheet, x="Legend.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw library sizes.", startRow=new_row, startCol=new_col + plot_cols + 1)
    openxlsx::writeData(wb, sheet=sheet, x="Non-zero genes.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_row <- new_row + 1
    new_plot <- metrics[["legend"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- new_col + plot_cols + 1
    new_plot <- metrics[["libsize"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    ## Same row, non-zero plot
    new_col <- new_col + plot_cols + 1
    new_plot <- metrics[["non_zero"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- new_col + plot_cols + 1

    ## Visualize distributions
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw data density plot.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw Boxplot.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_plot <- metrics[["density"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_plot <- metrics[["boxplot"]]
    new_col <- new_col + plot_cols + 1
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- 1

    ## Move down next set of rows, heatmaps
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw correlation heatmap.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw distance heatmap.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_plot <- metrics[["corheat"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_plot <- metrics[["disheat"]]
    new_col <- new_col + plot_cols + 1
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- 1

    ## SM plots
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw standard median correlation.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw standard distance correlation.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_plot <- metrics[["smc"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_plot <- metrics[["smd"]]
    new_col <- new_col + plot_cols + 1
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- 1

    ## PCA and qq_log
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw PCA.", startRow=new_row, startCol=new_col)
    openxlsx::writeData(wb, sheet=sheet, x="Raw QQ, log scale.", startRow=new_row, startCol=new_col + plot_cols + 1)
    new_plot <- metrics[["pcaplot"]]
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_plot <- metrics[["qqlog"]]
    new_col <- new_col + plot_cols + 1
    print(new_plot)
    openxlsx::insertPlot(wb, sheet=sheet, width=plot_dim, height=plot_dim,
                         startCol=plot_column, startRow=new_row, fileType="png", units="in")
    new_col <- 1

    ## PCA table
    new_row <- new_row + plot_rows + 1
    new_col <- 1
    openxlsx::writeData(wb, sheet=sheet, x="Raw PCA res.", startRow=new_row, startCol=new_col)
    xls_result <- write_xls(data=metrics[["pcares"]], wb=wb, sheet="raw_graphs")
    new_col <- xls_result[["end_col"]] + 2
    openxlsx::writeData(wb, sheet=sheet, x="Raw PCA table.", startRow=new_row, startCol=new_col)
    xls_result <- write_xls(data=metrics[["pcatable"]], wb=wb, sheet="raw_graphs")


    ## Give a median-by-factor accounting of the data
    message("Writing the median reads by factor.")
    sheet <- "median_data"
    new_col <- 1
    new_row <- 1
    median_data <- median_by_factor(Biobase::exprs(norm[["expressionset"]]), norm[["condition"]])
    median_data <- merge(median_data, info, by="row.names")
    xls_result <- write_xls(wb, data=median_data, start_row=new_row, start_col=new_col, sheet=sheet, title="Median Reads by factor.")


    ## Save the result
    save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite=TRUE))
    return(save_result)
}

## EOF
