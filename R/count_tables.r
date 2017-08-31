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
  ## arglist <- list(...)
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
    tmp_count <- try(read.table(files[table], header=header))
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
    message(paste0(files[table], " contains ", pre_merge,
                   " rows and merges to ", post_merge, " rows."))
  }
  count_table <- as.data.frame(count_dt, stringsAsFactors=FALSE)
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
    columns <- grep(pattern=type, fact)
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

## EOF
