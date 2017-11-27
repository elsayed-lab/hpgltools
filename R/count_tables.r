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

## EOF
