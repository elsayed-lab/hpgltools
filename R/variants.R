#' Gather snp information for an expt
#'
#' I made some pretty significant changes to the set of data which I
#' retain when using mpileup/freebayes.  As a result, this function
#' needs to be reworked.
#'
#' This function attempts to gather a set of variant positions using an extant
#' expressionset.  This therefore seeks to keep the sample metadata consistent
#' with the original data.  In its current iteration, it therefore makes some
#' potentially bad assumptions about the naming conventions for its input
#' files.  It furthermore assumes inputs from the variant calling methods in
#' cyoa.
#'
#' @param expt an expressionset from which to extract information.
#' @param annot_column Column in the metadata for getting the table of bcftools calls.
#' @param tolower Lowercase stuff like 'HPGL'?
#' @param snp_column Which column of the parsed bcf table contains our interesting material?
#' @return A new expt object
#' @seealso [Biobase] freebayes:DOI:10.48550/arXiv.1207.3907,
#'  mpileup:DOI:10.1093/gigascience/giab008
#' @examples
#'   \dontrun{
#'  expt <- create_expt(metadata, gene_information)
#'  snp_expt <- count_expt_snps(expt)
#'  ## This assumes that the metadata has a column named 'bcftable' with one file per
#'  ## cell.  These files in turn should have a column named 'diff_count' which will
#'  ## be the source of the numbers found when doing exprs(snp_expt).
#' }
#' @export
count_expt_snps <- function(expt, annot_column = "bcftable", tolower = TRUE,
                            snp_column = "diff_count") {
  samples <- rownames(pData(expt))
  if (isTRUE(tolower)) {
    samples <- tolower(samples)
  }
  file_lst <- pData(expt)[[annot_column]]
  if (is.null(file_lst)) {
    stop("This requires a set of bcf filenames, the column: ", annot_column, " does not have any.")
  }

  snp_dt <- read_snp_columns(samples, file_lst, column = snp_column)
  test_names <- snp_dt[["rownames"]]
  test <- grep(pattern = "^chr", x = test_names)
  if (length(test) == 0) {
    message("The rownames are missing the chromosome identifier,
they probably came from an older version of this method.")
    pat <- "^(.+)_(.+)_(.+)_(.+)$"
    new <-"chr_\\1_pos_\\2_ref_\\3_alt_\\4"
    new_names <- gsub(pattern = pat, replacement = new, x = test_names)
    snp_dt[["rownames"]] <- new_names
  }
  ## Make sure there are no underscores in the chromosome names.
  snp_dt[["rownames"]] <- gsub(pattern = "chr_(.*)_(.*)_pos",
                               replacement = "chr_\\1-\\2_pos",
                               x = snp_dt[["rownames"]])
  ## Get rid of underscores if they are in the chromosome name.

  snp_exprs <- as.data.frame(snp_dt)
  rownames(snp_exprs) <- snp_exprs[["rownames"]]
  snp_exprs[["rownames"]] <- NULL

  snp_features <- data.table::as.data.table(snp_dt)
  snp_features[, c("chr", "chromosome",
                   "pos", "position",
                   "ref", "original",
                   "alt", "new") :=
                   data.table::tstrsplit(snp_features[["rownames"]], "_", fixed = TRUE)]
  snp_features <- snp_features[, c("chromosome", "position", "original", "new")]
  snp_features <- as.data.frame(snp_features)
  rownames(snp_features) <- snp_dt[["rownames"]]

  snp_metadata <- pData(expt)
  snp_metadata <- new("AnnotatedDataFrame", snp_metadata)
  Biobase::sampleNames(snp_metadata) <- colnames(snp_exprs)

  snp_features <- new("AnnotatedDataFrame", snp_features)
  Biobase::featureNames(snp_features) <- rownames(snp_exprs)

  expressionset <- new("ExpressionSet", exprs = snp_exprs,
                       phenoData = snp_metadata, featureData = snp_features)

  new_expt <- expt
  new_expt[["expressionset"]] <- expressionset
  new_expt[["original_expressionset"]] <- expressionset
  ## Make a matrix where the existence of a variant is 1 vs. 0.
  mtrx <- exprs(expressionset)
  idx <- mtrx > 0
  mtrx[idx] <- 1
  new_expt[["variants"]] <- mtrx
  return(new_expt)
}

#' Extract the observed snps unique to individual categories in a snp set.
#'
#' The result of get_snp_sets provides sets of snps for all possible
#' categories.  This is cool and all, but most of the time we just want the
#' results of a single group in that rather large set (2^number of categories)
#'
#' @param retlist  The result from get_snp_sets().
get_individual_snps <- function(retlist) {
  individual_lst <- list()
  possibilities <- retlist[["possibilities"]]
  elements <- retlist[["elements"]]
  names <- retlist[["names"]]
  for (p in possibilities) {
    foundp <- p == names
    individual_lst[[p]] <- elements[[foundp]]
  }
  return(individual_lst)
}

#' Create all possible sets of variants by sample (types).
#'
#' I like this function.  It generates an exhaustive catalog of the snps by
#' chromosome for all the various categories as defined by factor.
#'
#' @param snp_expt Expressionset of variants.
#' @param factor Use this metadata factor to split the data.
#' @param stringency Allow for some wiggle room in the calls.
#' @param do_save Save the results to an rda fil.
#' @param savefile This is redundant with do_save.
#' @param proportion Used with stringency to finetune the calls.
#' @return A funky list by chromosome containing:  'medians', the median number
#'   of hits / position by sample type; 'possibilities', the;
#'  'intersections', the groupings as detected by Vennerable;
#'  'chr_data', the raw data; 'set_names', a character list of the actual names
#'  of the groupings; 'invert_names', the opposite of set_names which is to say
#'  the names of groups which do _not_ include samples x,y,z; 'density', a list
#'  of snp densities with respect to chromosomes.  Note that this last one is
#'  approximate as I just calculate with the largest chromosome position
#'  number, not the explicit number of nucleotides in the chromosome.
#' @seealso [medians_by_factor()]
#' @examples
#'   \dontrun{
#'  expt <- create_expt(metadata, gene_information)
#'  snp_expt <- count_expt_snps(expt)
#'  snp_sets <- get_snp_sets(snp_expt, factor = "condition")
#'  ## This assumes a column in the metadata for the expt named 'condition'.
#' }
#' @export
get_snp_sets <- function(snp_expt, factor = "pathogenstrain",
                         stringency = NULL, do_save = FALSE,
                         savefile = "variants.rda", proportion = 0.9) {
  if (is.null(pData(snp_expt)[[factor]])) {
    stop("The factor does not exist in the expt.")
  } else {
    message("The samples represent the following categories: ")
    print(table(pData(snp_expt)[[factor]]))
  }
  if (isTRUE(do_save) && file.exists(glue("{savefile}_{factor}.rda"))) {
    retlist <- new.env()
    loaded <- load(savefile, envir = retlist)
    retlist <- retlist[["retlist"]]
    return(retlist)
  }

  if (!is.null(proportion)) {
    message("Using a proportion of observed variants, converting the data to binary observations.")
    nonzero <- exprs(snp_expt)
    nonzero[nonzero > 0] <- 1
    stringency <- "proportion"
    nonzero <- as.matrix(nonzero)
    exprs(snp_expt) <- nonzero
  }

  ## This function should be renamed to summarize_by_factor.
  by_factor_data <- median_by_factor(snp_expt, fact = factor)
  values <- data.frame()
  observed <- data.frame()
  if (stringency == "proportion") {
    values <- by_factor_data[["sums"]]
    observed <- values
    spc <- by_factor_data[["samples_per_condition"]]
    for (f in 1:length(spc)) {
      num <- spc[f]
      fact <- names(spc)[f]
      observed[[fact]] <- observed[[fact]] / num
    }
    not_observed_idx <- observed < proportion
    observed_idx <- observed >= proportion
    observed[not_observed_idx] <- 0
    observed[observed_idx] <- 1
  } else if (stringency == "min") {
    values <- by_factor_data[["mins"]]
    observed <- values
    observed[observed > 0] <- 1
  } else if (stringency == "max") {
    values <- by_factor_data[["maxs"]]
    observed <- values
    observed[observed > 0] <- 1
  } else if (stringency == "median") {
    observed <- by_factor_data[["values"]]
  } else {
    stop("I do not understand this stringency parameter.")
  }

  ## I am going to split this by chromosome, as a run of 10,000 took 2 seconds,
  ## 100,000 took a minute, and 400,000 took an hour.
  ##chr <- gsub(pattern = "^.+_(.+)_.+_.+_.+$", replacement = "\\1", x = rownames(values))
  chr <- gsub(pattern = "^chr_(.+)_pos_.+_ref_.+_alt_.+$",
              replacement = "\\1", x = rownames(observed))
  observed[["chr"]] <- chr

  tt <- sm(requireNamespace("parallel"))
  tt <- sm(requireNamespace("doParallel"))
  tt <- sm(requireNamespace("iterators"))
  tt <- sm(requireNamespace("foreach"))
  tt <- sm(try(attachNamespace("foreach"), silent = TRUE))
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  num_levels <- length(levels(as.factor(chr)))

  show_progress <- interactive() && is.null(getOption("knitr.in.progress"))
  pb_opts <- list()
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(max = num_levels, style = 3)
    progress <- function(n) {
      setTxtProgressBar(bar, n)
    }
    pb_opts[["progress"]] <- progress
  }
  returns <- list()
  i <- 1
  res <- list()
  res <- foreach(i = 1:num_levels, .packages = c("hpgltools", "doParallel"),
                 .options.snow = pb_opts, .export = c("snp_by_chr")) %dopar% {

    chromosome_name <- levels(as.factor(chr))[i]
    returns[[chromosome_name]] <- snp_by_chr(observed, chr_name = chromosome_name)
  }
  if (isTRUE(show_progress)) {
    close(bar)
  }
  message("Finished iterating over the chromosomes.")
  parallel::stopCluster(cl)

  ## Unpack the res data structure (which probably can be simplified)
  data_by_chr <- list()
  possibilities <- c()
  set_names <- list()
  invert_names <- list()
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style = 3)
  }
  end <- length(res)
  mesg("Iterating over ", end, " elements.")
  for (element in seq_len(end)) {
    if (isTRUE(show_progress)) {
      pct_done <- element / length(res)
      utils::setTxtProgressBar(bar, pct_done)
    }
    datum <- res[[element]]
    chromosome <- datum[["chromosome"]]
    data_by_chr[[chromosome]] <- datum
    possibilities <- data_by_chr[[chromosome]][["possibilities"]]
    set_names <- data_by_chr[[chromosome]][["set_names"]]
    invert_names <- data_by_chr[[chromosome]][["invert_names"]]
  }
  if (isTRUE(show_progress)) {
    close(bar)
  }

  ## Calculate approximate snp densities by chromosome
  all_intersections <- list()
  density_by_chr <- list()
  for (chr in names(data_by_chr)) {
    ## snps <- rownames(data_by_chr[[chr]][["medians"]])
    snps <- rownames(data_by_chr[[chr]][["observations"]])
    num_snps <- length(snps)
    ##last_position <- max(as.numeric(gsub(pattern = "^.+_.+_(.+)_.+_.+$",
    ##                                     replacement = "\\1", x = snps)))
    last_position <- max(
        as.numeric(gsub(pattern = "^chr_.+_pos_(.+)_ref_.+_alt_.+$",
                        replacement = "\\1", x = snps)))
    snp_density <- num_snps / as.numeric(last_position)
    density_by_chr[[chr]] <- snp_density
    for (inter in names(data_by_chr[[chr]][["intersections"]])) {
      if (is.null(all_intersections[[inter]])) {
        all_intersections[[inter]] <- data_by_chr[[chr]][["intersections"]][[inter]]
      } else {
        all_intersections[[inter]] <- c(all_intersections[[inter]],
                                        data_by_chr[[chr]][["intersections"]][[inter]])
      }
    }
  }
  density_by_chr <- unlist(density_by_chr)

  retlist <- list(
    "factor" = factor,
    "values" = values,
    "observations" = observed,
    "possibilities" = possibilities,
    "intersections" = all_intersections,
    "chr_data" = data_by_chr,
    "set_names" = set_names,
    "invert_names" = invert_names,
    "density" = density_by_chr)
  if (isTRUE(do_save)) {
    saved <- save(list = "retlist", file = savefile)
  }
  class(retlist) <- "snp_sets"
  return(retlist)
}

#' Read the output from bcfutils into a count-table-esque
#'
#' Previously, I put all my bcfutils output files into one directory.  This
#' function would iterate through every file in that directory and add the
#' contents as columns to this growing data table.  Now it works by accepting a
#' list of filenames (presumably kept in the metadata for the experiment) and
#' reading them into the data table.  It is worth noting that it can accept
#' either a column name or index -- which when you think about it is pretty much
#' always true, but in this context is particularly interesting since I changed
#' the names of all the columns when I rewrote this functionality.
#'
#' @param samples Sample names to read.
#' @param file_lst Set of files to read.
#' @param column Column from the bcf file to read.
#' @param verbose Print information about the input data.
#' @seealso [readr]
#' @return A big honking data table.
read_snp_columns <- function(samples, file_lst, column = "diff_count", verbose = FALSE) {
  ## Read the first file
  first_sample <- samples[1]
  if (isTRUE(verbose)) {
    mesg("Reading sample: ", first_sample, ".")
  }
  first_file <- file_lst[1]
  first_read <- readr::read_tsv(first_file, show_col_types = FALSE)
  colnames(first_read) <- make.names(
      gsub(x = colnames(first_read), pattern = "\\.+\\d+$", replacement = ""),
      unique = TRUE)
  if (is.null(first_read[[column]])) {
    stop("The column: ", column, " does not appear to exist in the variant summary file.")
  }
  ## Create a simplified data table from it.
  input_column <- as.numeric(first_read[[column]])
  first_column <- data.table::as.data.table(input_column)
  first_column[["rownames"]] <- first_read[[1]]
  colnames(first_column) <- c(first_sample, "rownames")
  ## Copy that dt to the final data structure.
  snp_columns <- first_column

  show_progress <- interactive() && is.null(getOption("knitr.in.progress"))
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style = 3)
  }
  ## Foreach sample, do the same read of the data and merge it onto the end of the
  ## final data table.
  for (sample_num in seq(from = 2, to = length(samples))) {
    if (isTRUE(show_progress)) {
      pct_done <- sample_num / length(samples)
      setTxtProgressBar(bar, pct_done)
    }
    sample <- samples[sample_num]
    mesg("Reading sample: ", sample, ", number ", sample_num, " of ",
         length(samples), ".")
    file <- file_lst[sample_num]
    if (is.na(file)) {
      mesg("There is no file listed for ", sample, ".")
      next
    }
    if (!file.exists(file)) {
      mesg("Unable to find file: ", file, " for ", sample, ", skipping it.")
      next
    }
    new_table <- try(readr::read_tsv(file, show_col_types = FALSE))
    colnames(new_table) <- make.names(
        gsub(x = colnames(new_table), pattern = "\\.+\\d+$", replacement = ""),
        unique = TRUE)
    if (class(new_table)[1] == "try-error") {
      next
    }
    ## I am not sure why, but a recent sample came up as a character vector
    ## instead of numeric...
    input_column <- as.numeric(new_table[[column]])
    new_column <- data.table::as.data.table(input_column)
    new_column[["rownames"]] <- new_table[[1]]
    colnames(new_column) <- c(sample, "rownames")
    snp_columns <- merge(snp_columns, new_column, by = "rownames", all = TRUE)
  }
  if (isTRUE(show_progress)) {
    close(bar)
  }
  na_positions <- is.na(snp_columns)
  snp_columns[na_positions] <- 0
  return(snp_columns)
}

#' The real worker.  This extracts positions for a single chromosome and puts
#' them into a parallelizable data structure.
#'
#' @param observations A set of observations by position to look through
#' @param chr_name Chromosome name to search
#' @param limit Minimum number of median hits/position to count as a snp.
#' @return A list of variant positions where each element is one chromosome.
#' @seealso [Vennerable]
snp_by_chr <- function(observations, chr_name = "01", limit = 1) {
  set_names <- list()
  possibilities <- c()
  count <- 0
  kept_rows <- observations[["chr"]] == chr_name
  observations <- observations[kept_rows, ]
  kept_cols <- "chr" != colnames(observations)
  observations <- observations[, kept_cols]
  data_by_chr <- list()
  data_by_chr[["chromosome"]] <- chr_name
  data_by_chr[["observations"]] <- observations
  limit_true_false <- as.data.frame(observations >= limit)
  x_lst <- list()
  possibilities <- unique(c(possibilities, colnames(limit_true_false)))
  for (d in seq_along(possibilities)) {
    column_name <- colnames(limit_true_false)[d]
    column_data <- limit_true_false[[column_name]]
    included_snps <- rownames(limit_true_false)[column_data]
    x_lst[[column_name]] <- included_snps
  }
  data_by_chr[["lst"]] <- x_lst
  set_list <- Vennerable::Venn(Sets = x_lst)
  data_by_chr[["venn"]] <- set_list
  elements_per_set <- set_list@IntersectionSets
  data_by_chr[["intersections"]] <- elements_per_set
  set_names <- list()
  invert_names <- list()
  inner_count <- 0
  for (symbolic in names(elements_per_set)) {
    inner_count <- inner_count + 1
    symbols <- strsplit(as.character(symbolic), "")[[1]]
    name <- c()
    invert_name <- c()
    for (s in seq_along(symbols)) {
      symbol <- symbols[s]
      if (symbol == 1) {
        name <- c(name, possibilities[s])
      } else {
        invert_name <- c(invert_name, possibilities[s])
      }
    }
    name <- toString(name)
    invert_name <- toString(invert_name)
    set_names[[inner_count]] <- name
    invert_names[[inner_count]] <- invert_name
  } ## End for symbolic in names(elements_per_set)
  names(set_names) <- names(elements_per_set)
  names(invert_names) <- names(elements_per_set)
  data_by_chr[["set_names"]] <- set_names
  data_by_chr[["invert_names"]] <- invert_names
  data_by_chr[["possibilities"]] <- possibilities
  return(data_by_chr)
}

#' Cross reference observed variants against the transcriptome annotation.
#'
#' This function should provide counts of how many variant positions were
#' observed with respect to each chromosome and with respect to each annotated
#' sequence (currently this is limited to CDS, but that is negotiable).
#'
#' @param expt The original expressionset.  This provides the annotation data.
#' @param snp_result The result from get_snp_sets or count_expt_snps.
#' @param start_column Metadata column with the start position of each ORF.
#' @param end_column Metadata column with the end position of each ORF.
#' @param chr_column Column in the annotation with the chromosome names.
#' @return List containing the set of intersections in the conditions contained
#'  in snp_result, the summary of numbers of variants per chromosome, and
#   summary of numbers per gene.
#' @seealso [snps_vs_genes()] [GenomicRanges::makeGRangesFromDataFrame()]
#'  [IRanges::subsetByOverlaps()] [IRanges::countOverlaps()]
#' @examples
#'  \dontrun{
#'  expt <- create_expt(metadata, gene_information)
#'  snp_expt <- count_expt_snps(expt)
#'  snp_result <- get_snp_sets(snp_expt)
#'  intersections <- snps_vs_intersections(expt, snp_result)
#' }
#' @export
snps_intersections <- function(expt, snp_result,
                               start_column = "start", end_column = "end",
                               chr_column = "seqnames") {
  features <- fData(expt)
  if (is.null(features[[start_column]])) {
    stop("The column containing the feature starts is missing.")
  }
  if (is.null(features[[end_column]])) {
    stop("The column containing the feature ends is missing.")
  }
  if (is.null(features[[chr_column]])) {
    stop("The column containing the feature chromosomes is missing.")
  }

  features[["start"]] <- sm(as.numeric(features[[start_column]]))
  na_starts <- is.na(features[["start"]])
  features <- features[!na_starts, ]
  features[["end"]] <- as.numeric(features[[end_column]])
  ## Again, remember that I modified the seqnames of the snp_expt.
  features[["seqnames"]] <- gsub(pattern = "_",
                                 replacement = "-",
                                 x = features[[chr_column]])
  features[["seqnames"]] <- gsub(pattern = "^.+_(.+)$",
                                 replacement = "\\1", x = features[["seqnames"]])
  grange_input <- features[, c("start", "end", "seqnames")]
  grange_dup <- duplicated(grange_input)
  expt_granges <- GenomicRanges::makeGRangesFromDataFrame(grange_input)

  set_names <- snp_result[["set_names"]]
  chr_summaries <- list()
  gene_summaries <- list()
  inters <- list()
  for (inter in names(snp_result[["intersections"]])) {
    inter_name <- set_names[[inter]]
    inter_df <- data.frame(row.names = snp_result[["intersections"]][[inter]])
    inter_df[["seqnames"]] <- gsub(pattern = "^chr_(.+)_pos_.+_ref.+_alt.+$",
                                   replacement = "\\1",
                                   x = rownames(inter_df))
    inter_df[["start"]] <- as.numeric(gsub(pattern = "^chr_.+_pos_(.+)_ref.+_alt.+$",
                                           replacement = "\\1",
                                           x = rownames(inter_df)))
    inter_df[["end"]] <- inter_df[["start"]] + 1
    inter_df[["strand"]] <- "+"
    inter_granges <- GenomicRanges::makeGRangesFromDataFrame(inter_df)
    inter_by_gene <- suppressWarnings(IRanges::subsetByOverlaps(inter_granges,
                                                                expt_granges,
                                                                type = "within",
                                                                ignore.strand = TRUE))
    inters[[inter_name]] <- inter_by_gene
    summarized_by_chr <- data.table::as.data.table(inter_by_gene)
    ## Faking out r cmd check with a couple empty variables which will be used by data.table
    seqnames <- count <- NULL
    .N <- NULL  ## .N is a read-only symbol in data.table
    summarized_by_chr[, count := .N, by = list(seqnames)]
    summarized_by_chr <- unique(summarized_by_chr[, c("seqnames", "count"), with = FALSE])
    chr_summaries[[inter_name]] <- summarized_by_chr

    summarized_by_gene <- suppressWarnings(IRanges::countOverlaps(query = expt_granges,
                                                                  subject = inter_granges,
                                                                  type = "any", ignore.strand = TRUE))
    summarized_idx <- order(summarized_by_gene, decreasing = TRUE)
    summarized_by_gene <- summarized_by_gene[summarized_idx]
    gene_summaries[[inter_name]] <- summarized_by_gene
  }

  retlist <- list(
      "inters" = inters,
      "chr_summaries" = chr_summaries,
      "gene_summaries" = gene_summaries)
  class(retlist) <- "snp_intersections"
  return(retlist)
}

#' Look for only the variant positions in a subset of genes.
#'
#' This was written in response to a query from Nancy and Maria Adelaida who
#' wanted to look only at the variant positions in a few specific genes.
#'
#' @param expt Initial expressionset.
#' @param snp_expt Variant position expressionset.
#' @param start_col Metadata column with the start positions for each gene.
#' @param end_col Metadata column with the end of the genes.
#' @param expt_name_col Metadata column with the chromosome names.
#' @param snp_name_col Ditto for the snp_expressionset.
#' @param snp_start_col Metadata column containing the variant positions.
#' @param expt_gid_column ID column for the genes.
#' @param genes Set of genes to cross reference.
#' @return New expressionset with only the variants for the genes of interest.
#' @seealso [GenomicRanges::makeGRangesFromDataFrame()] [IRanges::subsetByOverlaps()]
#' @export
snp_subset_genes <- function(expt, snp_expt, start_col = "start", end_col = "end",
                             expt_name_col = "chromosome", snp_name_col = "chromosome",
                             snp_start_col = "position", expt_gid_column = "gid",
                             genes = c("LPAL13_120010900", "LPAL13_340013000", "LPAL13_000054100",
                                       "LPAL13_140006100", "LPAL13_180018500", "LPAL13_320022300")) {
  features <- fData(expt)
  if (is.null(features[[start_col]])) {
    stop("Unable to find the ", start_col, " column in the annotation data.")
  }
  if (is.null(features[[end_col]])) {
    stop("Unable to find the ", end_col, " column in the annotation data.")
  }
  features[[start_col]] <- sm(as.numeric(features[[start_col]]))
  na_starts <- is.na(features[[start_col]])
  features <- features[!na_starts, ]
  features[[end_col]] <- as.numeric(features[[end_col]])
  ## Keep in mind that when creating the snp_expt, I removed '_' from
  ## the chromosome names and replaced them with '-'.
  features[[expt_name_col]] <- gsub(x = features[[expt_name_col]],
                                    pattern = "_",
                                    replacement = "-")

  ## Therefore, in order to cross reference, I need to do the same here.
  ## In this invocation, I need the seqnames to be the chromosome of each gene.
  expt_granges <- GenomicRanges::makeGRangesFromDataFrame(features,
                                                          seqnames.field = expt_name_col,
                                                          start.field = start_col,
                                                          end.field = end_col)
  expt_subset <- expt_granges[genes, ]

  snp_annot <- fData(snp_expt)
  snp_annot[[snp_start_col]] <- as.numeric(snp_annot[[snp_start_col]])
  snp_annot[["end"]] <- snp_annot[[snp_start_col]] + 1
  snp_granges <- GenomicRanges::makeGRangesFromDataFrame(snp_annot,
                                                         seqnames.field = snp_name_col,
                                                         start.field = snp_start_col,
                                                         end.field = "end")

  snps_in_genes <- IRanges::subsetByOverlaps(snp_granges, expt_subset,
                                             type = "within", ignore.strand = TRUE)
  snp_expt_ids <- names(snps_in_genes)
  snp_subset <- exclude_genes_expt(snp_expt, ids = snp_expt_ids, method = "keep")

  snp_fdata <- as.data.frame(snps_in_genes)
  expt_fdata <- as.data.frame(expt_subset)
  expt_fdata[["gene"]] <- rownames(expt_fdata)
  expt_fdata <- expt_fdata[, c("seqnames", "gene")]
  snp_fdata <- merge(snp_fdata, expt_fdata, by = "seqnames", all.x = TRUE)
  snp_fdata <- merge(snp_fdata, features, by.x = "gene", by.y = expt_gid_column, all.x = TRUE)
  fData(snp_subset[["expressionset"]]) <- snp_fdata
  return(snp_subset)
}

#' Make a summary of the observed snps by gene ID.
#'
#' Instead of cross referencing variant positions against experimental
#' condition, one might be interested in seeing what variants are observed per
#' gene.  This function attempts to answer that question.
#'
#' @param expt The original expressionset.
#' @param snp_result The result from get_snp_sets().
#' @param start_col Which column provides the start of each gene?
#' @param end_col and the end column of each gene?
#' @param snp_name_col Name of the column in the metadata with the sequence names.
#' @param expt_name_col Name of the metadata column with the chromosome names.
#' @param ignore_strand Ignore strand information when returning?
#' @return List with some information by gene.
#' @seealso [GenomicRanges::makeGRangesFromDataFrame()] [IRanges::subsetByOverlaps()]
#'  [IRanges::mergeByOverlaps()] [IRanges::countOverlaps()]
#' @examples
#'  \dontrun{
#'  expt <- create_expt(metadata, gene_information)
#'  snp_expt <- count_expt_snps(expt)
#'  snp_result <- get_snp_sets(snp_expt)
#'  gene_intersections <- snps_vs_genes(expt, snp_result)
#' }
#' @export
#' @importFrom S4Vectors mcols mcols<-
snps_vs_genes <- function(expt, snp_result, start_col = "start", end_col = "end",
                          snp_name_col = "seqnames", observed_in = NULL,
                          expt_name_col = "chromosome", ignore_strand = TRUE) {
  features <- fData(expt)
  if (is.null(features[[start_col]])) {
    stop("Unable to find the ", start_col, " column in the annotation data.")
  }
  if (is.null(features[[end_col]])) {
    stop("Unable to find the ", end_col, " column in the annotation data.")
  }
  features[[start_col]] <- sm(as.numeric(features[[start_col]]))
  na_starts <- is.na(features[[start_col]])
  features <- features[!na_starts, ]
  features[[end_col]] <- as.numeric(features[[end_col]])
  ## Keep in mind that when creating the snp_expt, I removed '_' from
  ## the chromosome names and replaced them with '-'.
  ## Therefore, in order to cross reference, I need to do the same here.
  ## I don't quite want 5'/3' UTRs, I just want the coordinates starting with
  ## (either 1 or) the end of the last gene and ending with the beginning of the
  ## current gene with respect to the beginning of each chromosome.
  ## That is a weirdly difficult problem for creatures with more than 1 chromosome.
  ## inter_features <- features[, c("start", "end", "seqnames")]
  ## inter_features[["chr_start"]] <- paste0(inter_features[["seqnames"]], "_",
  ##                                         inter_features[["start"]])
  ## inter_feature_order <- order(inter_features[["chr_start"]])
  ## inter_features <- inter_features[inter_feature_order, ]

  ## In this invocation, I need the seqnames to be the chromosome of each gene.
  expt_granges <- GenomicRanges::makeGRangesFromDataFrame(
    features, seqnames.field = expt_name_col,
    start.field = start_col, end.field = end_col)
  ## keep.extra.columns = FALSE
  ## ignore.strand = FALSE
  ## seqinfo = NULL
  ## seqnames.field = c("seqnames","chromosome", "chr", "seqid")
  ## start.field = "start"
  ## end.field = "end"
  ## strand.field = "strand"

  snp_positions <- snp_result[["observations"]]
  observations <- data.frame()
  if (!is.null(observed_in)) {
    observed_idx <- snp_positions[[observed_in]] > 0
    message("variants were observed at ", sum(observed_idx),
            " positions in group ", observed_in, ".")
    observations <- data.frame(row.names = rownames(snp_positions))
    observations[[observed_in]] <- 0
    observations[observed_idx, observed_in] <- 1
  }
  snp_positions[[snp_name_col]] <- gsub(
    pattern = "^chr_(.+)_pos_.+_ref.+_alt.+$",
    replacement = "\\1", x = rownames(snp_positions))
  snp_positions[[start_col]] <- as.numeric(
      gsub(pattern = "^chr_.+_pos_(.+)_ref.+_alt.+$",
           replacement = "\\1", x = rownames(snp_positions)))
  snp_positions[[end_col]] <- snp_positions[[start_col]]
  snp_positions[["strand"]] <- "+"
  snp_positions <- snp_positions[, c(snp_name_col, start_col, end_col, "strand")]
  ## Keep in mind that when creating the snp_expt, I removed '_' from
  ## the chromosome names and replaced them with '-'.
  snp_positions[[snp_name_col]] <- gsub(pattern = "-", replacement = "_",
                                        x = snp_positions[[snp_name_col]])
  snp_granges <- GenomicRanges::makeGRangesFromDataFrame(
    snp_positions, seqnames.field = snp_name_col,
    start.field = start_col, end.field = end_col)

  ## Faking out r cmd check with a couple empty variables which will be used by data.table
  seqnames <- count <- NULL
  ## This is how one sets the metadata for a GRanges thing.
  ## When doing mergeByOverlaps, countOverlaps, etc, this is useful.
  ## mcols(object)$column_name <- some data column
  mcols(expt_granges)[, "gene_name"] <- names(expt_granges)

  ## Lets add metadata columns for each column for the medians table
  ## This will let us find the positions unique to a condition.
  mcols(snp_granges)[, "snp_name"] <- names(snp_granges)
  snp_columns <- colnames(snp_result[["observations"]])
  for (c in seq_along(snp_columns)) {
    colname <- snp_columns[c]
    mcols(snp_granges)[, colname] <- snp_result[["observations"]][[colname]]
  }
  message("The snp grange data has ", length(snp_granges), " elements.")
  if (!is.null(observed_in)) {
    observed_snp_idx <- mcols(snp_granges)[[observed_in]] > 0
    message("The set observed in ", observed_in, " comprises ",
            sum(observed_snp_idx), " elements.")
    snp_granges <- snp_granges[observed_snp_idx, ]
  }

  snps_by_chr <- IRanges::subsetByOverlaps(snp_granges, expt_granges,
                                           type = "within", ignore.strand = ignore_strand)
  message("There are ", length(snps_by_chr), " overlapping variants and genes.")

  summarized_by_chr <- data.table::as.data.table(snps_by_chr)
  .N <- NULL  ## .N is a read-only symbol in data.table
  summarized_by_chr[, count := .N, by = list(seqnames)]

  ## I think I can replace this data table invocation with countOverlaps...
  ## Ahh no, the following invocation merely counts which snps are found in name,
  ## which is sort of the opposite of what I want.
  ## test <- IRanges::countOverlaps(query = snp_granges, subject = expt_granges,
  ##                               type = "within", ignore.strand = TRUE)
  summarized_by_chr <- unique(summarized_by_chr[, c("seqnames", "count"), with = FALSE])
  ## The ignore.strand is super important for this task.
  merged_grange <- IRanges::mergeByOverlaps(query = snp_granges, subject = expt_granges,
                                            ignore.strand = ignore_strand)

  count_by_gene_irange <- IRanges::countOverlaps(query = expt_granges, subject = snp_granges,
                                                 type = "any", ignore.strand = ignore_strand)

  ## I am getting odd results using countOverlaps,
  ## lets get a second opinion using dplyr and tally()
  second_opinion <- data.frame("gene" = merged_grange[["gene_name"]],
                               "snp" = merged_grange[["snp_name"]])
  count_by_gene_dplyr <- second_opinion %>%
    group_by(.data[["gene"]]) %>%
    dplyr::tally()
  count_by_gene_dplyr_names <- count_by_gene_dplyr[["gene"]]
  count_by_gene_dplyr <- count_by_gene_dplyr[["n"]]
  names(count_by_gene_dplyr) <- count_by_gene_dplyr_names
  summarized_idx <- order(count_by_gene_irange, decreasing = TRUE)
  count_by_gene_irange <- count_by_gene_irange[summarized_idx]
  summarized_idx <- order(count_by_gene_dplyr, decreasing = TRUE)
  count_by_gene_dplyr <- count_by_gene_dplyr[summarized_idx]
  retlist <- list(
      "expt_granges" = expt_granges,
      "snp_granges" = snp_granges,
      "snps_by_chr" = snps_by_chr,
      "merged_by_gene" = merged_grange,
      "count_by_gene" = count_by_gene_irange,
      "count_by_gene_dplyr" = count_by_gene_dplyr,
      "summary" = summarized_by_chr)
  class(retlist) <- "snps_genes"
  return(retlist)
}

#' A copy of the above function with padding for species without defined UTRs
#'
#' @param expt The original expressionset.
#' @param snp_result The result from get_snp_sets().
#' @param start_col Which column provides the start of each gene?
#' @param end_col and the end column of each gene?
#' @param strand_col Define strands.
#' @param padding Add this amount to each CDS.
#' @param normalize Normalize the returns to the length of the putative CDS.
#' @param snp_name_col Name of the column in the metadata with the sequence names.
#' @param expt_name_col Name of the metadata column with the chromosome names.
#' @param observed_in Print some information about how many variants were observed.
#' @param ignore_strand Ignore the strand information when returning?
#' @return List with some information by gene.
#' @seealso [GenomicRanges::makeGRangesFromDataFrame()] [IRanges::subsetByOverlaps()]
#'  [IRanges::mergeByOverlaps()] [IRanges::countOverlaps()]
#' @examples
#'  \dontrun{
#'  expt <- create_expt(metadata, gene_information)
#'  snp_expt <- count_expt_snps(expt)
#'  snp_result <- get_snp_sets(snp_expt)
#'  gene_intersections <- snps_vs_genes(expt, snp_result)
#' }
#' @export
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom IRanges %over%
snps_vs_genes_padded <- function(expt, snp_result, start_col = "start", end_col = "end",
                          strand_col = "strand", padding = 200, normalize = TRUE,
                          snp_name_col = "seqnames", expt_name_col = "chromosome",
                          observed_in = NULL, ignore_strand = TRUE) {
  features <- fData(expt)
  if (is.null(features[[start_col]])) {
    stop("Unable to find the ", start_col, " column in the annotation data.")
  }
  if (is.null(features[[end_col]])) {
    stop("Unable to find the ", end_col, " column in the annotation data.")
  }
  if (is.null(features[[strand_col]])) {
    stop("Unable to find the ", strand_col, " column in the annotation data.")
  }
  features[[start_col]] <- sm(as.numeric(features[[start_col]]))
  na_starts <- is.na(features[[start_col]])
  features <- features[!na_starts, ]
  features[[end_col]] <- as.numeric(features[[end_col]])
  ## Keep in mind that when creating the snp_expt, I removed '_' from
  ## the chromosome names and replaced them with '-'.
  ## Therefore, in order to cross reference, I need to do the same here.
  ## I don't quite want 5'/3' UTRs, I just want the coordinates starting with
  ## (either 1 or) the end of the last gene and ending with the beginning of the
  ## current gene with respect to the beginning of each chromosome.
  ## That is a weirdly difficult problem for creatures with more than 1 chromosome.
  ## inter_features <- features[, c("start", "end", "seqnames")]
  ## inter_features[["chr_start"]] <- paste0(inter_features[["seqnames"]], "_",
  ##                                         inter_features[["start"]])
  ## inter_feature_order <- order(inter_features[["chr_start"]])
  ## inter_features <- inter_features[inter_feature_order, ]

  plus_idx <- features[[strand_col]] == "+" | features[[strand_col]] > 0 |
    features[[strand_col]] == "forward"
  minus_idx <- ! plus_idx
  plus_features <- features[plus_idx, ]
  minus_features <- features[minus_idx, ]
  plus_features[["5p_start"]] <- plus_features[[start_col]] - padding
  neg_idx <- plus_features[["5p_start"]] < 0
  if (sum(neg_idx) > 0) {
    message("There are ", sum(neg_idx),
            " genes with less than 200 nt. before the start of the chromosome on the plus strand")
    plus_features[neg_idx, "5p_start"] <- 0
  }
  plus_features[["3p_end"]] <- plus_features[[end_col]] + padding
  ## I will need to extract the chromosome lengths to boundary check this.
  minus_features[["5p_start"]] <- minus_features[[end_col]] + padding
  minus_features[["3p_end"]] <- minus_features[[start_col]] - padding
  neg_idx <- plus_features[["3p_end"]] < 0
  if (sum(neg_idx) > 0) {
    message("There are ", sum(neg_idx),
            " genes with less than 200 nt. before the start of the chromosome on the minus strand")
    minus_features[neg_idx, "3p_end"] <- 0
  }
  message("There are ", sum(plus_idx), " plus strand features and ", sum(minus_idx),
          " minus strand features.")

  plus_5p_granges <- GenomicRanges::makeGRangesFromDataFrame(
    plus_features, seqnames.field = expt_name_col,
    start.field = "5p_start", end.field = start_col)
  plus_3p_granges <- GenomicRanges::makeGRangesFromDataFrame(
    plus_features, seqnames.field = expt_name_col,
    start.field = end_col, end.field = "3p_end")
  minus_5p_granges <- GenomicRanges::makeGRangesFromDataFrame(
    minus_features, seqnames.field = expt_name_col,
    start.field = end_col, end.field = "5p_start")
  minus_3p_granges <- GenomicRanges::makeGRangesFromDataFrame(
    minus_features, seqnames.field = expt_name_col,
    start.field = "3p_end", end.field = end_col)

  fivep_granges <- c(plus_5p_granges, minus_5p_granges)
  threep_granges <- c(plus_3p_granges, minus_3p_granges)

  snp_positions <- snp_result[["observations"]]
  observations <- data.frame()
    if (!is.null(observed_in)) {
    observed_idx <- snp_positions[[observed_in]] > 0
    message("variants were observed at ", sum(observed_idx),
            " positions in group ", observed_in, ".")
    observations <- data.frame(row.names = rownames(snp_positions))
    observations[[observed_in]] <- 0
    observations[observed_idx, observed_in] <- 1
  }
  snp_positions[[snp_name_col]] <- gsub(
      pattern = "^chr_(.+)_pos_.+_ref.+_alt.+$",
      replacement = "\\1",
      x = rownames(snp_positions))
  snp_positions[[start_col]] <- as.numeric(
      gsub(pattern = "^chr_.+_pos_(.+)_ref.+_alt.+$",
           replacement = "\\1",
           x = rownames(snp_positions)))
  snp_positions[[end_col]] <- snp_positions[[start_col]]
  snp_positions[["strand"]] <- "+"
  snp_positions <- snp_positions[, c(snp_name_col, start_col, end_col, "strand")]
  ## Keep in mind that when creating the snp_expt, I removed '_' from
  ## the chromosome names and replaced them with '-'.
  snp_positions[[snp_name_col]] <- gsub(pattern = "-", replacement = "_",
                                        x = snp_positions[[snp_name_col]])
  snp_granges <- GenomicRanges::makeGRangesFromDataFrame(
    snp_positions, seqnames.field = snp_name_col,
    start.field = start_col, end.field = end_col)

  ## Faking out r cmd check with a couple empty variables which will be used by data.table
  seqnames <- count <- NULL
  ## This is how one sets the metadata for a GRanges thing.
  ## When doing mergeByOverlaps, countOverlaps, etc, this is useful.
  ## mcols(object)$column_name <- some data column
  mcols(fivep_granges)[, "gene_name"] <- names(fivep_granges)
  mcols(threep_granges)[, "gene_name"] <- names(threep_granges)

  ## Lets add metadata columns for each column for the medians table
  ## This will let us find the positions unique to a condition.
  mcols(snp_granges)[, "snp_name"] <- names(snp_granges)
  snp_columns <- colnames(snp_result[["observations"]])
  for (c in seq_along(snp_columns)) {
    colname <- snp_columns[c]
    mcols(snp_granges)[, colname] <- snp_result[["observations"]][[colname]]
  }
  message("The snp grange data has ", length(snp_granges), " elements.")
  if (!is.null(observed_in)) {
    observed_snp_idx <- mcols(snp_granges)[[observed_in]] > 0
    message("The set observed in ", observed_in, " comprises ",
            sum(observed_snp_idx), " elements.")
    snp_granges <- snp_granges[observed_snp_idx, ]
  }

  snps_by_fivep <- IRanges::subsetByOverlaps(snp_granges, fivep_granges,
                                             type = "within", ignore.strand = ignore_strand)
  message("There are ", length(snps_by_fivep), " overlapping variants and 5' padded UTRs.")
  summarized_fivep_by_chr <- data.table::as.data.table(snps_by_fivep)
  snps_by_threep <- IRanges::subsetByOverlaps(snp_granges, threep_granges,
                                              type = "within", ignore.strand = ignore_strand)
  message("There are ", length(snps_by_threep), " overlapping variants and 3' padded UTRs.")
  summarized_threep_by_chr <- data.table::as.data.table(snps_by_threep)

  .N <- NULL  ## .N is a read-only symbol in data.table
  summarized_fivep_by_chr[, count := .N, by = list(seqnames)]
  summarized_threep_by_chr[, count := .N, by = list(seqnames)]
  ## I think I can replace this data table invocation with countOverlaps...
  ## Ahh no, the following invocation merely counts which snps are found in name,
  ## which is sort of the opposite of what I want.
  ## test <- IRanges::countOverlaps(query = snp_granges, subject = expt_granges,
  ##                               type = "within", ignore.strand = TRUE)
  snp_ranges <- unique(snp_granges)
  summarized_fivep_by_chr <- unique(summarized_fivep_by_chr[, c("seqnames", "count"),
                                                            with = FALSE])
  summarized_threep_by_chr <- unique(summarized_threep_by_chr[, c("seqnames", "count"),
                                                              with = FALSE])
  merged_fivep_grange <- IRanges::mergeByOverlaps(query = snp_ranges, subject = fivep_granges,
                                                  ignore.strand = ignore_strand)
  merged_threep_grange <- IRanges::mergeByOverlaps(query = snp_ranges, subject = threep_granges,
                                                   ignore.strand = ignore_strand)

  summarized_fivep_by_gene <- IRanges::countOverlaps(
    query = fivep_granges, subject = snp_granges, type = "any", ignore.strand = ignore_strand)
  tt <- fivep_granges %over% snp_ranges
  summarized_threep_by_gene <- IRanges::countOverlaps(
    query = threep_granges, subject = snp_granges, type = "any", ignore.strand = ignore_strand)
  summarized_fivep_idx <- order(summarized_fivep_by_gene, decreasing = TRUE)
  summarized_fivep_by_gene <- summarized_fivep_by_gene[summarized_fivep_idx]
  summarized_threep_idx <- order(summarized_threep_by_gene, decreasing = TRUE)
  summarized_threep_by_gene <- summarized_threep_by_gene[summarized_threep_idx]

  ## Normalize in this context just means dividing the numbers by the padding
  ## so that we can directly compare the result to the normalized by genelength
  ## cds values and/or other padding lengths and/or a real UTR feature set.
  if (isTRUE(normalize)) {
    summarized_fivep_by_gene <- summarized_fivep_by_gene / padding
    summarized_threep_by_gene <- summarized_threep_by_gene / padding
    summarized_fivep_by_chr <- summarized_fivep_by_chr / padding
    summarized_threep_by_chr <- summarized_threep_by_chr / padding
  }
  retlist <- list(
    "fivep_granges" = fivep_granges,
    "threep_granges" = threep_granges,
    "snp_granges" = snp_granges,
    "snps_by_fivep" = snps_by_fivep,
    "snps_by_threep" = snps_by_threep,
    "merged_by_fivep" = merged_fivep_grange,
    "merged_by_threep" = merged_threep_grange,
    "count_fivep_by_gene" = summarized_fivep_by_gene,
    "count_threep_by_gene" = summarized_threep_by_gene,
    "summary_fivep" = summarized_fivep_by_chr,
    "summary_threep" = summarized_threep_by_chr)
  return(retlist)
}

#' Write a matrix of variants in an alignment-esque format.
#'
#' @param expt variant expressionset.
#' @param output_file File to write, presumably to be passed to something like phyML.
write_snps <- function(expt, output_file = "funky.aln") {
  start_mtrx <- exprs(expt)
  samples <- colnames(start_mtrx)
  aln_string <- ""
  for (r in seq_len(ncol(start_mtrx))) {
    aln_string <- glue::glue("{aln_string}{samples[r]}     {start_mtrx[[r]]}
")
  }

  write_output <- file(output_file, open = "a+")
  cat(aln_string, file = write_output, sep = "")
  close(write_output)
  return(output_file)
}

#' If I were smart I would use an I/GRanges for this.
#'
#' But I was asked to get the closest feature if it is not inside one.
#' I am not sure how to do that with a ranges. Sadly, I think it will
#' be easier for me to just iterate over the sequence_df and query
#' each feature on that chromosome/scaffold.
#'
#' @param sequence_df dataframe of sequence regions of interest.
#' @param gff gff annotations against which to hunt.
#' @param bin_width size of the regions of interest (e.g. the amplicon size)
#' @param feature_type What feature type to hunt for?
#' @param feature_start Column containing the starts.
#' @param feature_end Column containing the ends.
#' @param feature_strand Column containing strand information.
#' @param feature_chr Column containing the chromosome names.
#' @param feature_type_column Column containing the feature types.
#' @param feature_id Column with the IDs (coming from the gff tags).
#' @param feature_name Column with the descriptive name.
xref_regions <- function(sequence_df, gff, bin_width = 600,
                         feature_type = "protein_coding_gene", feature_start = "start",
                         feature_end = "end", feature_strand = "strand",
                         feature_chr = "seqnames", feature_type_column = "type",
                         feature_id = "ID", feature_name = "description") {
  annotation <- load_gff_annotations(gff)
  wanted_columns <- c(feature_id, feature_name, feature_chr,
                      feature_start, feature_end, feature_strand)
  wanted_rows <- annotation[[feature_type_column]] == feature_type
  annotation <- annotation[wanted_rows, wanted_columns]
  colnames(annotation) <- c("id", "name", "chr", "start", "end", "strand")
  mesg("Now the annotation has ", nrow(annotation),  " rows.")

  sequence_df[["overlap_gene_id"]] <- ""
  sequence_df[["overlap_gene_description"]] <- ""
  sequence_df[["overlap_gene_start"]] <- ""
  sequence_df[["overlap_gene_end"]] <- ""
  sequence_df[["closest_gene_before_id"]] <- ""
  sequence_df[["closest_gene_before_description"]] <- ""
  sequence_df[["closest_gene_before_start"]] <- ""
  sequence_df[["closest_gene_before_end"]] <- ""
  sequence_df[["closest_gene_after_id"]] <- ""
  sequence_df[["closest_gene_after_description"]] <- ""
  sequence_df[["closest_gene_after_start"]] <- ""
  sequence_df[["closest_gene_after_end"]] <- ""
  for (r in seq_along(nrow(sequence_df))) {
    entry <- sequence_df[r, ]
    amplicon_chr <- entry[["chr"]]
    amplicon_start <- as.numeric(entry[["bin_start"]])
    amplicon_end <- amplicon_start + bin_width
    annot_idx <- annotation[["chr"]] == amplicon_chr
    xref_features <- annotation[annot_idx, ]
    if (nrow(xref_features) == 0) {
      next
    }

    ## Amplicon:           ------->
    ##                  ############# fs < qs && fe > qe
    ## Features:  ##     ### query_start > feature_start &&
    ##                               query_start < feature end
    ##                       ##  qstart < fstart && fend > qend
    ##                            ### qend > fstart && qend < fend

    ## So there are 4 obvious scenarios we want to collect
    ## When the feature is surrounding the amplicon
    ## The feature overlaps the 5' of the amplicon rosalind
    ## The feature overlaps the 3' of the amplicon rosalind
    ## The feature is entirely inside the amplicon

    ## If all those fail, we want to find the closest feature.
    ## Amplicon:             ------->
    ## Get the left one:  ##                ##

    hits <- 0
    hit_df <- data.frame()
    primer_inside_feature <- xref_features[["start"]] <= amplicon_start &
      xref_features[["end"]] >= amplicon_end
    if (sum(primer_inside_feature) > 0) {
      hits <- hits + sum(primer_inside_feature)
      sequence_df[r, "overlap_gene_id"] <- xref_features[primer_inside_feature, "id"]
      sequence_df[r, "overlap_gene_description"] <- xref_features[primer_inside_feature, "name"]
      sequence_df[r, "overlap_gene_start"] <- xref_features[primer_inside_feature, "start"]
      sequence_df[r, "overlap_gene_end"] <- xref_features[primer_inside_feature, "end"]
    }

    ##     ### amplicon_start > feature_start &&
    primer_overlap_fivep <- xref_features[["start"]] <= amplicon_start &
      xref_features[["end"]] >= amplicon_start
    if (sum(primer_overlap_fivep) > 0) {
      hits <- hits + sum(primer_overlap_fivep)
      sequence_df[r, "overlap_gene_id"] <- xref_features[primer_overlap_fivep, "id"]
      sequence_df[r, "overlap_gene_description"] <- xref_features[primer_overlap_fivep, "name"]
      sequence_df[r, "overlap_gene_start"] <- xref_features[primer_overlap_fivep, "start"]
      sequence_df[r, "overlap_gene_end"] <- xref_features[primer_overlap_fivep, "end"]
    }

    ##  Amplicon:           ------->
    ##                           ##### qend > fstart && qend < fend
    primer_overlap_threep <- xref_features[["start"]] <= amplicon_end &
      xref_features[["end"]] >= amplicon_end
    if (sum(primer_overlap_threep) > 0) {
      hits <- hits + sum(primer_overlap_threep)
      sequence_df[r, "overlap_gene_id"] <- xref_features[primer_overlap_threep, "id"]
      sequence_df[r, "overlap_gene_description"] <- xref_features[primer_overlap_threep, "name"]
      sequence_df[r, "overlap_gene_start"] <- xref_features[primer_overlap_threep, "start"]
      sequence_df[r, "overlap_gene_end"] <- xref_features[primer_overlap_threep, "end"]
    }

    ##  Amplicon:           ------->
    ##                        ##   qstart < fstart & qend > fend
    primer_overlap_surround <- xref_features[["start"]] >= amplicon_end &
      xref_features[["end"]] <= amplicon_end
    if (sum(primer_overlap_surround) > 0) {
      hits <- hits + sum(primer_overlap_surround)
      sequence_df[r, "overlap_gene_id"] <- xref_features[primer_overlap_surround, "id"]
      sequence_df[r, "overlap_gene_description"] <- xref_features[primer_overlap_surround, "name"]
      sequence_df[r, "overlap_gene_start"] <- xref_features[primer_overlap_surround, "start"]
      sequence_df[r, "overlap_gene_end"] <- xref_features[primer_overlap_surround, "end"]
    }

    ## Now report the closest genes before/after the amplicon
    xref_features[["fivep_dist"]] <- amplicon_start - xref_features[["end"]]
    fivep_wanted_idx <- xref_features[["fivep_dist"]] > 0
    if (sum(fivep_wanted_idx) > 0) {
      fivep_candidates <- xref_features[fivep_wanted_idx, ]
      fivep_min_idx <- fivep_candidates[["fivep_dist"]] == min(fivep_candidates[["fivep_dist"]])
      fivep_candidate <- fivep_candidates[fivep_min_idx, ]
      sequence_df[r, "closest_gene_before_id"] <- fivep_candidate["id"]
      sequence_df[r, "closest_gene_before_description"] <- fivep_candidate["name"]
      sequence_df[r, "closest_gene_before_start"] <- fivep_candidate["start"]
      sequence_df[r, "closest_gene_before_end"] <- fivep_candidate["end"]
    }

    xref_features[["threep_dist"]] <- xref_features[["start"]] - amplicon_end
    threep_wanted_idx <- xref_features[["threep_dist"]] > 0
    if (sum(threep_wanted_idx) > 0) {
      threep_candidates <- xref_features[threep_wanted_idx, ]
      threep_min_idx <- threep_candidates[["threep_dist"]] == min(threep_candidates[["threep_dist"]])
      threep_candidate <- threep_candidates[threep_min_idx, ]
      sequence_df[r, "closest_gene_after_id"] <- threep_candidate["id"]
      sequence_df[r, "closest_gene_after_description"] <- threep_candidate["name"]
      sequence_df[r, "closest_gene_after_start"] <- threep_candidate["start"]
      sequence_df[r, "closest_gene_after_end"] <- threep_candidate["end"]
    }

  } ## End iterating over every row of the sequence df.
  return(sequence_df)
}

## EOF
