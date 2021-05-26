#' Gather snp information for an expt
#'
#' This function attempts to gather a set of variant positions using an extant
#' expressionset.  This therefore seeks to keep the sample metadata consistent
#' with the original data.  In its current iteration, it therefore makes some
#' potentially bad assumptions about the naming conventions for its input
#' files.  It furthermore assumes inputs from the variant calling methods in
#' cyoa.
#'
#' @param expt an expressionset from which to extract information.
#' @param type Use counts / samples or ratios?
#' @param annot_column Column in the metadata for getting the table of bcftools calls.
#' @param tolower Lowercase stuff like 'HPGL'?
#' @param snp_column Which column of the parsed bcf table contains our interesting material?
#' @return A new expt object
#' @seealso [Biobase]
#' @examples
#'   \dontrun{
#'  expt <- create_expt(metadata, gene_information)
#'  snp_expt <- count_expt_snps(expt)
#'  ## This assumes that the metadata has a column named 'bcftable' with one file per
#'  ## cell.  These files in turn should have a column named 'diff_count' which will
#'  ## be the source of the numbers found when doing exprs(snp_expt).
#' }
#' @export
count_expt_snps <- function(expt, type = "counts", annot_column = "bcftable", tolower = TRUE,
                            snp_column = "diff_count") {
  samples <- rownames(pData(expt))
  if (isTRUE(tolower)) {
    samples <- tolower(samples)
  }
  file_lst <- pData(expt)[[annot_column]]
  if (is.null(file_lst)) {
    stop("This requires a set of bcf filenames, the column: ", annot_column, " does not have any.")
  }
  ## Create a data table of snp columns
  chosen_column <- snp_column
  if (is.null(chosen_column)) {
    if (type == "percent") {
      message("Making a matrix of percentages.")
      chosen_column <- "pct"
    } else {
      message("Making a matrix of counts.")
      chosen_column <- "diff_count"
    }
  }
  snp_dt <- read_snp_columns(samples, file_lst, column = chosen_column)
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
#' @param snp_expt The result of count_expt_snps()
#' @param factor Experimental factor to use for cutting and splicing the data.
#' @param limit Minimum median number of hits / factor to define a position as
#'  a hit.
#' @param do_save Save the result?
#' @param savefile Prefix for a savefile if one chooses to save the result.
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
get_snp_sets <- function(snp_expt, factor = "pathogenstrain", limit = 1,
                         do_save = FALSE, savefile = "variants.rda") {
  if (is.null(snp_expt[["design"]][[factor]])) {
    stop("The factor does not exist in the expt.")
  }
  if (isTRUE(do_save) &  file.exists(glue("{savefile}_{factor}.rda"))) {
    retlist <- new.env()
    loaded <- load(savefile, envir = retlist)
    retlist <- retlist[["retlist"]]
    return(retlist)
  }

  medians <- median_by_factor(snp_expt, fact = factor)[["medians"]]
  ## I am going to split this by chromosome, as a run of 10,000 took 2 seconds,
  ## 100,000 took a minute, and 400,000 took an hour.
  ##chr <- gsub(pattern = "^.+_(.+)_.+_.+_.+$", replacement = "\\1", x = rownames(medians))
  chr <- gsub(pattern = "^chr_(.+)_pos_.+_ref_.+_alt_.+$", replacement = "\\1", x = rownames(medians))
  medians[["chr"]] <- chr
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
                   returns[[chromosome_name]] <- snp_by_chr(medians,
                                                            chromosome_name)
                 }
  if (isTRUE(show_progress)) {
    close(bar)
  }
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
  message("Iterating over ", end, " elements.")
  for (element in 1:end) {
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
    snps <- rownames(data_by_chr[[chr]][["medians"]])
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
    "medians" = medians,
    "possibilities" = possibilities,
    "intersections" = all_intersections,
    "chr_data" = data_by_chr,
    "set_names" = set_names,
    "invert_names" = invert_names,
    "density" = density_by_chr)
  if (isTRUE(do_save)) {
    saved <- save(list = "retlist", file = savefile)
  }
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
#' @seealso [readr]
#' @return A big honking data table.
read_snp_columns <- function(samples, file_lst, column = "diff_count") {
  ## Read the first file
  first_sample <- samples[1]
  first_file <- file_lst[1]
  first_read <- readr::read_tsv(first_file)
  ## Create a simplified data table from it.
  first_column <- data.table::as.data.table(first_read[[column]])
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
  for (sample_num in 2:length(samples)) {
    if (isTRUE(show_progress)) {
      pct_done <- sample_num / length(samples)
      setTxtProgressBar(bar, pct_done)
    }
    sample <- samples[sample_num]
    file <- file_lst[sample_num]
    new_table <- try(readr::read_tsv(file))
    if (class(new_table)[1] == "try-error") {
      next
    }
    new_column <- data.table::as.data.table(new_table[[column]])
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

#' Use Rsamtools to read alignments and get snp coverage.
#'
#' This is horrifyingly slow.  I think I might remove this function.
#'
#' @param expt Expressionset to analyze
#' @param type counts or percent?
#' @param input_dir Directory containing the samtools results.
#' @param annot_column Passed along to count_expt_snps()
#' @param tolower lowercase the sample names?
#' @param bam_suffix In case the data came from sam.
#' @return It is so slow I no longer know if it works.
#' @seealso [count_expt_snps()] [Rsamtools] [GenomicRanges]
samtools_snp_coverage <- function(expt, type = "counts", input_dir = "preprocessing/outputs",
                                  tolower = TRUE, bam_suffix = ".bam", annot_column = annot_column) {
  snp_counts <- count_expt_snps(expt, type = type, tolower = tolower)
  snp_counts <- fData(snp_counts)
  samples <- rownames(pData(expt))
  if (isTRUE(tolower)) {
    samples <- tolower(samples)
  }
  if (type == "counts") {
    file_suffix <- "_parsed_count.txt"
  } else {
    file_suffix <- "_parsed_ratio.txt"
  }
  bam_filenames <- c()
  for (sample_num in 1:length(samples)) {
    bam_filenames[sample_num] <- paste0(file.path(input_dir, sample), bam_suffix)
    if (!file.exists(bam_filenames[sample_num])) {
      warning(glue("The bam filename does not exist: {bam_filenames[sample_num]}"))
    }
  }

  ## Now I have a data table of rownames and percentages
  ## Next I need to cross reference these against the coverage by position, ergo
  ## I must split out the rownames
  pileup_info <- Rsamtools::PileupFiles(bam_filenames)
  ## Taken directly from the Rsamtools manual
  queries <- glue("{snp_counts[['species']]}_{snp_counts[['chromosome']]}:\\
                  {as.numeric(snp_counts[['position']]) - 1}-\\
                  {as.numeric(snp_counts[['position']]) + 1}")
  which <- GenomicRanges::GRanges(queries)
  which_list <- split(which, GenomicRanges::seqnames(which))
  returns <- list()
  res <- list()
  what <- "seq"

  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  tt <- sm(requireNamespace("parallel"))
  tt <- sm(requireNamespace("doParallel"))
  tt <- sm(requireNamespace("iterators"))
  tt <- sm(requireNamespace("foreach"))
  res <- foreach(c = 1:length(which_list), .packages = c("hpgltools", "Rsamtools")) %dopar% {
    chr_name <- names(which_list)[c]
    chr_param <- Rsamtools::ApplyPileupsParam(which = which_list[[c]], what = what)
    snp_calc_coverage <- function(x) {
      ## information at each pile-up position
      qme <- function(y) {
        y <- y[c("A", "C", "G", "T"), , drop = FALSE]
        y <- y + 1L
        result <- colSums(y)
        return(result)
      }
      info <- apply(x[["seq"]], 2, qme)
      retlist <- list(seqnames = x[["seqnames"]], pos = x[["pos"]], info = info)
      return(retlist)
    }
    coverage_result <- Rsamtools::applyPileups(
                                    pileup_info,
                                    snp_calc_coverage,
                                    param = chr_param)
    result_list[[c]] <- coverage_result
  }
  stopped <- parallel::stopCluster(cl)

  ## result is a list of n elements where n is the number of rows in snp_dt
  ## Each element of result is in turn a list containing the following slots:
  ##  seqnames (chromosome), pos (position(s)), info (coverage by file)
  ## The piece of information we want to put into snp_dt is therefore:
  ## coverage_list <- result[[snp_dt_row]][[info]][2, ]
  ## coverage_list is in turn a character list named by filename (which begins
  ## with the sample ID) We will therefore extract the hpglID from it and the
  ## coverage for every sample.

  ## I am not sure if the following line is correct, but I think it is -- either
  ## way, snp_dt is missing without it.  Somewhere along the way I forgot to
  ## make sure to keep the variable snp_dt alive, the most logical existence for
  ## it is as a recasting of the matrix resulting from exprs(snp_counts)
  ## For a reminder of why I say this, check out the first 15 lines of count_expt_snps()
  snp_dt <- data.table::as.data.table(exprs(snp_counts))
  names(coverage_result) <- snp_dt[["rownames"]]

  ## Now extract from the rather strange coverage_result data the coverage by position/sample
  new_dt <- NULL
  snp_extract_coverage <- function(element) {
    row <- element[["info"]][2, ]
    names(row) <- samples
    new_dt <<- rbind(new_dt, row)
  }
  ## unused_var <- try(lapply(coverage_result, snp_extract_coverage))
  unused_var <- try(parallel::mclapply(coverage_result, snp_extract_coverage))
  if (class(unused_var) == "try-error") {
    message("There was an error when creating the data table of coverage.")
  }
  new_dt <- data.table::as.data.table(new_dt)
  new_dt[["rownames"]] <- snp_dt[["rownames"]]

}

#' The real worker.  This extracts positions for a single chromosome and puts
#' them into a parallelizable data structure.
#'
#' @param medians A set of medians by position to look through
#' @param chr_name Chromosome name to search
#' @param limit Minimum number of median hits/position to count as a snp.
#' @return A list of variant positions where each element is one chromosome.
#' @seealso [Vennerable]
snp_by_chr <- function(medians, chr_name = "01", limit = 1) {
  set_names <- list()
  possibilities <- c()
  count <- 0
  kept_rows <- medians[["chr"]] == chr_name
  medians <- medians[kept_rows, ]
  kept_cols <- "chr" != colnames(medians)
  medians <- medians[, kept_cols]
  data_by_chr <- list()
  data_by_chr[["chromosome"]] <- chr_name
  data_by_chr[["medians"]] <- medians
  limit_true_false <- as.data.frame(medians >= limit)
  x_lst <- list()
  possibilities <- unique(c(possibilities, colnames(limit_true_false)))
  for (d in 1:length(possibilities)) {
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
    for (s in 1:length(symbols)) {
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
snps_vs_genes <- function(expt, snp_result, start_col = "start", end_col = "end",
                          snp_name_col = "seqnames", expt_name_col = "chromosome") {
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
  expt_granges <- GenomicRanges::makeGRangesFromDataFrame(features,
                                                          seqnames.field = expt_name_col,
                                                          start.field = start_col,
                                                          end.field = end_col)
  ## keep.extra.columns = FALSE
  ## ignore.strand = FALSE
  ## seqinfo = NULL
  ## seqnames.field = c("seqnames","chromosome", "chr", "seqid")
  ## start.field = "start"
  ## end.field = "end"
  ## strand.field = "strand"

  snp_positions <- snp_result[["medians"]]
  ##snp_positions <- as.data.frame(exprs(snp_result))
  ##snp_positions[["seqnames"]] <- gsub(pattern = "^(.+_.+)_.+_.+_.+$",
  ##                                    replacement = "\\1",
  ##                                    x = rownames(snp_positions))
  snp_positions[[snp_name_col]] <- gsub(
    pattern = "^chr_(.+)_pos_.+_ref.+_alt.+$",
    replacement = "\\1",
    x = rownames(snp_positions))
  snp_positions[[start_col]] <- as.numeric(
    gsub(pattern = "^chr_.+_pos_(.+)_ref.+_alt.+$",
         replacement = "\\1",
         x = rownames(snp_positions)))
  snp_positions[[end_col]] <- snp_positions[[start_col]] + 1
  snp_positions[["strand"]] <- "+"
  snp_positions <- snp_positions[, c(snp_name_col, start_col, end_col, "strand")]
  ## Keep in mind that when creating the snp_expt, I removed '_' from
  ## the chromosome names and replaced them with '-'.
  snp_positions[[snp_name_col]] <- gsub(pattern = "-", replacement = "_",
                                        x = snp_positions[[snp_name_col]])
  snp_granges <- GenomicRanges::makeGRangesFromDataFrame(snp_positions,
                                                         seqnames.field = snp_name_col,
                                                         start.field = start_col,
                                                         end.field = end_col)

  ## Faking out r cmd check with a couple empty variables which will be used by data.table
  seqnames <- count <- NULL
  ## This is how one sets the metadata for a GRanges thing.
  ## When doing mergeByOverlaps, countOverlaps, etc, this is useful.
  ## mcols(object)$column_name <- some data column
  S4Vectors::mcols(expt_granges)[, "gene_name"] <- names(expt_granges)

  ## Lets add metadata columns for each column for the medians table
  ## This will let us find the positions unique to a condition.
  S4Vectors::mcols(snp_granges)[, "snp_name"] <- names(snp_granges)
  snp_columns <- colnames(snp_result[["medians"]])
  for (c in 1:length(snp_columns)) {
    colname <- snp_columns[c]
    S4Vectors::mcols(snp_granges)[, colname] <- snp_result[["medians"]][[colname]]
  }

  snps_by_chr <- IRanges::subsetByOverlaps(snp_granges, expt_granges,
                                           type = "within", ignore.strand = TRUE)
  summarized_by_chr <- data.table::as.data.table(snps_by_chr)
  .N <- NULL  ## .N is a read-only symbol in data.table
  summarized_by_chr[, count := .N, by = list(seqnames)]
  ## I think I can replace this data table invocation with countOverlaps...
  ## Ahh no, the following invocation merely counts which snps are found in name,
  ## which is sort of the opposite of what I want.
  ## test <- IRanges::countOverlaps(query = snp_granges, subject = expt_granges,
  ##                               type = "within", ignore.strand = TRUE)
  summarized_by_chr <- unique(summarized_by_chr[, c("seqnames", "count"), with = FALSE])
  merged_grange <- IRanges::mergeByOverlaps(query = snp_granges, subject = expt_granges)

  summarized_by_gene <- IRanges::countOverlaps(query = expt_granges, subject = snp_granges,
                                               type = "any", ignore.strand = TRUE)
  summarized_idx <- order(summarized_by_gene, decreasing = TRUE)
  summarized_by_gene <- summarized_by_gene[summarized_idx]
  retlist <- list(
    "expt_granges" = expt_granges,
    "snp_granges" = snp_granges,
    "snps_by_chr" = snps_by_chr,
    "merged_by_gene" = merged_grange,
    "summary_by_gene" = summarized_by_gene,
    "summary" = summarized_by_chr)
  return(retlist)
}

write_snps <- function(expt, output_file = "funky.aln") {
  start_mtrx <- exprs(expt)
  samples <- colnames(start_mtrx)
  aln_string <- ""
  for (r in 1:ncol(start_mtrx)) {
    aln_string <- glue::glue("{aln_string}{samples[r]}     {start_mtrx[[r]]}
")
  }

  write_output <- file(output_file, open = "a+")
  cat(aln_string, file = write_output, sep = "")
  close(write_output)
  return(output_file)
}

## EOF
