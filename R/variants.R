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
  ## Create a data table of snp columns
  if (is.null(snp_column)) {
    stop("This requires a tag column to extract from the freebayes output.")
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
  mesg("Iterating over ", end, " elements.")
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
  first_read <- readr::read_tsv(first_file, show_col_types = FALSE)
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
    new_table <- try(readr::read_tsv(file, show_col_types = FALSE))
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
samtools_snp_coverage <- function(expt, input_dir = "preprocessing/outputs",
                                  tolower = TRUE, bam_suffix = ".bam", annot_column = annot_column) {
  snp_counts <- count_expt_snps(expt, tolower = tolower)
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

#' Create a density function given a variant output and some metadata
#'
#' It is hoped that this will point out regions of a genome which
#' might prove useful when designing PCR primers for a specific
#' condition in a dataset of variants.
#'
#' @param snp_count Result from count_expt_snps()
#' @param pdata_column Metadata column containing the condition of
#'  interest.
#' @param condition Chosen condition to search for variants.
#' @param cutoff Minimum number of variants in a region.
#' @param bin_width Bin size/region of genome to consider.
#' @param divide Normalize by bin width?
#' @param topn Keep only this number of candidates.
#' @param target_temp Try to get primers with this Tm.
#' @param max_primer_length Keep primers at or less than this length.
#' @param bsgenome Genome package containing the sequence of interest.
#' @param gff GFF to define regions of interest.
#' @param feature_type GFF feature type to search against.
#' @param feature_start GFF column with the starts (needed?)
#' @param feature_end GFF column with the ends (needed?)
#' @param feature_strand GFF column with strand information (needed?)
#' @param feature_chr GFF column with chromosome information.
#' @param feature_type_column GFF column with type information.
#' @param feature_id GFF tag with the ID information.
#' @param feature_name GFF tag with the names.
#' @param truncate Truncate the results to just the columns I think
#'  are useful.
#' @export
snp_density_primers <- function(snp_count, pdata_column = "condition",
                                condition = "z2.3", cutoff = 20, bin_width = 600,
                                divide = FALSE, topn = 400,
                                target_temp = 53, max_primer_length = 50,
                                bsgenome = "BSGenome.Leishmania.panamensis.MHOMCOL81L13.v52",
                                gff = "reference/lpanamensis_col_v46.gff",
                                feature_type = "protein_coding_gene", feature_start = "start",
                                feature_end = "end", feature_strand = "strand",
                                feature_chr = "seqnames", feature_type_column = "type",
                                feature_id = "ID", feature_name = "description",
                                truncate = TRUE) {

  ## Start out by loading the bsgenome data
  genome <- NULL
  if (!is.null(bsgenome)) {
    library(bsgenome, character.only=TRUE)
    genome <- get0(bsgenome)
  }

  samples_by_condition <- pData(snp_count)[[pdata_column]]
  ## Keep only those samples of the condition of interest for now,
  ## maybe make it a loop to iterate over conditions later.
  interest_idx <- samples_by_condition == condition
  snp_table <- exprs(snp_count)[, interest_idx]
  start_rows <- nrow(snp_table)
  interest_rows <- rowMeans(snp_table) >= cutoff
  snp_table <- snp_table[interest_rows, ]
  mesg("Started with ", start_rows, ", after cutoff there are ",
          sum(interest_rows), " rows left.")
  position_table <- data.frame(row.names = rownames(snp_table))
  position_table[["chromosome"]] <- gsub(x = rownames(position_table),
                                         pattern = "chr_(.*)_pos_.*",
                                         replacement = "\\1")
  position_table[["position"]] <- gsub(x = rownames(position_table),
                                       pattern = ".*_pos_(\\d+)_.*",
                                       replacement = "\\1")
  ## In the case of lpanamensis at least, chromosomes have '-' in the name, this is not allowed.
  position_table[["chromosome"]] <- gsub(pattern = "-", replacement = "_",
                                         x = position_table[["chromosome"]])
  position_table[["ref"]] <- gsub(x = rownames(position_table),
                                  pattern = ".*_pos_\\d+_ref_(\\w).*",
                                  replacement = "\\1")
  position_table[["alt"]] <- gsub(x = rownames(position_table),
                                  pattern = ".*_pos_\\d+_ref_.+_alt_(\\w)$",
                                  replacement = "\\1")
  ## Lets make a dataframe of chromosomes, their lengths, and number of variants/chromosome
  chromosomes <- levels(as.factor(position_table[["chromosome"]]))

  chromosome_df <- data.frame(row.names = chromosomes)
  chromosome_df[["length"]] <- 0
  chromosome_df[["variants"]] <- 0
  density_lst <- list()
  variant_lst <- list()
  ## This was written in an attempt to make it work if one does or does not
  ## have a BSgenome from which to get the actual chromosome lengths.
  for (ch in 1:length(chromosomes)) {
    chr <- chromosomes[ch]
    mesg("Starting chromosome: ", chr, ".")
    chr_idx <- position_table[["chromosome"]] == chr
    chr_data <- position_table[chr_idx, ]
    if (is.null(genome)) {
      chromosome_df[chr, "length"] <- max(chr_data[["position"]])
    } else {
      this_length <- length(genome[[chr]])
      chromosome_df[chr, "length"] <- this_length
    }
    chromosome_df[chr, "variants"] <- nrow(chr_data)

    ## Stop scientific notation for this operation
    current_options <- options(scipen = 999)
    density_vector <- rep(0, ceiling(this_length / bin_width))
    variant_vector <- rep('', ceiling(this_length / bin_width))
    density_name <- 0
    for (i in 1:length(density_vector)) {
      range_min <- density_name
      range_max <- (density_name + bin_width) - 1
      ## We want to subtract the maximum primer length from the position
      ## and have each extension reaction start near/at the bin.
      chr_data[["position"]] <- as.numeric(chr_data[["position"]])
      density_idx <- chr_data[["position"]] >= range_min &
        chr_data[["position"]] <= range_max
      ## Note that if I want to completely correct in how I do this
      ## it should not be blindly bin_width, but should be the smaller of
      ## bin_width or the remainder of the chromosome.  Leaving it at bin_width
      ## will under-represent the end of each chromosome/contig, but I don't think
      ## I mind that at all.
      if (isTRUE(divide)) {
        density_vector[i] <- sum(density_idx) / bin_width
      } else {
        density_vector[i] <- sum(density_idx)
      }
      names(density_vector)[i] <- as.character(density_name)
      ## Add the observed variants to the variant_vector
      ## variant_set <- chr_data[density_idx, "mutation"]
      variant_relative_positions <- chr_data[density_idx, "position"] - range_min
      variant_set <- paste0(chr_data[density_idx, "ref"],
                            variant_relative_positions,
                            chr_data[density_idx, "alt"])
      variant_vector[i] <- toString(variant_set)
      names(variant_vector)[i] <- as.character(density_name)

      density_name <- density_name + bin_width
    }
    density_lst[[chr]] <- density_vector
    variant_lst[[chr]] <- variant_vector
  }
  new_options <- options(current_options)

  long_density_vector <- vector()
  long_variant_vector <- vector()
  for (ch in 1:length(density_lst)) {
    chr <- names(density_lst)[ch]
    density_vec <- density_lst[[chr]]
    variant_vec <- variant_lst[[chr]]
    vecnames <- paste0(chr, "_start_", names(density_vec))
    names(density_vec) <- vecnames
    names(variant_vec) <- vecnames
    long_density_vector <- c(long_density_vector, density_vec)
    long_variant_vector <- c(long_variant_vector, variant_vec)
  }
  most_idx <- order(long_density_vector, decreasing = TRUE)
  long_density_vector <- long_density_vector[most_idx]
  long_variant_vector <- long_variant_vector[most_idx]

  mesg("Extracting primer regions.")
  sequence_df <- sm(choose_sequence_regions(long_variant_vector, topn = topn,
                                            max_primer_length = max_primer_length,
                                            genome = genome))
  mesg("Searching for overlapping/closest genes.")
  sequence_df <- xref_regions(sequence_df, gff, bin_width = bin_width,
                              feature_type = feature_type, feature_start = feature_start,
                              feature_end = feature_end, feature_strand = feature_strand,
                              feature_chr = feature_chr, feature_type_column = feature_type_column,
                              feature_id = feature_id, feature_name = feature_name)
  ## Get rid of any regions with 'N' in them

  dropme <- grepl(x = sequence_df[["fivep_superprimer"]], pattern = "N")
  mesg("Dropped ", sum(dropme), " regions with Ns in the 5' region.")
  sequence_df <- sequence_df[!dropme, ]
  dropme <- grepl(x = sequence_df[["threep_superprimer"]], pattern = "N")
  mesg("Dropped ", sum(dropme), " regions with Ns in the 3' region.")
  sequence_df <- sequence_df[!dropme, ]

  keepers <- c("variants", "threep_variants", "fivep_superprimer", "threep_superprimer",
               "fivep_primer", "threep_primer", "overlap_gene_id", "overlap_gene_description",
               "closest_gene_before_id", "closest_gene_before_description")
  if (isTRUE(truncate)) {
    sequence_df <- sequence_df[, keepers]
  }

  retlist <- list(
      "density_vector" = long_density_vector,
      "variant_vector" = long_variant_vector,
      "favorites" = sequence_df)
  return(retlist)

  ## TODO:
  ## 1. Reindex based on primer start
  ## 2. Report closest CDS regions and be able to filter on multicopies
  ## 3. Report if in CDS or not
  ## 4. Extra credit: report polyN runs in the putative amplicon

}

#' Given a named vector of fun regions, make a dataframe which
#' includes putative primers and the spec strings for expected
#' variants.
#'
#' This function came out of our TMRC2 work and seeks to provide an
#' initial set of potential PCR primers which are able to distinguish
#' between different aspects of the data.  In the actual data, we were
#' looking for differences between the zymodemes 2.2 and 2.3.
#'
#' @param vector variant-based set of putative regions with variants
#'  between conditions of interest.
#' @param max_primer_length given this length as a start, whittle down
#'  to a hopefully usable primer size.
#' @param topn Choose this number of variant regions from the rather
#'  larger set of possibilities..
#' @param bin_size Separate the genome into chunks of this size when
#'  hunting for primers, this size will therefore be the approximate
#'  PCR amplicon length.
#' @param genome (BS)Genome to search.
#' @param target_temp PCR temperature to attempt to match.
choose_sequence_regions <- function(vector, max_primer_length = 45,
                                    topn = 200, bin_size = 600,
                                    genome = NULL, target_temp = 58) {
  ## Now get the nucleotides of the first 30 nt of each window
  sequence_df <- as.data.frame(head(vector, n = topn))
  colnames(sequence_df) <- "variants"
  ## At this point we have one column which contains variant specs which start
  ## at the beginning of the region of interest.  So let us choose
  ## primers which end 1 nt before that region.
  sequence_df[["chr"]] <- gsub(x = rownames(sequence_df), pattern = "^(.*)_start_.*$",
                               replacement = "\\1")
  sequence_df[["bin_start"]] <- as.numeric(gsub(x = rownames(sequence_df), pattern = "^(.*)_start_(.*)$",
                                                replacement = "\\2"))
  ## I think it would be nice to have spec strings which count from right->left to make
  ## reading off a 3' primer easier.
  sequence_df[["fivep_references"]] <- gsub(x = sequence_df[["variants"]],
                                            pattern = "([[:alpha:]])(\\d+)([[:alpha:]])",
                                            replacement = "\\1")
  sequence_df[["threep_references"]] <- chartr("ATGC", "TACG", sequence_df[["fivep_references"]])

  sequence_df[["fivep_positions"]] <- gsub(x = sequence_df[["variants"]],
                                           pattern = "([[:alpha:]])(\\d+)([[:alpha:]])",
                                           replacement = "\\2")
  ## Getting the relative 3' position cannot be done with a simple regex/tr operation.
  ## but instead is binwidth - position
  ## which therefore requires a little apply() shenanigans
  recount_positions <- function(position_string, bin_width = 600) {
    positions <- strsplit(position_string, ",")[[1]]
    positions <- gsub(x = positions, pattern = " ", replacement = "")
    positions <- as.numeric(positions)
    new_positions <- toString(as.character(bin_width - positions))
    return(new_positions)
  }
  sequence_df[["threep_positions"]] <- mapply(FUN = recount_positions,
                                              sequence_df[["fivep_positions"]])

  sequence_df[["fivep_alternates"]] <- gsub(x = sequence_df[["variants"]],
                                            pattern = "([[:alpha:]])(\\d+)([[:alpha:]])",
                                            replacement = "\\3")
  sequence_df[["threep_alternates"]] <- chartr("ATGC", "TACG", sequence_df[["fivep_alternates"]])
  ## Now lets rewrite the variants as the combinations of threep_reference/position/alternate
  threep_variants <- function(ref_string, position_string, alt_string) {
    refs <- gsub(pattern = " ", replacement = "", x = strsplit(ref_string, ",")[[1]])
    pos <- gsub(pattern = " ", replacement = "", x = strsplit(position_string, ",")[[1]])
    alts <- gsub(pattern = " ", replacement = "", x = strsplit(alt_string, ",")[[1]])
    new <- toString(paste0(refs, pos, alts))
    return(new)
  }
  sequence_df[["threep_variants"]] <- mapply(FUN = threep_variants,
                                             sequence_df[["threep_references"]],
                                             sequence_df[["threep_positions"]],
                                             sequence_df[["threep_alternates"]])

  ## I am making explicit columns for these so I can look manually
  ## and make sure I am not trying to get sequences which run off the chromosome.
  sequence_df[["fivep_superprimer_start"]] <- sequence_df[["bin_start"]] - max_primer_length - 1
  sequence_df[["threep_superprimer_start"]] <- sequence_df[["bin_start"]] + bin_size + max_primer_length + 1
  sequence_df[["fivep_superprimer_end"]] <- sequence_df[["bin_start"]] - 1
  sequence_df[["threep_superprimer_end"]] <- sequence_df[["bin_start"]] + bin_size
  ## The super_primers are the entire region of length == max_primer_length, I will
  ## iteratively remove bases from the 5' until the target is reached.
  sequence_df[["fivep_superprimer"]] <- ""
  sequence_df[["threep_superprimer"]] <- ""
  ## The fivep_primer threep_primer columns will hold the final primers.
  sequence_df[["fivep_primer"]] <- ""
  sequence_df[["threep_primer"]] <- ""

  ## For me, this is a little too confusing to do in an apply.
  ## I want to fill in the fivep/threep_superprimer columns with the longer
  ## sequence/revcomp regions of interest.
  ## Then I want to chip away at the beginning/end of them to get the actual fivep/threep primer.
  for (i in 1:nrow(sequence_df)) {
    chr <- sequence_df[i, "chr"]

    silence = TRUE
    fivep_superprimer <- try(Biostrings::subseq(genome[[chr]],
                                                sequence_df[i, "fivep_superprimer_start"],
                                                sequence_df[i, "fivep_superprimer_end"]),
                             silent = silence)
    if ("try-error" %in% class(fivep_superprimer)) {
      sequence_df[i, "fivep_superprimer"] <- "Ran over chromosome end"
    } else if (sequence_df[i, "fivep_superprimer_start"] < 0) {
      ## Added this because subseq will assume a negative value means a circular chromosome.
      sequence_df[i, "fivep_superprimer"] <- "Ran over chromosome end"
    } else {
      fivep_superprimer <- as.character(fivep_superprimer)
      sequence_df[i, "fivep_superprimer"] <- fivep_superprimer
    }
    threep_superprimer <- try(Biostrings::subseq(genome[[chr]],
                                                 sequence_df[i, "threep_superprimer_end"],
                                                 sequence_df[i, "threep_superprimer_start"]),
                              silent = silence)
    if ("try-error" %in% class(threep_superprimer)) {
      sequence_df[i, "threep_superprimer"] <- "Ran over chromosome end"
    } else {
      threep_superprimer <- toupper(as.character(spgs::reverseComplement(threep_superprimer)))
      sequence_df[i, "threep_superprimer"] <- threep_superprimer
    }

    fivep_primer <- try(find_subseq_target_temp(fivep_superprimer,
                                                target = target_temp,
                                                direction = "forward"),
                        silent = silence)
    if ("try-error" %in% class(fivep_primer)) {
      fivep_primer <- "bad sequence for priming"
    }
    sequence_df[i, "fivep_primer"] <- fivep_primer

    threep_primer <- try(find_subseq_target_temp(threep_superprimer,
                                                 target = target_temp,
                                                 direction = "reverse"),
                         silent = silence)
    if ("try-error" %in% class(threep_primer)) {
      threep_primer <- "bad sequence for priming"
    }
    sequence_df[i, "threep_primer"] <- threep_primer
  } ## End iterating over sequence_df

  ## Drop anything which is definitely useless
  dropme <- grepl(x = sequence_df[["fivep_superprimer"]], pattern = "Ran")
  sequence_df <- sequence_df[!dropme, ]
  dropme <- grepl(x = sequence_df[["threep_superprimer"]], pattern = "Ran")
  sequence_df <- sequence_df[!dropme, ]

  wanted_columns <- c("variants", "threep_variants", "fivep_superprimer", "threep_superprimer",
                      "fivep_primer", "threep_primer", "overlap_gene_description",
                      "overlap_gene_start", "overlap_gene_end", "closest_gene_before_description",
                      "closest_gene_before_start", "closest_gene_before_end",
                      "closest_gene_after_id", "closest_gene_after_description",
                      "closest_gene_after_start", "closest_gene_after_end",
                      "chr", "bin_start", "fivep_references", "fivep_positions",
                      "fivep_alternates", "threep_references", "threep_positions",
                      "threep_alternates")
  ## sequence_df <- sequence_df[, wanted_columns]
  return(sequence_df)
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
  for (r in 1:nrow(sequence_df)) {
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

find_subseq_target_temp <- function(sequence, target=53, direction="forward", verbose=FALSE) {
  cheapo <- cheap_tm(sequence)
  if (isTRUE(verbose)) {
    mesg("Starting with ", sequence, " and ", cheapo)
  }
  if (nchar(sequence) < 1) {
    stop("Never hit the target tm, something is wrong.")
  }
  if (cheapo[["tm"]] > target) {
    if (direction == "forward") {
      ## Then drop the first character
      subseq <- stringr::str_sub(sequence, 2, nchar(sequence))
    } else {
      ## Then drop the last character
      subseq <- stringr::str_sub(sequence, 1, nchar(sequence) - 1)
    }
    result <- find_subseq_target_temp(subseq, target=target, direction=direction)
  } else {
    mesg("Final tm is: ", cheapo[["tm"]], " from sequence: ", sequence, ".")
    return(sequence)
  }
}


cheap_tm <- function(sequence) {
  ## Taken from: https://www.biostars.org/p/58437/
  n <- nchar(sequence)
  ng <- stringr::str_count(toupper(sequence), "G")
  nc <- stringr::str_count(toupper(sequence), "C")
  n_gc <- ng + nc
  n_at <- n - n_gc
  tm <- 0
  if (n_gc < 18) {
    method <- "2*AT + 4*GC"
    tm <- ((2 * n_at) + (4 * n_gc)) - 7
  } else {
    method <- "64.9 + (41*(GC-16.4)/n)"
    tm <- 64.9 + (41 * ((n_gc - 16.4) / n))
  }
  ret <- list(
      "method" = method,
      "tm" = tm)
  return(ret)
}

## EOF
