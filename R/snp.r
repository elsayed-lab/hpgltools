#' Gather snp information for an expt
#'
#' I have some initial code for working with snps, but it seems that it will be getting more use, so
#' make it testable etc.
#'
#' @param expt an expressionset from which to extract information.
#' @param input_dir  Directory to scan for snps output files.
#' @param file_suffix  What to add on the end of the files for the resulting output.
#' @param bam_suffix  How do we find the bam files?
#' @param tolower  Lowercase stuff like 'HPGL'?
#' @return  A new expt object
#' @export
count_expt_snps <- function(expt, type="counts", input_dir="preprocessing/outputs", tolower=TRUE) {
  samples <- rownames(pData(expt))
  if (isTRUE(tolower)) {
    samples <- tolower(samples)
  }
  if (type == "counts") {
    file_suffix <- "_parsed_count.txt"
  } else {
    file_suffix <- "_parsed_ratio.txt"
  }
  ## Create a data table of snp columns
  snp_dt <- read_snp_columns(samples, input_dir, file_suffix)

  snp_exprs <- as.data.frame(snp_dt)
  rownames(snp_exprs) <- snp_exprs[["rownames"]]
  snp_exprs <- as.matrix(snp_exprs[, -1])

  snp_features <- data.table::as.data.table(snp_dt)
  snp_features[, c("species","chromosome","position","original","new") :=
                   data.table::tstrsplit(snp_features[["rownames"]], "_", fixed=TRUE)]
  snp_features <- snp_features[, c("species", "chromosome", "position", "original", "new")]
  snp_features <- as.data.frame(snp_features)
  rownames(snp_features) <- snp_dt[["rownames"]]

  snp_metadata <- pData(expt)
  snp_metadata <- new("AnnotatedDataFrame", snp_metadata)
  Biobase::sampleNames(snp_metadata) <- colnames(snp_exprs)

  snp_features <- new("AnnotatedDataFrame", snp_features)
  Biobase::featureNames(snp_features) <- rownames(snp_exprs)

  expressionset <- new("ExpressionSet", exprs=snp_exprs,
                       phenoData=snp_metadata, featureData=snp_features)

  new_expt <- expt
  new_expt[["expressionset"]] <- expressionset
  new_expt[["original_expressionset"]] <- expressionset
  return(new_expt)
}

#' Extract the observed snps unique to individual categories in a snp set.
#'
#' The result of get_snp_sets provides sets of snps for all possible
#' categories.  This is cool and all, but most of the time we just want the
#' results of a single group in that rather large set (2^number of categories)
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
#' @param snp_expt  The result of count_expt_snps()
#' @param factor  Experimental factor to use for cutting and splicing the data.
#' @param limit  Minimum median number of hits / factor to define a position as
#'   a hit.
#' @param do_save  Save the result?
#' @param savefile  Prefix for a savefile if one chooses to save the result.
#' @return A funky list by chromosome containing:  'medians', the median number
#'   of hits / position by sample type; 'possibilities', the;
#'  'intersections', the groupings as detected by Vennerable;
#'  'chr_data', the raw data; 'set_names', a character list of the actual names
#'  of the groupings; 'invert_names', the opposite of set_names which is to say
#'  the names of groups which do _not_ include samples x,y,z; 'density', a list
#'  of snp densities with respect to chromosomes.  Note that this last one is
#'  approximate as I just calculate with the largest chromosome position
#'  number, not the explicit number of nucleotides in the chromosome.
#' @export
get_snp_sets <- function(snp_expt, factor="pathogenstrain", limit=1,
                         do_save=FALSE, savefile="variants") {
  if (is.null(snp_expt[["design"]][[factor]])) {
    stop("The factor does not exist in the expt.")
  }
  if (isTRUE(do_save) &  file.exists(paste0(savefile, "_", factor, ".rda"))) {
    retlist <- new.env()
    loaded <- load(savefile, envir=retlist)
    retlist <- retlist[["retlist"]]
    return(retlist)
  }

  medians <- median_by_factor(snp_expt, fact=factor)
  ## I am going to split this by chromosome, as a run of 10,000 took 2 seconds,
  ## 100,000 took a minute, and 400,000 took an hour.
  chr <- gsub(pattern="^.+_(.+)_.+_.+_.+$", replacement="\\1", x=rownames(medians))
  medians[["chr"]] <- chr
  tt <- sm(requireNamespace("parallel"))
  tt <- sm(requireNamespace("doParallel"))
  tt <- sm(requireNamespace("iterators"))
  tt <- sm(requireNamespace("foreach"))
  tt <- sm(try(attachNamespace("foreach"), silent=TRUE))
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  num_levels <- length(levels(as.factor(chr)))
  bar <- utils::txtProgressBar(max=num_levels, style=3)
  progress <- function(n) {
    setTxtProgressBar(bar, n)
  }
  pb_opts <- list(progress=progress)
  returns <- list()
  res <- list()
  res <- foreach(i=1:num_levels, .packages=c("hpgltools", "doParallel"), .options.snow=pb_opts, .export=c("snp_by_chr")) %dopar% {
    chromosome_name <- levels(as.factor(chr))[i]
    returns[[chromosome_name]] <- snp_by_chr(medians,
                                             chromosome_name)
  }
  close(bar)
  parallel::stopCluster(cl)

  ## Unpack the res data structure (which probably can be simplified)
  data_by_chr <- list()
  possibilities <- c()
  set_names <- list()
  invert_names <- list()
  for (element in 1:length(res)) {
    datum <- res[[element]]
    chromosome <- datum[["chromosome"]]
    data_by_chr[[chromosome]] <- datum
    possibilities <- data_by_chr[[chromosome]][["possibilities"]]
    set_names <- data_by_chr[[chromosome]][["set_names"]]
    invert_names <- data_by_chr[[chromosome]][["invert_names"]]
  }

  ## Calculate approximate snp densities by chromosome
  all_intersections <- list()
  density_by_chr <- list()
  for (chr in names(data_by_chr)) {
    snps <- rownames(data_by_chr[[chr]][["medians"]])
    num_snps <- length(snps)
    last_position <- max(as.numeric(gsub(pattern="^.+_.+_(.+)_.+_.+$", replacement="\\1", x=snps)))
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
    saved <- save(list="retlist", file=savefile)
  }
  return(retlist)
}

#' Read the output from bcfutils into a count-table-esque
#'
#' I put all my bcfutils output files into one directory, so hunt them down and
#' read them into a data table.
#' @param samples  Sample names to read
#' @param input_dir  Directory from which to read them.
#' @param file_suffix  The suffix of my output files.
#' @return A big honking data table.
read_snp_columns <- function(samples,
                             input_dir="preprocessing/outputs",
                             file_suffix="_parsed_ratio.txt") {
  sample <- samples[1]
  filename_prefix <- file.path(input_dir, sample)
  count_filename <- paste0(filename_prefix, file_suffix)
  column <- read.table(count_filename)
  snp_columns <- data.table::as.data.table(column)
  rownames(snp_columns) <- snp_columns[["V1"]]
  colnames(snp_columns) <- c("rownames", sample)

  bar <- utils::txtProgressBar(style=3)
  for (sample_num in 2:length(samples)) {
    pct_done <- sample_num / length(samples)
    setTxtProgressBar(bar, pct_done)
    sample <- samples[sample_num]
    filename_prefix <- file.path(input_dir, samples[sample_num])
    count_filename <- paste0(filename_prefix, file_suffix)
    column <- read.table(count_filename)
    snp_column <- data.table::as.data.table(column)
    rownames(snp_column) <- snp_column[["V1"]]
    colnames(snp_column) <- c("rownames", sample)
    snp_columns <- merge(snp_columns, snp_column, by="rownames", all=TRUE)
  }
  close(bar)
  na_positions <- is.na(snp_columns)
  snp_columns[na_positions] <- 0
  return(snp_columns)
}

#' Use Rsamtools to read alignments and get snp coverage.
#'
#' This is horrifyingly slow.
#'
#' @param expt  Expressionset to analyze
#' @param type  counts or percent?
#' @param input_dir  Directory containing the samtools results.
#' @param tolower lowercase the sample names?
#' @return  It is so slow I no longer know if it works.
samtools_snp_coverage <- function(expt,
                                  type="counts",
                                  input_dir="preprocessing/outputs",
                                  tolower=TRUE) {
  snp_counts <- count_expt_snps(expt, type=type, input_dir=input_dir, tolower=tolower)
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
      warn(paste0("The bam filename does not exist: ", bam_filenames[sample_num]))
    }
  }

  ## Now I have a data table of rownames and percentages
  ## Next I need to cross reference these against the coverage by position, ergo I must split out the rownames
  pileup_info <- Rsamtools::PileupFiles(bam_filenames)
  ## Taken directly from the Rsamtools manual
  queries <- paste0(snp_counts[["species"]], "_",
                    snp_counts[["chromosome"]], ":",
                    as.numeric(snp_counts[["position"]]) - 1, "-",
                    as.numeric(snp_counts[["position"]]) + 1)
  which <- GenomicRanges::GRanges(queries)
  which_list <- split(which, seqnames(which))
  returns <- list()
  res <- list()
  what <- "seq"

  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  tt <- sm(requireNamespace("parallel"))
  tt <- sm(requireNamespace("doParallel"))
  tt <- sm(requireNamespace("iterators"))
  tt <- sm(requireNamespace("foreach"))
  res <- foreach(c=1:length(which_list), .packages=c("hpgltools", "Rsamtools")) %dopar% {
    chr_name <- names(which_list)[c]
    chr_param <- Rsamtools::ApplyPileupsParam(which=which_list[[c]], what=what)
    snp_calc_coverage <- function(x) {
      ## information at each pile-up position
      qme <- function(y) {
        y <- y[c("A", "C", "G", "T"), , drop=FALSE]
        y <- y + 1L
        result <- colSums(y)
        return(result)
      }
      info <- apply(x[["seq"]], 2, qme)
      retlist <- list(seqnames=x[["seqnames"]], pos=x[["pos"]], info=info)
      return(retlist)
    }
    coverage_result <- Rsamtools::applyPileups(
                                    pileup_info,
                                    snp_calc_coverage,
                                    param=chr_param)
    result_list[[c]] <- coverage_result
  }
  stopped <- parallel::stopCluster(cl)

  ## result is a list of n elements where n is the number of rows in snp_dt
  ## Each element of result is in turn a list containing the following slots:
  ##  seqnames (chromosome), pos (position(s)), info (coverage by file)
  ## The piece of information we want to put into snp_dt is therefore:
  ## coverage_list <- result[[snp_dt_row]][[info]][2, ]
  ## coverage_list is in turn a character list named by filename (which begins with the sample ID)
  ## We will therefore extract the hpglID from it and the coverage for every sample
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
    message("There was an error when creating the data table of coverage, hopefully the data is salvageable.")
  }
  new_dt <- as.data.table(new_dt)
  new_dt[["rownames"]] <- snp_dt[["rownames"]]

}

#' The real worker.  This extracts positions for a single chromosome and puts
#' them into a parallelizable data structure.
#'
#' @param medians  A set of medians by position to look through
#' @param chr_name  Chromosome name to search
#' @param limit  Minimum number of median hits/position to count as a snp.
#' @return  A fun list by chromosome!
snp_by_chr <- function(medians, chr_name="01", limit=1) {
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
  set_list <- Vennerable::Venn(Sets=x_lst)
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

snps_vs_intersections <- function(expt, snp_result) {
  features <- fData(expt)
  features[["start"]] <- sm(as.numeric(features[["start"]]))
  na_starts <- is.na(features[["start"]])
  features <- features[!na_starts, ]
  features[["end"]] <- as.numeric(features[["end"]])
  features[["seqnames"]] <- gsub(pattern="^.+_(.+)$", replacement="\\1", x=features[["seqnames"]])
  expt_granges <- GenomicRanges::makeGRangesFromDataFrame(features)

  set_names <- snp_result[["set_names"]]
  summaries <- list()
  inters <- list()
  for (inter in names(snp_result[["intersections"]])) {
    inter_name <- set_names[[inter]]
    inter_df <- data.frame(row.names=snp_result[["intersections"]][[inter]])
    inter_df[["seqnames"]] <- gsub(pattern="^.+_(.+)_.+_.+_.+$",
                                   replacement="\\1",
                                   x=rownames(inter_df))
    inter_df[["start"]] <- as.numeric(gsub(pattern="^.+_.+_(.+)_.+_.+$",
                                           replacement="\\1",
                                           x=rownames(inter_df)))
    inter_df[["end"]] <- inter_df[["start"]] + 1
    inter_df[["strand"]] <- "+"
    inter_granges <- GenomicRanges::makeGRangesFromDataFrame(inter_df)
    inter_by_gene <- IRanges::subsetByOverlaps(inter_granges,
                                               expt_granges,
                                               type="within",
                                               ignore.strand=TRUE)
    inters[[inter_name]] <- inter_by_gene
    summarized_by_gene <- data.table::as.data.table(inter_by_gene)
    summarized_by_gene[, count := .N, by = list(seqnames)]
    summarized_by_gene <- unique(summarized_by_gene[, c("seqnames", "count"), with=FALSE])
    summaries[[inter_name]] <- summarized_by_gene
  }

  retlist <- list(
    "inters" = inters,
    "summaries" = summaries)
  return(retlist)
}

#' Make a summary of the observed snps/gene
#'
#' @param expt  The original expressionset
#' @param snp_result  The result from get_snp_sets()
#' @return a fun list with some information by gene.
#' @export
snps_vs_genes <- function(expt, snp_result) {
  features <- fData(expt)
  features[["start"]] <- sm(as.numeric(features[["start"]]))
  na_starts <- is.na(features[["start"]])
  features <- features[!na_starts, ]
  features[["end"]] <- as.numeric(features[["end"]])

  ## I don't quite want 5'/3' UTRs, I just want the coordinates starting with
  ## (either 1 or) the end of the last gene and ending with the beginning of the
  ## current gene with respect to the beginning of each chromosome.
  ## That is a weirdly difficult problem for creatures with more than 1 chromosome.
  ## inter_features <- features[, c("start", "end", "seqnames")]
  ## inter_features[["chr_start"]] <- paste0(inter_features[["seqnames"]], "_",
  ##                                         inter_features[["start"]])
  ## inter_feature_order <- order(inter_features[["chr_start"]])
  ## inter_features <- inter_features[inter_feature_order, ]

  expt_granges <- GenomicRanges::makeGRangesFromDataFrame(features)

  snp_positions <- snp_result[["medians"]]
  snp_positions[["seqnames"]] <- gsub(pattern="^(.+_.+)_.+_.+_.+$",
                                      replacement="\\1",
                                      x=rownames(snp_positions))
  snp_positions[["start"]] <- as.numeric(gsub(pattern="^.+_.+_(.+)_.+_.+$",
                                              replacement="\\1",
                                              x=rownames(snp_positions)))
  snp_positions[["end"]] <- snp_positions[["start"]] + 1
  snp_positions[["strand"]] <- "+"
  snp_positions <- snp_positions[, c("seqnames", "start", "end", "strand")]
  snp_granges <- GenomicRanges::makeGRangesFromDataFrame(snp_positions)

  snps_by_gene <- IRanges::subsetByOverlaps(snp_granges, expt_granges, type="within", ignore.strand=TRUE)
  summarized_by_gene <- data.table::as.data.table(snps_by_gene)
  summarized_by_gene[, count := .N, by = list(seqnames)]
  summarized_by_gene <- unique(summarized_by_gene[, c("seqnames", "count"), with=FALSE])
  retlist <- list(
    "expt_granges" = expt_granges,
    "snp_granges" = snp_granges,
    "snps_by_gene" = snps_by_gene,
    "summary" = summarized_by_gene)
  return(retlist)
}

## EOF
