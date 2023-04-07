#' Make a plot and some simple numbers about tnseq saturation
#'
#' This function takes as input a tab separated file from essentiality_tas.pl
#' This is a perl script written to read a bam alignment of tnseq reads against
#' a genome and count how many hits were observed on every TA in the given
#' genome. It furthermore has some logic to tell the difference between reads
#' which were observed on the forward vs. reverse strand as well as reads which
#' appear to be on both strands (eg. they start and end with 'TA').
#'
#' @param data data to plot
#' @param column which column to use for plotting
#' @param ylimit Define the y axis?
#' @param adjust Prettification parameter from ggplot2.
#' @return A plot and some numbers:
#'  \enumerate{
#'   \item maximum_reads = The maximum number of reads observed in a single position.
#'   \item hits_by_position = The full table of hits / position
#'   \item num_hit_table = A table of how many times every number of hits was observed.
#'   \item eq_0 = How many times were 0 hits observed?
#'   \item gt_1 = How many positions have > 1 hit?
#'   \item gt_2 = How many positions have > 2 hits?
#'   \item gt_4 = How many positions have > 4 hits?
#'   \item gt_8 = How many positions have > 8 hits?
#'   \item gt_16 = How many positions have > 16 hits?
#'   \item gt_32 = How many positions have > 32 hits?
#'   \item ratios = Character vector of the ratios of each number of hits vs. 0 hits.
#'   \item hit_positions = 2 column data frame of positions and the number of observed hits.
#'   \item hits_summary = summary(hit_positions)
#'   \item plot = Histogram of the number of hits observed.
#'  }
#' @seealso [ggplot2]
#' @examples
#'  \dontrun{
#'  input <- "preprocessing/hpgl0837/essentiality/hpgl0837-trimmed_ca_ta-v0M1.wig"
#'  saturation <- tnseq_saturation(file = input)
#' }
#' @export
tnseq_saturation <- function(data, column = "Reads", ylimit = 100, adjust = 2) {
  table <- NULL
  if (class(data) == "character") {
    table <- read.table(file = data, header = 1, comment.char = "")
    if (colnames(table)[1] == "X.Start") {
      colnames(table)[1] <- "Start"
    }
  } else {
    table <- data
  }

  ## wig_comment removal
  variablestep_idx <- grepl(pattern = "^variable", x = table[["Start"]])
  table <- table[-variablestep_idx, ]
  table[[column]] <- as.numeric(table[[column]])
  table[["Start"]] <- as.numeric(table[["Start"]])
  max_reads <- max(table[[column]], na.rm = TRUE)
  table[["l2"]] <- log2(table[[column]] + 1)
  table[["sample"]] <- "sample"
  count <- NULL
  density_plot <- ggplot(data = table, mapping = aes(x = .data[["l2"]])) +
    ggplot2::geom_density(y = ggplot2::after_stat(count),
                          position = "identity", adjust = adjust) +
    ggplot2::scale_y_continuous(limits = c(0, 0.25)) +
    ggplot2::labs(x = "log2(Number of reads observed)", y = "Number of TAs")

  violin_plot <- ggplot(data = table, mapping = aes(x = .data[["sample"]], y = .data[["l2"]])) +
    ggplot2::geom_violin()

  data_list <- as.numeric(table[, column])
  max_reads <- max(data_list, na.rm = TRUE)
  log2_data_list <- as.numeric(log2(data_list + 1))
  data_plot <- plot_histogram(log2_data_list, bins = 300, adjust = adjust)
  data_plot <- data_plot +
    ggplot2::scale_x_continuous(limits = c(0, 6)) +
    ggplot2::scale_y_continuous(limits = c(0, 2))

  raw <- table(unlist(data_list))
  num_zeros <- raw[as.numeric(names(raw)) == 0]
  num_gt_ones <- raw[as.numeric(names(raw)) >= 1]
  num_gt_one <- sum(num_gt_ones)
  num_gt_twos <- raw[as.numeric(names(raw)) >= 2]
  num_gt_two <- sum(num_gt_twos)
  num_gt_fours <- raw[as.numeric(names(raw)) >= 4]
  num_gt_four <- sum(num_gt_fours)
  num_gt_eights <- raw[as.numeric(names(raw)) >= 8]
  num_gt_eight <- sum(num_gt_eights)
  num_gt_sixteens <- raw[as.numeric(names(raw)) >= 16]
  num_gt_sixteen <- sum(num_gt_sixteens)
  num_gt_thirtytwos <- raw[as.numeric(names(raw)) >= 32]
  num_gt_thirtytwo <- sum(num_gt_thirtytwos)

  saturation_ratios <- vector(length = 6)
  saturation_ratios[1] <- sprintf("%.8f", num_gt_one / num_zeros)
  saturation_ratios[2] <- sprintf("%.8f", num_gt_two / num_zeros)
  saturation_ratios[3] <- sprintf("%.8f", num_gt_four / num_zeros)
  saturation_ratios[4] <- sprintf("%.8f", num_gt_eight / num_zeros)
  saturation_ratios[5] <- sprintf("%.8f", num_gt_sixteen / num_zeros)
  saturation_ratios[6] <- sprintf("%.8f", num_gt_thirtytwo / num_zeros)
  names(saturation_ratios) <- c(1, 2, 4, 8, 16, 32)

  hit_positions <- table
  ## Recast
  hit_positions[[1]] <- as.numeric(hit_positions[[1]])
  hit_positions[[2]] <- as.numeric(hit_positions[[2]])
  hit_idx <- hit_positions[[column]] != 0
  hit_positions <- hit_positions[hit_idx, 1]
  hit_averages <- vector(length = length(hit_positions - 1))
  for (position in seq(from = 2, to = length(hit_positions))) {
    last_position <- position - 1
    hit_averages[last_position] <- hit_positions[position] - hit_positions[last_position]
  }
  hits_summary <- summary(hit_averages)

  retlist <- list(
    "maximum_reads" = max_reads,
    "hits_by_position" = table,
    "num_hit_table" = raw,
    "eq_0" = num_zeros,
    "gt_1" = num_gt_one,
    "gt_2" = num_gt_two,
    "gt_4" = num_gt_four,
    "gt_8" = num_gt_eight,
    "gt_16" = num_gt_sixteen,
    "gt_32" = num_gt_thirtytwo,
    "ratios" = saturation_ratios,
    "hit_positions" = hit_positions,
    "hits_summary" = hits_summary,
    "density" = density_plot,
    "plot" = data_plot
  )

  return(retlist)
}

#' Plot the essentiality of a library as per DeJesus et al.
#'
#' This provides a plot of the essentiality metrics 'zbar' with respect to
#' gene.  In my pipeline, I use their stand alone mh_ess and tn_hmm packages.
#' The result files produced are named
#' mh_ess-{sequence_prefix}-{mapping_parameters}_gene_tas_{m_parameter}.csv
#' where {sequence_prefix} is the basename() of the input sequence file,
#' {mapping_parameters} are a string describing the bowtie mapping used,
#' and {m_parameter} is usually one of 1,2,4,8,16,32 and defines the lower limit
#' of read depth to be considered useful by the mh_ess package.  Thus, before
#' using this, one may want to look at the result from tnseq_saturation() to see
#' if there is a most-appropriate m_parameter.  I think I should figure out a
#' heuristic to choose the m, but I am not sure what it would be, perhaps the
#' median of the hits summary?
#'
#' @param file Result from the DeJesus essentiality package.  I think this has
#'   been effectively replaced by their TRANSIT package.
#' @param order_by What column to use when ordering the data?
#' @param keep_esses Keep entries in the data which are 'S' meaning insufficient
#'   evidence.
#' @param min_sig Minimal value below which a gene is deemed non-essential and
#'   above which it is uncertain.
#' @param max_sig Maximum value above which a gene is deemed essential and below
#'   which it is uncertain.
#' @return A couple of plots
#' @seealso [ggplot2]
#' @export
plot_essentiality <- function(file, order_by = "posterior_zbar", keep_esses = FALSE,
                              min_sig = 0.0371, max_sig = 0.9902) {
  ess <- readr::read_tsv(file = file, comment = "#",
                         col_names = c("gene", "orf_hits", "orf_tas", "max_run",
                                       "max_run_span", "posterior_zbar", "call"),
                         col_types = c("gene"="c", "orf_hits"="i", "orf_tas"="i",
                                       "max_run"="i", "max_run_span"="i",
                                       "posterior_zbar"="d", "call"="f"))
  if (!is.null(order_by)) {
    order_idx <- order(ess[[order_by]], decreasing = FALSE)
    ess <- ess[order_idx, ]
  }
  dropped_esses <- 0
  num_essential <- sum(ess[["posterior_zbar"]] >= max_sig)
  num_uncertain <- sum(ess[["posterior_zbar"]] < max_sig &
                         ess[["posterior_zbar"]] > min_sig)
  num_insig <- sum(ess[["posterior_zbar"]] <= min_sig)
  insig_border <- num_insig
  esses_idx <- ess[["posterior_zbar"]] > -1
  dropped_esses <- sum(!esses_idx)
  if (!isTRUE(keep_esses)) {
    ess <- ess[esses_idx, ]
    insig_border <- insig_border - dropped_esses
  }
  sig_border <- insig_border + num_uncertain

  ess[["rank"]] <- 1:nrow(ess)
  zbar_plot <- ggplot(
    data = ess,
    aes(x = .data[["rank"]], y = .data[["posterior_zbar"]], colour = .data[["call"]])) +
    ggplot2::geom_point(stat = "identity", size = 2) +
    ggplot2::geom_hline(color = "grey", yintercept = min_sig) +
    ggplot2::geom_hline(color = "grey", yintercept = max_sig) +
    ggplot2::geom_vline(color = "grey", xintercept = insig_border) +
    ggplot2::geom_vline(color = "grey", xintercept = sig_border) +
    ggplot2::labs(caption = glue("Insufficient evidence: {dropped_esses}
Essential genes: {num_essential}
Uncertain genes: {num_uncertain}
Non-Essential genes: {num_insig}")) +
  ggplot2::theme_bw()

  span_df <- ess[, c("max_run", "max_run_span")]
  span <- plot_linear_scatter(as.data.frame(span_df))
  retlist <- list(
    "zbar" = zbar_plot,
    "span_plot" = span[["scatter"]],
    "span_cor" = span[["correlation"]],
    "span_hist" = span[["both_histogram"]],
    "span_model" = span[["lm_model"]])
  return(retlist)
}

#' A scoring function for the mh_ess TNSeq method.
#'
#' I dunno, I might delete this function, I am not sure if it will ever get use.
#'
#' @param expt Input expressionset with a metadata column with the ess output files.
#' @param ess_column Metadata column containing the mh_ess output files.
#' @return List containing the scores along with the genes which have changed using it.
score_mhess <- function(expt, ess_column = "essm1") {
  expr <- expt[["expressionset"]]
  design <- pData(expt)
  file_lst <- design[[ess_column]]
  counts <- exprs(expt)
  scores <- data.frame(stringsAsFactors = FALSE)
  for (f in seq_along(file_lst)) {
    sample <- colnames(counts)[f]
    file <- file_lst[f]
    ess <- readr::read_tsv(file = file, comment = "#",
                           col_names = c("gene", "orf_hits", "orf_tas", "max_run",
                                         "max_run_span", "posterior_zbar", "call"),
                           col_types = c("gene"="c", "orf_hits"="i", "orf_tas"="i",
                                         "max_run"="i", "max_run_span"="i",
                                         "posterior_zbar"="d", "call"="f"))
    ess_df <- as.data.frame(ess[, c("gene", "call")])
    rownames(ess_df) <- gsub(x = ess_df[["gene"]], pattern = "^cds_",
                             replacement = "")
    dropped_idx <- rownames(ess_df) == "unfound"
    ess_df <- ess_df[!dropped_idx, ]
    if (f == 1) {
      scores <- ess_df
      scores[["gene"]] <- NULL
    } else {
      scores <- cbind(scores, ess_df[["call"]])
    }
    colnames(scores)[f] <- sample
    scores[[f]] <- as.character(scores[[f]])
  }

  ## S gets a score of 0, NE gets 1, U gets 10, E gets 100.
  zero_idx <- scores == "S"
  scores[zero_idx] <- 0
  one_idx <- scores == "NE"
  scores[one_idx] <- 1
  ten_idx <- scores == "U"
  scores[ten_idx] <- 10
  hun_idx <- scores == "E"
  scores[hun_idx] <- 100
  for (c in colnames(scores)) {
    scores[[c]] <- as.numeric(scores[[c]])
  }
  exprs(expt[["expressionset"]]) <- as.matrix(scores)
  cond_scores <- median_by_factor(expt, fun = "mean")[["medians"]]

  mscores <- rowMeans(cond_scores)
  changed_idx <- mscores != cond_scores[[1]]
  changed_genes <- rownames(cond_scores)[changed_idx]
  changed_df <- cond_scores[changed_genes, ]
  zero_idx <- changed_df == 0
  one_idx <- changed_df == 1
  ten_idx <- changed_df == 10
  hun_idx <- changed_df == 100
  changed_df[zero_idx] <- "S"
  changed_df[one_idx] <- "NE"
  changed_df[ten_idx] <- "U"
  changed_df[hun_idx] <- "E"

  retlist <- list(
    "score_df" = scores,
    "changed_genes" = changed_genes,
    "changed_state" = changed_df)
  return(retlist)
}

#' Plot the saturation of multiple libraries simultaneously.
#'
#' @param meta Experimental metadata
#' @param meta_column Metadata column containing the filenames to query.
#' @param ylimit Maximum y axis
#' @param column Data file column to use for density calculation.
#' @param adjust Density adjustment.
#' @param ggstatsplot Include pretty ggstatsplot plot?
#' @return a plot and table of the saturation for all samples.
#' @export
tnseq_multi_saturation <- function(meta, meta_column, ylimit = 100,
                                   column = "Reads", adjust = 1, ggstatsplot = FALSE) {
  table <- NULL
  filenames <- meta[[meta_column]]
  for (f in seq_along(filenames)) {
    file <- filenames[f]
    sample <- rownames(meta)[f]
    column_data <- read.table(file = file, header = 1, comment.char = "#")
    colnames(column_data) <- c("start", sample)
    column_data <- data.table::as.data.table(column_data)
    if (f == 1) {
      table <- column_data
    } else {
      table <- merge(table, column_data, by.x = "start", by.y = "start")
    }
  }

  melted <- reshape2::melt(data = table, id.vars = "start",
                           value.name = "reads", varnames = "sample")
  colnames(melted) <- c("start", "sample", "reads")
  melted[["log2"]] <- log2(melted[["reads"]] + 1)
  count <- NULL
  plt <- ggplot(data = melted, mapping = aes(x = .data[["log2"]], fill = .data[["sample"]])) +
    ggplot2::geom_density(mapping = aes(y = ggplot2::after_stat(count)), position = "identity",
                          adjust = adjust, alpha = 0.3) +
    ggplot2::scale_y_continuous(limits = c(0, ylimit)) +
    ggplot2::labs(x = "log2(Number of reads observed)", y = "Number of TAs")

  violin <- ggplot(data = melted, mapping = aes(x = .data[["sample"]], y = .data[["log2"]])) +
    ggplot2::geom_violin()

  ggstats <- NULL
  if (isTRUE(ggstatsplot)) {
    ggstats <- ggstatsplot::ggbetweenstats(data = melted, x = sample, y = log2,
                                           notch = TRUE, mean.ci = TRUE, k = 3, outlier.tagging = FALSE,
                                           ggtheme = ggthemes::theme_fivethirtyeight(),
                                           messages = TRUE)
  }

  retlist <- list(
    "table" = table,
    "plot" = plt,
    "ggstats" = ggstats)
  return(retlist)
}

## EOF
