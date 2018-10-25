#' Make a plot and some simple numbers about tnseq saturation
#'
#' This function takes as input a tab separated file from essentiality_tas.pl
#' This is a perl script written to read a bam alignment of tnseq reads against a genome
#' and count how many hits were observed on every TA in the given genome.
#' It furthermore has some logic to tell the difference between reads which were
#' observed on the forward vs. reverse strand as well as reads which appear to be on both strands
#' (eg. they start and end with 'TA').
#'
#' @param data  data to plot
#' @param column which column to use for plotting
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
#' @seealso \pkg{ggplot2}
#' @examples
#'  \dontrun{
#'  input <- "preprocessing/hpgl0837/essentiality/hpgl0837-trimmed_ca_ta-v0M1.wig"
#'  saturation <- tnseq_saturation(file=input)
#' }
#' @export
tnseq_saturation <- function(data, column="Reads") {
  table <- NULL
  if (class(data) == "character") {
    table <- read.table(file=data, header=1, comment.char="")
    if (colnames(table)[1] == "X.Start") {
      colnames(table)[1] <- "Start"
    }
  } else {
    table <- data
  }

  ## wig_comment removal
  variablestep_idx <- grepl(pattern="^variable", x=table[["Start"]])
  table <- table[-variablestep_idx, ]

  data_list <- as.numeric(table[, column])
  max_reads <- max(data_list, na.rm=TRUE)
  log2_data_list <- as.numeric(log2(data_list + 1))
  data_plot <- plot_histogram(log2_data_list, bins=300)
  data_plot <- data_plot + ggplot2::scale_x_continuous(limits=c(0, 6)) +
    ggplot2::scale_y_continuous(limits=c(0, 2))

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

  saturation_ratios <- vector(length=6)
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
  hit_averages <- vector(length=length(hit_positions - 1))
  for (position in 2:length(hit_positions)) {
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
    "plot" = data_plot
  )

  return(retlist)
}

#' Plot the essentiality of a library as per DeJesus et al.
#'
#' This provides a plot of the essentiality metrics 'zbar' with respect to gene.
#'
#' @param file  a file created using the perl script 'essentiality_tas.pl'
#' @return A couple of plots
#' @seealso \pkg{ggplot2}
#' @export
plot_essentiality <- function(file) {
  ess <- read.csv(file=file, comment.char="#", sep="\t", header=FALSE)
  colnames(ess) <- c("gene", "orf_hits", "orf_tas", "max_run",
                     "max_run_span", "posterior_zbar", "call")
  ess <- ess[with(ess, order("posterior_zbar")), ]
  ##ess <- subset(ess, posterior_zbar > -1)
  ess <- ess[ess[["posterior_zbar"]] >= -1, ]
  ess <- transform(ess, rank=ave("posterior_zbar", FUN=function(x) order(x, decreasing=FALSE)))
  zbar_plot <- ggplot2::ggplot(data=ess, ggplot2::aes_string(x="rank", y="posterior_zbar")) +
    ggplot2::geom_point(stat="identity", size=2) +
    ## What the crab apples are these static numbers (I think to define calling boundaries)
    ggplot2::geom_hline(color="grey", yintercept=0.0371) +
    ggplot2::geom_hline(color="grey", yintercept=0.9902) +
    ggplot2::theme_bw()
  span_df <- ess[, c("max_run", "max_run_span")]
  span_plot <- plot_linear_scatter(span_df)
  returns <- list("zbar" = zbar_plot, "scatter" = span_plot[["scatter"]])
  return(returns)
}

## EOF
