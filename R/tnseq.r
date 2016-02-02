## Time-stamp: <Tue Feb  2 15:49:25 2016 Ashton Trey Belew (abelew@gmail.com)>

#' tnseq_saturation()  Make a plot and some simple numbers about tnseq saturation
#'
#' This function takes as input a tab separated file from essentiality_tas.pl
#' This is a perl script written to read a bam alignment of tnseq reads against a genome
#' and count how many hits were observed on every TA in the given genome.
#' It furthermore has some logic to tell the difference between reads which were
#' observed on the forward vs. reverse strand as well as reads which appear to be on both strands
#' (eg. they start and end with 'TA').
#' @param file  a file created using the perl script 'essentiality_tas.pl'
#'
#' @return A plot and some numbers
#' @export
tnseq_saturation <- function(file) {
    ## Test arguments
    ##    file = "data/essentiality/mh_ess/5448/tas_t0v0.txt"
    ## End test arguments
    table <- read.table(file=file, header=1)
    second_file <- paste(file, "-inter", sep="")
    second_table <- read.table(file=second_file, header=1)
    table <- rbind(table, second_table)
    data_list <- as.numeric(table$reads)
    max_reads <- max(data_list, na.rm=TRUE)
    log2_data_list <- as.numeric(log2(table$reads + 1))
    data_plot <- hpgl_histogram(log2_data_list, bins=500)
    data_plot <- data_plot + ggplot2::scale_x_continuous(limits=c(0,6)) +
        ggplot2::scale_y_continuous(limits=c(0,2))
    print(sprintf("The maximum value is: %s ", max_reads))
    print(summary(data_list))
    raw <- table(unlist(data_list))
    num_zeros <- raw[as.numeric(names(raw)) == 0]
    num_zeros
    num_gt_ones <- raw[as.numeric(names(raw)) >= 1]
    num_gt_one <- sum(num_gt_ones)
    num_gt_one
    num_gt_twos <- raw[as.numeric(names(raw)) >= 2]
    num_gt_two <- sum(num_gt_twos)
    num_gt_two
    num_gt_fours <- raw[as.numeric(names(raw)) >= 4]
    num_gt_four <- sum(num_gt_fours)
    num_gt_four
    num_gt_eights <- raw[as.numeric(names(raw)) >= 8]
    num_gt_eight <- sum(num_gt_eights)
    num_gt_eight
    num_gt_sixteens <- raw[as.numeric(names(raw)) >= 16]
    num_gt_sixteen <- sum(num_gt_sixteens)
    num_gt_sixteen
    num_gt_thirtytwos <- raw[as.numeric(names(raw)) >= 32]
    num_gt_thirtytwo <- sum(num_gt_thirtytwos)
    num_gt_thirtytwo
    saturation_ratio_1 <- num_gt_one / num_zeros
    print(sprintf("Saturation ratio of 1s to 0s is: %s", saturation_ratio_1))
    saturation_ratio_8 <- num_gt_eight / num_zeros
    print(sprintf("Saturation ratio of 8 to 0s is: %s", saturation_ratio_8))
    saturation_ratio_32 <- num_gt_thirtytwo / num_zeros
    print(sprintf("Saturation ratio of 32 to 0s is: %s", saturation_ratio_32))
    print(sprintf("The number of zeros, ones, twos, fours, eights, sixteens, thirtytwos, saturation1, saturation8, saturation32 are: %s %s %s %s %s %s %s %s %s %s",
                  num_zeros, num_gt_one, num_gt_two, num_gt_four, num_gt_eight, num_gt_sixteen,
                  num_gt_thirtytwo, saturation_ratio_1, saturation_ratio_8, saturation_ratio_32))
    return(data_plot)
}

#' plot_essentiality()  Plot the essentiality of a library as per DeJesus et al.
#'
#' @param file  a file created using the perl script 'essentiality_tas.pl'
#'
#' @return A couple of plots
plot_essentiality <- function(file) {
    ess <- read.csv(file=file, comment.char="#", sep="\t", header=FALSE)
    colnames(ess) <- c("gene","orf_hits","orf_tas","max_run","max_run_span","posterior_zbar","call")
    ess <- ess[with(ess, order("posterior_zbar")), ]
    ##ess <- subset(ess, posterior_zbar > -1)
    ess <- ess[ess$posterior_zbar >= -1, ]
    ess <- transform(ess, rank=ave("posterior_zbar", FUN=function(x) order(x, decreasing=FALSE)))
    zbar_plot <- ggplot2::ggplot(data=ess, ggplot2::aes_string(x="rank", y="posterior_zbar")) +
        ggplot2::geom_point(stat="identity", size=2) +
        ggplot2::geom_hline(color="grey", yintercept=0.0371) +
        ggplot2::geom_hline(color="grey", yintercept=0.9902) +
        ggplot2::theme_bw()
    span_df <- ess[,c("max_run","max_run_span")]
    span_plot <- hpgl_linear_scatter(span_df)
    returns <- list("zbar" = zbar_plot, "scatter" = span_plot$scatter)
    return(returns)
}

## EOF
