#' Plot the average mass and expected intensity of a set of sequences given an enzyme.
#'
#' This uses the cleaver package to generate a plot of expected intensities vs. weight
#' for a list of protein sequences.
#'
#' @param pep_sequences  Set of protein sequences.
#' @param enzyme  One of the allowed enzymes for cleaver.
#' @param start  Limit the set of fragments from this point
#' @param end  to this point.
#' @return  List containing the distribution of weights and the associated plot.
#' @export
plot_cleaved <- function(pep_sequences, enzyme="trypsin", start=600, end=1500) {
  products <- cleaver::cleave(pep_sequences, enzym=enzyme)
  ## old_par <- par(pin=c(6,3))
  pep_sizes <- data.frame()
  plot(NA, xlim=c(start, end), ylim=c(0, 1),
       xlab="mass in Daltons", ylab="relative intensity",
       main=paste0("Digested sequences with: ", enzyme))
  for (pep in 1:length(products)) {
    seq <- names(pep_sequences)[[pep]]
    prod <- products[[pep]]
    for (i in 1:length(prod)) {
      atoms <- try(BRAIN::getAtomsFromSeq(prod[[i]]), silent=TRUE)
      if (class(atoms) != "try-error") {
        d <- BRAIN::useBRAIN(atoms)
        avg_mass <- d[["avgMass"]]
        most_likely <- max(d[["isoDistr"]])
        lines(avg_mass, most_likely, type="h", col=2)
        row <- c(seq, avg_mass, most_likely)
        pep_sizes <- rbind(pep_sizes, row)
      }
    }
  }
  colnames(pep_sizes) <- c("cds", "average_mass", "highest_likelihood")
  plot <- grDevices::recordPlot()
  retlist <- list("plot" = plot, "sizes" = pep_sizes)
  return(retlist)
}

#' Make a histogram of how many peptides are expected at every integer dalton
#' from a given start to end size for a given enzyme digestion.
#'
#' This is very similar to plot_cleaved() above, but tries to be a little bit smarter.
#'
#' @param pep_sequences  Protein sequences as per plot_cleaved().
#' @param enzyme  Compatible enzyme name from cleaver.
#' @param start  Print histogram from here
#' @param end to here.
#' @param color  Make the bars this color.
#' @return  List containing the plot and size distribution.
#' @export
cleavage_histogram <- function(pep_sequences, enzyme="trypsin",
                               start=600, end=1500, color="black") {
  products <- cleaver::cleave(pep_sequences, enzym=enzyme)
  prod_df <- as.data.frame(products)
  prod_df <- dplyr::as.tbl(prod_df[, c("group_name", "value")])
  colnames(prod_df) <- c("group_name", "sequence")

  gather_masses <- function(sequence) {
    atoms <- try(BRAIN::getAtomsFromSeq(sequence), silent=TRUE)
    if (class(atoms) != "try-error") {
      d <- BRAIN::useBRAIN(atoms)
      ret <- round(d[["avgMass"]])
    } else {
      ret <- 0
    }
    return(ret)
  }

  new_df <- prod_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(mass=gather_masses(sequence))

  plot <- ggplot2::ggplot(data=new_df, ggplot2::aes_string(x="mass")) +
    ggplot2::geom_histogram(binwidth=1, colour=color) +
    ggplot2::scale_x_continuous(limits=c(start, end))

  retlist <- list(
    "plot" = plot,
    "masses" = new_df)
  return(retlist)
}

## EOF
