% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_shared.R
\name{compare_significant_contrasts}
\alias{compare_significant_contrasts}
\title{Implement a cleaner version of 'subset_significants' from analyses with Maria
Adelaida.}
\usage{
compare_significant_contrasts(
  sig_tables,
  second_sig_tables = NULL,
  compare_by = "deseq",
  weights = FALSE,
  contrasts = c(1, 2, 3)
)
}
\arguments{
\item{sig_tables}{Set of significance tables to poke at.}

\item{second_sig_tables}{Separate set of significant results, intra vs. inter comparisons.}

\item{compare_by}{Use which program for the comparisons?}

\item{weights}{When printing venn diagrams, weight them?}

\item{contrasts}{List of contrasts to compare.}
}
\value{
List containing the intersections of the contrasts and plots describing them.
}
\description{
This should provide nice venn diagrams and some statistics to compare 2 or 3
contrasts in a differential expression analysis.
}
\seealso{
[Vennerable]
}
