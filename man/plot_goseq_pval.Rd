% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ontology_plots.R
\name{plot_goseq_pval}
\alias{plot_goseq_pval}
\title{Make a pvalue plot from goseq data.}
\usage{
plot_goseq_pval(
  goterms,
  wrapped_width = 30,
  cutoff = 0.1,
  x_column = "score",
  order_by = "score",
  decreasing = FALSE,
  n = 30,
  mincat = 5,
  level = NULL,
  ...
)
}
\arguments{
\item{goterms}{Some data from goseq!}

\item{wrapped_width}{Number of characters before wrapping to help legibility.}

\item{cutoff}{Pvalue cutoff for the plot.}

\item{x_column}{Choose the data column to put on the x-axis of the plot.}

\item{order_by}{Choose the data column for ordering the bars.}

\item{decreasing}{When ordering the bars, go up or down?}

\item{n}{How many groups to include?}

\item{mincat}{Minimum size of the category for inclusion.}

\item{level}{Levels of the ontology tree to use.}

\item{...}{Arguments passed from simple_goseq()}
}
\value{
Plots!
}
\description{
With minor changes, it is possible to push the goseq results into a
clusterProfiler-ish pvalue plot.  This handles those changes and returns the
ggplot results.
}
\seealso{
[ggplot2]
}
