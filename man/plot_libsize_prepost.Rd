% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_bar.R
\name{plot_libsize_prepost}
\alias{plot_libsize_prepost}
\title{Visualize genes observed before/after filtering.}
\usage{
plot_libsize_prepost(expt, low_limit = 2, filter = TRUE, ...)
}
\arguments{
\item{expt}{Input expressionset.}

\item{low_limit}{Threshold to define 'low-representation.'}

\item{filter}{Method used to low-count filter the data.}

\item{...}{Extra arbitrary arguments to pass to normalize_expt()}
}
\value{
Bar plot showing the number of genes below the low_limit before and
 after filtering the data.
}
\description{
Thanks to Sandra Correia for this!  This function attempts to represent the
change in the number of genes which are well/poorly represented in the data
before and after performing a low-count filter.
}
\seealso{
[plot_libsize()] [filter_counts()]
}
