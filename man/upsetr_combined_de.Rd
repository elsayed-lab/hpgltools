% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_plots.R
\name{upsetr_combined_de}
\alias{upsetr_combined_de}
\title{Make an upset plot of all up/down genes in a set of contrasts.}
\usage{
upsetr_combined_de(
  combined,
  according_to = "deseq",
  lfc = 1,
  adjp = 0.05,
  text_scale = 2,
  color_by = NULL,
  desired_contrasts = NULL
)
}
\arguments{
\item{combined}{Result from combine_de_tables.}

\item{according_to}{Choose the lfc column to use.}

\item{lfc}{Choose the logFC}

\item{adjp}{and the p-value.}

\item{desired_contrasts}{Use factors from a few contrasts.}
}
\description{
This is intended to give a quick and dirty view of the genes
observed in a series of de comparisons.
}
