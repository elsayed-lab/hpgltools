% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proteomics.R
\name{impute_expt}
\alias{impute_expt}
\title{Impute missing values using code from DEP reworked for expressionsets.}
\usage{
impute_expt(
  expt,
  filter = TRUE,
  p = 0.5,
  fun = c("bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero", "mixed",
    "nbavg"),
  ...
)
}
\arguments{
\item{expt}{An ExpressionSet (well, expt), I think it is assumed that this should have
been normalized and filtered for features which have no values across 'most' samples.}

\item{filter}{Use normalize_expt() to filter the data?}

\item{p}{When filtering with pofa, use this p parameter.}

\item{fun}{"bpca", "knn", "QRILC", "MLE", "MinDet",
"MinProb", "man", "min", "zero", "mixed" or "nbavg",
Function used for data imputation based on
[MSnbase::impute-methods()]}

\item{...}{Additional arguments for imputation functions.}
}
\value{
An imputed expressionset.
}
\description{
[impute_expt()] imputes missing values in a proteomics dataset.
}
\seealso{
[MSnbase]
}
