% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_noiseq.R
\name{noiseq_pairwise}
\alias{noiseq_pairwise}
\title{Perform pairwise comparisons using noiseq.}
\usage{
noiseq_pairwise(
  input = NULL,
  conditions = NULL,
  batches = NULL,
  model_cond = TRUE,
  model_batch = TRUE,
  annot_df = NULL,
  k = 0.5,
  norm = "rpkm",
  factor = "condition",
  lc = 1,
  r = 20,
  adj = 1.5,
  a0per = 0.9,
  filter = 1,
  keepers = NULL,
  ...
)
}
\arguments{
\item{input}{Expressionset to compare.}

\item{conditions}{Set of conditions to query}

\item{batches}{known batches in the data, or a surrogate estimator.}

\item{model_cond}{Add condition to the model?}

\item{model_batch}{Add batch to the model, noiseq has its own combat-like method,
so maybe not necessary?}

\item{annot_df}{Extra annotations.}

\item{k}{Taken from the noiseq docs.}

\item{norm}{Normalization method (noiseq oddly defaults to rpkm).}

\item{factor}{Metadata factor over which to iterate.}

\item{lc}{taken from the noiseq docs.}

\item{r}{taken from the noiseq docs.}

\item{adj}{taken from the noiseq docs.}

\item{a0per}{taken from the noiseq docs.}

\item{filter}{Filter the data?}

\item{keepers}{Perform the comparison only over these specific contrasts instead of all.}

\item{...}{Extra arguments.}
}
\value{
List similar to deseq_pairwise/edger_pairwise/etc.
}
\description{
Perform pairwise comparisons using noiseq.
}
\seealso{
DOI:10.1093/nar/gkv711
}
