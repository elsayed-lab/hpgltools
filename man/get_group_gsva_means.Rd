% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsva.r
\name{get_group_gsva_means}
\alias{get_group_gsva_means}
\title{Create dataframe which gets the maximum within group mean gsva score for each gene set}
\usage{
get_group_gsva_means(gsva_scores, groups, keep_single = TRUE, method = "mean")
}
\arguments{
\item{groups}{list of groups for which to calculate the means}

\item{keep_single}{Keep categories with only 1 element.}

\item{method}{mean or median?}

\item{gsva_result}{Result from simple_gsva()}
}
\value{
dataframe containing max_gsva_score, and within group means for gsva scores
}
\description{
Create dataframe which gets the maximum within group mean gsva score for each gene set
}
\seealso{
[simple_gsva()]
}
