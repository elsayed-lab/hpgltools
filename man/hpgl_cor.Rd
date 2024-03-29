% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_misc.R
\name{hpgl_cor}
\alias{hpgl_cor}
\title{Wrap cor() to include robust correlations.}
\usage{
hpgl_cor(df, method = "pearson", ...)
}
\arguments{
\item{df}{Data frame to test.}

\item{method}{Correlation method to use. Includes pearson, spearman, kendal, robust.}

\item{...}{Other options to pass to stats::cor().}
}
\value{
Some fun correlation statistics.
}
\description{
Take covRob's robust correlation coefficient and add it to the set of
correlations available when one calls cor().  I should reimplement this using
S4.
}
\examples{
\dontrun{
 hpgl_cor(df = df)
 hpgl_cor(df = df, method = "robust")
}
}
\seealso{
[robust]
}
