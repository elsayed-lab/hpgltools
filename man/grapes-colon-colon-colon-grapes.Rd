% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_hpgltools.R
\name{\%:::\%}
\alias{\%:::\%}
\title{R CMD check is super annoying about :::.}
\usage{
pkg \%:::\% fun
}
\arguments{
\item{pkg}{on the left hand side}

\item{fun}{on the right hand side}
}
\description{
In a fit of pique, I did a google search to see if anyone else has been
annoyed in the same way as was I.  Yihui Xie was, and in his email to r-devel
in 2013 he proposed a game of hide-and-seek; which I am repeating here.
}
\details{
This just implements ::: as an infix operator that will not trip check.
}
