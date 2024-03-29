% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ontology_kegg.R
\name{pct_all_kegg}
\alias{pct_all_kegg}
\title{Extract the percent differentially expressed genes for all KEGG pathways.}
\usage{
pct_all_kegg(
  all_ids,
  sig_ids,
  organism = "dme",
  pathways = "all",
  pathdir = "kegg_pathways",
  verbose = FALSE,
  ...
)
}
\arguments{
\item{all_ids}{Set of all gene IDs in a given analysis.}

\item{sig_ids}{Set of significant gene IDs.}

\item{organism}{KEGG organism identifier.}

\item{pathways}{What pathways to look at?}

\item{pathdir}{Directory into which to copy downloaded pathway files.}

\item{verbose}{Talky talky?}

\item{...}{Options I might pass from other functions are dropped into arglist.}
}
\value{
Dataframe including the filenames, percentages, nodes included, and
 differential nodes.
}
\description{
KEGGgraph provides some interesting functionality for mapping KEGGids and
examining the pieces. This attempts to use that in order to evaluate how many
'significant' genes are in a given pathway.
}
\seealso{
[KEGGgraph] [KEGGREST]
}
