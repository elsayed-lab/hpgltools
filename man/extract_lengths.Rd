% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ontology_shared.R
\name{extract_lengths}
\alias{extract_lengths}
\title{Take gene/exon lengths from a suitable data source (gff/TxDb/OrganismDbi)}
\usage{
extract_lengths(
  db = NULL,
  gene_list = NULL,
  type = "GenomicFeatures::transcripts",
  id = "TXID",
  possible_types = c("GenomicFeatures::genes", "GenomicFeatures::cds",
    "GenomicFeatures::transcripts"),
  ...
)
}
\arguments{
\item{db}{Object containing data, if it is a string then a filename is
assumed to a gff file.}

\item{gene_list}{Set of genes to query.}

\item{type}{Function name used for extracting data from TxDb objects.}

\item{id}{Column from the resulting data structure to extract gene IDs.}

\item{possible_types}{Character list of types I have previously used.}

\item{...}{More arguments are passed to arglist.}
}
\value{
Dataframe containing 2 columns: ID, length
}
\description{
Primarily goseq, but also other tools on occasion require a set of gene IDs
and lengths. This function is resposible for pulling that data from either a
gff, or TxDb/OrganismDbi.
}
\seealso{
[GenomicFeatures]
}
