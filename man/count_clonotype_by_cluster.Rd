% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/singlecell.R
\name{count_clonotype_by_cluster}
\alias{count_clonotype_by_cluster}
\title{Add a df of clonotype observations by cell cluster to @misc of a
Seurat datastructure.}
\usage{
count_clonotype_by_cluster(
  scd,
  column = "res0p2_clusters",
  clono_column = "raw_clonotype_id",
  add_sum = TRUE
)
}
\arguments{
\item{scd}{Seurat single cell datastructure.}

\item{column}{Cluster column in the metadata.}

\item{clono_column}{Column containing VDJ annotations.}

\item{add_sum}{Add sums of the clusters to the metadata?}
}
\value{
The scd with some new metadata.
}
\description{
This seeks to count up and provide a couple of metric of how many
B/T cells are in each cluster of a VDJ single cell dataset.
}
