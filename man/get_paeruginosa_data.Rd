% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_hpgltools.r
\name{get_paeruginosa_data}
\alias{get_paeruginosa_data}
\title{Pseudomonas areuginosa strain PA14 data files.}
\usage{
get_paeruginosa_data()
}
\description{
The sample sheet from an experiment with another ESKAPE pathogen,
Pesudomonas!  This is emblematic of how I like to organize samples.  The most
relevant columns for creating an expressionset with create_expt() include:
'Sample ID', 'Condition', 'Batch', and 'file'.  This actually provides a
subset of an experiment in which we were looking simultaneously at the
'large' and 'small' RNA populations in two PA strains, one of which is
deficient in an oligonucleotide degradation enzyme 'orn'.  We were also
seeking to find changes from exponential growth to stationary.  The portions
of the experiment included in this sample sheet are only 3 replicates of the
large RNA samples.  The gff and fasta correspond to the genome and
annotations used when mapping and may be recreated with the accompanying
genbank flat file.  Finally, 'counts' contains the count directory from the
Pseudomonas experiment subset, which is an archive file containing the raw
tables created via bowtie2 -> samtools -> htseq-count.
}
