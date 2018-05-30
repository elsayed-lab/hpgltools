hpgltools
---------

Status: Travis CI [![Build Status](https://travis-ci.org/abelew/hpgltools.svg?branch=master)]
(https://travis-ci.org/abelew/hpgltools),

## Overview

A bunch of R functions to make playing with high throughput data easier,
developed for applications at the UMD
[Host-Pathogen Genomics Laboratory (HPGL)](http://www.najibelsayed.org/research.aspx).

Although the functionality in this packaged was developed primarily to answer
questions related to host-pathogen genomics and transcriptomics, much of the
functionality is general enough that it may be useful for many other kinds of
analyses.

## Installation

This package contains a mix of bioconductor and cran packages.  This is annoying because
install_github and friends will fail as a result.

There are two ways around this problem:

* Using bioconductor, devtools, and remotes

From the fresh R installation:

> source("http://bioconductor.org/biocLite.R")
> biocLite("devtools")
> devtools::install_github("mangothecat/remotes")
> remotes::install_github("abelew/hpgltools", dependencies=TRUE)

* Otherwise, using make and bioconductor

Download the package via 'git pull' or a zip or whatever, go
into the hpgltools/ directory and:

> make install

One important caveat for newer versions of hpgltools: it now uses packrat to
keep a database of the package versions which I used.  If one wishes to attempt
installing my versions of these packages:

> make packrat_install

A few packages appear to have fallen out of bioconductor/CRAN and so sometimes
packrat installation has annoying errors.

If you wish to run some tests and (re)build the documentation:

> make

There are a bunch of other potential targets available in the Makefile which
might be interesting. Simply running make rebuilds the documentation strings,
runs check, build, create vignettes, and runs the tests.  make install does what
it says on the tin.

Instead, one may perform:

> make prereq

or

> make build

to have it regenerate the vignettes and check for (new) problems.

## Exploring

The easiest way to poke at this and see what it can do is:

R> library(hpgltools)
R> browseVignettes("hpgltools")

As of last count, there were a couple examples using the data(fission)
set, pasilla, and a bacterial data set.

## Functionality

Listed by function then filename in R/.

* annotations: Handles loading annotation data from sources like OrganismDb.
    - **biomart**:  Uses the biomaRt interface to query ensembl.
    - **genbank**:  Uses Rentrez to query genbank.
    - **microbesonline**:  Uses some XML/MySQL to query microbesonline.org --
      though I apparently queried them too often in testing from my workstation,
      so I got blacklisted :(
    - **orgdb**:  Uses the DBI interface to query OrganismDbi/TxDb/OrgDb
      instances.
    - **eupathdb**:  Uses a mix of http/xml/text queries to query up
      eupathdb.org.
    - **uniprot**:  Extract data from either uniprot text files or their
      webservices.
    - **gff**: Read gff files and take annotation data from them.
    - **text**: Extract putative annotation data from various text formats,
      primarily the csv files from trinity.
    - **kegg**: Use the KEGG webservices api.
* ExpressionSets:  Attempts to simplify gathering together
  counts/annotations/designs.
    - **tximport**:  Uses tximport when given filenames ending in supprted extensions.
    - **read metadata**:  Able to read metadata from xlsx/csv/etc and extract
      count table filenames from them.
    - **write expressionsets**:  Writes rdata files of the expressionsets along
      with xlsx files with the raw data, normalized data, and plots.
    - **expt**:  A light S3 object wrapper around ExpressionSets to hopefully make them easier to play with.
* Helpers:  Miscellaneous functions which don't fit elsewhere
    - **installation**:  please_install() attempts to install packages from
      various sources, including biocLite(), install.packages(), install_github(), etc for
      given packages.  Also provides a function to install _all_ of bioconductor
      for the truly self-flagellating.
    - **miscellaneous**:  Arbitrary (hopefully) helpful functions including sm()
      to silence most (all?) functions and pp() to fill in default values when
      making plots.
* Visualization and Sample Metrics:  Plotters and distribution visualizers before doing expression analysis
    - **bar**: Some re-used bar plots like library size.
    - **circos**: Make using circos less obnoxious.
    - **distribution**: Plot some distributions/density plots.
    - **dotplot**: Plot some common dot plots.
    - **genplot**: Attempt to make genoplotR less annoying.
    - **gvis**: Small wrappers for gvis.
    - **heatmap**: Various heatmap shortcuts.
    - **hist**: Reused histograms.
    - **misc**: Some misc plots which have no other home.
    - **point**: Creates various reused scatter/point plots.
    - **shared**: Invokes a series of plots on new data.
* Model Tests / Batch Evaluation:  Play with models, visualize the changes they evoke.
    - **pca**: Functions to simplify PCA analyses.
    - **TSNE**: Functions to simplify TSNE analyses.
    - **surrogates**: Different ways of evaluating surrogate variables(often batch) as per Leek et al.
    - **model testing**: Query the rank of various models.
    - **varpartition**: Use variancePartition to evaluate the contribution of experimental factors in data sets.
* Normalization:  Various data normalizers under a single roof
    - **batch**: Batch correction via limma, sva, ruv.
    - **convert**: Perform cpm/rpkm/patternkm normalizations.
    - **filter**: Low-count filtering via cbcb and genefilter.
    - **norm**: Normalizations taken from cbcb, DESeq, and edgeR.
    - **transform**: Performs logn transformations of data.
    - **shared**: Calls functions in the norm_* files.
* Differential Expression:  Functions to simplify invoking the various differential expression tools
    - **basic**: Completely naive differential expression, provides baseline for comparison.
    - **deseq**: Invokes DESeq2 for all pairwise comparisons of a data set.
    - **edger**: Invokes edgeR for all pairwise comparisons of a data set.
    - **limma**: Invokes limma for all pairwise comparisons of a data set.
    - **shared**: Shared functions across all tools.
* Ontology, KEGG, and friends:  Query if given sets of genes are significantly over represented in various annotation schemes.
    - **clusterprofiler**: Invokes the clusterProfiler family of ontology searches.
    - **goseq**: Invokes goseq for ontology searches, helps graph results.
    - **gostats**: Invokes gostats for ontology searches, helps graph results.
    - **gprofiler**: Invokes gProfiler for ontology searches.
    - **kegg**: Functions to poke at KEGG enrichment.
    - **reactome**: Functions to poke at reactome enrichment.
    - **topgo**: Invokes topGO for ontology searches, graphs results.
    - **shared**: Calls functions in the ontology_* files.
* Random bits:
    - **snp**: Some functions to gather snp data into hopefully easier data
               structures.
    - **proteomics**: Parsers for mzXML data, plotters for the data they contain.
    - **kmeans**: A couple functions for kmeans clustering and visualization from Keith and Ginger.
    - **motif**: Some motif analysis wrappers.
    - **nmer**: Simplifying functions when working with nmers.
    - **xlsx**: Writing files for Excel is terrible, this hopefully makes it less so.

# Notes to self

* Just use data.table.
* Single and double quotes are actually interchangeable in R, I have been assuming there was a
  difference with respect to character/vector or interpolation or escaping and just that I hadn't
  hit a case where it mattered.  This is untrue, in fact the only arbiter I have seen suggests that
  double quotes are preferred style because they remind programmers of other languages that they are
  not characters.  This works for me.
* I have been subsetting poorly.  It turns out that R has two subsetting 'modes', one which
  preserves the type of the original data, and one which does not.  For consistency one should
  probably be careful to use the preserving method: (http://adv-r.had.co.nz/Subsetting.html#subsetting-operators)
  Here are the preserving methods for the first and third elements of each type:

    * vec <- as.vector(c(1,2,3,4,5,6))   :: vec[c(1,3)]
    * lst <- as.list(c(1,2,3,4,5,6))     :: lst[c(1,3)]
    * fac <- as.factor(c(1,2,3,4,5,6))   :: fac[c(1,3), drop=FALSE]

* FALSE is default, which is important if you want to get rid of unused levels, at least for my work
  drop=TRUE is helpful

    * arr <- as.array(c(1,2,3,4,5,6))    :: arr[c(1,3), drop=FALSE]
             Hadley's example doesn't work for me x[, 1, drop=FALSE]
    * df <- data.frame(id=c('a','b','c','d'), data=c(1,2,3,4))  :: df[, 'id', drop=FALSE]

* In contrast, the simplifying methods coerce the results into a different data type, which may be
  either awesome or obnoxious depending on context.  An example of awesome would include times when
  you want to subset a list and just get the character representation of what is left after the
  subset.  If you take the above lst[c(1,3)] you get back the names and values, but if you instead do
  lst[[1]] you get just '1' -- trying lst[[c(1,2)]] ends badly, I am not sure why.

    * vec[[1]]
    * lst[[1]]
    * fac[1:2, drop=TRUE]
    * array[c(1,2), drop=TRUE]  TRUE is the default
    * df[['id']] or df[, 'id']  drop=TRUE is default
