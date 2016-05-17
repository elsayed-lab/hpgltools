hpgltools
---------

Status: Travis CI [![Build Status](https://travis-ci.org/abelew/hpgltools.svg?branch=master)](https://travis-ci.org/abelew/hpgltools),

## Overview

A bunch of R functions to make playing with high throughput data easier,
developed for applications at the UMD
[Host-Pathogen Genomics Laboratory (HPGL)](http://www.najibelsayed.org/research.aspx).

Although the functionality in this packaged was developed primarily to answer
questions related to host-pathogen genomics and transcriptomics, much of the
functionality is general enough that it may be useful for many other kinds of
analyses.

## Installation

There are two likely methods of installing this package:

* If one has Hadley's devtools installed, then run the following from R:

> install_github("elsayed-lab/hpgltools")

Upon completion, if one wishes to update/install the various
packages used in this, perform:

> autoload_all(update=TRUE)

* Otherwise:

Download the package via 'git pull' or a zip or whatever, go
into the hpgltools/ directory and perform:

> make
> make install

The make command will rebuild the documentation strings, run check, build, create vignettes, and run
the tests.  make install does what it says on the tin.

Instead, one may perform:

> make prereq to install some datasets and helpers

or

> make build

to have it regenerate the vignettes and check for (new) problems.

## Exploring

The easiest way to poke at this and see what it can do is:

R> library(hpgltools)
R> browseVignettes("hpgltools")

As of last count, there were a couple examples using the data(fission)
set and some made up data.

## Functionality

Listed by filename in R/.

- **annotations**:     Handles loading annotation data from sources like OrganismDb.
- **autoloads**:       Invokes biocLite() and install_github() for dependencies.
- **biomart**:         Sample biomart searches.
- **circos**:          Creates syntatically valid circos configuration stanzas/files.
- **count_tables**:    Collates count tables and puts them into ExpressionSets.
- **de_basic**:        Completely naive differential expression, provides baseline for comparison.
- **de_deseq**:        Invokes DESeq2 for all pairwise comparisons of a data set.
- **de_edger**:        Invokes edgeR for all pairwise comparisons of a data set.
- **de_limma**:        Invokes limma for all pairwise comparisons of a data set.
- **gvis_plots**:      Create googleVis plots with click-able dots!
- **hpgltools**:       Placeholder for this namespace.
- **kegg**:            Uses KEGGREST and pathview to look at KEGG representation.
- **kmeans**:          Simplistic and incomplete kmeans clustering analyses.
- **misc_functions**:  Functions with no other logical home.
- **models**:          Test models before passing them to things like limma/combat/etc.  Incomplete.
- **motif**:           Perform some simple motif analyses.  Incomplete.
- **norm_batch**:      Batch correction via limma, sva, ruv.
- **norm_convert**:    Perform cpm/rpkm/patternkm normalizations.
- **norm_filter**:     Low-count filtering via cbcb and genefilter.
- **norm_norm**:       Normalizations taken from cbcb, DESeq, and edgeR.
- **norm_shared**:     Calls functions in the norm_* files.
- **norm_transform**:  Performs logn transformations of data.
- **ont_cluster**:     Invokes the clusterProfiler family of ontology searches.
- **ont_goseq**:       Invokes goseq for ontology searches, helps graph results.
- **ont_gostats**:     Invokes gostats for ontology searches, helps graph results.
- **ont_gprofiler**:   Invokes gProfiler for ontology searches.
- **ont_topgo**:       Invokes topGO for ontology searches, graphs results.
- **ont_shared**:      Calls functions in the ontology_* files.
- **pca**:             Simple metrics for pca plotting and interpretation.
- **plot_bar**:        Some re-used bar plots like library size.
- **plot_distribution**: Plot some distributions/density plots.
- **plot_dotplot**:    Plot some common dot plots.
- **plot_heatmap**:    Various heatmap shortcuts.
- **plot_hist**:       Reused histograms.
- **plot_misc**:       Some misc plots which have no other home.
- **plot_point**:      Creates various reused scatter/point plots.
- **plot_shared**:     Invokes a series of plots on new data.
- **reactome**:        Use reactome!
- **surrogate_estimators**: Use sva/ruvg as per Jeff Leek's examples.
- **tnseq**:           Some reused functions for tnseq data.
- **tritrypdb**:       Functions for manipulating/loading data from http://tritrypdb.org/

# Notes to self

* Single and double quotes are actually interchangeable in R, I have been assuming there was a
  difference with respect to character/vector or interpolation or escaping and just that I hadn't
  hit a case where it mattered.  This is untrue, in fact the only arbiter I have seen suggests that
  double quotes are preferred style because they remind programmers of other languages that they are
  not characters.  This works for me.
* I have been subsetting poorly.  It turns out that R has two subsetting 'modes', one which
  preserves the type of the original data, and one which does not.  For consistency one should
  probably be careful to use the preserving method: (http://adv-r.had.co.nz/Subsetting.html#subsetting-operators)
  Here are the preserving methods for the first and third elements of each type:
  ** vec <- as.vector(c(1,2,3,4,5,6))   :: vec[c(1,3)]
  ** lst <- as.list(c(1,2,3,4,5,6))     :: lst[c(1,3)]
  ** fac <- as.factor(c(1,2,3,4,5,6))   :: fac[c(1,3), drop=FALSE]
  FALSE is default, which is important if you want to get rid of unused levels, at least for my work
  drop=TRUE is helpful
  ** arr <- as.array(c(1,2,3,4,5,6))    :: arr[c(1,3), drop=FALSE]
  Hadley's example doesn't work for me x[, 1, drop=FALSE]
  ** df <- data.frame(id=c('a','b','c','d'), data=c(1,2,3,4))  :: df[, 'id', drop=FALSE]

In contrast, the simplifying methods coerce the results into a different data type, which may be
either awesome or obnoxious depending on context.  An example of awesome would include times when
you want to subset a list and just get the character representation of what is left after the
subset.  If you take the above lst[c(1,3)] you get back the names and values, but if you instead do
lst[[1]] you get just '1' -- trying lst[[c(1,2)]] ends badly, I am not sure why.
vec[[1]]
lst[[1]]
fac[1:2, drop=TRUE]
array[c(1,2), drop=TRUE]  TRUE is the default
df[['id']] or df[, 'id']  drop=TRUE is default

