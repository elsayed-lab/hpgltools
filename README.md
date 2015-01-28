hpgltools
---------

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

one may perform

> make prereq to install some datasets and helpers

or

> make build

to have it regenerate the vignettes and check for (new) problems.

## Functionality

### autoloads.R

- **require.auto**: Handles autoloading of packages.

### count_tables.R

Handle simple count table operations

- **expt_subset**:   Pulls pieces out of my expt superclass of ExpressionSet
- **filter_counts**: Low-count filtering of count tables
- **my_read_files**: Read in a set of count tables using some filenames and
                     sample IDs

### misc_functions.R

- **hpgl_cor**: A correlation wrapper
- **write_xls**: Writing excel sheets is easier this way
- **annotate_data**: Ramzi's annotation function from his RNAseq analyses
- **hpgl_voom**: A couple hacks on voom so that it handles normalized data, prints pretty graphs
- **sillydist**: A silly distance function which is a slightly modified euclidean distance from both axes/medians of a graph
- **simple_comparison**: An example function showing a simplified voom/limma workflow
- **simple_pairwise**: A wrapper for simple_comparison which fills in default values for most of the parameters
- **Beta.NA**: A quick solver for voom added by Kwame

### normalization.R

Some simple normalization functions.

- **hpgl_rpkm**: Calculate RPKM using a gff annotation and data frame
- **divide_seq**: Calculate RPKseqM -- reads per Kilo-sequence Million reads, useful for TNSeq and normalizing by # of TAs
- **hpgl_norm**: Wrapper for the many possible normalization schemes, standardizes them all into a DGEList

### ontology.R

Functions to help simplify/combine gene ontology analyses.

- **extract_go**: Make a dataframe of ontology terms using GOTERM
- **dirty_go**: Combine some GO analyses
- **hpgl_enrichGO**: A slightly more robust version of the enrichment function in clusterProfiler
- **hpgl_enrich.internal**: A slightly more robust version of the enrichment function in clusterProfiler

