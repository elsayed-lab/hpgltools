# hpgltools

A bunch of R functions() to make playing with high throughput data
easier.

# A map of files and functions

## autoloads.R

require.auto handles autoloading of packages.

## count_tables.R

Handle simple count table operations

expt_subset():  Pulls pieces out of my expt superclass of ExpressionSet
filter_counts(): Low-count filtering of count tables
my_read_files(): Read in a set of count tables using some filenames and sample IDs

## misc_functions.R

hpgl_cor(): a correlation wrapper
write_xls(): writing excel sheets is easier this way
annotate_data(): Ramzi's annotation function from his RNAseq analyses
hpgl_voom(): A couple hacks on voom() so that it handles normalized data, prints pretty graphs
sillydist(): A silly distance function which is a slightly modified euclidean distance from both axes/medians of a graph
simple_comparison(): An example function showing a simplified voom/limma workflow
simple_pairwise(): A wrapper for simple_comparison() which fills in default values for most of the parameters
Beta.NA(): A quick solver for voom() added by Kwame

## normalization.R

Some simple normalization functions

hpgl_rpkm(): Calculate RPKM using a gff annotation and data frame
divide_seq(): Calculate RPKseqM -- reads per Kilo-sequence Million reads, useful for TNSeq and normalizing by # of TAs
hpgl_norm(): Wrapper for the many possible normalization schemes, standardizes them all into a DGEList

## ontology.R

Functions to help simplify/combine gene ontology analyses

extract_go(): Make a dataframe of ontology terms using GOTERM()
dirty_go(): Combine some GO analyses
hpgl_enrichGO(): A slightly more robust version of the enrichment function in clusterProfiler
hpgl_enrich.internal(): A slightly more robust version of the enrichment function in clusterProfiler
