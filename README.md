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

## Exploring

The easiest way to poke at this and see what it can do is:

R> library(hpgltools)
R> browseVignettes("hpgltools")

As of last count, there were a couple examples using the data(fission)
set and some made up data.

## Functionality

### autoloads.R

- **require.auto**: Handles autoloading/updating of packages (stolen and modified from Ramzi Temanni).
- **autoloads_ontology**: Loads a set of useful ontology packages.
- **autoloads_genome**: Loads a set of packages used in manipulating genomic ranges or downloading genomes.
- **autoloads_deseq**: Loads tools explicitly for differential expression analysis.
- **autoloads_graphs**: Loads extensions to R's plotting capabilities.
- **autoloads_helpers**: Loads development tools, reporting tools, and generally useful toys.
- **autoloads_stats**: Loads extensions to R's statistics framework.
- **autoloads_misc**: Miscellaneous packages which I play with on occasion.
- **autoloads_all**: Runs the above and sets some options for knitr and fonts.

### count_tables.R

Handle simple count table operations

- **create_expt**:   Creates an expt superclass using a csv file of sample definitions.
This includes extraneous metadata like sample batches, colors, etc.
- **create_experiment**:  Does the actual work for create_expt().  Looks for count tables to read; sets colors, batches, conditions; gathers metadata; creates the experiment and returns it.
- **expt_subset**:   Extracts portions of an expt class and returns them.
- **hpgl_read_files**: Read in a set of count tables using some filenames and
sample IDs

### differential_expression.R

Generic methods for working with differential expression data

- **write_limma**: Takes a large # of toptables(), adds qvalues, puts them into a list, and optionally writes them as csv/excel spreadsheets.
- **hpgl_voom**: A minor hack on limma's voom() function to give a prettier plot.
- **unbalanced_pairwise**: Sets up model matrices and contrasts as per a discussion with Kwame in the case when there are irreconcilable batches/conditions, then performs every possible pairwise comparison in limma.
- **balanced_pairwise**: Sets up model matrices and contrasts when one has balanced batch/conditions, then performs every possible pairwise comparison in limma.
- **simple_comparison**: An example function of a changed vs. control comparison, this provides an example implementation of tasks that may be performed with voom/limma/combat/etc.
- **make_exampledata**: A hack of limma's exampleData() to allow for arbitrary data set sizes.

### misc_functions.R

- **hpgl_cor**: A correlation wrapper
- **write_xls**: Writing excel sheets is easier this way
- **sillydist**: A very simplistic distance function in order to get weights from arbitrary axes.
- **Beta.NA**: A quick solver for voom added by Kwame

### normalization.R

Some simple normalization functions.

- **hpgl_rpkm**: Calculate RPKM using a gff annotation and data frame
- **hpgl_log2cpm**: Convert count matrix to log2 counts-per-million reads
- **divide_seq**: Calculate RPKseqM -- reads per Kilo-sequence Million reads, useful for TNSeq and normalizing by # of TAs
- **filter_counts**: Kwame's suggested filter_counts() from one of our lab meetings.
- **normalize_expt**: Takes an expt class, performs a normalization, and returns the class with the data normalized and backed up.
- **hpgl_norm**: Wrapper for the many possible normalization schemes, standardizes them all into a DGEList

### ontology.R

Functions to help simplify/combine gene ontology analyses.

- **goterm**: Pull go terms from GO.db
- **gosyn**: Pull go synonyms from GO.db
- **gosyn**: Pull go secondary terms from GO.db
- **godef**: Pull go definitions from GO.db
- **goont**: Pull go top-level ontologies from GO.db
- **golev**: Pull approximate go levels from GO.db
- **gotest**: Test whether a given GOid is valid in GO.db (some TriTrypDB ones are not.)
- **goseq_table**: Add some annotation information to the tables generated by goseq() so that they become human readable
- **simple_goseq**: Make performing a goseq() analysis simpler and hopefully more robust.  Attempts to handle using non-standard organisms.
- **topgo_pval_plots**: Create clusterProfiler() p-value barplots using topGO data.
- **goseq_pval_plots**: Create clusterProfiler() p-value barplots using goseq data.
- **pval_plot**: Makes a ggplot2 pvalue plot from IDs, scores, and p-values as per clusterProfiler().
- **limma_ontology**: Given a list of limma toptables() (provided by write_limma()), perform ontology searches of the top/bottom n or z genes.
- **simple_topgo**: Make performing a topGO() analysis simpler and hopefully more robust.  Attempts to handle using non-standard organisms.
- **topgo_tables**: Given the difficult results from a topGO() search, return a consistent table.
- **topgo_trees**: Simplify the generation and printing of ontology trees from topGO()
- **simple_clusterprofiler**: Make performing a clusterProfiler() analysis simpler and hopefully handle non-standard organisms.  Also attempts to get around some error conditions often thrown by clusterProfiler.
- **make_id2gomap**: Creates topGO() formatted mappings of gene ID to go ID.
- **goseq_trees**: Generates topGO-like trees from goseq() data.
- **cluster_trees**: Generates topGO-like trees from clusterProfiler() data.
- **hpgl_pathview**: Attempts to color every available KEGG pathway given some DE data.  Is moderately non-stupid in handling mappings from known gene-names to KEGG gene-names.
- **hpgl_enrichGO**: A small hack of clusterProfiler to avoid some errors.
- **hpgl_enrich.internal**: Continues the hack from enrichGO(), makes it useful for sparser organisms like Leishmania.
- **hpgl_GOplot**: A hack on the topGO GOplot() function to make it possible to change line widths etc as per Najib's preferences.
- **Gff2GeneTable**: A hack on clusterProfiler()'s function to make geneID to GOID mappings.  The original implementation did not attempt to limit the gff file in any way and therefore would easily run any computer out of memory.

### plots.R

- **graph_metrics**: Graph all the various metrics of samples in one function call.  Normalize the data and do it again.  removeBatchEffect() and do it again.
- **graph_nonzero**: A minor hack of Ramzi's nonzero_plot() function.  Plots the # of non-zero genes with respect to CPM by library.
- **hpgl_libsize**: A pretty bargraph of the library size of every sample in an experiment.
- **hpgl_boxplot**: A prettified boxplot of the samples in an experiment.
- **hpgl_qq_all**: Performs a qq plot of every sample in an experiment vs. the mean of all samples.
- **hpgl_qq_plot**: Makes a qqplot (log and ratio) between two samples.
- **multiplot**: A grid layout of lots of plots.
- **hpgl_corheat**: A heatmap.3 correlation heatmap.
- **hpgl_disheat**: A heatmap.3 distance heatmap.
- **hpgl_heatmap**: Does the work for the above-two functions' heatmaps.
- **hpgl_sample_heatmap**: A heatmap.3 of all the genes of a set of samples.
- **hpgl_smc**: Plot the standard median correlation of a set of samples.
- **hpgl_smd**: Plot the standard median distance of a set of samples.
- **hpgl_pca**: Standardizes PCA plots of samples.
- **hpgl_ma_plot**: Performs an MA plot, may give a googleVis version.
- **hpgl_dist_scatter**: Performs a scatter plot and colors dots with weights, may give a googleVis version.
- **hpgl_scatter**: Performs a scatter plot of two samples, may give a googleVis version.
- **hpgl_linear_scatter**: Performs a scatter plot, estimates the linear relationship between the two samples, and colors dots accordingly, may give a googleVis version.
- **hpgl_histogram**: ggplot2() histograms with some defaults set to make it less annoying.
- **hpgl_multihistogram**: ggplot2() histograms of multiple distributions with some defaults set.
- **hpgl_density_plot**: Plot a density plot of the samples in an experiment.  It makes a rainbow!
- **hpgl_volcano_plot**: Plot a volcano plot of some toptable() data, may give a googleVis version.
- **heatmap.3**: A hack of heatmap.2 to allow for changes to line widths as per Najib's perference.
- **hpgl_plot_bcv**: An implementation of EdgeR's plotBCV() (essentially voom's plot)

