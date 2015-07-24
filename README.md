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
- **autoloads_elsayedlab**: Loads packages not in CRAN and specific for the El-Sayed lab.
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

- **all_pairwise**: Performs a full pairwise analysis of an expt using limma, deseq2, and edger by calling limma_pairwise(), deseq2_pairwise(), and edger_pairwise().
- **coefficient_scatter**: Perform arbitrary scatter plots of columns across toptable() tables.
- **compare_tables**: Perform simple comparisons of the toptable-like outputs from limma/deseq/edger.
- **deseq2_pairwise**: Sets up model matrices and all pairwise contrasts, then performs the comparisons in limma.
- **edger_pairwise**: Sets up model matrices and all pairwise contrasts, then performs the comparisons in EdgeR.
- **hpgl_voom**: A minor hack on limma's voom() function to give a prettier plot.
- **limma_pairwise**: Sets up model matrices and all pairwise contrasts, then performs the comparisons in limma.
- **limma_scatter**: Scatterplot arbitrary data across limma toptables.
- **limma_subset**: Grabs the top/bottom n genes from toptable(), or takes the genes outside a given z-score from the median.
- **make_exampledata**: A hack of limma's exampledata to have arbitrarily sized matrices returned.
- **make_pairwise_contrasts**: Make all possible pairwise contrasts in formats suitable for limma/deseq/edger.
- **makeSVD**: Calls fast.svd and returns the v,u,d elements, not particularly interesting.
- **remove_batch_effect**: An alternate batch removal method resulting from discussions with the Corrada-Bravo lab.
- **simple_comparison**: An example function of a changed vs. control comparison, this provides an example implementation of tasks that may be performed with voom/limma/combat/etc.
- **write_deseq2**: Currently incomplete, intended to take the tables from DESeq2 and write them to csv/excel spreadsheets.
- **write_edger**: Currently incomplete, intended to take the tables from EdgeR and write them to csv/excel spreadsheets.
- **write_limma**: Takes a large # of toptables(), adds qvalues, puts them into a list, and optionally writes them as csv/excel spreadsheets.

### misc_functions.R

- **Beta.NA**: A quick solver for voom added by Kwame
- **hpgl_cor**: A correlation wrapper
- **sillydist**: A very simplistic distance function in order to get weights from arbitrary axes.
- **write_xls**: Writing excel sheets is easier this way

### normalization.R

Some simple normalization functions.

- **batch_counts**: Try out different batch/covariant correction/removal methods. Steals ideas/code from cbcbSEQ/sva/ruvg/limma.
- **cbcb_filter_counts**: A stolen copy of cbcbSEQ's filter_counts() function with some extra talky talky.
- **convert_counts**: Perform cpm/rpkm/etc conversions.
- **divide_seq**: Calculate RPKseqM -- reads per Kilo-sequence Million reads, useful for TNSeq and normalizing by # of TAs.
- **genefilter_pofa_counts**: Use genefilter's PofA method to low-count filter.
- **genefilter_cv_counts**: Use genefilter's cv method to low-count filter.
- **genefilter_kofa_counts**: Use genefilter's kOverA method to low-count filter.
- **hpgl_combatMod**: A stolen copy of cbcbSEQ's combatMod with some hacks to make it not die on odd data.
- **hpgl_log2cpm**: Log2CPM stolen from cbcbSEQ.
- **hpgl_norm**: Wrapper for the many possible normalization schemes, standardizes them all into a DGEList
- **hpgl_qshrink**: A stolen/hacked copy of Kwame's qshrink/qstats functions
- **hpgl_qstats**: A stolen/hacked copy of Kwame's qshrink/qstats functions
- **hpgl_rpkm**: Calculate RPKM using a gff annotation and data frame
- **lowfilter_counts**: The low-count filter stolen from Kwame/cbcbSEQ with extra talkytalky.
- **normalize_counts**: Perform sf|quant|qsmooth|rle|upperquartile|tmm|whatever normalization.
- **normalize_expt**: Takes an expt class, performs a normalization, and returns the class with the data normalized and backed up.
- **transform_counts**: Perform logn transformation (or not)

### gvis_plots.R

Fun click-able googleVis html/javascript plots.

- **hpgl_gvis_ma_plot**: A clicky ma plot
- **hpgl_gvis_volcano_plot**: A clicky volcano plot
- **hpgl_gvis_scatter**: A clicky scatter

### kegg.R

Some functions to help deal with KEGG data

- **hpgl_base_pathview**: A hacked copy of the pathview function.  There was some dumb error in this, but I don't remember what it was anymore.
- **hpgl_pathview**: Iterate through every pathway and fill in gene for a given species with pretty colors.
- **gostats_kegg**: Use gostats' hyperGTest on kegg pathways
- **kegg_get_orgn**: Figure out the KEGG species identifier given a species string like 'Leishmania'

### ontology_clusterprofiler.R

Working with clusterProfiler doesn't have to be obnoxious.

- **simple_clusterprofiler**: Perform a clusterprofiler search given a set of de genes and goids.
- **cluster_trees**: Use topgo's showSigOfNodes() to make pretty trees from clusterprofiler data.
- **hpgl_enrichGO**: A hack of clusterprofiler's enrichGO to make it not fail on corner-cases.
- **hpgl_enrich.internal**: A hack of clusterprofiler's enrich.internal to make it not fail on corner-cases.

### ontology_goseq.R

Working with goseq doesn't have to be irksome.

- **goseq_table**: Prettyify GO tables from goseq.
- **simple_goseq**: Perform a goseq search given a set of de genes and goids.
- **goseq_trees**: Use topgo's showSigOfNodes() to make pretty trees from goseq data.
- **goseq_pval_plots**: Make pretty clusterProfiler-ish barplots of goseq data.

### ontology_gostats.R

Working with gostats doesn't have to be tiresome.

- **simple_gostats**: Perform a gostats search given a set of de genes and goids.
- **gostats_trees**: Use topgo's showSigOfNodes() to make pretty trees from gostats data.
- **gostats_pval_plots**: Make pretty clusterProfiler-ish barplots of gostats data.

### ontology_topgo.R

Working with goseq doesn't have to be irksome.

- **topgo_tables**: Prettyify GO tables from topgo.
- **simple_topgo**: Perform a topgo search given a set of de genes and goids.
- **topgo_trees**: Use topgo's showSigOfNodes() to make pretty trees.
- **make_id2gomap**: Make a id2go.map file in the topgo preferred format given a df of goids in the goseq format.
- **hpgl_topdiffgenes**: A copy of the topdiffgenes function, with the hopes of making it more interesting
- **topgo_pval_plots**: Make pretty clusterProfiler-ish barplots of topgo data.
- **getEdgeWeights**: A hack of the DAG graphing functions to make the lines less bold
- **hpgl_GOplot**: A hack of the DAG graphing functions to make the lines less bold


### ontology_shared.R

Functions to help simplify/combine gene ontology analyses, shared among gostats/goseq/clusterprofiler/topgo

- **deparse_go_value**: The GOTERM function returns complex data structures which don't print well.  Clean it up.
- **goterm**: Pull go terms from GO.db
- **gosyn**: Pull go synonyms from GO.db
- **gosec**: Pull go secondary terms from GO.db
- **godef**: Pull go definitions from GO.db
- **goont**: Pull go top-level ontologies from GO.db
- **golev**: Pull approximate go levels from GO.db
- **gotest**: Test whether a given GOid is valid in GO.db (some TriTrypDB ones are not.)
- **pval_plot**: A generic p-value plotter of GO data.
- **get_genelengths**: Grab gene lengths from a gff file, I should put this with the rpkm function.
- **limma_ontology**: Pass a limma toptable(s) straight to goseq/clusterprofiler/topgo/gostats
- **golevel_df**: Use clusterprofiler's getgolevel() to quickly extract ontology levels to make more succinct plots/trees/tables.

### plots.R

Plots!

- **graph_metrics**: Graph all the various metrics of samples in one function call.  Normalize the data and do it again.  removeBatchEffect() and do it again.
- **hpgl_bcv_plot**: EdgeR's plotBCV as a ggplot2 plot.
- **hpgl_boxplot**: A prettified boxplot of the samples in an experiment.
- **hpgl_density_plot**: A prettified density plotter.
- **hpgl_dist_scatter**: Performs a scatter plot and colors dots with weights, may give a googleVis version.
- **hpgl_corheat**: A heatmap.3 correlation heatmap.
- **hpgl_disheat**: A heatmap.3 distance heatmap.
- **hpgl_heatmap**: Does the work for the above-two functions' heatmaps.
- **hpgl_histogram**: ggplot2() histograms with some defaults set to make it less annoying.
- **hpgl_libsize**: A pretty bargraph of the library size of every sample in an experiment.
- **hpgl_linear_scatter**: Performs a scatter plot, estimates the linear relationship between the two samples, and colors dots accordingly, may give a googleVis version.
- **hpgl_ma_plot**: Performs an MA plot, may give a googleVis version.
- **hpgl_multihistogram**: ggplot2() histograms of multiple distributions with some defaults set.
- **hpgl_nonzero**: A minor hack of Ramzi's nonzero_plot() function.  Plots the # of non-zero genes with respect to CPM by library.
- **hpgl_pairwise_ma**: Plot all the ma relationships between samples of an experiment.
- **hpgl_qq_all**: Performs a qq plot of every sample in an experiment vs. the mean of all samples.
- **hpgl_qq_plot**: Makes a qqplot (log and ratio) between two samples.
- **hpgl_sample_heatmap**: A heatmap.3 of all the genes of a set of samples.
- **hpgl_scatter**: Performs a scatter plot of two samples, may give a googleVis version.
- **hpgl_smc**: Plot the standard median correlation of a set of samples.
- **hpgl_smd**: Plot the standard median distance of a set of samples.
- **hpgl_volcano_plot**: Plot a volcano plot of some toptable() data, may give a googleVis version.
- **multiplot**: A grid layout of lots of plots.
- **heatmap.3**: A hack of heatmap.2 to allow for changes to line widths as per Najib's perference.

### pca.R

PCA is weird.

- **hpgl_pca**: A PCA plotter mutated from cbcbSEQ's
- **plot_pcs**: A Principle Component plotter.
- **pca_plot_largebatch**: A PCA for experiments with lots of batches.
- **pca_plot_smallbatch**: A PCA for experiments with > 6 batches. (prettier)
- **factor_rsquared**: Get r^2 values between a given experimental factor and the data's SVD.
- **u_plot**: Plot the rank-order slope change of svd's u field.
- **pca_information**: Play with the various data returned from svd.
- **pca_highscores**: Get the highest scoring genes for each PC.
