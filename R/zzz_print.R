## I do not fully understand how R collates the files in a package and uses that order
## to define what functions/methods get examined in what order.
## It appears that when I had the `[.expt` and print.expt in expt.R, they were not
## getting evaluated on expt objects, but instead the list print/subsets were being used.
## When I moved them here, they started working.  I am therefore
## thinking to move all of my various methods here.

#' Print a summary of a set of abundant genes.
#'
#' @param x Abundant gene list comprising an element 'high' and 'low'
#'  for the most and least abundant genes observed.
#' @param ... Other args to match the generic.
#' @export
print.abundant_genes <- function(x, ...) {
  message("A set of genes deemed abundant.")
  return(invisible(x))
}

#' Print the result of all_gprofiler()
#'
#' @param x List composed of simple_gprofiler() results for every
#'  up/down set provided by extract_significant_genes().
#' @param ... Other args to match the generic.
#' @export
print.all_gprofiler <- function(x, ...) {
  summary_df <- data.frame()
  for (gp in x) {
    summary_df <- rbind(summary_df, gp[["num_hits"]])
  }
  rownames(summary_df) <- names(x)
  colnames(summary_df) <- names(x[[1]][["num_hits"]])
  message("Running gProfiler on every set of significant genes found:")
  print(summary_df)
  return(invisible(x))
}

#' Print function for a pairwise differential expression result.
#'
#' @param x List containing basic, deseq, edger, ebseq, limma, and
#'  noiseq pairwise comparisons as well as some information about the
#'  contrasts and models used.
#' @param ... Other args to match the generic.
#' @export
print.all_pairwise <- function(x, ...) {
  included <- c()
  for (b in c("basic", "deseq", "ebseq", "edger", "limma", "noiseq")) {
    if (!is.null(x[[b]])) {
      included <- c(included, b)
    }
  }
  message("A pairwise differential expression with results from: ",
          toString(included), ".")
  message("This used a surrogate/batch estimate from: ", x[["batch_type"]], ".")
  message("The primary analysis performed ", nrow(x[["comparison"]][["comp"]]), " comparisons.")
  if (is.null(x[["comparison"]][["heat"]])) {
    message("The logFC agreement among the methods follows:")
    print(x[["comparison"]][["comp"]])
  } else {
    print(x[["comparison"]][["heat"]])
  }
  return(invisible(x))
}

#' Print function for a set of annotations downloaded from biomart.
#'
#' @param x List containing the relevant information gathered from ensembl's biomart.
#' @param ... Other args to match the generic.
#' @export
print.annotations_biomart <- function(x, ...) {
  num_genes <- nrow(x[["annotation"]])
  num_annot <- ncol(x[["annotation"]])
  summary_string <- glue("{num_annot} annotation types for {prettyNum(num_genes, big.mark = ',')} \\
genes/transcripts downloaded from {x[['host']]}.")
  message(summary_string)
  return(invisible(x))
}

#' Print the result from calculate_aucc().
#'
#' @param x List containing the AUCC, correlation, and a plot
#'  describing the AUCC.
#' @param ... Other args to match the generic.
#' @export
print.aucc_info <- function(x, ...) {
  summary_string <- glue("These two tables have an aucc value of: \\
{x[['aucc']]} and correlation: ")
  message(summary_string)
  print(x[["cor"]])
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print the results of load_biomart_go().
#'
#' @param x List from load_biomart_go() containing the table of data,
#'  mart used, host used, name of the mart, and attributes.
#' @param ... Other args to match the generic.
#' @export
print.biomart_go <- function(x, ...) {
  summary_string <- glue("The GO annotations from biomart host {x[['host']]} \\
for species {x[['species']]} provided {prettyNum(nrow(x[['go']]), big.mark = ',')} rows.")
  message(summary_string)
  return(invisible(x))
}

#' Print the results of snp_cds_primers().
#'
#' @param x List currently only containing the dataframe of putative primers.
#' @param ... Other args to match the generic.
#' @export
print.cds_variant_primers <- function(x, ...) {
  summary_string <- glue("A set of potential primers explicitly overlapping known CDS.")
  message(summary_string)
  return(invisible(x))
}

#' Print the results of classify_snps().
#'
#' @param x List containing some fun stats of variants observed.
#' @param ... Other args to match the generic.
#' @export
print.classified_mutations <- function(x, ...) {
  message("The set of classified mutations observed, including: ")
  print(x[["mutations_by_sample"]])
  pheatmap::pheatmap(x[["sample_mutations_norm"]])
  return(invisible(x))
}

#' Print the result from self_evaluate_model().
#'
#' @param x List showing AUC/ROC curves of the test performed, summary
#'  thereof, the confusion matrix, and vector of incorrectly called samples.
#' @param ... Other args to match the generic.
#' @export
print.classifier_evaluation <- function(x, ...) {
  message("The summary of the (in)correct calls is: ")
  print(x[["self_summary"]])
  message("The missed samples are: ")
  print(x[["wrong_samples"]])
  message("The confusion matrix is:")
  print(x[["confusion_mtrx"]])
  message("The ROC AUC is: ", x[["auc"]], ".")
  print(x[["roc_plot"]])
  return(invisible(x))
}

#' Print a clusterprofiler over representation search.
#'
#' @param x Monstrous list of the various results, including but not
#'  limited to plots, go-gene mappings, enrichmed, kegg, david, GO
#'  analyses.
#' @param ... Other args to match the generic.
#' @export
print.clusterprofiler_result <- function(x, ...) {
  message("A set of ontologies produced by clusterprofiler.")
  return(invisible(x))
}

#' Print a combined differential expression analysis.
#'
#' @param x List containing the dataframes for each contrast, the
#'  various plots, the set of wanted contrasts, models used, and
#'  summaries of the data.
#' @param ... Other args to match the generic.
#' @export
print.combined_de <- function(x, ...) {
  message("A set of combined differential expression results.")
  summary_table <- x[["de_summary"]][, c("table", "deseq_sigup", "deseq_sigdown",
                                         "edger_sigup", "edger_sigdown",
                                         "limma_sigup", "limma_sigdown")]
  print(summary_table)
  print(upsetr_combined_de(x))
  return(invisible(x))
}

#' Print a single combined DE result.
#'
#' @param x Data table of combined differential expression results.
#' @param ... Other args to match the generic.
#' @export
print.combined_table <- function(x, ...) {
  message("A combined differential expression table.")
  message("Comprising ", nrow(x[["comb"]]), " genes and including ", x[["includes"]], ".")
  return(invisible(x))
}

#' Print the result of plot_corheat().
#'
#' @param x List containing the correlations observed and a recorded heatmap.3().
#' @param ... Other args to match the generic.
#' @export
print.correlation_heatmap <- function(x, ...) {
  min_cor <- min(x[["data"]])
  non_one_idx <- x[["data"]] == 1
  non_one <- as.data.frame(x[["data"]])
  non_one[non_one_idx] <- 0
  max_cor <- max(non_one)
  summary_string <- glue("A heatmap of pairwise sample correlations ranging from: \
{min_cor} to {max_cor}.")
  message(summary_string)
  print(x[["plot"]])
  return(invisible(x))
}

#' Print the result of plot_disheat().
#'
#' @param x List containing the distances observed and a recorded heatmap.3().
#' @param ... Other args to match the generic.
#' @export
print.distance_heatmap <- function(x, ...) {
  max_distance <- max(x[["data"]])
  non_zero_idx <- x[["data"]] == 0
  non_zero <- as.data.frame(x[["data"]])
  non_zero[non_zero_idx] <- Inf
  min_distance <- min(non_zero)
  summary_string <- glue("A heatmap of pairwise sample distances ranging from: \
{min_distance} to {max_distance}.")
  message(summary_string)
  print(x[["plot"]])
  return(invisible(x))
}

#' Print a representation of compare_de_tables().
#' Note I think I want to have that function return slightly different types
#' depending on how the function call was set up.
#'
#' @param x List provided by plot_linear_scatter() containing
#'  correlations, plots, linear model.
#' @param ... Other args to match the generic.
#' @export
print.cross_table_comparison <- function(x, ...) {
  summary_string <- glue("Comparison of the {x[['fcx_column']]} vs {x[['fcy_column']]}
of two DE tables.")
  message(summary_string)
  plot(x[["scatter"]])
  return(invisible(x))
}

#' Print the result of running deseq_lrt().
#'
#' @param x List containing the DESeq2 result, the associated table,
#'  clusters from degPatterns, list of associated genes, and dataframes
#'  of the most significant genes.
#' @param ... Other args to match the generic.
#' @export
print.deseq_lrt <- function(x, ...) {
  summary_string <- glue("The result from deseq_lrt().")
}

#' Print the result from plot_density().
#'
#' @param x List containing the plot, summary of the
#'  conditions/batches/samples, and the melted table of reads/gene.
#' @param ... Other args to match the generic.
#' @export
print.density_plot <- function(x, ...) {
  summary_string <- glue("Density plot describing {nrow(x[['sample_summary']])} samples.")
  message(summary_string)
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a summary of putative PCR primers based on variant density.
#'
#' @param x List from snp_density_primers()
#' @param ... Other args for the generic.
#' @export
print.density_primers <- function(x, ...) {
  summary_string <- glue("A set of putative PCR primers to distinguish strain variants.")
  message(summary_string)
  return(invisible(x))
}

#' Modified print function for an expt.
#'
#' I am trying to understand how R collates functions.
#' @param x List from create_expt containing the expressionSet,
#'  annotation data, batches, conditions, colors, libsizes, etc.
#' @param ... Other args to match the generic.
#' @export
print.expt <- function(x, ...) {
  summary_string <- glue("An expressionSet containing experiment with {nrow(exprs(x))}
{x[['feature_type']]} and {ncol(exprs(x))} samples. There are {ncol(pData(x))} metadata columns and
{ncol(fData(x))} annotation columns; the primary condition is comprised of:
{toString(levels(as.factor(pData(x)[['condition']])))}.
Its current state is: {what_happened(x)}.")
  message(summary_string)
  return(invisible(x))
}

#' Print a goseq over representation search.
#'
#' @param x List containing the various goseq results, plots,
#'  significant subsets, enrichResult coercions, etc.
#' @param ... Other args to match the generic.
#' @export
print.goseq_result <- function(x, ...) {
  hit_string <- glue("Ontologies observed by goseq using {x[['num_genes']]} genes
with significance cutoff {x[['threshold']]}.
There are {nrow(x[['mf_interesting']])} MF hits, {nrow(x[['bp_interesting']])}, \\
BP hits, and {nrow(x[['cc_interesting']])} CC hits.")
  message(hit_string)
  most <- NULL
  most_num <- 0
  tries <- c("mf", "bp", "cc")
  for (t in tries) {
    key <- paste0(t, "p_plot_over")
    if (is.null(x[["pvalue_plots"]][[key]])) {
      next
    }
    this_num <- nrow(x[["pvalue_plots"]][[key]][["data"]])
    if (this_num > most_num) {
      most_num <- this_num
      most <- key
    }
  }
  if (is.null(most)) {
    message("No categories were deemed worth plotting.")
  } else {
    message("Category ", most, " is the most populated with ", most_num, " hits.")
    plot(x[["pvalue_plots"]][[most]])
  }
  return(invisible(x))
}

#' Print a gostats over representation search.
#'
#' @param x List containing the various gostats results, plots,
#'  significant subsets, enrichResult coercions, etc.
#' @param ... Other args to match the generic.
#' @export
print.gostats_result <- function(x, ...) {
  bp_entries <- nrow(x[["tables"]][["bp_over_enriched"]])
  mf_entries <- nrow(x[["tables"]][["mf_over_enriched"]])
  cc_entries <- nrow(x[["tables"]][["cc_over_enriched"]])
  summary_string <- glue("topgo found {bp_entries} BP categories, {mf_entries} MF categories, and \\
{cc_entries} CC categories.")
  message(summary_string)
  enrichplot::dotplot(x[["enrich_results"]][["bp"]])
  return(invisible(x))
}

#' Print a gprofiler over representation search.
#'
#' @param x List from gProfiler2 containing its various plots, tables
#'  of significant categories for GO, reactome, KEGG, miRNA,
#'  transcription factors, CORUM, wiki pathways, etc; along with the
#'  coerced enrichResult versions.
#' @param ... Other args to match the generic.
#' @export
print.gprofiler_result <- function(x, ...) {
  hit_string <- glue("A set of ontologies produced by gprofiler using {x[['num_genes']]}
genes against the {x[['species']]} annotations and significance cutoff {x[['threshold']]}.
There are {nrow(x[['GO']])} GO hits, {nrow(x[['KEGG']])}, \\
KEGG hits, {nrow(x[['REAC']])} reactome hits, \\
{nrow(x[['WP']])} wikipathway hits, {nrow(x[['TF']])} transcription factor hits, \\
{nrow(x[['MIRNA']])} miRNA hits, {nrow(x[['HPA']])} HPA hits, \\
{nrow(x[['HP']])} HP hits, and {nrow(x[['CORUM']])} CORUM hits.")
  message(hit_string)
  most <- NULL
  most_num <- 0
  tries <- c("MF", "BP", "CC", "KEGG", "REAC", "WP", "TF",
             "MIRNA", "HPA", "HP", "CORUM")
  for (t in tries) {
    if (is.null(x[["pvalue_plots"]][[t]])) {
      next
    }
    this_num <- nrow(x[["pvalue_plots"]][[t]][["data"]])
    if (this_num > most_num) {
      most_num <- this_num
      most <- t
    }
  }
  if (is.null(most)) {
    message("No categories were deemed worth plotting.")
  } else {
    message("Category ", most, " is the most populated with ", most_num, " hits.")
    plot(x[["pvalue_plots"]][[most]])
  }
  return(invisible(x))
}

#' Print a message about the results from graph_metrics().
#'
#' @param x List containing a large number of plots and some tables.
#' @param ... Other args to match the generic.
#' @export
print.graphed_metrics <- function(x, ...) {
  summary_string <- glue("A large number of plots produced by graph_metrics(), \
here are the elements:")
  message(summary_string)
  print(names(x))
  return(invisible(x))
}

#' Print a gsva category search.
#'
#' @param x List containing signature annotations, the result from
#'  GSVA, a modified expressionset, the signatures used, and method.
#' @param ... Other args to match the generic.
#' @export
print.gsva_result <- function(x, ...) {
  summary_string <- glue("GSVA result using method: {x[['method']]} against the \\
{x[['signature_category']]} dataset.
Scores range from: {prettyNum(min(exprs(x[['expt']])))} \\
to: {prettyNum(max(exprs(x[['expt']])))}.")
  message(summary_string)
  return(invisible(x))
}

#' Print gsva categories deemed 'significant'.
#'
#' @param x List of scored GSVA results, including some plots,
#'  likelihood tables, subsets of significant categories, etc.
#' @param ... Other args to match the generic.
#' @export
print.gsva_sig <- function(x, ...) {
  summary_string <- glue("The set of GSVA categories deemed significantly higher than the
distribution of all scores.  It comprises {nrow(x[['subset_table']])} gene sets.")
  message(summary_string)
  print(x[["subset_plot"]])
  return(invisible(x))
}

#' Print a legend of an expressionset.
#'
#' @param x List containing the condition factor, colors used, and plot.
#' @param ... Other args to match the generic.
#' @export
print.legend_plot <- function(x, ...) {
  summary_string <- glue("The colors used in the expressionset are: \\
{toString(x[['color_fact']])}.")
  message(summary_string)
  print(x[["plot"]])
  return(invisible(x))
}

#' Print the library sizes from an experiment.
#'
#' @param x List containing a summary of the library sizes, the plot,
#'  and table.
#' @param ... Other args to match the generic.
#' @export
print.libsize_plot <- function(x, ...) {
  min_value <- min(x[["table"]][["sum"]])
  max_value <- max(x[["table"]][["sum"]])
  message("Library sizes of ", nrow(x[["table"]]), " samples, \
ranging from ", prettyNum(min_value, big.mark = ","),
" to ", prettyNum(max_value, big.mark = ","), ".")
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a set of mapped keepers from combine_de_tables()
#'
#' @param x List full of kept information.
#' @param ... Other args to match the generic.
#' @export
print.mapped_keepers <- function(x, ...) {
  message("The set of mappings between the wanted data and available data in a pairwise comparison.")
  print(names(x))
  return(invisible(x))
}

#' Print a metadata sankey plot.
#'
#' @param x List containing the table of connected nodes and a ggplot2
#'  sankey.
#' @param ... Other args to match the generic.
#' @export
print.meta_sankey <- function(x, ...) {
  summary_string <- glue("A sankey plot describing the metadata of {nrow(x[['design']])} samples,
including {length(x[['observed_nodes']])} out of {length(x[['permutations']])} nodes \\
and traversing metadata factors:
{toString(x[['factors']])}.")
  message(summary_string)
  plot(x[["ggplot"]])
  return(invisible(x))
}

#' Print a nonzero plot.
#'
#' @param x List containing the plot and table describing the data.
#' @param ... Other args to match the generic.
#' @export
print.nonzero_plot <- function(x, ...) {
  summary_string <- glue("A non-zero genes plot of {nrow(x[['table']])} samples.
These samples have an average {prettyNum(mean(x[['table']][['cpm']]))} CPM coverage and \\
{as.integer(mean(x[['table']][['nonzero_genes']]))} genes observed, ranging from \\
{as.integer(min(x[['table']][['nonzero_genes']]))} to
{as.integer(max(x[['table']][['nonzero_genes']]))}.")
  message(summary_string)
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print something useful about the result of create_partitions()
#'
#' @param x List containing the n sets of partitioned data test/train.
#' @param ... Other args to match the generic.
#' @export
print.partitioned_data <- function(x, ...) {
  train_sets <- list()
  count <- 0
  for (tr in x[["trainers"]]) {
    count <- count + 1
    name <- names(x[["trainers"]])[[count]]
    train_sets[[name]] <- rownames(x[["trainers"]][[name]])
  }
  upset_input <- UpSetR::fromList(train_sets)
  upset_plot <- UpSetR::upset(upset_input)
  print(upset_plot)
  summary_string <- glue("A series of {x[['times']]} data partitions with a {x[['p']]} proportion of train/test.")
  message(summary_string)
  return(invisible(x))
}

#' Print some information about a pattern counted genome
#'
#' @param x Dataframe containing how many instances of the pattern
#'  were observed in every gene.
#' @param ... Other args to match the generic.
#' @export
print.pattern_counted <- function(x, ...) {
  summary_string <-
    glue("The pattern was observed {sum(x[['number']])} times ober {nrow(x)} genes.")
  message(summary_string)
  return(invisible(x))
}

#' Print the result from one of the various dimension reductions.
#'
#' @param x List comprised of the residuals, variance summary, tables,
#'  the PCA-esque plot, experimental design, etc.
#' @param ... Other args to match the generic.
#' @export
print.pca_result <- function(x, ...) {
  cond_column <- x[["cond_column"]]
  batch_column <- x[["batch_column"]]

  color_levels <- toString(levels(as.factor(x[["design"]][[cond_column]])))
  batch_levels <- toString(levels(as.factor(x[["design"]][[batch_column]])))
  message("The result of performing a ", x[["pc_method"]], " dimension reduction.
The x-axis is PC", x[["x_pc"]], " and the y-axis is PC", x[["y_pc"]], "
Colors are defined by ", color_levels, "
Shapes are defined by ", batch_levels, ".")
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a representation of the pre vs. post filtered data.
#'
#' @param x List containing the information before/after filtering,
#'  the plots, and summary information.
#' @param ... Other args to match the generic.
#' @export
print.prepost_filter <- function(x, ...) {
  na_idx <- is.na(x[["table"]][["sub_low"]])
  x[["table"]][na_idx, "sub_low"] <- 0
  changed_idx <- x[["table"]][["sub_low"]] > 0
  changed <- x[["table"]][changed_idx, "sub_low"]
  min_range <- min(changed)
  max_range <- max(changed)
  summary_string <- glue("A comparison of the counts before and after filtering.
The number of genes with low coverage changes by {min_range}-{max_range} genes.")
  message(summary_string)
  plot(x[["lowgene_plot"]])
  return(invisible(x))
}

#' Print the result from gather_preprocessing_metadata().
#'
#' @param x List composed of the xlsx output file, new columns added
#'  to it, and copies of the metadata before/after modification.
#' @param ... Other args to match the generic.
#' @export
print.preprocessing_metadata <- function(x, ...) {
  cond_column <- x[["cond_column"]]
  batch_column <- x[["batch_column"]]

  color_levels <- toString(levels(as.factor(x[["design"]][[cond_column]])))
  batch_levels <- toString(levels(as.factor(x[["design"]][[batch_column]])))
  message("The result of performing a ", x[["pc_method"]], " dimension reduction.
The x-axis is PC", x[["x_pc"]], " and the y-axis is PC", x[["y_pc"]], "
Colors are defined by ", color_levels, "
Shapes are defined by ", batch_levels, ".")
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print the result from simple_proper().
#'
#' @param x List including the various plots from PROPER, the
#'  associated tables, simulation options, and example text for a
#'  paper/grant.
#' @param ... Other args to match the generic.
#' @export
print.proper_estimate <- function(x, ...) {
  message(x[[1]][["interpolated_text"]])
  plot(x[[1]][["power_plot"]])
  return(invisible(x))
}

#' Print the result of a reordered variance partition analysis.
#'
#' @param x List of a resorted variance partition analysis and its plot.
#' @param ... Other args to match the generic.
#' @export
print.reordered_varpart <- function(x, ...) {
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print some significantly differentially expressed genes.
#'
#' @param x List containing the parameters used, gene subset tables,
#'  plots, xlsx output file, etc.
#' @param ... Other args to match the generic.
#' @export
print.sig_genes <- function(x, ...) {
  message("A set of genes deemed significant according to ", toString(x[["according"]]), ".")
  params <- ""
  if (!is.null(x[["lfc"]])) {
    params <- glue("LFC cutoff: {x[['lfc']]}")
  }
  if (!is.null(x[["p"]])) {
    params <- glue("{params} {x[['p_type']]} P cutoff: {x[['p']]}")
  }
  if (!is.null(x[["z"]])) {
    params <- glue("{params} Z cutoff: {x[['z']]}")
  }
  if (!is.null(x[["n"]])) {
    params <- glue("{params} topN: {x[['n']]}")
  }
  message("The parameters defining significant were:")
  message(params)
  print(x[["summary_df"]])
  if (sum(rowSums(x[["summary_df"]])) > 0) {
    plot(x[["sig_bar_plots"]][["deseq"]])
  }
  return(invisible(x))
}

#' Print the intersection of significant genes from multiple analyses.
#'
#' @param x List containing some venn diagrams, summaries of
#'  intersections, subsets of the intersections, etc.
#' @param ... Other args to match the generic.
#' @export
print.sig_intersect <- function(x, ...) {
  message("A set of genes deemed significant by multiple methods.")
  return(invisible(x))
}

#' Print some information about the result of snp_intersections().
#'
#' @param x List containing a datatable of intersections, summaries by
#'  chromosome and gene.
#' @param ... Other args to match the generic.
#' @export
print.snp_intersections <- function(x, ...) {
  summary_string <- glue("The combinations of variants, \\
chromosomes, and genes which are unique to every factor
and combination of factors in the data.")
  ## TODO: Decide on something useful to provide here.
  message(summary_string)
  return(invisible(x))
}

#' Print the result of get_snp_sets().
#'
#' @param x List containing the cross references of variants by
#'  factor, the set of observed variants, the possible combinations of
#'  the factor, etc.
#' @param ... Other args to match the generic.
#' @export
print.snp_sets <- function(x, ...) {
  summary_string <- glue("A set of variants observed when cross referencing all variants against
the samples associated with each metadata factor: {x[['factor']]}.  {ncol(x[['values']])}
categories and {nrow(x[['values']])} variants were observed with {length(x[['intersections']])}
combinations among them.  {length(x[['chr_data']])} chromosomes/scaffolds were observed with a
density of variants ranging from {min(x[['density']])} to {max(x[['density']])}.")
  message(summary_string)
  return(invisible(x))
}

#' Print the result of snps_vs_genes().
#'
#' @param x List containing granges of variants, variants observed by
#'  chromosome, gene, and summaries of the result.
#' @param ... Other args to match the generic.
#' @export
print.snps_genes <- function(x, ...) {
  gt_zero <- sum(x[["count_by_gene"]] > 0)
  most_num <- max(x[["count_by_gene"]])
  most_idx <- x[["count_by_gene"]] == most_num
  most_name <- names(x[["count_by_gene"]])[most_idx]
  summary_string <- glue("When the variants observed were cross referenced against annotated genes,
{gt_zero} genes were observed with at least 1 variant.
{most_name} had the most variants, with {most_num}.")
  message(summary_string)
  return(invisible(x))
}

#' Print the result of plot_sm()
#'
#' @param x List containing the pairwise distances/correlations, median/mean values,
#'  quartiles, and the standard median plot.
#' @param ... Other args to match the generic.
#' @export
print.standardmedian_plot <- function(x, ...) {
  min_comp <- min(x[["measurement"]])
  max_comp <- max(x[["measurement"]])
  first_quart <- x[["quantile"]][1]
  third_quart <- x[["quantile"]][2]
  summary_string <- glue("When the standard median metric was plotted, the values observed range
from {min_comp} to {max_comp} with quartiles at {first_quart} and {third_quart}.")
  message(summary_string)
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a topgo over representation search.
#'
#' @param x List of the various over/under representation analyses
#'  provided by topGO, the associated plots, and coerced enrichResults.
#' @param ... Other args to match the generic.
#' @export
print.topgo_result <- function(x, ...) {
  bp_entries <- nrow(x[["tables"]][["bp_over_enriched"]])
  mf_entries <- nrow(x[["tables"]][["mf_over_enriched"]])
  cc_entries <- nrow(x[["tables"]][["cc_over_enriched"]])
  summary_string <- glue("topgo found {bp_entries} BP categories, {mf_entries} MF categories, and \\
{cc_entries} CC categories.")
  message(summary_string)
  enrichplot::dotplot(x[["enrich_results"]][["bp"]])
  return(invisible(x))
}

#' Print a result from plot_topn().
#'
#' @param x List with the topn plot and summary table.
#' @param ... Other args to match the generic.
#' @export
print.topn_plot <- function(x, ...) {
  summary_string <- glue("Plot describing the top-n genes from every sample of a dataset.")
  message(summary_string)
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a result from plot_variance_coefficients().
#'
#' @param x List containing the coefficient of variance plot and summary.
#' @param ... Other args to match the generic.
#' @export
print.varcoef_plot <- function(x, ...) {
  summary_string <- glue("Plot describing the observed variance coefficients on a per-gene basis.")
  message(summary_string)
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print variance partition results.
#'
#' @param x List of results from variancePartition including the model
#'  information, percent/partition plots, dataframes of the
#'  fitted/sorted data by variance, etc.
#' @param ... Other args to match the generic.
#' @export
print.varpart <- function(x, ...) {
  summary_string <- glue("The result of using variancePartition with the model:
{x[['model_string']]}")
  message(summary_string)
  plot(x[["partition_plot"]])
  return(invisible(x))
}

#' Print the result from write_expt.
#'
#' @param x List containing all the many plots, the dataframes, etc.
#' @param ... Other args to match the generic.
#' @export
print.written_expt <- function(x, ...) {
  result_string <- glue("The result from write_expt() sent to:
{x[['excel']]}")
  message(result_string)
  return(invisible(x))
}

#' Print the result from write_xlsx.
#'
#' @param x List containing some information about the xlsx file.
#' @param ... Other args for the generic.
#' @export
print.written_xlsx <- function(x, ...) {
  result_string <- glue("write_xlsx() wrote {x[['file']]}.
 The cursor is on sheet {x[['sheet']]}, row: {x[['end_row']]} column: {x[['end_col']]}.")
  message(result_string)
  return(invisible(x))
}

## EOF
