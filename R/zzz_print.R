## I do not fully understand how R collates the files in a package and uses that order
## to define what functions/methods get examined in what order.
## It appears that when I had the `[.expt` and print.expt in expt.R, they were not
## getting evaluated on expt objects, but instead the list print/subsets were being used.
## When I moved them here, they started working.  I am therefore
## thinking to move all of my various methods here.

#' Print a summary of a set of abundant genes.
#'
#' @param x Abundant gene list.
#' @export
print.abundant_genes <- function(x) {
  message("A set of genes deemed abundant.")
  return(invisible(x))
}

#' Print the result of all_gprofiler()
#'
#' @param x The set of results produced by all_gprofiler().
#' @export
print.all_gprofiler <- function(x) {
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
#' @param x pairwise comparisons.
#' @export
print.all_pairwise <- function(x) {
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

#' Print a summary of the result from calculate_aucc().
#'
#' @param x list of results from calculate_aucc().
#' @export
print.aucc_info <- function(x) {
  summary_string <- glue("These two tables have an aucc value of: \\
{x[['aucc']]} and correlation: ")
  message(summary_string)
  print(x[["cor"]])
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a summary of the result from self_evaluate_model().
#'
#' @param x list of results from self_evaluate_model().
#' @export
print.classifier_evaluation <- function(x) {
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

#' Print a summary of a clusterprofiler over representation search.
#'
#' @param x list of results from clusterprofiler.
#' @export
print.clusterprofiler_result <- function(x) {
  message("A set of ontologies produced by clusterprofiler.")
  return(invisible(x))
}

#' Print a summary of a combined differential expression analysis.
#'
#' @param x Set of DE results in one big pile.
#' @export
print.combined_de <- function(x) {
  message("A set of combined differential expression results.")
  summary_table <- x[["de_summary"]][, c("table", "deseq_sigup", "deseq_sigdown",
                                         "edger_sigup", "edger_sigdown",
                                         "limma_sigup", "limma_sigdown")]
  print(summary_table)
  print(upsetr_all(x))
  return(invisible(x))
}

#' Print a summary of a single combined DE result.
#'
#' @param x Big table and some summaries.
#' @export
print.combined_table <- function(x) {
  message("A combined differential expression table.")
  message("Comprising ", nrow(x[["comb"]]), " genes and including ", includes, ".")
  return(invisible(x))
}

#' Print the result of plot_corheat().
#'
#' @param x plot_corheat()
#' @export
print.correlation_heatmap <- function(x) {
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
#' @param x plot_disheat()
#' @export
print.distance_heatmap <- function(x) {
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
#' @param x Result from compare_de_tables().
#' @export
print.cross_table_comparison <- function(x) {
  summary_string <- glue("Comparison of the {x[['fcx_column']]} vs {x[['fcy_column']]}
of two DE tables.")
  message(summary_string)
  plot(x[["scatter"]])
  return(invisible(x))
}

#' Print the result of running deseq_lrt().
#'
#' @param x The list provided by deseq_lrt().
#' @export
print.deseq_lrt <- function(x) {
  summary_string <- glue("The result from deseq_lrt().")
}

#' Print a summary of the result from plot_density().
#'
#' @param x The density plot and associated data in question.
#' @export
print.density_plot <- function(x) {
  summary_string <- glue("Density plot describing {nrow(x[['sample_summary']])} samples.")
  message(summary_string)
  plot(x[["plot"]])
  return(invisible(x))
}

#' Modified print function for an expt.
#'
#' I am trying to understand how R collates functions.
#' @param x expt to print
#' @export
print.expt <- function(x) {
  summary_string <- glue("An expressionSet containing experiment with {nrow(exprs(x))}
genes and {ncol(exprs(x))} samples. There are {ncol(pData(x))} metadata columns and
{ncol(fData(x))} annotation columns; the primary condition is comprised of:
{toString(levels(as.factor(pData(x)[['condition']])))}.
Its current state is: {what_happened(x)}.")
  message(summary_string)
  return(invisible(x))
}

#' Print a summary of a goseq over representation search.
#'
#' @param x list of results from goseq.
#' @export
print.goseq_result <- function(x) {
  message("A set of ontologies produced by goseq.")
  return(invisible(x))
}

#' Print a summary of a gostats over representation search.
#'
#' @param x list of results from gostats.
#' @export
print.gostats_result <- function(x) {
  bp_entries <- nrow(x[["tables"]][["bp_over_enriched"]])
  mf_entries <- nrow(x[["tables"]][["mf_over_enriched"]])
  cc_entries <- nrow(x[["tables"]][["cc_over_enriched"]])
  summary_string <- glue("topgo found {bp_entries} BP categories, {mf_entries} MF categories, and \\
{cc_entries} CC categories.")
  message(summary_string)
  enrichplot::dotplot(x[["enrich_results"]][["bp"]])
  return(invisible(x))
}

#' Print a summary of a gprofiler over representation search.
#'
#' @param x list of results from gprofiler.
#' @export
print.gprofiler_result <- function(x) {
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

#' Print a summary of a gsva category search.
#'
#' @param x list of results from gprofiler.
#' @export
print.gsva_result <- function(x) {
  summary_string <- glue("GSVA result using method: {x[['method']]} against the \\
{x[['signature_category']]} dataset.
Scores range from: {prettyNum(min(exprs(x[['expt']])))} \\
to: {prettyNum(max(exprs(x[['expt']])))}.")
  message(summary_string)
  return(invisible(x))
}

#' Print a summary of gsva categories deemed 'significant'.
#'
#' @param x list of results from gprofiler.
#' @export
print.gsva_sig <- function(x) {
  summary_string <- glue("The set of GSVA categories deemed significantly higher than the
distribution of all scores.  It comprises {nrow(x[['subset_table']])} gene sets.")
  message(summary_string)
  print(x[["subset_plot"]])
  return(invisible(x))
}

#' Print a legend of an expressionset.
#'
#' @param x list of results from plot_legend().
#' @export
print.legend_plot <- function(x) {
  summary_string <- glue("The colors used in the expressionset are: \\
{toString(x[['color_fact']])}.")
  message(summary_string)
  print(x[["plot"]])
  return(invisible(x))
}

#' Print a summary of the library sizes from an experiment.
#'
#' @param x list of results from plot_libsize().
#' @export
print.libsize_plot <- function(x) {
  min_value <- min(x[["table"]][["sum"]])
  max_value <- max(x[["table"]][["sum"]])
  message("Library sizes of ", nrow(x[["table"]]), " samples, \
ranging from ", prettyNum(min_value, big.mark = ","),
" to ", prettyNum(max_value, big.mark = ","), ".")
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a summary of a metadata sankey plot.
#'
#' @param x list of results from plot_meta_sankey().
#' @export
print.meta_sankey <- function(x) {
  summary_string <- glue("A sankey plot describing the metadata of {nrow(x[['design']])} samples,
including {length(x[['observed_nodes']])} out of {length(x[['permutations']])} nodes \\
and traversing metadata factors:
{toString(x[['factors']])}.")
  message(summary_string)
  plot(x[["ggplot"]])
  return(invisible(x))
}

#' Print a summary of a nonzero plot.
#'
#' @param x list of results from plot_libsize().
#' @export
print.nonzero_plot <- function(x) {
  summary_string <- glue("A non-zero genes plot of {nrow(x[['table']])} samples.
These samples have an average {prettyNum(mean(x[['table']][['cpm']]))} CPM coverage and \\
{as.integer(mean(x[['table']][['nonzero_genes']]))} genes observed, ranging from \\
{as.integer(min(x[['table']][['nonzero_genes']]))} to
{as.integer(max(x[['table']][['nonzero_genes']]))}.")
  message(summary_string)
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a summary of the result from one of the various dimension reductions.
#'
#' @param x Result from plot_pca()
#' @export
print.pca_result <- function(x) {
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
#' @param x Result from plot_libsize_prepost().
#' @export
print.prepost_filter <- function(x) {
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

#' Print a summary of the result from gather_preprocessing_metadata().
#'
#' @param x Result from gather_preprocessing_metadata()
#' @export
print.preprocessing_metadata <- function(x) {
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

#' Print a summary of the result from simple_proper().
#'
#' @param x list of results from simple_proper().
#' @export
print.proper_estimate <- function(x) {
  message(x[[1]][["interpolated_text"]])
  plot(x[[1]][["power_plot"]])
  return(invisible(x))
}

#' Print a summary of some significantly differentially expressed genes.
#'
#' @param x Significant gene data structure.
#' @export
print.sig_genes <- function(x) {
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

#' Print a summary of a set of intersection significant genes from multiple analyses.
#'
#' @param x Intersection gene lists.
#' @export
print.sig_intersect <- function(x) {
  message("A set of genes deemed significant by multiple methods.")
  return(invisible(x))
}

#' Print some information about the result of snp_intersections().
#'
#' @param x list of results from snp_intersections().
#' @export
print.snp_intersections <- function(x) {
  summary_string <- glue("The combinations of variants, \\
chromosomes, and genes which are unique to every factor
and combination of factors in the data.")
  ## TODO: Decide on something useful to provide here.
  message(summary_string)
  return(invisible(x))
}

#' Print the result of get_snp_sets().
#'
#' @param x list of results from get_snp_sets().
#' @export
print.snp_sets <- function(x) {
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
#' @param x list of results from snps_vs_genes().
#' @export
print.snps_genes <- function(x) {
  gt_zero <- sum(x[["count_by_gene"]] > 0)
  most_num <- max(x[["count_by_gene"]])
  most_idx <- x[["count_by_gene"]] == most_num
  most_name <- names(x[["count_by_gene"]])[most_idx]
  summary_string <- glue("When the variants observed were cross referenced against annotated genes,
{gt_zero} genes were observed with at least 1 variant.  {most_name} had the most variants, with
{most_num}.")
  message(summary_string)
  return(invisible(x))
}

#' Print a summary of a topgo over representation search.
#'
#' @param x list of results from topgo.
#' @export
print.topgo_result <- function(x) {
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
#' @param x list of results from plot_topn().
#' @export
print.topn_plot <- function(x) {
  summary_string <- glue("Plot describing the top-n genes from every sample of a dataset.")
  message(summary_string)
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a result from plot_variance_coefficients().
#'
#' @param x list of results from plot_variance_coefficients().
#' @export
print.varcoef_plot <- function(x) {
  summary_string <- glue("Plot describing the observed variance coefficients on a per-gene basis.")
  message(summary_string)
  plot(x[["plot"]])
  return(invisible(x))
}

#' Print a summary of variance partition results.
#'
#' @param x list of results from simple_varpart().
#' @export
print.varpart <- function(x) {
  summary_string <- glue("The result of using variancePartition with the model: \\
x[['model_string']]")
  message(summary_string)
  plot(x[["partition_plot"]])
  return(invisible(x))
}

#' Print a summary from write_expt.
#'
#' @param x result of write_expt()
#' @export
print.written_expt <- function(x) {
  result_string <- glue("The result from write_expt() sent to:
{x[['excel']]}")
  return(invisible(x))
}

## EOF
