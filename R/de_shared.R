## de_shared.r: Merge differential expression methods into single invocations.
## The functionsin this file seek to make it possible to treat all the DE
## methods identically and/or simultaneously.

#' Perform limma, DESeq2, EdgeR pairwise analyses.
#'
#' This takes an expt object, collects the set of all possible pairwise
#' comparisons, sets up experimental models appropriate for the differential
#' expression analyses, and performs them.
#'
#' Tested in test_29de_shared.R
#' This runs limma_pairwise(), deseq_pairwise(), edger_pairwise(),
#' basic_pairwise() each in turn. It collects the results and does some simple
#' comparisons among them.
#'
#' @param input Dataframe/vector or expt class containing count tables,
#'  normalization state, etc.
#' @param conditions Factor of conditions in the experiment.
#' @param batches Factor of batches in the experiment.
#' @param model_cond Include condition in the model?  This is likely always true.
#' @param modify_p Depending on how it is used, sva may require a modification
#'  of the p-values.
#' @param model_batch Include batch in the model?  This may be true/false/"sva"
#'  or other methods supported by all_adjusters().
#' @param filter Added because I am tired of needing to filter the data before
#'  invoking all_pairwise().
#' @param model_intercept Use an intercept model instead of cell means?
#' @param extra_contrasts Optional extra contrasts beyone the pairwise
#'  comparisons.  This can be pretty neat, lets say one has conditions
#'  A,B,C,D,E and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like:
#'   "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A, de_vs_cb = (E-D)-(C-B)".
#' @param alt_model Alternate model to use rather than just condition/batch.
#' @param libsize Library size of the original data to help voom().
#' @param test_pca Perform some tests of the data before/after applying a given
#'  batch effect.
#' @param annot_df Annotations to add to the result tables.
#' @param parallel Use dopar to run limma, deseq, edger, and basic simultaneously.
#' @param do_basic Perform a basic analysis?
#' @param do_deseq Perform DESeq2 pairwise?
#' @param do_ebseq Perform EBSeq (caveat, this is NULL as opposed to TRUE/FALSE
#'  so it can choose).
#' @param do_edger Perform EdgeR?
#' @param do_limma Perform limma?
#' @param convert Modify the data with a 'conversion' method for PCA?
#' @param norm Modify the data with a 'normalization' method for PCA?
#' @param verbose Print extra information while running?
#' @param surrogates Either a number of surrogates or method to estimate it.
#' @param ...  Picks up extra arguments into arglist, currently only passed to
#'  write_limma().
#' @return A list of limma, deseq, edger results.
#' @seealso [limma_pairwise()] [edger_pairwise()] [deseq_pairwise()] [ebseq_pairwise()]
#'  [basic_pairwise()]
#' @examples
#' \dontrun{
#'  lotsodata <- all_pairwise(input = expt, model_batch = "svaseq")
#'  summary(lotsodata)
#'  ## limma, edger, deseq, basic results; plots; and summaries.
#' }
#' @export
all_pairwise <- function(input = NULL, conditions = NULL,
                         batches = NULL, model_cond = TRUE,
                         modify_p = FALSE, model_batch = TRUE, filter = NULL,
                         model_intercept = FALSE, extra_contrasts = NULL,
                         alt_model = NULL, libsize = NULL, test_pca = TRUE,
                         annot_df = NULL, parallel = TRUE,
                         do_basic = TRUE, do_deseq = TRUE, do_ebseq = NULL,
                         do_edger = TRUE, do_limma = TRUE,
                         convert = "cpm", norm = "quant", verbose = TRUE,
                         surrogates = "be", ...) {
  arglist <- list(...)
  if (is.null(model_cond)) {
    model_cond <- TRUE
  }
  if (is.null(model_batch)) {
    model_batch <- FALSE
  }
  if (is.null(model_intercept)) {
    model_intercept <- FALSE
  }

  if (isTRUE(model_cond)) {
    message("This DE analysis will perform all pairwise comparisons among:")
    print(table(pData(input)[["condition"]]))
    if (isTRUE(model_batch)) {
      message("This analysis will include a batch factor in the model comprised of:")
      print(table(pData(input)[["batch"]]))
    } else if ("character" %in% class(model_batch)) {
      message("This analysis will include surrogate estimates from: ", model_batch, ".")
    } else if ("matrix" %in% class(model_batch)) {
      message("This analysis will include a matrix of surrogate estimates.")
    }
    if (!is.null(filter)) {
      message("This will pre-filter the input data using normalize_expt's: ",
              filter, " argument.")
    }
  } else {
    message("This analysis is not using the condition factor from the data.")
  }

  if (!is.null(filter)) {
    input <- sm(normalize_expt(input, filter = filter))
  }
  null_model <- NULL
  sv_model <- NULL
  model_type <- model_batch
  if (class(model_batch)[1] == "character") {
    model_params <- all_adjusters(input, estimate_type = model_batch,
                                  surrogates = surrogates)
    model_batch <- model_params[["model_adjust"]]
    null_model <- model_params[["null_model"]]
    sv_model <- model_batch
  }

  ## Add a little logic to do a before/after batch PCA plot.
  pre_pca <- NULL
  post_pca <- NULL
  if (isTRUE(test_pca)) {
    pre_batch <- sm(normalize_expt(input, filter = TRUE, batch = FALSE,
                                   transform = "log2", convert = convert, norm = norm))
    mesg("Plotting a PCA before surrogate/batch inclusion.")
    pre_pca <- plot_pca(pre_batch, plot_labels = FALSE,
                        ...)
    post_batch <- pre_batch
    if (isTRUE(model_type)) {
      model_type <- "batch in model/limma"
      if (isTRUE(verbose)) {
        message("Using limma's removeBatchEffect to visualize with(out) batch inclusion.")
      }
      post_batch <- sm(normalize_expt(input, filter = TRUE, batch = TRUE, transform = "log2"))
    } else if (class(model_type)[1] == "character") {
      mesg("Using ", model_type, " to visualize before/after batch inclusion.")
      test_norm <- "quant"
      if (model_type != "TRUE" & model_type != FALSE) {
        ## Then it is probably some sort of sva which will have a hard time with quantile.
        test_norm <- "raw"
      }
      mesg("Performing a test normalization with: ", test_norm)
      if (!isFALSE(model_batch)) {
        post_batch <- try(normalize_expt(input, filter = TRUE, batch = model_type,
                                         transform = "log2", convert = "cpm",
                                         norm = test_norm))
      } else {
        post_batch <- NULL
      }
    } else {
      model_type <- "none"
      mesg("Assuming no batch in model for testing pca.")
    }
    post_pca <- plot_pca(post_batch, plot_labels = FALSE,
                         ...)
  }

  ## do_ebseq defaults to NULL, this is so that we can query the number of
  ## conditions and choose accordingly. EBSeq is very slow, so if there are more
  ## than 3 or 4 conditions I do not think I want to wait for it.
  num_conditions <- 0
  if (is.null(conditions)) {
    num_conditions <- length(levels(as.factor(input[["conditions"]])))
  } else {
    num_conditions <- length(levels(as.factor(conditions)))
  }
  if (is.null(do_ebseq)) {
    if (num_conditions > 4) {
      do_ebseq <- FALSE
    } else if (num_conditions < 2) {
      stop("Unable to find the number of conditions in the data.")
    } else {
      do_ebseq <- TRUE
    }
  }

  ## Put a series of empty lists in the final results data structure
  ## so that later I will know to perform each of these analyses without
  ## having to query do_method.
  results <- list()
  if (isTRUE(do_basic)) {
    results[["basic"]] <- list()
  }
  if (isTRUE(do_deseq)) {
    results[["deseq"]] <- list()
  }
  if (isTRUE(do_ebseq)) {
    results[["ebseq"]] <- list()
  }
  if (isTRUE(do_edger)) {
    results[["edger"]] <- list()
  }
  if (isTRUE(do_limma)) {
    results[["limma"]] <- list()
  }

  res <- NULL
  if (isTRUE(parallel)) {
    cl <- parallel::makeCluster(4)
    doParallel::registerDoParallel(cl)
    tt <- sm(requireNamespace("parallel"))
    tt <- sm(requireNamespace("doParallel"))
    tt <- sm(requireNamespace("iterators"))
    tt <- sm(requireNamespace("foreach"))
    res <- foreach(c = 1:length(names(results)), .packages = c("hpgltools")) %dopar% {
      type <- names(results)[c]
      results[[type]] <- do_pairwise(
          type, input = input, conditions = conditions, batches = batches,
          model_cond = model_cond, model_batch = model_batch, model_intercept = model_intercept,
          extra_contrasts = extra_contrasts, alt_model = alt_model, libsize = libsize,
          annot_df = annot_df,
          ...)
    } ## End foreach() %dopar% { }
    parallel::stopCluster(cl)
    if (isTRUE(verbose)) {
      message("Finished running DE analyses, collecting outputs.")
    }
    ## foreach returns the results in no particular order
    ## Therefore, I will reorder the results now and ensure that they are happy.
    for (r in seq_along(res)) {
      a_result <- res[[r]]
      type <- a_result[["type"]]
      results[[type]] <- a_result
    }
    rm(res)
    ## End performing parallel comparisons
  } else {
    for (type in names(results)) {
      if (isTRUE(verbose)) {
        message("Starting ", type, "_pairwise().")
      }
      results[[type]] <- do_pairwise(
          type, input = input, conditions = conditions, batches = batches,
          model_cond = model_cond, model_batch = model_batch, model_intercept = model_intercept,
          extra_contrasts = extra_contrasts, alt_model = alt_model, libsize = libsize,
          annot_df = annot_df,
          ...)
    }
  } ## End performing a serial comparison

  original_pvalues <- NULL
  ## Add in a little work to re-adjust the p-values in the situation where sva
  ## was used Only perform this f adjustment if you modify the data without
  ## making limma/deseq/edger aware of the modified model.  Ergo, if we feed
  ## sv_model to this function, then by definition, we do not want to use this
  ## function.
  ## Thus we will use modified_data to note the data was modified by sva.
  modified_data <- FALSE
  if (is.null(sv_model) & isTRUE(modified_data)) {
    ret <- sva_modify_pvalues(results)
    results <- ret[["results"]]
    original_pvalues <- ret[["original_pvalues"]]
  }

  result_comparison <- correlate_de_tables(results, annot_df = annot_df,
                                           extra_contrasts = extra_contrasts)
  ## The first few elements of this list are being passed through into the return
  ## So that if I use combine_tables() I can report in the resulting tables
  ## some information about what was performed.
  ret <- list(
      "basic" = results[["basic"]],
      "deseq" = results[["deseq"]],
      "ebseq" = results[["ebseq"]],
      "edger" = results[["edger"]],
      "limma" = results[["limma"]],
      "batch_type" = model_type,
      "comparison" = result_comparison,
      "extra_contrasts" = extra_contrasts,
      "input" = input,
      "model_cond" = model_cond,
      "model_batch" = model_batch,
      "original_pvalues" = original_pvalues,
      "pre_batch" = pre_pca,
      "post_batch" = post_pca)
  class(ret) <- c("all_pairwise", "list")

  if (!is.null(arglist[["combined_excel"]])) {
    if (isTRUE(verbose)) {
      message("Invoking combine_de_tables().")
    }
    combined <- combine_de_tables(ret, excel = arglist[["combined_excel"]], ...)
    ret[["combined"]] <- combined
  }
  return(ret)
}

#' Calculate the Area under the Concordance Curve.
#'
#' This is taken verbatim from a recent paper sent to me by Julie
#' Cridland.  I will put the link in shortly, I need to go.
#'
#' @param tbl DE table
#' @param tbl2 Second table
#' @param px first set of p-values column
#' @param py second set
#' @param lx first set of logFCs column
#' @param ly second set
#' @param topn Number of genes to consider (or percentage of the
#'  whole).
#' @export
calculate_aucc <- function(tbl, tbl2 = NULL, px = "deseq_adjp", py = "edger_adjp",
                           lx = "deseq_logfc", ly = "edger_logfc",
                           topn = 0.1) {
  ## If the topn argument is an integer, the just ask for that number.
  ## If it is a floating point 0<x<1, then set topn to that proportion
  ## of the number of genes.
  if (topn <= 0) {
    stop("topn need to be either a float from 0-1 or an interger bigger than 100.")
  } else if (topn > 1 & topn <= 100) {
    stop("topn need to be either a float from 0-1 or an interger bigger than 100.")
  } else if (topn < 1) {
    topn <- ceiling(nrow(tbl) * topn)
  }

  ## By default, assume all the relevant columns are in the tbl variable.
  ## However, if a second tbl is provided, then do not assume that.
  if (is.null(tbl2)) {
    message("A second table was not providing, performing comparison among columns of the first.")
  } else {
    tmp_tbl <- merge(tbl, tbl2, by = "row.names")
    rownames(tmp_tbl) <- tmp_tbl[["Row.names"]]
    tmp_tbl[["Row.names"]] <- NULL
    tbl <- tmp_tbl
    rm("tmp_tbl")
    if (px == py) {
      px <- paste0(px, ".x")
      py <- paste0(py, ".y")
      lx <- paste0(lx, ".x")
      ly <- paste0(ly, ".y")
    }
  }

  x_df <- tbl[, c(px, lx)]
  x_order <- rownames(x_df)
  y_df <- tbl[x_order, c(py, ly)]

  ## Do a simple correlation test just to have it as a comparison point
  simple_cor <- cor.test(tbl[[lx]], tbl[[ly]])

  ## curve (AUCC), we ranked genes in both the single-cell and bulk datasets in
  ## descending order by the statistical significance of their differential expression.
  x_idx <- order(x_df[[px]], decreasing = FALSE)
  x_df <- x_df[x_idx, ]
  y_idx <- order(y_df[[py]], decreasing = FALSE)
  y_df <- y_df[y_idx, ]

  ## Then, we created lists of the top-ranked genes in each dataset of
  ## matching size, up to some maximum size k. For each of these lists
  ## (that is, for the top-1 genes, top-2 genes, top-3 genes, and so
  ## on), we computed the size of the intersection between the
  ## single-cell and bulk DE genes. This procedure yielded a curve
  ## relating the number of shared genes between datasets to the
  ## number of top-ranked genes considered.

  intersections <- rep(0, topn)
  for (i in seq_len(topn)) {
    if (i == 1) {
      intersections[i] <- rownames(x_df)[i] == rownames(y_df)[i]
    } else {
      x_set <- rownames(x_df)[1:i]
      y_set <- rownames(y_df)[1:i]
      intersections[i] <- sum(x_set %in% y_set)
    }
  }

  ## The area under this curve was computed by summing the size of all
  ## intersections, and normalized to the range [0, 1] by dividing it
  ## by its maximum possible value, k × ( k +1) / 2. To evaluate the
  ## concordance of DE analysis, we used k =500 except where otherwise
  ## noted, but found our results were insensitive to the
  sumint <- sum(intersections)
  norm <- (topn * (topn + 1)) /2
  aucc <- sumint / norm

  intersection_df <- data.frame(x = 1:topn, y = intersections)
  intersection_lm <- lm(intersection_df, formula = y ~ x)
  inter <- as.numeric(coef(intersection_lm))
  inter_inter <- inter[1]
  inter_slope <- inter[2]
  intersection_plot <- ggplot(data = intersection_df,
                              aes_string(x = "x", y = "y")) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(colour = "grey", slope = 1, intercept = 0) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, topn)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::geom_abline(colour = "blue", slope = inter_slope, intercept = inter_inter) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x) +
    ggplot2::annotate("text", x = topn / 1.5, y = topn / 10,
                      label = glue::glue("AUCC: {signif(aucc, 3)},
 AUCC slope: {signif(inter_slope, 3)}"))

  retlist <- list(
      "aucc" = aucc,
      "cor" = simple_cor,
      "plot" = intersection_plot)
 return(retlist)
}

#' Use sva's f.pvalue to adjust p-values for data adjusted by combat.
#'
#' This is from section 5 of the sva manual:  "Adjusting for surrogate
#' values using the f.pvalue function." The following chunk of code is longer
#' and more complex than I would like. This is because f.pvalue() assumes a
#' pairwise comparison of a data set containing only two experimental
#' factors. As a way to provide an example of _how_ to calculate
#' appropriately corrected p-values for surrogate factor adjusted models,
#' this is great; but when dealing with actual data, it falls a bit short.
#'
#' @param results Table of differential expression results.
#' @seealso [sva]
sva_modify_pvalues <- function(results) {
  input <- results[["input"]]
  original_pvalues <- data.table::data.table(
                                      rownames = rownames(results[["edger"]][["all_tables"]][[1]]))
  if (isTRUE(verbose)) {
    message("Using f.pvalue() to modify the returned p-values of deseq/limma/edger.")
  }
  for (it in seq_along(results[["edger"]][["all_tables"]])) {
    name <- names(results[["edger"]][["all_tables"]])[it]
    if (isTRUE(verbose)) {
      message("Readjusting the p-values for comparison: ", name)
    }
    namelst <- strsplit(x = name, split = "_vs_")
    ## something like 'mutant'
    first <- namelst[[1]][[1]]
    ## something like 'wildtype', ergo the contrast was "mutant_vs_wildtype"
    second <- namelst[[1]][[2]]
    ## The comments that follow will use mutant and wildtype as examples

    ## I am going to need to extract the set of data for the samples in
    ## 'first' and 'second'. I will need to also extract the surrogates for
    ## those samples from sv_model. Then I rewrite null_model as the
    ## subset(null_model, samples included) and rewrite sv_model as
    ## subset(sv_model, samples_included) in that rewrite, there will just be
    ## conditions a, b where a and b are the subsets for first and
    ## second. Then the sv_model will be a for the first samples and b for the
    ## second. With that information, I should e able to feed sva::f.pvalue
    ## the appropriate information for it to run properly. The resulting
    ## pvalues will then be appropriate for backfilling the various tables
    ## from edger/limma/deseq.

    ## Get the samples from the limma comparison which are condition 'mutant'
    samples_first_idx <- results[["limma"]][["conditions"]] == first
    num_first <- sum(samples_first_idx)
    ## Subset the expressionset and make a new matrix of only the 'mutant' samples
    samples_first <- exprs(input)[, samples_first_idx]
    ## Repeat for the 'wildtype' samples, when finished the m columns for
    ## samples_first 'mutant' and the n samples of samples_second will be
    ## 'wildtype'
    samples_second_idx <- results[["limma"]][["conditions"]] == second
    num_second <- sum(samples_second_idx)
    samples_second <- exprs(results[["input"]])[, samples_second_idx]
    ## Concatenate the 'mutant' and 'wildtype' samples by column
    included_samples <- cbind(samples_first, samples_second)
    ## Arbitrarily call them 'first' and 'second'
    colnames(included_samples) <- c(rep("first", times = num_first),
                                    rep("second", times = num_second))
    ## Do the same thing, but using the rows of the sva model adjustment
    first_sva <- results[["limma"]][["sv_model"]][samples_first_idx, ]
    second_sva <- results[["limma"]][["sv_model"]][samples_second_idx, ]
    ## But instead of just appending them, I need a matrix of 0s and 1s
    ## identifying which sv rows correspond to 'wildtype' and 'mutant'
    first_model <- append(rep(1, num_first), rep(0, num_second))
    second_model <- append(rep(0, num_first), rep(1, num_second))
    ## Once I have them, make the subset model matrix with append and cbind
    new_sv_model <- append(first_sva, second_sva)
    new_model <- cbind(first_model, second_model, new_sv_model)
    colnames(new_model) <- c("first", "second", "sv")
    ## The sva f.pvalue requires a null model of the appropriate size, create
    ## that here.
    new_null <- cbind(rep(1, (num_first + num_second)), new_sv_model)
    ## And give its columns suitable names
    colnames(new_null) <- c("null", "sv")
    ## Now all the pieces are in place, call f.pvalue().
    new_pvalues <- try(sva::f.pvalue(included_samples, new_model, new_null), silent = TRUE)
    ## For some things, f.pvalue returns NA, this is unacceptable.
    na_pvalues_idx <- is.na(new_pvalues)
    ## For the NaN pvalues, just set it to 1 under the assumption that
    ## something is fubar.
    new_pvalues[na_pvalues_idx] <- 2.2E-16
    ## For strange non-pairwise contrasts, the f.pvalue() should fail.
    if (class(new_pvalues)[1] == "try-error") {
      new_pvalues <- NULL
      warning("Unable to adjust pvalues for: ", name)
      warning("If this was not for an extra contrast, then this is a serious problem.")
    } else {
      ## Most of the time it should succeed, so do a BH adjustment of the new values.
      new_adjp <- p.adjust(new_pvalues, method = "BH")
      ## Now I need to fill in the tables with these new values.
      ## This section is a little complex.  In brief, it pulls the appropriate
      ## columns from each of the limma, edger, and deseq tables
      ## copies them to a new data.table named 'original_pvalues', and
      ## then replaces them with the just-calculated (adj)p-values.

      ## Start with limma, make no assumptions about table return-order
      limma_table_order <- rownames(results[["limma"]][["all_tables"]][[name]])
      reordered_pvalues <- new_pvalues[limma_table_order]
      reordered_adjp <- new_adjp[limma_table_order]
      ## Create a temporary dt with the old p-values and merge it into original_pvalues.
      tmpdt <- data.table::data.table(
                               results[["limma"]][["all_tables"]][[name]][["P.Value"]])
      tmpdt[["rownames"]] <- rownames(results[["limma"]][["all_tables"]][[name]])
      ## Change the column name of the new data to reflect that it is from limma.
      colnames(tmpdt) <- c(glue("limma_{name}"), "rownames")
      original_pvalues <- merge(original_pvalues, tmpdt, by = "rownames")
      ## Swap out the p-values and adjusted p-values.
      results[["limma"]][["all_tables"]][[name]][["P.Value"]] <- reordered_pvalues
      results[["limma"]][["all_tables"]][[name]][["adj.P.Val"]] <- reordered_adjp

      ## Repeat the above verbatim, but for edger
      edger_table_order <- rownames(results[["edger"]][["all_tables"]][[name]])
      reordered_pvalues <- new_pvalues[edger_table_order]
      reordered_adjp <- new_adjp[edger_table_order]
      tmpdt <- data.table::data.table(
                               results[["edger"]][["all_tables"]][[name]][["PValue"]])
      tmpdt[["rownames"]] <- rownames(results[["edger"]][["all_tables"]][[name]])
      colnames(tmpdt) <- c(glue("edger_{name}"), "rownames")
      original_pvalues <- merge(original_pvalues, tmpdt, by = "rownames")
      results[["edger"]][["all_tables"]][[name]][["PValue"]] <- reordered_pvalues
      results[["edger"]][["all_tables"]][[name]][["FDR"]] <- reordered_adjp

      ## Ibid.
      deseq_table_order <- rownames(results[["deseq"]][["all_tables"]][[name]])
      tmpdt <- data.table::data.table(
                               results[["deseq"]][["all_tables"]][[name]][["P.Value"]])
      tmpdt[["rownames"]] <- rownames(results[["deseq"]][["all_tables"]][[name]])
      colnames(tmpdt) <- c(glue("deseq_{name}"), "rownames")
      original_pvalues <- merge(original_pvalues, tmpdt, by = "rownames")
      results[["deseq"]][["all_tables"]][[name]][["P.Value"]] <- reordered_pvalues
      results[["deseq"]][["all_tables"]][[name]][["adj.P.Val"]] <- reordered_adjp
    } ## End checking that f.pvalue worked.
  }  ## End foreach table
  original_pvalues <- as.data.frame(original_pvalues)
  retlist <- list(
      "original" = original_pvalues,
      "results" = results)
  return(retlist)
}

#' A sanity check that a given set of data is suitable for methods which assume
#' a negative binomial distribution of input.
#'
#' Take an expt and poke at it to ensure that it will not result in troubled
#' results.
#'
#' Invoked by deseq_pairwise() and edger_pairwise().
#'
#' @param input Expressionset containing expt object.
#' @param verbose Print some information about what is happening?
#' @param force Ignore every warning and just use this data.
#' @param ... Extra arguments passed to arglist.
#' @return dataset suitable for limma analysis
#' @seealso [DESeq2] [edgeR] [choose_basic_dataset()] [choose_limma_dataset()]
choose_binom_dataset <- function(input, verbose = TRUE, force = FALSE, ...) {
  ## arglist <- list(...)
  input_class <- class(input)
  ## I think I would like to make this function smarter so that it will remove
  ## the log2 from transformed data.
  data <- NULL
  warn_user <- 0
  libsize <- NULL
  if ("expt" %in% input_class) {
    conditions <- input[["conditions"]]
    batches <- input[["batches"]]
    data <- as.data.frame(exprs(input))
    ## As I understand it, EdgeR fits a binomial distribution
    ## and expects data as integer counts, not floating point nor a log2
    ## transformation. Thus, having the 'normalization' state set to something
    ## other than 'raw' is a likely violation of its stated preferred/demanded
    ## input.  There are of course ways around this but one should not take them
    ## lightly, or ever.
    tran_state <- input[["state"]][["transform"]]
    if (is.null(tran_state)) {
      tran_state <- "raw"
    }
    conv_state <- input[["state"]][["conversion"]]
    if (is.null(conv_state)) {
      conv_state <- "raw"
    }
    norm_state <- input[["state"]][["normalization"]]
    if (is.null(norm_state)) {
      norm_state <- "raw"
    }
    filt_state <- input[["state"]][["filter"]]
    if (is.null(filt_state)) {
      filt_state <- "raw"
    }
    if (norm_state == "round") {
      norm_state <- "raw"
    }
    libsize <- input[["libsize"]]

    if (isTRUE(force)) {
      ## Setting force to TRUE allows one to round the data to fool edger/deseq
      ## into accepting it. This is a pretty terrible thing to do
      if (isTRUE(verbose)) {
        message("About to round the data, this is a pretty terrible thing to do. ",
                "But if you, like me, want to see what happens when you put ",
                "non-standard data into deseq, then here you go.")
      }
      data <- round(data)
      less_than <- data < 0
      data[less_than] <- 0
      na_idx <- is.na(data)
      data[na_idx] <- 0
      warn_user <- 1
    } else if (norm_state != "raw" & tran_state != "raw" & conv_state != "raw") {
      ## These if statements may be insufficient to check for the appropriate
      ## input for deseq.
      data <- exprs(input[["original_expressionset"]])
    } else if (norm_state != "raw" | tran_state != "raw") {
      ## This makes use of the fact that the order of operations in the
      ## normalization function is
      ## static. filter->normalization->convert->batch->transform. Thus, if the
      ## normalized state is not raw, we can look back either to the filtered or
      ## original data. The same is true for the transformation state.
      if (isTRUE(verbose)) {
        message("EdgeR/DESeq expect raw data as input, reverting to count filtered data.")
      }
      data <- input[["normalized"]][["intermediate_counts"]][["filter"]][["count_table"]]
      if (is.null(data)) {
        data <- input[["normalized"]][["intermediate_counts"]][["original"]]
      }
    } else {
      if (isTRUE(verbose)) {
        message("The data should be suitable for EdgeR/DESeq/EBSeq. ",
                "If they freak out, check the state of the count table ",
                "and ensure that it is in integer counts.")
      }
    }
    ## End testing if normalization has been performed
  } else {
    data <- as.data.frame(input)
    libsize <- colSums(data)
  }

  retlist <- list(
      "libsize" = libsize,
      "conditions" = conditions,
      "batches" = batches,
      "data" = data)
  if (warn_user == 1) {
    warning("This data was inappropriately forced into integers.")
  }
  return(retlist)
}

#' Choose a suitable data set for Edger/DESeq
#'
#' The _pairwise family of functions all demand data in specific formats.
#' This tries to make that consistent.
#'
#' Invoked by _pairwise().
#'
#' @param input Expt input.
#' @param choose_for One of limma, deseq, edger, or basic.  Defines the
#'  requested data state.
#' @param force Force non-standard data?
#' @param verbose Print some information about what is happening?
#' @param ... More options for future expansion.
#' @return List the data, conditions, and batches in the data.
#' @seealso [choose_binom_dataset()] [choose_limma_dataset()] [choose_basic_dataset()]
#' @examples
#' \dontrun{
#'  starting_data <- create_expt(metadata)
#'  modified_data <- normalize_expt(starting_data, transform = "log2", norm = "quant")
#'  a_dataset <- choose_dataset(modified_data, choose_for = "deseq")
#'  ## choose_dataset should see that log2 data is inappropriate for DESeq2 and
#'  ## return it to a base10 state.
#' }
#' @export
choose_dataset <- function(input, choose_for = "limma", force = FALSE, verbose = TRUE, ...) {
  ## arglist <- list(...)
  result <- NULL
  if (choose_for == "limma") {
    result <- choose_limma_dataset(input, force = force, ...)
  } else if (choose_for == "basic") {
    result <- choose_basic_dataset(input, force = force, ...)
  } else if (choose_for == "edger") {
    result <- choose_binom_dataset(input, force = force, ...)
  } else if (choose_for == "deseq") {
    result <- choose_binom_dataset(input, force = force, ...)
  } else {
    if (isTRUE(verbose)) {
      message("Unknown tool for which to choose a data set.")
    }
  }
  return(result)
}

#' A sanity check that a given set of data is suitable for analysis by limma.
#'
#' Take an expt and poke at it to ensure that it will not result in troubled
#' limma results.
#'
#' @param input Expressionset containing expt object.
#' @param force Ingore warnings and use the provided data asis.
#' @param which_voom Choose between limma'svoom, voomWithQualityWeights, or the
#'  hpgl equivalents.
#' @param verbose Print some information about what is happening?
#' @param ... Extra arguments passed to arglist.
#' @return dataset suitable for limma analysis
#' @seealso [limma] [choose_dataset()]
choose_limma_dataset <- function(input, force = FALSE, which_voom = "limma", verbose = TRUE, ...) {
  ## arglist <- list(...)
  input_class <- class(input)
  data <- NULL
  warn_user <- 0
  libsize <- NULL
  ## It turns out, a more careful examination of how normalization affects the
  ## results, the above seems only to be true if the following are true:
  ## 1.  There are >2-3k features(genes/transcripts) with a full range of count values.
  ## 2.  One does not attempt to use sva, or at least one uses sva before
  ##     messing with the normalization state.
  ## 2a. #2 primarily applies if one is using quantile normalization, it looks
  ##     like tmm/rle does not have so profound an effect, and this effect is
  ##     tightly bound with the state of #1 -- in other words, if one has nice
  ##     dense data with low->high counts in an even distribution, then
  ##     quantile+sva might be ok. But if that is not true, then one should
  ##     expect a poor result.
  ## For these reasons I am telling this function to revert to non-normalized
  ## data unless force is on, just like I do for edger/deseq.  I think to help
  ## this, I will add a parameter which allows one to to turn on/off
  ## normalization at the voom() step.

  if ("expt" %in% input_class) {
    conditions <- input[["conditions"]]
    batches <- input[["batches"]]
    libsize <- input[["libsize"]]
    data <- as.data.frame(exprs(input))

    tran_state <- input[["state"]][["transform"]]
    ## Note that voom will take care of this for us.
    if (is.null(tran_state)) {
      tran_state <- "raw"
    }
    conv_state <- input[["state"]][["conversion"]]
    ## Note that voom takes care of this for us.
    if (is.null(conv_state)) {
      conv_state <- "raw"
    }
    norm_state <- input[["state"]][["normalization"]]
    if (is.null(norm_state)) {
      norm_state <- "raw"
    }
    filt_state <- input[["state"]][["filter"]]
    if (is.null(filt_state)) {
      filt_state <- "raw"
    }

    ## ready <- input
    data <- exprs(input)
    if (isTRUE(force)) {
      if (isTRUE(verbose)) {
        message("Leaving the data alone, regardless of normalization state.")
      }
      retlist <- list(
          "libsize" = libsize,
          "conditions" = conditions,
          "batches" = batches,
          "data" = data)
      return(retlist)
    }

    ## If we are using limma::voom*, then make sure we do things the limma way.
    ## If we use the hpgltools::hpgl_voom*, let the freak flags fly.
    if (grepl(pattern = "limma", x = which_voom)) {
      ## Limma's voom requires we return log2(cpm()) to base 10.
      ## Otherwise it should accept pretty much anything.
      if (tran_state == "log2") {
        if (isTRUE(verbose)) {
          message("Using limma's voom, returning to base 10.")
        }
        data <- (2 ^ data) - 1
      }
    }
  } else {
    data <- as.data.frame(input)
    libsize <- colSums(data)
  }
  retlist <- list(
      "libsize" = libsize,
      "conditions" = conditions,
      "batches" = batches,
      "data" = data)
  return(retlist)
}

#' Try out a few experimental models and return a likely working option.
#'
#' The _pairwise family of functions all demand an experimental model.  This
#' tries to choose a consistent and useful model for all for them.  This does
#' not try to do multi-factor, interacting, nor dependent variable models, if
#' you want those do them yourself and pass them off as alt_model.
#'
#' Invoked by the _pairwise() functions.
#'
#' @param input Input data used to make the model.
#' @param conditions Factor of conditions in the putative model.
#' @param batches Factor of batches in the putative model.
#' @param model_batch Try to include batch in the model?
#' @param model_cond Try to include condition in the model? (Yes!)
#' @param model_intercept Use an intercept model instead of cell-means?
#' @param alt_model Use your own model.
#' @param alt_string String describing an alternate model.
#' @param intercept Choose an intercept for the model as opposed to 0.
#' @param reverse Reverse condition/batch in the model?  This shouldn't/doesn't
#'  matter but I wanted to test.
#' @param contr List of contrasts.arg possibilities.
#' @param surrogates Number of or method used to choose the number of surrogate
#'  variables.
#' @param verbose Print some information about what is happening?
#' @param ... Further options are passed to arglist.
#' @return List including a model matrix and strings describing cell-means and
#'  intercept models.
#' @seealso [stats::model.matrix()]
#' @examples
#' \dontrun{
#'  a_model <- choose_model(expt, model_batch = TRUE, model_intercept = FALSE)
#'  a_model$chosen_model
#'  ## ~ 0 + condition + batch
#' }
#' @export
choose_model <- function(input, conditions = NULL, batches = NULL, model_batch = TRUE,
                         model_cond = TRUE, model_intercept = FALSE,
                         alt_model = NULL, alt_string = NULL,
                         intercept = 0, reverse = FALSE, contr = NULL,
                         surrogates = "be", verbose = TRUE, ...) {
  arglist <- list(...)
  design <- NULL
  if (class(input)[1] != "matrix" & class(input)[1] != "data.frame") {
    design <- pData(input)
  }
  if (is.null(design)) {
    conditions <- as.factor(conditions)
    batches <- as.factor(batches)
    design <- data.frame("condition" = conditions,
                         "batch" = batches,
                         stringsAsFactors = TRUE)
  }
  ## Make a model matrix which will have one entry for
  ## each of the condition/batches
  ## It would be much smarter to generate the models in the following if() {} blocks
  ## But I have it in my head to eventually compare results using different models.

  ## The previous iteration of this had an explicit contrasts.arg set, like this:
  ## contrasts.arg = list(condition = "contr.treatment"))
  ## Which looked like this for a full invocation:
  ## condbatch_int_model <- try(stats::model.matrix(~ 0 + conditions + batches,
  ##                                   contrasts.arg = list(condition = "contr.treatment",
  ##                                                      batch = "contr.treatment")),
  ## The contrasts.arg has been removed because it seems to result in the same model.

  clist <- list("condition" = "contr.treatment")
  blist <- list("batch" = "contr.treatment")
  cblist <- list("condition" = "contr.treatment", "batch" = "contr.treatment")
  if (!is.null(contr)) {
    if (!is.null(contr[["condition"]]) & !is.null(contr[["batch"]])) {
      cblist <- list("condition" = contr[["condition"]], "batch" = contr[["batch"]])
    } else if (!is.null(contr[["condition"]])) {
      clist <- list("condition" = contr[["condition"]])
      cblist[["condition"]] <- contr[["condition"]]
    } else if (!is.null(contr[["batch"]])) {
      blist <- list("batch" = contr[["batch"]])
      cblist[["batch"]] <- contr[["batch"]]
    }
  }

  cond_noint_string <- "~ 0 + condition"
  cond_noint_model <- try(stats::model.matrix(
                                     object = as.formula(cond_noint_string),
                                     contrasts.arg = clist,
                                     data = design), silent = TRUE)

  batch_noint_string <- "~ 0 + batch"
  batch_noint_model <- try(stats::model.matrix(
                                      object = as.formula(batch_noint_string),
                                      contrasts.arg = blist,
                                      data = design), silent = TRUE)

  condbatch_noint_string <- "~ 0 + condition + batch"
  condbatch_noint_model <- try(stats::model.matrix(
                                          object = as.formula(condbatch_noint_string),
                                          contrasts.arg = cblist,
                                          data = design), silent = TRUE)

  batchcond_noint_string <- "~ 0 + batch + condition"
  batchcond_noint_model <- try(stats::model.matrix(
                                          object = as.formula(batchcond_noint_string),
                                          contrasts.arg = cblist,
                                          data = design), silent = TRUE)

  cond_int_string <- "~ condition"
  cond_int_model <- try(stats::model.matrix(
                                   object = as.formula(cond_int_string),
                                   contrasts.arg = clist,
                                   data = design), silent = TRUE)

  batch_int_string <- "~ batch"
  batch_int_model <- try(stats::model.matrix(
                                    object = as.formula(batch_int_string),
                                    contrasts.arg = blist,
                                    data = design), silent = TRUE)

  condbatch_int_string <- "~ condition + batch"
  condbatch_int_model <- try(stats::model.matrix(
                                        object = as.formula(condbatch_int_string),
                                        contrasts.arg = cblist,
                                        data = design), silent = TRUE)

  batchcond_int_string <- "~ batch + condition"
  batchcond_int_model <- try(stats::model.matrix(
                                        object = as.formula(batchcond_int_string),
                                        contrasts.arg = cblist,
                                        data = design), silent = TRUE)

  noint_model <- NULL
  int_model <- NULL
  noint_string <- NULL
  int_string <- NULL
  including <- NULL
  if (!is.null(alt_model)) {
    chosen_model <- stats::model.matrix(object = as.formula(alt_model),
                                        data = design)
    if (!is.null(contr)) {
      chosen_model <- stats::model.matrix(object = as.formula(alt_model),
                                          contrasts.arg = contr,
                                          data = design)
    }
    int_model <- chosen_model
    noint_model <- chosen_model
    int_string <- alt_model
    notint_string <- alt_model
    including <- "alt"
  } else if (is.null(model_batch)) {
    int_model <- cond_int_model
    noint_model <- cond_noint_model
    int_string <- cond_int_string
    noint_string <- cond_noint_string
    including <- "condition"
  } else if (isTRUE(model_cond) & isTRUE(model_batch)) {
    if (class(condbatch_int_model)[1] == "try-error") {
      if (isTRUE(verbose)) {
        message("The condition+batch model failed. ",
                "Does your experimental design support both condition and batch? ",
                "Using only a conditional model.")
      }
      int_model <- cond_int_model
      noint_model <- cond_noint_model
      int_string <- cond_int_string
      noint_string <- cond_noint_string
      including <- "condition"
    } else if (isTRUE(reverse)) {
      int_model <- batchcond_int_model
      noint_model <- batchcond_noint_model
      int_string <- batchcond_int_string
      noint_string <- batchcond_noint_string
      including <- "batch+condition"
    } else {
      int_model <- condbatch_int_model
      noint_model <- condbatch_noint_model
      int_string <- condbatch_int_string
      noint_string <- condbatch_noint_string
      including <- "condition+batch"
    }
  } else if (class(model_batch)[1] == "character") {
    ## Then calculate the estimates using all_adjusters
    if (isTRUE(verbose)) {
      message("Extracting surrogate estimates from ", model_batch,
              " and adding them to the model.")
    }
    model_batch_info <- all_adjusters(input, estimate_type = model_batch,
                                      surrogates = surrogates)
    ## Changing model_batch from 'sva' to the resulting matrix.
    ## Hopefully this will simplify things later for me.
    model_batch <- model_batch_info[["model_adjust"]]
    int_model <- stats::model.matrix(~ condition + model_batch,
                                     contrasts.arg = clist,
                                     data = design)
    noint_model <- stats::model.matrix(~ 0 + condition + model_batch,
                                       contrasts.arg = clist,
                                       data = design)
    sv_names <- glue("SV{1:ncol(model_batch)}")
    noint_string <- cond_noint_string
    int_string <- cond_int_string
    sv_string <- ""
    for (sv in sv_names) {
      sv_string <- glue("{sv_string} + {sv}")
    }
    noint_string <- glue("{noint_string}{sv_string}")
    int_string <- glue("{int_string}{sv_string}")
    rownames(model_batch) <- rownames(int_model)
    including <- glue("condition{sv_string}")
  } else if (class(model_batch)[1] == "numeric" | class(model_batch)[1] == "matrix") {
    if (isTRUE(verbose)) {
      message("Including batch estimates from sva/ruv/pca in the model.")
    }
    int_model <- stats::model.matrix(~ condition + model_batch,
                                     contrasts.arg = clist,
                                     data = design)
    noint_model <- stats::model.matrix(~ 0 + condition + model_batch,
                                       contrasts.arg = clist,
                                       data = design)
    sv_names <- glue("SV{1:ncol(model_batch)}")
    int_string <- cond_int_string
    noint_string <- cond_noint_string
    sv_string <- ""
    for (sv in sv_names) {
      sv_string <- glue("{sv_string} + {sv}")
    }
    int_string <- glue("{int_string}{sv_string}")
    noint_string <- glue("{noint_string}{sv_string}")
    rownames(model_batch) <- rownames(int_model)
    including <- glue("condition{sv_string}")
  } else if (isTRUE(model_cond)) {
    int_model <- cond_int_model
    noint_model <- cond_noint_model
    int_string <- cond_int_string
    noint_string <- cond_noint_string
    including <- "condition"
  } else if (isTRUE(model_batch)) {
    int_model <- batch_int_model
    noint_model <- batch_noint_model
    int_string <- batch_int_string
    noint_string <- batch_noint_string
    including <- "batch"
  } else {
    ## Default to the conditional model
    int_model <- cond_int_model
    noint_model <- cond_noint_model
    int_string <- cond_int_string
    noint_string <- cond_noint_string
    including <- "condition"
  }

  tmpnames <- colnames(int_model)
  tmpnames <- gsub(pattern = "model_batch", replacement = "SV1", x = tmpnames)
  tmpnames <- gsub(pattern = "data[[:punct:]]", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "-", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "+", replacement = "", x = tmpnames)
  ## The next lines ensure that conditions/batches which are all numeric will
  ## not cause weird errors for contrasts. Ergo, if a condition is something
  ## like '111', now it will be 'c111' Similarly, a batch '01' will be 'b01'
  tmpnames <- gsub(pattern = "^condition(\\d+)$", replacement = "c\\1", x = tmpnames)
  tmpnames <- gsub(pattern = "^batch(\\d+)$", replacement = "b\\1", x = tmpnames)
  tmpnames <- gsub(pattern = "condition", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "batch", replacement = "", x = tmpnames)
  colnames(int_model) <- tmpnames

  tmpnames <- colnames(noint_model)
  tmpnames <- gsub(pattern = "model_batch", replacement = "SV1", x = tmpnames)
  tmpnames <- gsub(pattern = "data[[:punct:]]", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "-", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "+", replacement = "", x = tmpnames)
  ## The next lines ensure that conditions/batches which are all numeric will
  ## not cause weird errors for contrasts. Ergo, if a condition is something
  ## like '111', now it will be 'c111' Similarly, a batch '01' will be 'b01'
  tmpnames <- gsub(pattern = "condition^(\\d+)$", replacement = "c\\1", x = tmpnames)
  tmpnames <- gsub(pattern = "batch^(\\d+)$", replacement = "b\\1", x = tmpnames)
  tmpnames <- gsub(pattern = "condition", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "batch", replacement = "", x = tmpnames)
  colnames(noint_model) <- tmpnames

  chosen_model <- NULL
  chosen_string <- NULL
  if (isTRUE(model_intercept)) {
    if (isTRUE(verbose)) {
      message("Choosing the intercept containing model.")
    }
    chosen_model <- int_model
    chosen_string <- int_string
  } else {
    if (isTRUE(verbose)) {
      message("Choosing the non-intercept containing model.")
    }
    chosen_model <- noint_model
    chosen_string <- noint_string
  }

  retlist <- list(
      "int_model" = int_model,
      "noint_model" = noint_model,
      "int_string" = int_string,
      "noint_string" = noint_string,
      "chosen_model" = chosen_model,
      "chosen_string" = chosen_string,
      "model_batch" = model_batch,
      "including" = including)
  return(retlist)
}

#' Compare the results of separate all_pairwise() invocations.
#'
#' Where compare_led_tables looks for changes between limma and friends, this
#' function looks for differences/similarities across the models/surrogates/etc
#' across invocations of limma/deseq/edger.
#'
#' Tested in 29de_shared.R
#'
#' @param first One invocation of combine_de_tables to examine.
#' @param second A second invocation of combine_de_tables to examine.
#' @param cor_method Method to use for cor.test().
#' @param try_methods List of methods to attempt comparing.
#' @return A list of compared columns, tables, and methods.
#' @seealso [all_pairwise()]
#' @examples
#' \dontrun{
#'  first <- all_pairwise(expt, model_batch = FALSE, excel = "first.xlsx")
#'  second <- all_pairwise(expt, model_batch = "svaseq", excel = "second.xlsx")
#'  comparison <- compare_de_results(first$combined, second$combined)
#' }
#' @export
compare_de_results <- function(first, second, cor_method = "pearson",
                               try_methods = c("limma", "deseq", "edger")) {

  result <- list()
  logfc_result <- list()
  p_result <- list()
  adjp_result <- list()
  comparisons <- c("logfc", "p", "adjp")
  ## First make sure we can collect the data for each differential expression method.
  methods <- c()
  for (m in seq_along(try_methods)) {
    method <- try_methods[m]
    message("Testing method: ", method, ".")
    test_column <- glue("{method}_logfc")
    if (is.null(first[["data"]][[1]][[test_column]])) {
      message("The first datum is missing method: ", method, ".")
    } else if (is.null(second[["data"]][[1]][[test_column]])) {
      message("The second datum is missing method: ", method, ".")
    } else {
      message("Adding method: ", method, " to the set.")
      methods <- c(methods, method)
    }
  }

  ## Now start building tables containing the correlations between the methods/contrasts.
  for (m in seq_along(methods)) {
    method <- methods[m]
    result[[method]] <- list()
    tables <- names(first[["data"]])
    for (t in seq_along(tables)) {
      table <- tables[t]
      message(" Starting method ", method, ", table ", table, ".")
      result[[method]][[table]] <- list()
      for (c in seq_along(comparisons)) {
        comparison <- comparisons[c]
        column_name <- glue("{method}_{comparison}")
        f_column <- as.vector(as.numeric(first[["data"]][[table]][[column_name]]))
        s_column <- as.vector(as.numeric(second[["data"]][[table]][[column_name]]))
        if (length(f_column) == 0 | length(s_column) == 0) {
          ## The column of data does not exist.
          break
        }
        names(f_column) <- rownames(first[["data"]][[table]])
        names(s_column) <- rownames(second[["data"]][[table]])
        fs <- merge(f_column, s_column, by = "row.names")
        comp <- cor(x = fs[["x"]], y = fs[["y"]], method = cor_method)
        result[[method]][[table]][[comparison]] <- comp
        if (comparison == "logfc") {
          logfc_result[table] <- comp
          tmp <- as.vector(as.numeric(logfc_result))
          names(tmp) <- names(logfc_result)
          logfc_result <- tmp
        } else if (comparison == "p") {
          p_result[table] <- comp
          tmp <- as.vector(as.numeric(p_result))
          names(tmp) <- names(p_result)
          p_result <- tmp
        } else if (comparison == "adjp") {
          adjp_result[table] <- comp
          tmp <- as.vector(as.numeric(adjp_result))
          names(tmp) <- names(adjp_result)
          adjp_result <- tmp
        }
      }
    }
  }
  comp_df <- data.frame(row.names = names(result[[1]]))
  p_df <- data.frame(row.names = names(result[[1]]))
  adjp_df <- data.frame(row.names = names(result[[1]]))
  cols <- names(result)
  rows <- names(result[[1]])
  for (i in seq_along(rows)) {
    row <- rows[i]
    for (j in seq_along(cols)) {
      col <- cols[j]
      if (col == "basic") {
        next
      }
      element <- result[[col]][[row]][["logfc"]]
      comp_df[row, col] <- element
      element <- result[[col]][[row]][["p"]]
      p_df[row, col] <- element
      element <- result[[col]][[row]][["adjp"]]
      adjp_df[row, col] <- element
    }
  }

  original <- par(mar = c(7, 4, 4, 2) + 0.1)
  text_size <- 1.0
  heat_colors <- grDevices::colorRampPalette(c("white", "darkblue"))
  lfc_heatmap <- try(heatmap.3(as.matrix(comp_df), scale = "none",
                               trace = "none", keysize = 1.5, linewidth = 0.5,
                               margins = c(12, 8), cexRow = text_size, cexCol = text_size,
                               col = heat_colors, dendrogram = "none",
                               Rowv = FALSE, Colv = FALSE,
                               main = "Compare lFC results"), silent = TRUE)
  lfc_heat <- NULL
  if (! "try-error" %in% class(lfc_heatmap)) {
    lfc_heat <- recordPlot()
  }
  heat_colors <- grDevices::colorRampPalette(c("white", "darkred"))
  p_heatmap <- try(heatmap.3(as.matrix(p_df), scale = "none",
                             trace = "none", keysize = 1.5, linewidth = 0.5,
                             margins = c(12, 8), cexRow = text_size, cexCol = text_size,
                             col = heat_colors, dendrogram = "none",
                             Rowv = FALSE, Colv = FALSE,
                             main = "Compare p-values"), silent = TRUE)
  p_heat <- NULL
  if (! "try-error" %in% class(p_heatmap)) {
    p_heat <- recordPlot()
  }
  heat_colors <- grDevices::colorRampPalette(c("white", "darkgreen"))
  adjp_heatmap <- try(heatmap.3(as.matrix(adjp_df), scale = "none",
                                trace = "none", keysize = 1.5, linewidth = 0.5,
                                margins = c(12, 8), cexRow = text_size, cexCol = text_size,
                                col = heat_colors, dendrogram = "none",
                                Rowv = FALSE, Colv = FALSE,
                                main = "Compare adjp-values"), silent = TRUE)
  adjp_heat <- NULL
  if (! "try-error" %in% class(adjp_heatmap)) {
    adjp_heat <- recordPlot()
  }
  new <- par(original)

  retlist <- list(
      "result" = result,
      "logfc" = logfc_result,
      "p" = p_result,
      "adjp" = adjp_result,
      "lfc_heat" = lfc_heat,
      "p_heat" = p_heat,
      "adjp_heat" = adjp_heat)
  return(retlist)
}

compare_de_tables <- function(first, second, fcx = "deseq_logfc", px = "deseq_adjp",
                              fcy = "deseq_logfc", py = "deseq_adjp",
                              first_table = NULL, second_table = NULL) {
  if (!is.null(first_table)) {
    ## Then assume this is a combine_de_tables() result.
    first <- first[["data"]][[first_table]]
  }
  if (!is.null(second_table)) {
    second <- second[["data"]][[second_table]]
  }
  merged <- merge(first, second, by = "row.names")
  rownames(merged) <- merged[["Row.names"]]
  merged[["Row.names"]] <- NULL
  kept_columns <- c(fcx, px, fcy, py)
  if (fcx == fcy) {
    kept_columns <- c(paste0(fcx, ".x"), paste0(px, ".x"),
                      paste0(fcy, ".y"), paste0(fcy, ".y"))
  }
  merged <- merged[, kept_columns]
  colnames(merged) <- c("first_lfc", "first_p", "second_lfc", "second_p")
  scatter <- plot_linear_scatter(merged[, c("first_lfc", "second_lfc")])
  return(scatter)
}

#' See how similar are results from limma/deseq/edger/ebseq.
#'
#' limma, DEseq2, and EdgeR all make somewhat different assumptions.
#' and choices about what makes a meaningful set of differentially.
#' expressed genes.  This seeks to provide a quick and dirty metric
#' describing the degree to which they (dis)agree.
#'
#' Invoked by all_pairwise().
#'
#' @param results Data from do_pairwise()
#' @param annot_df Include annotation data?
#' @return Heatmap showing how similar they are along with some
#'  correlations betwee the three players.
#' @param extra_contrasts include some extra contrasts when comparing results.
#' @seealso [limma_pairwise()] [edger_pairwise()] [deseq_pairwise()]
#' @examples
#' \dontrun{
#'  l = limma_pairwise(expt)
#'  d = deseq_pairwise(expt)
#'  e = edger_pairwise(expt)
#'  fun = compare_led_tables(limma = l, deseq = d, edger = e)
#' }
#' @export
correlate_de_tables <- function(results, annot_df = NULL, extra_contrasts = NULL) {
  ## Fill each column/row of these with the correlation between tools for one
  ## contrast performed
  retlst <- list()
  methods <- c()
  if (class(results[["limma"]])[1] == "limma_result") {
    retlst[["limma"]] <- results[["limma"]][["all_tables"]]
    methods <- c(methods, "limma")
  }
  if (class(results[["deseq"]])[1] == "deseq_result") {
    retlst[["deseq"]] <- results[["deseq"]][["all_tables"]]
    methods <- c(methods, "deseq")
  }
  if (class(results[["edger"]])[1] == "edger_result") {
    retlst[["edger"]] <- results[["edger"]][["all_tables"]]
    methods <- c(methods, "edger")
  }
  if (class(results[["ebseq"]])[1] == "ebseq_result") {
    retlst[["ebseq"]] <- results[["ebseq"]][["all_tables"]]
    methods <- c(methods, "ebseq")
  }
  if (class(results[["basic"]])[1] == "basic_result") {
    retlst[["basic"]] <- results[["basic"]][["all_tables"]]
    methods <- c(methods, "basic")
  }

  extra_eval_names <- NULL
  if (!is.null(extra_contrasts)) {
    extra_eval_strings <- strsplit(extra_contrasts, ",")[[1]]
    extra_eval_names <- extra_eval_strings
    extra_eval_names <- stringi::stri_replace_all_regex(extra_eval_strings,
                                                        "^(\\s*)(\\w+)\\s*=\\s*.*$", "$2")
    extra_eval_names <- gsub(pattern = "^\\s+", replacement = "", x = extra_eval_names, perl = TRUE)
  }

  complst <- list()
  plotlst <- list()
  comparison_df <- data.frame()
  lenminus <- length(methods) - 1
  message("Comparing analyses.")
  meth <- methods[1]
  len <- length(names(retlst[[meth]]))
  total_comparisons <- lenminus * (length(methods) - 1) * len
  progress_count <- 0
  if (isTRUE(verbose)) {
    bar <- utils::txtProgressBar(style = 3)
  }
  for (c in seq_len(lenminus)) {
    c_name <- methods[c]
    nextc <- c + 1
    for (d in seq(from = nextc, to = length(methods))) {
      d_name <- methods[d]
      method_comp_name <- glue("{c_name}_vs_{d_name}")
      contrast_name_list <- c()
      for (l in seq_len(len)) {
        progress_count <- progress_count + 1
        if (isTRUE(verbose)) {
          pct_done <- progress_count / total_comparisons
          utils::setTxtProgressBar(bar, pct_done)
        }
        contr <- names(retlst[[c_name]])[l]
        if (contr %in% extra_eval_names) {
          next
        }

        ## assume all three have the same names() -- note that limma has more
        ## than the other two though
        num_den_names <- strsplit(x = contr, split = "_vs_")[[1]]
        num_name <- num_den_names[1]
        den_name <- num_den_names[2]
        rev_contr <- glue("{den_name}_vs_{num_name}")
        num_reversed <- 0
        fst <- retlst[[c_name]][[contr]]
        scd <- retlst[[d_name]][[contr]]
        if (is.null(fst)) {
          fst <- retlst[[c_name]][[rev_contr]]
          fst[["logFC"]] <- fst[["logFC"]] * -1
          message("Used reverse contrast for ", c_name, ".")
          num_reversed <- num_reversed + 1
        }
        if (is.null(scd)) {
          scd <- retlst[[d_name]][[rev_contr]]
          scd[["logFC"]] <- scd[["logFC"]] * -1
          message("Used reverse contrast for ", d_name, ".")
          num_reversed <- num_reversed + 1
        }
        ## An extra test condition in case of extra contrasts not performed by all methods.
        if (is.na(contr)) {
          next
        }
        contrast_name_list <- c(contr, contrast_name_list)
        fs <- merge(fst, scd, by = "row.names")
        if (nrow(fs) == 0) {
          warning("The merge of ", c_name, ", ", contr, " and ",
                  d_name, ", ", contr, " failed.")
          next
        }
        fs <- fs[, c("logFC.x", "logFC.y")]
        colnames(fs) <- c(glue("{c_name} logFC"), glue("{d_name} logFC"))
        fs_cor <- stats::cor.test(x = fs[, 1], y = fs[, 2])[["estimate"]]
        comparison_df[method_comp_name, contr] <- fs_cor
        fs_plt <- plot_scatter(fs) +
          ggplot2::labs(title = glue("{contr}: {c_name} vs. {d_name}.")) +
          ggplot2::geom_abline(intercept = 0.0, slope = 1.0, colour = "blue")
        complst[[method_comp_name]] <- fs_cor
        plotlst[[method_comp_name]] <- fs_plt
      } ## End iterating through the contrasts
    } ## End the second method loop
  } ## End the first method loop
  if (isTRUE(verbose)) {
    close(bar)
  }

  comparison_df <- as.matrix(comparison_df)
  ## I think this next line is a likely source of errors because
  ## of differences when using extra_contrasts.
  ## colnames(comparison_df) <- names(retlst[["deseq"]])
  colnames(comparison_df) <- contrast_name_list

  heat_colors <- grDevices::colorRampPalette(c("white", "black"))
  original <- par(mar = c(7, 4, 4, 2) + 0.1)
  comparison_heatmap <- try(heatmap.3(comparison_df, scale = "none",
                                      trace = "none", keysize = 1.5,
                                      cexCol = 1.0, cexRow = 1.0,
                                      linewidth = 0.5, margins = c(12, 8),
                                      col = heat_colors, dendrogram = "none",
                                      Rowv = FALSE, Colv = FALSE,
                                      main = "Compare DE tools"), silent = TRUE)
  new <- par(original)
  heat <- NULL
  if (! "try-error" %in% class(comparison_heatmap)) {
    heat <- recordPlot()
  }
  ret <- append(complst, plotlst)
  ret[["comp"]] <- comparison_df
  ret[["heat"]] <- heat
  return(ret)
}

#' Compare logFC values from limma and friends
#'
#' There are some peculiar discrepencies among these tools, what is up with that?
#'
#' Invoked by combine_de_tables() in order to compare the results.
#'
#' @param combined_tables The combined tables from limma et al.
#' @return Some plots
#' @seealso [plot_linear_scatter()]
#' @examples
#' \dontrun{
#'  limma_vs_deseq_vs_edger <- compare_logfc_plots(combined)
#'  ## Get a list of plots of logFC by contrast of LvD, LvE, DvE
#'  ## It provides comparisons against the basic analysis, but who cares about that.
#' }
#' @export
compare_logfc_plots <- function(combined_tables) {
  plots <- list()
  data <- NULL
  if (!is.null(combined_tables[["data"]])) {
    data <- combined_tables[["data"]]
  } else {
    data <- combined_tables
  }
  tnames <- names(data)
  retlist <- list()
  for (c in seq_along(tnames)) {
    tname <- tnames[c]
    tab <- data[[tname]]
    if (!is.null(tab[["limma_logfc"]]) & !is.null(tab[["edger_logfc"]])) {
      le_data <- tab[, c("limma_logfc", "edger_logfc", "limma_adjp", "edger_adjp")]
      le <- sm(plot_linear_scatter(le_data, pretty_colors = FALSE)[["scatter"]])
    } else {
      le <- NULL
    }
    if (!is.null(tab[["limma_logfc"]]) & !is.null(tab[["deseq_logfc"]])) {
      ld_data <- tab[, c("limma_logfc", "deseq_logfc", "limma_adjp", "deseq_adjp")]
      ld <- sm(plot_linear_scatter(ld_data, pretty_colors = FALSE)[["scatter"]])
    } else {
      ld <- NULL
    }
    if (!is.null(tab[["deseq_logfc"]]) & !is.null(tab[["edger_logfc"]])) {
      de_data <- tab[, c("deseq_logfc", "edger_logfc", "deseq_adjp", "edger_adjp")]
      de <- sm(plot_linear_scatter(de_data, pretty_colors = FALSE)[["scatter"]])
    } else {
      de <- NULL
    }
    if (!is.null(tab[["limma_logfc"]]) & !is.null(tab[["basic_logfc"]])) {
      lb_data <- tab[, c("limma_logfc", "basic_logfc", "limma_adjp", "basic_p")]
      lb <- sm(plot_linear_scatter(lb_data, pretty_colors = FALSE)[["scatter"]])
    } else {
      lb <- NULL
    }
    if (!is.null(tab[["deseq_logfc"]]) & !is.null(tab[["basic_logfc"]])) {
      db_data <- tab[, c("deseq_logfc", "basic_logfc", "deseq_adjp", "basic_p")]
      db <- sm(plot_linear_scatter(db_data, pretty_colors = FALSE)[["scatter"]])
    } else {
      db <- NULL
    }
    if (!is.null(tab[["edger_logfc"]]) & !is.null(tab[["basic_logfc"]])) {
      eb_data <- tab[, c("edger_logfc", "basic_logfc", "edger_adjp", "basic_p")]
      eb <- sm(plot_linear_scatter(eb_data, pretty_colors = FALSE)[["scatter"]])
    } else {
      db <- NULL
    }
    compared <- list(
        "name" = tname,
        "le" = le, "ld" = ld, "de" = de,
        "lb" = lb, "db" = db, "eb" = eb)
    retlist[[tname]] <- compared
  }
  return(retlist)
}

#' Implement a cleaner version of 'subset_significants' from analyses with Maria
#' Adelaida.
#'
#' This should provide nice venn diagrams and some statistics to compare 2 or 3
#' contrasts in a differential expression analysis.
#'
#' @param sig_tables Set of significance tables to poke at.
#' @param compare_by Use which program for the comparisons?
#' @param weights When printing venn diagrams, weight them?
#' @param contrasts List of contrasts to compare.
#' @return List containing the intersections of the contrasts and plots describing them.
#' @seealso [Vennerable]
#' @export
compare_significant_contrasts <- function(sig_tables, compare_by = "deseq",
                                          weights = FALSE, contrasts = c(1, 2, 3)) {
  retlist <- NULL
  contrast_names <- names(sig_tables[[compare_by]][["ups"]])
  for (i in seq_along(contrasts)) {
    contr <- contrasts[i]
    if (is.numeric(contr)) {
      contrasts[i] <- contrast_names[i]
    } else if (!is.na(sm(as.numeric(contr)))) {
      ## If one changes one number to character, then they all get recast, ergo this foolishness.
      contrasts[i] <- contrast_names[i]
    }
  }
  up_lst <- list()
  down_lst <- list()
  for (c in seq_along(contrasts)) {
    contr <- contrasts[c]
    up_lst[[contr]] <- rownames(sig_tables[[compare_by]][["ups"]][[contr]])
    down_lst[[contr]] <- rownames(sig_tables[[compare_by]][["downs"]][[contr]])
  }

  up_venn <- Vennerable::Venn(Sets = up_lst)
  up_intersect <- rename_vennerable_intersections(up_venn, up_lst)
  down_venn <- Vennerable::Venn(Sets = down_lst)
  down_intersect <- rename_vennerable_intersections(down_venn, down_lst)

  up_tables <- get_vennerable_rows(sig_tables[[compare_by]][["ups"]], up_intersect)
  down_tables <- get_vennerable_rows(sig_tables[[compare_by]][["downs"]], down_intersect)

  up_plot <- Vennerable::plot(up_venn, doWeights = weights)
  up_plot <- grDevices::recordPlot(up_plot)
  down_plot <- Vennerable::plot(down_venn, doWeights = weights)
  down_plot <- grDevices::recordPlot(down_plot)

  retlst <- list(
      "up_intersections" = up_intersect,
      "down_intersections" = down_intersect,
      "up_venn" = up_venn,
      "down_venn" = down_venn,
      "up_tables" = up_tables,
      "down_tables" = down_tables,
      "up_plot" = up_plot,
      "down_plot" = down_plot)
  return(retlst)
}

#' Test for infected/control/beads -- a placebo effect?
#'
#' This was a function I copied out of Keith/Hector/Laura/Cecilia's paper in
#' which they sought to discriminate the effect of inert beads on macrophages
#' vs. the effect of parasites.  The simpler way of expressing it is: take the
#' worst p-value observed for the pair of contrasts, infected/uninfected and
#' beads/uninfected.
#'
#' The goal is therefore to find responses different than beads
#' The null hypothesis is (H0): (infected == uninfected) | (infected == beads)
#' The alt hypothesis is (HA): (infected != uninfected) & (infected != beads)
#'
#' @param contrast_fit Result of lmFit.
#' @param cellmeans_fit Result of a cellmeans fit.
#' @param conj_contrasts Result from the makeContrasts of the first set.
#' @param disj_contrast Result of the makeContrasts of the second set.
disjunct_pvalues <- function(contrast_fit, cellmeans_fit, conj_contrasts, disj_contrast) {
  ## arglist <- list(...)
  contr_level_counts <- rowSums(
      contrast_fit[["contrasts"]][, c(conj_contrasts, disj_contrast)] != 0)
  ## Define the condition levels involved in the tests
  levels_to_use <- names(contr_level_counts)[contr_level_counts > 0]
  ## Extract the average counts for each, make into table
  ave_expression_mat <- cellmeans_fit[["coef"]][, levels_to_use]
  exp_table <- data.frame("ID" = rownames(ave_expression_mat), stringsAsFactors = FALSE)
  exp_table <- cbind(exp_table, as.data.frame(ave_expression_mat))
  names(exp_table)[-1] <- paste("AveExpr", gsub("condition", "", levels_to_use), sep = ":")

  stat <- BiocGenerics::pmin(abs(contrast_fit[["t"]][, conj_contrasts]))
  pval <- BiocGenerics::pmax(contrast_fit[["p.value"]][, conj_contrasts])

  adj.pval <- p.adjust(pval, method = "BH")
  fcs <- as.data.frame(contrast_fit[["coef"]][, conj_contrasts])
  names(fcs) <- paste("logFC", names(fcs), sep = ":")
  conj_pvals <- as.data.frame(
      apply(contrast_fit[["p.value"]][, conj_contrasts], 2, p.adjust, method = "BH"))
  names(conj_pvals) <- paste("adj.P.Val", names(conj_pvals), sep = ":")
  conj_table <- data.frame("ID" = rownames(contrast_fit), stringsAsFactors = FALSE)
  conj_table <- cbind(conj_table, fcs, conj_pvals, stat = stat, adj.P.Value = adj.pval)
  names(conj_table)[seq(2 + 2 * length(conj_contrasts), ncol(conj_table))] <- paste(
      c("stat", "adj.P.Value"), paste(conj_contrasts, collapse = ":"), sep = ":")

  ## Make the table for the 'other' test
  disj_table <- data.frame(
      "ID" = rownames(contrast_fit),
      "logFC" = contrast_fit[["coef"]][, disj_contrast],
      "adj.P.Value" = p.adjust(contrast_fit[["p.value"]][, disj_contrast], method = "BH"),
      stringsAsFactors = FALSE)
  names(disj_table)[-1] <- paste(c("logFC", "adj.P.Value"), disj_contrast, sep = ":")
  ## Combine tables, making sure all tables are in the same order
  stopifnot(all(exp_table$ID == conj_table[["ID"]] & exp_table[["ID"]] == disj_table[["ID"]]))
  out_table <- cbind(exp_table, conj_table[, -1], disj_table[, -1])
  ## order output table by the statistic in the disjunctive test
  o <- order(-stat)
  out_table <- out_table[o, ]
  return(out_table)
}

#' Generalize pairwise comparisons
#'
#' I want to multithread my pairwise comparisons, this is the first step in
#' doing so.
#'
#' Used to make parallel operations easier.
#'
#' @param type Which type of pairwise comparison to perform
#' @param ... Set of arguments intended for limma_pairwise(),
#'  edger_pairwise(), and friends.
#' @return Result from limma/deseq/edger/basic
#' @seealso [all_pairwise()]
#' @export
do_pairwise <- function(type, ...) {
  res <- NULL
  if (type == "limma") {
    res <- try(limma_pairwise(...))
  } else if (type == "ebseq") {
    res <- try(ebseq_pairwise(...))
  } else if (type == "edger") {
    res <- try(edger_pairwise(...))
  } else if (type == "deseq") {
    res <- try(deseq2_pairwise(...))
  } else if (type == "basic") {
    res <- try(basic_pairwise(...))
  }
  res[["type"]] <- type
  return(res)
}

#' Find the set of most/least abundant genes according to limma and friends
#' following a differential expression analysis.
#'
#' Given a data set provided by limma, deseq, edger, etc; one might want to know
#' what are the most and least abundant genes, much like get_sig_genes() does to
#' find the most significantly different genes for each contrast.
#'
#' @param datum Output from the _pairwise() functions.
#' @param type Extract abundant genes according to what?
#' @param n Perhaps take just the top/bottom n genes.
#' @param z Or take genes past a given z-score.
#' @param unique Unimplemented: take only the genes unique among the conditions surveyed.
#' @return List of data frames containing the genes of interest.
#' @seealso [get_sig_genes()]
#' @examples
#' \dontrun{
#'  abundant <- get_abundant_genes(all_pairwise_output, type = "deseq", n = 100)
#'  ## Top 100 most abundant genes from deseq
#'  least <- get_abundant_genes(all_pairwise_output, type = "deseq", n = 100, least = TRUE)
#'  ## Top 100 least abundant genes from deseq
#'  abundant <- get_abundant_genes(all_pairwise_output, type = "edger", z = 1.5)
#'  ## Get the genes more than 1.5 standard deviations from the mean.
#' }
#' @export
get_abundant_genes <- function(datum, type = "limma", n = NULL, z = NULL,
                               unique = FALSE) {
  if (is.null(z) & is.null(n)) {
    n <- 100
  }
  if (!is.null(datum[["limma"]])) {
    if (type == "basic") {
      datum <- datum[["basic"]]
    } else if (type == "edger") {
      datum <- datum[["edger"]]
    } else if (type == "deseq") {
      datum <- datum[["deseq"]]
    } else {
      datum <- datum[["limma"]]
    }
  }

  ## Extract the coefficent df
  if (type == "edger") {
    ## In this case of edger, this can be a little tricky.
    ## I should probably therefore improve the returns from edger_pairwise()
    ## FIXME: I think this is the wrong way to handle this.
    coefficient_df <- datum[["lrt"]][[1]][["coefficients"]]
    if (max(coefficient_df) <= 0) {
      coefficient_df <- coefficient_df * -1.0
    }
    ## There are a couple of extraneous columns in this table.
    removers <- c("b", "z")
    keepers <- !(colnames(coefficient_df) %in% removers)
    coefficient_df <- coefficient_df[, keepers]
  } else if (type == "limma") {
    coefficient_df <- datum[["identity_comparisons"]][["coefficients"]]
    all_coefficients <- colnames(coefficient_df)
    keepers <- !grepl(pattern = "_vs_", x = all_coefficients)
    coefficient_df <- coefficient_df[, keepers]
  } else if (type == "deseq") {
    coefficient_df <- datum[["coefficients"]]
  } else if (type == "basic") {
    coefficient_df <- datum[["medians"]]
  }

  abundant_list <- list(
      "high" = list(),
      "low" = list())
  coefficient_df <- as.data.frame(coefficient_df)
  coefficients <- colnames(coefficient_df)
  coefficient_rows <- rownames(coefficient_df)
  coef_ordered <- NULL
  for (coef in coefficients) {
    new_order <- order(coefficient_df[[coef]], decreasing = TRUE)
    coef_ordered <- coefficient_df[new_order, ][[coef]]
    names(coef_ordered) <- coefficient_rows
    kept_rows <- NULL
    if (is.null(n)) {
      ## Then do it on a z-score
      tmp_summary <- summary(coef_ordered)
      tmp_mad <- stats::mad(as.numeric(coef_ordered, na.rm = TRUE))
      tmp_up_median_dist <- tmp_summary["Median"] + (tmp_mad * z)
      tmp_down_median_dist <- tmp_summary["Median"] - (tmp_mad * z)
      high_idx <- coef_ordered >= tmp_up_median_dist
      low_idx <- coef_ordered <= tmp_down_median_dist
      high_rows <- coef_ordered[high_idx]
      low_rows <- coef_ordered[low_idx]
      abundant_list[["high"]][[coef]] <- high_rows
      abundant_list[["low"]][[coef]] <- low_rows
    } else {
      ## Then do it in a number of rows
      abundant_list[["high"]][[coef]] <- head(coef_ordered, n = n)
      abundant_list[["low"]][[coef]] <- tail(coef_ordered, n = n)
    }
  }
  class(abundant_list) <- c("abundant_genes", "list")
  return(abundant_list)
}

#' A companion function for get_abundant_genes()
#'
#' Instead of pulling to top/bottom abundant genes, get all abundances and
#' variances or stderr.
#'
#' @param datum Output from _pairwise() functions.
#' @param type According to deseq/limma/ed ger/basic?
#' @param excel Print this to an excel file?
#' @return List containing the expression values and some metrics of
#'  variance/error.
#' @seealso [get_abundant_genes()]
#' @examples
#' \dontrun{
#'  abundance_excel <- get_pairwise_gene_abundances(combined, excel = "abundances.xlsx")
#'  ## This should provide a set of abundances after voom by condition.
#' }
#' @export
get_pairwise_gene_abundances <- function(datum, type = "limma", excel = NULL) {
  if (type == "limma") {
    ## Make certain we have a consistent set of column and row orders for the
    ## future operations
    conds <- names(datum[["limma"]][["identity_tables"]])
    row_order <- rownames(datum[["limma"]][["identity_tables"]][[1]])
    nconds <- length(conds)
    num_rows <- length(row_order)
    ## Pre-allocate and set the row/colnames
    expression_mtrx <- matrix(ncol = nconds, nrow = num_rows)
    stdev_mtrx <- matrix(ncol = nconds, nrow = num_rows)
    sigma_mtrx <- matrix(ncol = nconds, nrow = num_rows)
    s2post_mtrx <- matrix(ncol = nconds, nrow = num_rows)
    t_mtrx <- matrix(ncol = nconds, nrow = num_rows)
    colnames(expression_mtrx) <- conds
    colnames(stdev_mtrx) <- conds
    colnames(sigma_mtrx) <- conds
    colnames(s2post_mtrx) <- conds
    colnames(t_mtrx) <- conds
    rownames(expression_mtrx) <- row_order
    rownames(stdev_mtrx) <- row_order
    rownames(sigma_mtrx) <- row_order
    rownames(s2post_mtrx) <- row_order
    rownames(t_mtrx) <- row_order
    ## Recast these as data frames because they are easier syntactically.
    ## E.g. I like setting columns with df[[thingie]]
    ## I am explicitly setting the row order because it seems that the output
    ## from limma is not set from one table to the next.
    for (cond in conds) {
      expression_mtrx[, cond] <- datum[["limma"]][["identity_tables"]][[cond]][row_order, "logFC"]
      t_mtrx[, cond] <- datum[["limma"]][["identity_tables"]][[cond]][row_order, "t"]

      ## When attempting to extract something suitable for error bars:
      ## https://support.bioconductor.org/p/79349/
      ## The discussion seems to me to have settled on:
      ## std.error =  fit$stdev.unscaled * fit$sigma
      stdev_mtrx[, cond] <- as.numeric(
          datum[["limma"]][["identity_comparisons"]][["stdev.unscaled"]][row_order, cond])
      sigma_mtrx[, cond] <- as.numeric(
          datum[["limma"]][["identity_comparisons"]][["sigma"]])
      s2post_mtrx[, cond] <- as.numeric(
          datum[["limma"]][["identity_comparisons"]][["s2.post"]])
      std_error <- stdev_mtrx * sigma_mtrx
      another_error <- stdev_mtrx * s2post_mtrx
      ## another_error <- stdev_table * s2post_table
    }
  } ## End if type is limma

  retlist <- list(
      "expression_values" = expression_mtrx,
      "t_values" = t_mtrx,
      "error_values" = std_error,
      "another_error" = another_error,
      "stdev_values" = stdev_mtrx)
  if (!is.null(excel)) {
    annotations <- fData(datum[["input"]])
    expressions <- retlist[["expression_values"]]
    colnames(expressions) <- glue("expr_{colnames(expressions)}")
    errors <- retlist[["error_values"]]
    colnames(errors) <- glue("err_{colnames(errors)}")
    expression_table <- merge(annotations, expressions, by = "row.names")
    rownames(expression_table) <- expression_table[["Row.names"]]
    expression_table <- expression_table[, -1]
    expression_table <- merge(expression_table, errors, by = "row.names")
    rownames(expression_table) <- expression_table[["Row.names"]]
    expression_table <- expression_table[, -1]
    expression_written <- write_xlsx(
        data = expression_table,
        sheet = "expression_values",
        title = "Values comprising the logFCs and errors (expression / t-statistic)")
    write_result <- openxlsx::saveWorkbook(wb = expression_written[["workbook"]],
                                           file = excel, overwrite = TRUE)
  }
  return(retlist)
}

#' Make sure the outputs from limma and friends are in a format suitable for IHW.
#'
#' IHW seems like an excellent way to improve the confidence in the p-values
#' provided by the various DE methods.  It expects inputs fairly specific to
#' DESeq2, however, it is trivial to convert other methods to this, ergo this
#' function.
#'
#' https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html
#'
#' @param de_result Table which should have the 2 types of requisite columns:
#'  mean value of counts and p-value.
#' @param pvalue_column Name of the column of p-values.
#' @param type If specified, this will explicitly perform the calculation for
#'  the given type of differential expression analysis: limma, edger, deseq,
#'  etc.
#' @param mean_column Name of the column of mean values.
#' @param significance IHW uses this parameter, I don't know why.
#' @return weight adjusted p-values.
#' @seealso [IHW]
ihw_adjust <- function(de_result, pvalue_column = "pvalue", type = NULL,
                       mean_column = "baseMean", significance = 0.05) {
  ## We need to know the method used, because the values returned are not
  ## necessarily in the scale expected by IHW.
  if (is.null(type)) {
    if (grepl(pattern = "^limma", x = pvalue_column)) {
      type <- "limma"
    } else if (grepl(pattern = "^deseq", x = pvalue_column)) {
      type <- "deseq"
    } else if (grepl(pattern = "^edger", x = pvalue_column)) {
      type <- "edger"
    } else if (grepl(pattern = "^basic", x = pvalue_column)) {
      type <- "basic"
    } else if (grepl(pattern = "^ebseq", x = pvalue_column)) {
      type <- "ebseq"
    } else {
      stop("Unable to determine the type of pvalues in this data.")
    }
  }

  ## Now that we know the method used, coerce the result from method x into something
  ## which makes sense to IHW.
  tmp_table <- de_result
  if (type == "limma") {
    tmp_table[["base10_mean"]] <- 2 ^ tmp_table[[mean_column]]
    mean_column <- "base10_mean"
  } else if (type == "edger") {
    tmp_table[["base10_mean"]] <- 2 ^ tmp_table[[mean_column]]
    mean_column <- "base10_mean"
  } else if (type == "basic") {
    tmp_table[["base10_median"]] <- 2 ^ ((tmp_table[["basic_nummed"]] + tmp_table[["basic_denmed"]]) / 2)
    mean_column <- "base10_median"
  }

  ## Add a hopefully unneccessary check that everything is numeric (which is currently not true
  ## because I have invoked format() on some numbers, a task I subsequently
  ## removed since it was bad and dumb.
  tmp_table[[pvalue_column]] <- as.numeric(tmp_table[[pvalue_column]])
  tmp_table[[mean_column]] <- as.numeric(tmp_table[[mean_column]])

  ## Finally, invoke IHW and get its interpretation of adjusted p-values.
  formula <- as.formula(glue::glue("{pvalue_column} ~ {mean_column}"))
  ihw_result <- sm(IHW::ihw(formula, data = tmp_table, alpha = significance))
  adjusted_p_values <- IHW::adj_pvalues(ihw_result)
  return(adjusted_p_values)
}

#' Get a set of up/down differentially expressed genes.
#'
#' Take one or more criteria (fold change, rank order, (adj)p-value,
#' z-score from median FC) and use them to extract the set of genes
#' which are defined as 'differentially expressed.'  If no criteria
#' are provided, it arbitrarily chooses all genes outside of 1-z.
#'
#' Tested in test_29de_shared.R
#'
#' @param table Table from limma/edger/deseq.
#' @param n Rank-order top/bottom number of genes to take.
#' @param z Number of z-scores >/< the median to take.
#' @param lfc Fold-change cutoff.
#' @param p P-value cutoff.
#' @param fold Identifier reminding how to get the bottom portion of a
#'  fold-change (plusminus says to get the negative of the
#'  positive, otherwise 1/positive is taken).  This effectively
#'  tells me if this is a log fold change or not.
#' @param column Table's column used to distinguish top vs. bottom.
#' @param p_column Table's column containing (adjusted or not)p-values.
#' @return Subset of the up/down genes given the provided criteria.
#' @seealso [extract_significant_genes()] [get_abundant_genes()]
#' @examples
#' \dontrun{
#'  sig_table <- get_sig_genes(table, lfc = 1)
#' }
#' @export
get_sig_genes <- function(table, n = NULL, z = NULL, lfc = NULL, p = NULL,
                          column = "logFC", fold = "plusminus", p_column = "adj.P.Val") {
  if (is.null(z) & is.null(n) & is.null(lfc) & is.null(p)) {
    message("No n, z, p, nor lfc provided, setting p to 0.05 and lfc to 1.0.")
    p <- 0.05
    lfc <- 1.0
  }
  up_genes <- table
  down_genes <- table

  if (is.null(table[[column]])) {
    message("There is no ", column, " column in the table.")
    message("The columns are: ", toString(colnames(table)))
    stop("There is no ", column, " column in the table.")
  }

  if (!is.null(p)) {
    up_idx <- as.numeric(up_genes[[p_column]]) <= p
    ## Remember we have these reformatted as scientific
    up_genes <- up_genes[up_idx, ]
    down_idx <- as.numeric(down_genes[[p_column]]) <= p
    down_genes <- down_genes[down_idx, ]
    ## Going to add logic in case one does not ask for fold change
    ## In that case, a p-value assertion should still know the difference
    ## between up and down. But it should also still know the difference between
    ## ratio and log changes.
    if (fold == "plusminus" | fold == "log") {
      ## up_idx <- up_genes[, column] > 0.0
      up_idx <- as.numeric(up_genes[[column]]) > 0.0
      up_genes <- up_genes[up_idx, ]
      down_idx <- as.numeric(down_genes[[column]]) < 0.0
      down_genes <- down_genes[down_idx, ]
    } else {
      ## plusminus refers to a positive/negative number of logfold changes from
      ## a logFC(1) = 0
      up_idx <- as.numeric(up_genes[[column]]) >= 1.0
      up_genes <- up_genes[up_idx, ]
      down_idx <- as.numeric(down_genes[[column]]) <= -1.0
      down_genes <- down_genes[down_idx, ]
    }
  }

  if (!is.null(lfc)) {
    up_idx <- as.numeric(up_genes[[column]]) >= lfc
    up_genes <- up_genes[up_idx, ]
    if (fold == "plusminus" | fold == "log") {
      ## plusminus refers to a positive/negative number of logfold changes from
      ## a logFC(1) = 0
      down_idx <- as.numeric(down_genes[[column]]) <= (lfc * -1.0)
      down_genes <- down_genes[down_idx, ]
    } else {
      ## If it isn't log fold change, then values go from 0..x where 1 is unchanged
      down_idx <- as.numeric(down_genes[[column]]) <= (1.0 / lfc)
      down_genes <- down_genes[down_idx, ]
    }
  }

  if (!is.null(z)) {
    ## Take an arbitrary number which are >= and <= a value which is z zscores
    ## from the median.
    message("Getting the genes >= ", z, " z scores away from the median of all.")
    ## Use the entire table for the summary
    out_summary <- summary(as.numeric(table[[column]]))
    out_mad <- stats::mad(as.numeric(table[[column]]), na.rm = TRUE)
    up_median_dist <- out_summary["Median"] + (out_mad * z)
    down_median_dist <- out_summary["Median"] - (out_mad * z)
    ## But use the (potentially already trimmed) up/down tables for indexing
    up_idx <- as.numeric(up_genes[[column]]) >= up_median_dist
    up_genes <- up_genes[up_idx, ]
    down_idx <- as.numeric(down_genes[[column]]) <= down_median_dist
    down_genes <- down_genes[down_idx, ]
  }

  if (!is.null(n)) {
    ## Take a specific number of genes at the top/bottom of the rank ordered list.
    message("Getting the top and bottom ", n, " genes.")
    upranked <- up_genes[order(as.numeric(up_genes[[column]]), decreasing = TRUE), ]
    up_genes <- head(upranked, n = n)
    downranked <- down_genes[order(as.numeric(down_genes[[column]])), ]
    down_genes <- head(downranked, n = n)
  }
  up_genes <- up_genes[order(as.numeric(up_genes[[column]]), decreasing = TRUE), ]
  down_genes <- down_genes[order(as.numeric(down_genes[[column]]), decreasing = FALSE), ]
  ret <- list(
      "up_genes" = up_genes,
      "down_genes" = down_genes)
  return(ret)
}

#' Run makeContrasts() with all pairwise comparisons.
#'
#' In order to have uniformly consistent pairwise contrasts, I decided
#' to avoid potential human erors(sic) by having a function generate
#' all contrasts.
#'
#' Invoked by the _pairwise() functions.
#'
#' @param model Describe the conditions/batches/etc in the experiment.
#' @param conditions Factor of conditions in the experiment.
#' @param do_identities Include all the identity strings? Limma can
#'  use this information while edgeR can not.
#' @param do_extras Include extra contrasts?  This seems redundant with extra_contrasts
#'  below, but there is a reason for it.
#' @param do_pairwise Include all pairwise strings? This shouldn't
#'  need to be set to FALSE, but just in case.
#' @param extra_contrasts Optional string of extra contrasts to include.
#' @param ... Extra arguments passed here are caught by arglist.
#' @return List including the following information:
#' \enumerate{
#'  \item all_pairwise_contrasts = the result from makeContrasts(...)
#'  \item identities = the string identifying each condition alone
#'  \item all_pairwise = the string identifying each pairwise comparison alone
#'  \item contrast_string = the string passed to R to call makeContrasts(...)
#'  \item names = the names given to the identities/contrasts
#' }
#' @seealso [limma::makeContrasts()]
#' @examples
#' \dontrun{
#'  pretend <- make_pairwise_contrasts(model, conditions)
#' }
#' @export
make_pairwise_contrasts <- function(model, conditions, do_identities = FALSE,
                                    do_extras = TRUE, do_pairwise = TRUE, extra_contrasts = NULL, ...) {
  arglist <- list(...)
  tmpnames <- colnames(model)
  tmpnames <- gsub(pattern = "data[[:punct:]]", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "-", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "+", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "conditions^(\\d+)$", replacement = "c\\1", x = tmpnames)
  tmpnames <- gsub(pattern = "conditions", replacement = "", x = tmpnames)
  colnames(model) <- tmpnames
  conditions <- gsub(pattern = "^(\\d+)$", replacement = "c\\1", x = conditions)
  condition_table <- table(conditions)
  identities <- list()
  contrast_string <- ""
  eval_strings <- list()
  for (c in seq_along(condition_table)) {
    identity_name <- names(condition_table[c])
    identity_string <- paste(identity_name, " = ", identity_name, ",", sep = "")
    identities[identity_name] <- identity_string
  }
  ## If I also create a sample condition 'alice', and also perform a subtraction
  ## of 'alice' from 'bob', then the full makeContrasts() will be:
  ## makeContrasts(bob = bob, alice = alice, bob_vs_alice=(bob)-(alice), levels = design)
  ## The parentheses in this case are superfluous, but they remind me that I am
  ## finally doing some match, and also they remind me that we can do much more
  ## complex things like:
  ## ((bob-alice)-(jane-alice)) or whatever....
  all_pairwise <- list()
  identity_names <- names(identities)
  lenminus <- length(identities) - 1
  for (c in seq_len(lenminus)) {
    c_name <- names(identities[c])
    nextc <- c + 1
    for (d in seq(from = nextc, to = length(identities))) {
      d_name <- names(identities[d])
      minus_string <- paste(d_name, "_vs_", c_name, sep = "")
      exprs_string <- paste(minus_string, "=", d_name, "-", c_name, ",", sep = "")
      all_pairwise[minus_string] <- exprs_string
    }
  }
  ## At this point, I have strings which represent the definition of every
  ## sample condition as well as strings which represent every possible
  ## B-A where B comes somewhere after A in the model matrix.
  ## The goal now is to create the variables in the R environment
  ## and add them to makeContrasts()
  if (isTRUE(do_identities)) {
    eval_strings <- append(eval_strings, identities)
  }
  if (isTRUE(do_pairwise)) {
    eval_strings <- append(eval_strings, all_pairwise)
  }
  eval_names <- names(eval_strings)

  if (!is.null(extra_contrasts) & isTRUE(do_extras)) {
    extra_eval_strings <- strsplit(extra_contrasts, ",")[[1]]
    extra_eval_names <- extra_eval_strings
    extra_eval_names <- stringi::stri_replace_all_regex(extra_eval_strings,
                                                        "^(\\s*)(\\w+)\\s*=\\s*.*$", "$2")
    extra_eval_names <- gsub(pattern = "^\\s+", replacement = "", x = extra_eval_names, perl = TRUE)
    for (i in seq_along(extra_eval_strings)) {
      new_name <- extra_eval_names[[i]]
      extra_eval_string <- extra_eval_strings[[i]]
      extra_eval_string <- gsub(pattern = "^\\s+", replacement = "", x = extra_eval_string, perl = TRUE)
      extra_contrast <- glue("{extra_eval_string}, ")
      eval_strings <- append(eval_strings, extra_contrast)
      eval_names <- append(eval_names, new_name)
      all_pairwise[new_name] <- extra_contrast
    }
    names(eval_strings) <- eval_names
  }
  ## Add them to makeContrasts()
  contrast_string <- "all_pairwise_contrasts = mymakeContrasts("
  for (f in seq_along(eval_strings)) {
    eval_name = names(eval_strings[f])
    ## Get a little defensive to make sure I do not have contrasts which start with
    ## silly things like numbers of punctuation.
    if (grepl(x=eval_name, pattern="^([[:digit:]]|[[:punct:]])")) {
      stop("This function requires contrast names to start with a letter.")
    }
    eval_string <- eval_strings[f]

    contrast_string <- glue("{contrast_string} {eval_string}")
  }
  ## The final element of makeContrasts() is the design from voom()
  contrast_string <- glue("{contrast_string} levels = model)")
  eval(parse(text = contrast_string))
  ## I like to change the column names of the contrasts because by default
  ## they are kind of obnoxious and too long to type

  colnames(all_pairwise_contrasts) <- eval_names
  result <- list(
      "all_pairwise_contrasts" = all_pairwise_contrasts,
      "identities" = identities,
      "identity_names" = identity_names,
      "all_pairwise" = all_pairwise,
      "contrast_string" = contrast_string,
      "names" = eval_names)
  return(result)
}

#' A copy of limma::makeContrasts() with special sauce.
#'
#' This is a copy of limma::makeContrasts without the test of make.names()
#' Because I want to be able to use it with interaction models potentially
#' and if a model has first:second, make.names() turns the ':' to a '.'
#' and then the equivalence test fails, causing makeContrasts() to error
#' spuriously (I think).
#'
#' @param ... Conditions used to make the contrasts.
#' @param contrasts Actual contrast names.
#' @param levels contrast levels used.
#' @return Same contrasts as used in makeContrasts, but with unique names.
#' @seealso [limma::makeContrasts()]
mymakeContrasts <- function(..., contrasts = NULL, levels) {
  e <- substitute(list(...))
  if (is.factor(levels)) {
    levels <- levels(levels)
  }
  if (!is.character(levels)) {
    levels <- colnames(levels)
  }
  if (levels[1] == "(Intercept)") {
    levels[1] <- "Intercept"
    warning("Renaming (Intercept) to Intercept")
  }
  n <- length(levels)
  if (n < 1)
    stop("No levels to construct contrasts from")
  indicator <- function(i, n) {
    out <- rep(0, n)
    out[i] <- 1
    out
  }
  levelsenv <- new.env()
  for (i in seq_len(n)) {
    assign(levels[i], indicator(i, n), pos = levelsenv)
  }
  if (!is.null(contrasts)) {
    if (length(e) > 1)
      stop("Can't specify both ... and contrasts")
    e <- as.character(contrasts)
    ne <- length(e)
    cm <- matrix(0, n, ne, dimnames = list(Levels = levels,
                                           Contrasts = e))
    if (ne == 0)
      return(cm)
    for (j in seq_len(ne)) {
      ej <- parse(text = e[j])
      cm[, j] <- eval(ej, envir = levelsenv)
    }
    return(cm)
  }
  ne <- length(e)
  enames <- names(e)[2:ne]
  easchar <- as.character(e)[2:ne]
  if (is.null(enames))
    cn <- easchar
  else cn <- ifelse(enames == "", easchar, enames)
  cm <- matrix(0, n, ne - 1, dimnames = list(Levels = levels,
                                             Contrasts = cn))
  if (ne < 2)
    return(cm)
  for (j in seq(from = 1, to = ne - 1)) {
    ej <- e[[j + 1]]
    if (is.character(ej))
      ej <- parse(text = ej)
    ej <- eval(ej, envir = levelsenv)
    if (!is.numeric(ej)) {
      colnames(cm)[j] <- as.character(ej)[1]
      if (is.character(ej))
        ej <- parse(text = ej)
      ej <- eval(ej, envir = levelsenv)
    }
    cm[, j] <- ej
  }
  cm
}

#' Get rid of characters which will mess up contrast making and such before
#' playing with an expt.
#'
#' @param expt An expt object to clean.
sanitize_expt <- function(expt) {
  design <- pData(expt)
  conditions <- gsub(
      pattern = "^(\\d+)$", replacement = "c\\1", x = as.character(design[["condition"]]))
  batches <- gsub(
      pattern = "^(\\d+)$", replacement = "b\\1", x = as.character(design[["batch"]]))
  ## To be honest, there is absolutely no way I would have thought of this
  ## regular expression:
  ## https://stackoverflow.com/questions/30945993
  ## In theory I am pretty good with regexes, but this is devious to me!
  conditions <- gsub(pattern = "[^\\PP_]", replacement = "", x = conditions, perl = TRUE)
  batches <- gsub(pattern = "[^\\PP_]", replacement = "", x = batches, perl = TRUE)
  conditions <- gsub(pattern = "[[:blank:]]", replacement = "", x = conditions)
  batches <- gsub(pattern = "[[:blank:]]", replacement = "", x = batches)
  conditions <- gsub(pattern="[[:punct:]]", replacement = "", x = conditions)
  batches <- gsub(pattern="[[:punct:]]", replacement = "", x = batches)
  conditions <- as.factor(conditions)
  batches <- as.factor(batches)
  expressionset <- expt[["expressionset"]]
  Biobase::pData(expressionset)[["condition"]] <- conditions
  Biobase::pData(expressionset)[["batch"]] <- batches
  expt[["expressionset"]] <- expressionset
  expt[["conditions"]] <- conditions
  expt[["batches"]] <- batches
  return(expt)
}

#' Extract multicopy genes from up/down gene expression lists.
#'
#' The function semantic_copynumber_filter() is the inverse of this.
#'
#' Currently untested, used for Trypanosome analyses primarily, thus the default
#' strings.
#'
#' @param ... Arguments for semantic_copynumber_filter()
#' @export
semantic_copynumber_extract <- function(...) {
  ret <- semantic_copynumber_filter(..., invert = TRUE)
  return(ret)
}

#' Remove multicopy genes from up/down gene expression lists.
#'
#' In our parasite data, there are a few gene types which are
#' consistently obnoxious.  Multi-gene families primarily where the
#' coding sequences are divergent, but the UTRs nearly identical.  For
#' these genes, our sequence based removal methods fail and so this
#' just excludes them by name.
#'
#' Currently untested, used for Trypanosome analyses primarily, thus the default
#' strings.
#'
#' @param input List of sets of genes deemed significantly
#'  up/down with a column expressing approximate count numbers.
#' @param max_copies Keep only those genes with <= n putative
#'  copies.
#' @param use_files Use a set of sequence alignments to define the copy numbers?
#' @param invert Keep these genes rather than drop them?
#' @param semantic Set of strings with gene names to exclude.
#' @param semantic_column Column in the DE table used to find the
#'  semantic strings for removal.
#' @return Smaller list of up/down genes.
#' @seealso [semantic_copynumber_extract()]
#' @examples
#' \dontrun{
#'  pruned <- semantic_copynumber_filter(table, semantic = c("ribosomal"))
#'  ## Get rid of all genes with 'ribosomal' in the annotations.
#' }
#' @export
semantic_copynumber_filter <- function(input, max_copies = 2, use_files = FALSE, invert = TRUE,
                                       semantic = c("mucin", "sialidase", "RHS",
                                                    "MASP", "DGF", "GP63"),
                                       semantic_column = "product") {
  if ("expt" %in% class(input)) {
    result <- semantic_expt_filter(input, invert = invert, semantic = semantic,
                                   semantic_column = semantic_column)
    return(result)
  }

  table_type <- "significance"
  if (!is.null(input[["data"]])) {
    table_type <- "combined"
  }
  type <- "Kept"

  table_list <- NULL
  if (table_type == "combined") {
    table_list <- input[["data"]]
  } else {
    ## The set of significance tables will be 2x the number of contrasts
    ## Therefore, when we get to > 1x the number of contrasts, all the tables
    ## will be 'down'
    table_list <- c(input[["ups"]], input[["downs"]])
    upnames <- glue("up_{names(input[['ups']])}")
    downnames <- glue("down_{names(input[['downs']])}")
    names(table_list) <- c(upnames, downnames)
    up_to_down <- length(input[["ups"]])
  }

  numbers_removed <- list()
  for (count in seq_along(table_list)) {
    tab <- table_list[[count]]
    table_name <- names(table_list)[[count]]
    numbers_removed[[table_name]] <- list()
    message("Working on ", table_name)
    if (isTRUE(use_files)) {
      file <- ""
      if (table_type == "combined") {
        file <- file.path("singletons", "gene_counts",
                          glue("{table_name}.fasta.out.count"))
      } else {
        file <- file.path("singletons", "gene_counts",
                          glue("up_{table_name}.fasta.out.count"))
      }
      tmpdf <- try(read.table(file), silent = TRUE)
      if (class(tmpdf)[1] == "data.frame") {
        colnames(tmpdf) <- c("ID", "members")
        tab <- merge(tab, tmpdf, by.x = "row.names", by.y = "ID")
        rownames(tab) <- tab[["Row.names"]]
        tab <- tab[, -1, drop = FALSE]
        tab <- tab[count <= max_copies, ]
      }
    }  ## End using empirically defined groups of multi-gene families.

    ## Now remove genes by name.
    kept_list <- new_table <- NULL
    for (string in semantic) {
      pre_remove_size <- nrow(tab)
      idx <- NULL
      if (semantic_column == "rownames") {
        idx <- grepl(pattern = string, x = rownames(tab))
      } else {
        idx <- grepl(pattern = string, x = tab[, semantic_column])
      }
      type <- "Removed"
      if (!isTRUE(invert)) {
        idx <- !idx
        type <- "Kept"
      }
      num_removed <- sum(idx)
      numbers_removed[[table_name]][[string]] <- num_removed
      if (num_removed > 0) {
        tab <- tab[idx, ]
        table_list[[count]] <- tab
        message("Table started with: ", pre_remove_size, ". ", type,
                " entries with string ", string,
                ", found ", num_removed, "; table has ",
                nrow(tab),  " rows left.")
      } else {
        message("Found no entries of type ", string, ".")
      }
    } ## End of the foreach semantic thing to remove

    ## Now recreate the original table lists as either de tables or significance.
    if (table_type == "combined") {
      for (count in seq_along(table_list)) {
        input[["data"]][[count]] <- table_list[[count]]
      }
    } else {
      ## Then it is a set of significance tables.
      if (count <= up_to_down) {
        table_name <- names(table_list)[count]
        if (grep(pattern = "^up_", table_name)) {
          newname <- gsub(pattern = "^up_", replacement = "", x = table_name)
          input[["ups"]][[newname]] <- table_list[[count]]
        } else {
          newname <- gsub(pattern = "^down_", replacement = "", x = table_name)
          input[["downs"]][[newname]] <- table_list[[count]]
        }
      }
    }
  }
  ## Now the tables should be reconstructed.
  input[["numbers_removed"]] <- numbers_removed
  input[["type"]] <- type
  return(input)
}

## EOF
