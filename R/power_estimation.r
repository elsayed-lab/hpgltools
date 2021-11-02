#' Invoke PROPER and replace its default data set with data of interest.
#'
#' Recent reviewers of Najib's grants have taken an increased interest in
#' knowing the statistical power of the various experiments.  He queried
#' Dr. Corrada-Bravo who suggested PROPER.  I spent some time looking through it
#' and, with some revervations, modified its workflow to (at least in theory) be
#' able to examine any dataset.  The workflow in question is particularly odd
#' and warrants further discussion/analysis.  This function invokes PROPER
#' exactly as it was performed in their paper.
#'
#' @param de_tables A set of differential expression results, presumably from
#'  EdgeR or DESeq2.
#' @param p Cutoff
#' @param experiment The default data set in PROPER is entitled 'cheung'.
#' @param nsims Number of simulations to perform.
#' @param reps Simulate these number of experimental replicates.
#' @param de_method There are a couple choices here for tools which are pretty
#'  old, my version of this only accepts deseq or edger.
#' @param alpha_type I assume p-adjust type.
#' @param alpha Accepted fdr rate.
#' @param stratify There are a few options here, I don't fully understand them.
#' @param target Cutoff.
#' @param filter Apply a filter?
#' @param delta Not epsilon! (E.g. I forget what this does).
#' @return List containing the various results and plots from proper.
#' @seealso [PROPER]
default_proper <- function(de_tables, p = 0.05, experiment = "cheung", nsims = 20,
                           reps = c(3, 5, 7, 10), de_method = "edger", alpha_type = "fdr",
                           alpha = 0.1, stratify = "expr", target = "lfc",
                           filter = "none", delta = 0.5) {
  DEmethod <- "edgeR"
  if (de_method == "edger") {
    DEmethod <- "edgeR"
  } else if (de_method == "deseq") {
    DEmethod <- "DESeq2"
  } else {
    stop("This accepts only 'deseq' or 'edger'.")
  }
  genes <- nrow(de_tables[["data"]][[1]])
  attached <- attachNamespace("PROPER")
  ds <- c(experiment)
  loaded <- data(list = ds, package = "PROPER")
  simulation_options <- PROPER::RNAseq.SimOptions.2grp(
                                    ngenes = genes,
                                    p.DE = p,
                                    lOD = experiment,
                                    lBaselineExpr = experiment)
  simulation_result <- PROPER::runSims(Nreps = reps,
                                       sim.opts = simulation_options,
                                       DEmethod = DEmethod,
                                       nsims = nsims)
  ## alpha_types : fdr, pval
  ## stratify types : expr, dispersion
  ## filter types : none, expr
  ## target types : lfc, effectsize
  powers <- PROPER::comparePower(simulation_result,
                                 alpha.type = alpha_type,
                                 alpha.nominal = alpha,
                                 stratify.by = stratify,
                                 filter.by = filter,
                                 target.by = target,
                                 delta = delta)

  tmp_file <- tempfile(pattern = "power", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  PROPER::plotPower(powers)
  power_plot <- grDevices::recordPlot()
  dev.off()
  file.remove(tmp_file)

  tmp_file <- tempfile(pattern = "power", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  PROPER::plotPowerTD(powers)
  powertd_plot <- grDevices::recordPlot()
  dev.off()
  file.remove(tmp_file)

  tmp_file <- tempfile(pattern = "power", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  PROPER::plotPowerFD(powers)
  powerfd_plot <- grDevices::recordPlot()
  dev.off()
  file.remove(tmp_file)

  tmp_file <- tempfile(pattern = "power", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  PROPER::plotFDcost(powers)
  fdcost_plot <- grDevices::recordPlot()
  dev.off()
  file.remove(tmp_file)

  tmp_file <- tempfile(pattern = "power", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  PROPER::plotPowerHist(powerOutput = powers, simResult = simulation_result)
  powerhist_plot <- grDevices::recordPlot()
  dev.off()
  file.remove(tmp_file)

  tmp_file <- tempfile(pattern = "power", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  PROPER::plotPowerAlpha(powers)
  poweralpha_plot <- grDevices::recordPlot()
  dev.off()
  file.remove(tmp_file)

  power_table <- PROPER::power.seqDepth(simResult = simulation_result, powerOutput = powers)
  retlist <- list(
      "options" = simulation_options,
      "simulation" = simulation_result,
      "powers" = powers,
      "power_plot" = power_plot,
      "powertd_plot" = powertd_plot,
      "powerfd_plot" = powerfd_plot,
      "fdcost_plot" = fdcost_plot,
      "powerhist_plot" = powerhist_plot,
      "poweralpha_plot" = poweralpha_plot,
      "power_table" = power_table)
  return(retlist)
}

#' Invoke PROPER and replace its default data set with data of interest.
#'
#' Recent reviewers of Najib's grants have taken an increased interest in
#' knowing the statistical power of the various experiments.  He queried
#' Dr. Corrada-Bravo who suggested PROPER.  I spent some time looking through it
#' and, with some revervations, modified its workflow to (at least in theory) be
#' able to examine any dataset.  The workflow in question is particularly odd
#' and warrants further discussion/analysis.  This function is a modified
#' version of 'default_proper()' above and invokes PROPER after re-formatting a
#' given dataset in the way expected by PROPER.
#'
#' @param de_tables A set of differential expression results, presumably from
#'  EdgeR or DESeq2.
#' @param p Cutoff
#' @param experiment The default data set in PROPER is entitled 'cheung'.
#' @param nsims Number of simulations to perform.
#' @param reps Simulate these number of experimental replicates.
#' @param de_method There are a couple choices here for tools which are pretty
#'  old, my version of this only accepts deseq or edger.
#' @param alpha_type I assume p-adjust type.
#' @param alpha Accepted fdr rate.
#' @param stratify There are a few options here, I don't fully understand them.
#' @param target Cutoff.
#' @param mean_or_median Use mean or median values?
#' @param filter Apply a filter?
#' @param delta Not epsilon! (E.g. I forget what this does).
#' @return List containin the various tables and plots returned by PROPER.
#' @seealso [PROPER]
#' @export
simple_proper <- function(de_tables, p = 0.05, experiment = "cheung", nsims = 20,
                          reps = c(3, 5, 7, 10), de_method = "edger", alpha_type = "fdr",
                          alpha = 0.1, stratify = "expr", target = "lfc", mean_or_median = "mean",
                          filter = "none", delta = 0.5) {
  DEmethod <- "edgeR"
  if (de_method == "edger") {
    DEmethod <- "edgeR"
  } else if (de_method == "deseq") {
    DEmethod <- "DESeq2"
  } else {
    stop("This accepts only 'edger' or 'deseq'.")
  }
  exprs_mtrx <- exprs(de_tables[["input"]][["input"]])
  genes <- nrow(de_tables[["data"]][[1]])
  contrasts <- names(de_tables[["input"]][["edger"]][["all_tables"]])
  result_list <- list()
  ## For the moment this will hard-assume edger, but should be trivially changed for the others.
  count <- 0
  for (c in contrasts) {
    datum <- de_tables[["input"]][["edger"]]
    count <- count + 1
    samples <- datum[["contrast_list"]][[c]]
    model <- datum[["model"]]
    used_levels_idx <- samples != 0
    used_levels <- rownames(samples)[used_levels_idx]
    all_samples <- rownames(model)
    used_samples_idx <- rowSums(model[, used_levels]) > 0
    used_samples <- all_samples[used_samples_idx]
    used_mtrx <- exprs_mtrx[, used_samples]
    sample_coverages <- NULL
    all_coverage <- NULL
    if (mean_or_median == "mean") {
      sample_coverages <- colMeans(used_mtrx)
      all_coverage <- mean(sample_coverages)
    } else {
      sample_coverages <- miscTools::colMedians(used_mtrx)
      all_coverage <- median(sample_coverages)
    }
    x_intercept <- 0
    ## The points on the x-axis are: 0, 10, 20, 40, 80, 160, 320, 640, 1280, inf
    cutoffs <- c(10, 20, 40, 80, 160, 320, 640, 1280, 2560)
    clen <- length(cutoffs)
    for (d in 1:clen) {
      cutoff <- cutoffs[d]
      adder <- d - 1
      if (all_coverage < cutoff) {
        x_intercept <- adder + ((all_coverage - cutoffs[adder]) / cutoff)
      }
    }

    message("Working on contrast ", count, ": ", c, ".")
    simulation_options <- list(
        "ngenes" = genes,
        "p.DE" = p,
        "lBaselineExpr" = datum[["lrt"]][[c]][["AveLogCPM"]],
        "lOD" = log2(datum[["lrt"]][[c]][["dispersion"]]),
        "lfc" = datum[["all_tables"]][[c]][["logFC"]],
        "sim.seed" = 1,
        "design" = "2grp")
    simulation_result <- my_runsims(Nreps = reps,
                                    sim.opts = simulation_options,
                                    DEmethod = DEmethod,
                                    nsims = nsims)
    powers <- PROPER::comparePower(simulation_result,
                                   alpha.type = alpha_type,
                                   alpha.nominal = alpha,
                                   stratify.by = stratify,
                                   filter.by = filter,
                                   target.by = target,
                                   delta = delta)

    tmp_file <- tempfile(pattern = "power", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    PROPER::plotPower(powers)
    abline(v = x_intercept)
    power_plot <- grDevices::recordPlot()
    dev.off()
    file.remove(tmp_file)

    tmp_file <- tempfile(pattern = "power", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    PROPER::plotPowerTD(powers)
    abline(v = x_intercept)
    powertd_plot <- grDevices::recordPlot()
    dev.off()
    file.remove(tmp_file)

    tmp_file <- tempfile(pattern = "power", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    PROPER::plotPowerFD(powers)
    abline(v = x_intercept)
    powerfd_plot <- grDevices::recordPlot()
    dev.off()
    file.remove(tmp_file)

    tmp_file <- tempfile(pattern = "power", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    PROPER::plotFDcost(powers)
    abline(v = x_intercept)
    fdcost_plot <- grDevices::recordPlot()
    dev.off()
    file.remove(tmp_file)

    tmp_file <- tempfile(pattern = "power", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    PROPER::plotPowerHist(powerOutput = powers, simResult = simulation_result)
    abline(v = x_intercept)
    powerhist_plot <- grDevices::recordPlot()
    dev.off()
    file.remove(tmp_file)

    tmp_file <- tempfile(pattern = "power", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    PROPER::plotPowerAlpha(powers)
    abline(v = x_intercept)
    poweralpha_plot <- grDevices::recordPlot()
    dev.off()
    file.remove(tmp_file)

    power_table <- PROPER::power.seqDepth(simResult = simulation_result, powerOutput = powers)
    retlist <- list(
        "options" = simulation_options,
        "simulation" = simulation_result,
        "powers" = powers,
        "power_plot" = power_plot,
        "powertd_plot" = powertd_plot,
        "powerfd_plot" = powerfd_plot,
        "fdcost_plot" = fdcost_plot,
        "powerhist_plot" = powerhist_plot,
        "poweralpha_plot" = poweralpha_plot,
        "power_table" = power_table)
    result_list[[c]] <- retlist
  }
  return(result_list)
}

## Use my cheater infix ::: to handle some unexported stuff from PROPER.
update.RNAseq.SimOptions.2grp <- "PROPER" %:::% "update.RNAseq.SimOptions.2grp"
run.edgeR <- "PROPER" %:::% "run.edgeR"
run.DSS <- "PROPER" %:::% "run.DSS"
run.DESeq2 <- "PROPER" %:::% "run.DESeq2"

#' A version of PROPER:::runsims which is (hopefully) a little more robust.
#'
#' When I was testing PROPER, it fell down mysteriously on a few occasions.  The
#' source ended up being in runsims(), ergo this function.
#'
#' @param Nreps Vector of numbers of replicates to simulate.
#' @param Nreps2 Second vector of replicates.
#' @param nsims How many simulations to perform?
#' @param sim.opts Options provided in a list which include information about the expression,
#'  numbers of genes, logFC values, etc.
#' @param DEmethod I suggest using only either edgeR or DESeq2.
#' @param verbose Print some information along the way?
#' @seealso [PROPER]
my_runsims <- function (Nreps = c(3, 5, 7, 10), Nreps2, nsims = 100, sim.opts,
                        DEmethod = c("edgeR", "DSS", "DESeq", "DESeq2"), verbose = TRUE) {
  DEmethod <- match.arg(DEmethod)
  if (missing(Nreps2)) {
    Nreps2 <- Nreps
  } else {
    if (length(Nreps2) != length(Nreps)) {
      stop("Nreps and Nreps2 must be vectors of the same length")
    }
  }
  n1 <- max(Nreps)
  n2 <- max(Nreps2)
  old_seed <- set.seed(sim.opts[["sim.seed"]])
  pvalue <- fdrs <- xbar <- array(NA,
                                  dim = c(sim.opts[["ngenes"]],
                                        length(Nreps), nsims))
  DEids <- lfcs <- NULL
  for (i in 1:nsims) {
    sim.opts[["sim.seed"]] <- sim.opts[["sim.seed"]] + 1
    sim.opts <- update.RNAseq.SimOptions.2grp(sim.opts)
    dat.sim.big <- PROPER::simRNAseq(sim.opts, n1, n2)
    DEids[[i]] <- dat.sim.big[["DEid"]]
    lfcs[[i]] <- dat.sim.big[["simOptions"]][["lfc"]]
    for (j in 1:length(Nreps)) {
      nn1 <- Nreps[j]
      nn2 <- Nreps2[j]
      idx <- c(1:nn1, n1 + (1:nn2))
      this.design <- dat.sim.big[["designs"]][idx]
      this.X <- dat.sim.big[["counts"]][, idx]
      this.simOpts <- sim.opts
      ss <- rowSums(this.X)
      ix_na <- is.na(ss)
      if (sum(ix_na) > 0) {
        ss[ix_na] <- 0
      }
      ix.valid <- ss > 0
      this.X.valid <- this.X[ix.valid, , drop = FALSE]
      data0 <- list(counts = this.X.valid, designs = this.design)
      if (DEmethod == "edgeR") {
        res1 <- run.edgeR(data0)
      }
      if (DEmethod == "DSS") {
        res1 <- run.DSS(data0)
      }
      if (DEmethod == "DESeq2") {
        res1 <- run.DESeq2(data0)
      }
      pval <- fdr <- rep(1, nrow(this.X))
      X.bar1 <- rep(0, nrow(this.X))
      pval[ix.valid] <- res1[, "pval"]
      fdr[ix.valid] <- res1[, "fdr"]
      sizeF <- colSums(data0[["counts"]])
      sizeF <- sizeF / median(sizeF)
      X.bar1[ix.valid] <- rowMeans(sweep(data0[["counts"]], 2,
                                         sizeF, FUN = "/"))
      pvalue[, j, i] <- pval
      fdrs[, j, i] <- fdr
      xbar[, j, i] <- X.bar1
    }
  }
  retlist <- list(
      "pvalue" = pvalue,
      "fdrs" = fdrs,
      "xbar" = xbar,
      "DEid" = DEids,
      "lfcs" = lfcs,
      "Nreps1" = Nreps,
      "Nreps2" = Nreps2,
      "sim.opts" = sim.opts)
  return(retlist)
}
