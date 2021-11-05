## normalize_batch.r: Bring together various surrogate/batch
## estimation/correction methods.  I want to be able to compare/contrast
## surrogate estimators and batch correction methods.  These functions attempt
## to make that possible.

#' Combine all surrogate estimators and batch correctors into one function.
#'
#' For a long time, I have mostly kept my surrogate estimators and batch correctors
#' separate.  However, that separation was not complete, and it really did not make
#' sense.  This function brings them together. This now contains all the logic from
#' the freshly deprecated get_model_adjust().
#'
#' This applies the methodologies very nicely explained by Jeff Leek at
#' https://github.com/jtleek/svaseq/blob/master/recount.Rmd
#' and attempts to use them to acquire estimates which may be applied to an
#' experimental model by either EdgeR, DESeq2, or limma.  In addition, it
#' modifies the count tables using these estimates so that one may play with the
#' modified counts and view the changes (with PCA or heatmaps or whatever).
#' Finally, it prints a couple of the plots shown by Leek in his document. In
#' other words, this is entirely derivative of someone much smarter than me.
#'
#' @param input Dataframe or expt or whatever as the data to analyze/modify.
#' @param design If the data is not an expt, then put the design here.
#' @param estimate_type Name of the estimator.
#' @param batch1 Column in the experimental design for the first known batch.
#' @param batch2 Only used by the limma method, a second batch column.
#' @param surrogates Either a number of surrogates or a method to
#'  search for them.
#' @param low_to_zero Move elements which are <0 to 0?
#' @param cpus Use parallel and split intensive operations?
#' @param na_to_zero Set any NA entries to 0?
#' @param expt_state If this is not an expt, provide the state of the data here.
#' @param confounders List of confounded factors for smartSVA/iSVA.
#' @param chosen_surrogates Somewhat redundant with surrogates above,
#'  but provides a second place to enter because of the way I use
#'  ... in normalize_expt().
#' @param adjust_method Choose the method for applying the estimates
#'  to the data.
#' @param filter Filter the data?
#' @param thresh If filtering, use this threshold.
#' @param noscale If using combat, scale the data?
#' @param prior_plots Plot the priors?
#' @return List containing surrogate estimates, new counts, the models, and
#'   some plots, as available.
#' @seealso [all_adjuster()] [isva] [sva] [limma::removeBatchEffect()]
#'  [corpcor] [edgeR] [RUVSeq] [SmartSVA] [variancePartition] [counts_from_surrogates()]
#' @export
all_adjusters <- function(input, design = NULL, estimate_type = "sva", batch1 = "batch",
                          batch2 = NULL, surrogates = "be", low_to_zero = FALSE, cpus = 4,
                          na_to_zero = TRUE, expt_state = NULL, confounders = NULL,
                          chosen_surrogates = NULL, adjust_method = "ruv",
                          filter = "raw", thresh = 1, noscale = FALSE, prior_plots = FALSE) {
  my_design <- NULL
  my_data <- NULL

  ## Gather all the likely pieces we can use
  ## Without the following requireNamespace(ruv)
  ## we get an error 'unable to find an inherited method for function RUVr'
  ## ruv_loaded <- try(require(package = "ruv", quietly = TRUE))
  lib_result <- sm(requireNamespace("ruv"))
  att_result <- sm(try(attachNamespace("ruv"), silent = TRUE))
  ## In one test, this seems to have been enough, but in another, perhaps not.

  if ("expt" %in% class(input)) {
    ## Gather all the likely pieces we can use
    my_design <- input[["design"]]
    my_data <- exprs(input)
    expt_state <- input[["state"]]
  } else {  ## This is not an expt
    if (is.null(design)) {
      stop("If an expt is not passed, then design _must_ be.")
    }
    my_design <- design
    my_data <- input
    if (is.null(expt_state)) {
      message("Not able to discern the state of the data.")
      message("Going to use a simplistic metric to guess if it is log scale.")
      if (max(input) > 100) {
        expt_state[["transform"]] <- "raw"
      } else {
        expt_state[["transform"]] <- "log2"
      }
    }
  } ## Ending the tests of the input and its state.

  ## Once again this is a place where my new stance vis a vis NAs is relevant.
  if (isTRUE(na_to_zero)) {
    na_idx <- is.na(my_data)
    my_data[na_idx] <- 0
  }

  ## Different tools expect different inputs
  linear_mtrx <- NULL
  ## So here I will be creating separate matrices on the linear and log scales.
  log2_mtrx <- NULL
  ## I also want to record this along with the type of output returned by each method.
  input_scale <- NULL
  output_scale <- "linear"
  if (expt_state[["transform"]] == "linear" || expt_state[["transform"]] == "raw") {
    linear_mtrx <- as.matrix(my_data)
    log2_mtrx <- as.matrix(log2(linear_mtrx + 1))
    input_scale <- "linear"
  } else if (expt_state[["transform"]] == "log2") {
    log2_mtrx <- as.matrix(my_data)
    linear_mtrx <- as.matrix((2 ^ my_data) - 1)
    input_scale <- "log2"
  } else if (expt_state[["transform"]] == "log") {
    warning("Unexpected call for base e data.")
    log2_mtrx <- as.matrix(my_data / log(2))
    linear_mtrx <- as.matrix(exp(my_data) - 1)
    input_scale <- "loge"
  }

  zero_rows <- rowSums(linear_mtrx, na.rm = TRUE) == 0
  num_zero_rows <- sum(zero_rows, na.rm = TRUE)
  if (num_zero_rows > 0) {
    warning("batch_counts: Before batch/surrogate estimation, ",
            zero_rows, " rows are 0.")
  }

  elements <- nrow(linear_mtrx) * ncol(linear_mtrx)
  num_normal <- sum(linear_mtrx > 1, na.rm = TRUE)
  normal_pct <- scales::percent(num_normal / elements)
  mesg("batch_counts: Before batch/surrogate estimation, ",
       num_normal, " entries are x>1: ", normal_pct, ".")
  num_zero <- sum(linear_mtrx == 0, na.rm = TRUE)
  zero_pct <- scales::percent(num_zero / elements)
  message("batch_counts: Before batch/surrogate estimation, ",
          num_zero, " entries are x==0: ", zero_pct, ".")
  num_low <- sum(linear_mtrx < 1 & linear_mtrx > 0, na.rm = TRUE)
  low_pct <- scales::percent(num_low / elements)
  if (is.null(num_low)) {
    num_low <- 0
  }
  if (num_low > 0) {
    message("batch_counts: Before batch/surrogate estimation, ",
            num_low, " entries are 0<x<1: ", low_pct, ".")
  }
  num_neg <- sum(linear_mtrx < 0, na.rm = TRUE)
  if (num_neg > 0) {
    neg_pct <- scales::percent(num_neg / elements)
    message("batch_counts: Before batch/surrogate estimation, ",
            num_zero, " entries are x<0: ", neg_pct, ".")
  }

  conditions <- droplevels(as.factor(my_design[["condition"]]))
  batches <- droplevels(as.factor(my_design[[batch1]]))
  conditional_model <- model.matrix(~ conditions, data = my_design)
  sample_names <- colnames(input)
  null_model <- conditional_model[, 1]

  if (is.null(chosen_surrogates) & is.null(surrogates)) {
    chosen_surrogates <- 1
  } else if (is.null(surrogates)) {
    surrogates <- chosen_surrogates
  }

  chosen_surrogates <- 1
  if (is.null(surrogates)) {
    message("No estimate nor method to find surrogates was provided. ",
            "Assuming you want 1 surrogate variable.")
  } else {
    if (class(surrogates) == "character") {
      ## num.sv assumes the log scale.
      if (surrogates == "smartsva") {
        lm_rslt <- lm(t(linear_mtrx) ~ condition, data = my_design)
        sv_estimate_data <- t(resid(lm_rslt))
        chosen_surrogates <- isva::EstDimRMT(sv_estimate_data, FALSE)[["dim"]] + 1
      } else if (surrogates == "isva") {
        chosen_surrogates <- isva::EstDimRMT(log2_mtrx)
      } else if (surrogates != "be" & surrogates != "leek") {
        message("A string was provided, but it was neither 'be' nor 'leek', assuming 'be'.")
        chosen_surrogates <- sm(sva::num.sv(dat = log2_mtrx, mod = conditional_model))
      } else {
        chosen_surrogates <- sm(sva::num.sv(dat = log2_mtrx,
                                            mod = conditional_model, method = surrogates))
      }
      vword <- "variable"
      if (chosen_surrogates > 1) {
        vword <- "variables"
      }
      mesg("The ", surrogates, " method chose ",
           chosen_surrogates, " surrogate ", vword, ".")
    } else if (class(surrogates) == "numeric") {
      mesg("A specific number of surrogate variables was chosen: ", surrogates, ".")
      chosen_surrogates <- surrogates
    }
    if (chosen_surrogates < 1) {
      warning("One must have greater than 0 surrogates, setting chosen_surrogates to 1.")
      chosen_surrogates <- 1
    }
  }

  if (is.null(cpus)) {
    cpus <- parallel::detectCores() - 2
  }
  if (cpus < 0) {
    cpus <- 1
  }

  if (isTRUE(prior_plots)) {
    message("When using ComBat, using prior.plots may result in an error due to infinite ylim.")
  }

  ## empirical controls can take either log or base 10 scale depending on 'control_type'
  control_type <- "norm"
  control_likelihoods <- NULL
  if (control_type == "norm") {
    control_likelihoods <- try(sm(sva::empirical.controls(dat = log2_mtrx,
                                                          mod = conditional_model,
                                                          mod0=null_model,
                                                          n.sv = chosen_surrogates,
                                                          type = control_type)))
  } else {
    control_likelihoods <- try(sm(sva::empirical.controls(dat = linear_mtrx,
                                                          mod = conditional_model,
                                                          mod0=null_model,
                                                          n.sv = chosen_surrogates,
                                                          type = control_type)))
  }
  if (class(control_likelihoods) == "try-error") {
    mesg("The most likely error in sva::empirical.controls() ",
         "is a call to density in irwsva.build. ",
         "Setting control_likelihoods to zero and using unsupervised sva.")
    warning("It is highly likely that the underlying reason for this ",
            "error is too many 0's in the dataset, ",
            "please try doing a filtering of the data and retry.")
    control_likelihoods <- 0
  }
  if (sum(control_likelihoods) == 0) {
    if (estimate_type == "sva_supervised") {
      message("Unable to perform supervised estimations, changing to unsupervised_sva.")
      estimate_type <- "sva_unsupervised"
    } else if (estimate_type == "ruv_supervised") {
      message("Unable to perform supervised estimations, changing to empirical_ruv.")
      estimate_type <- "ruv_empirical"
    }
  }

  if (is.null(estimate_type)) {
    estimate_type <- "limma"
  }
  ## I use 'sva' as shorthand fairly often
  if (estimate_type == "sva") {
    estimate_type <- "sva_unsupervised"
    mesg("Estimate type 'sva' is shorthand for 'sva_unsupervised'.")
    mesg("Other sva options include: sva_supervised and svaseq.")
  }
  if (estimate_type == "ruv") {
    estimate_type <- "ruv_empirical"
    mesg("Estimate type 'ruv' is shorthand for 'ruv_empirical'.")
    mesg("Other ruv options include: ruv_residual and ruv_supervised.")
  }

  surrogate_result <- NULL
  model_adjust <- NULL
  adjusted_counts <- NULL
  type_color <- NULL
  source_counts <- NULL
  new_counts <- NULL
  matrx_scale <- "linear"
  sword <- "surrogate"
  if (chosen_surrogates > 1) {
    sword <- "surrogates"
  }
  ## Just an aside, calling this a base 10 matrix is stupid.  Just because something is
  ## put on a log scale does not suddently make it octal or binary or imaginary!
  surrogate_input <- linear_mtrx
  switchret <- switch(
      estimate_type,
      "combat" = {
        ## This peculiar syntax should match combat and combat_noscale
        ## to the same result.
        mesg("batch_counts: Using combat with a prior, no scaling, and a null model.")
        new_counts <- sm(sva::ComBat(linear_mtrx, batches, mod = NULL,
                                     par.prior = TRUE, prior.plots = prior_plots,
                                     mean.only = TRUE))
      },
      "combat_noprior" = {
        mesg("batch_counts: Using combat without a prior and no scaling.")
        mesg("This takes a long time!")
        new_counts <- sm(sva::ComBat(linear_mtrx, batches, mod = conditions,
                                     par.prior = FALSE, prior.plots = prior_plots,
                                     mean.only = TRUE))
      },
      "combat_noprior_scale" = {
        mesg("batch_counts: Using combat without a prior and with scaling.")
        new_counts <- sm(sva::ComBat(linear_mtrx, batches, mod = conditions,
                                     par.prior = FALSE, prior.plots = prior_plots,
                                     mean.only = FALSE))
      },
      "combat_notnull" = {
        mesg("batch_counts: Using combat with a prior, no scaling, and a conditional model.")
        new_counts <- sm(sva::ComBat(linear_mtrx, batches, mod = conditions,
                                     par.prior = TRUE, prior.plots = prior_plots,
                                     mean.only = TRUE))
      },
      "combat_scale" = {
        mesg("batch_counts: Using combat with a prior and with scaling.")
        new_counts <- sm(sva::ComBat(linear_mtrx, batches, mod = conditions,
                                     par.prior = TRUE, prior.plots = prior_plots,
                                     mean.only = FALSE))
      },
      "combatmod" = {
        mesg("batch_counts: Using a modified cbcbSEQ combatMod for batch correction.")
        new_counts <- cbcb_combat(dat = linear_mtrx, batch = batches,
                                  mod = conditions, noscale = noscale)
      },
      "fsva" = {
        ## Ok, I have a question:
        ## If we perform fsva using log2(data) and get back SVs on a scale of ~ -1
        ## to 1, then why are these valid for changing and visualizing the linear
        ## data.  That does not really make sense to me.
        mesg("Attempting fsva surrogate estimation with ",
             chosen_surrogates, " ", sword, ".")
        type_color <- "darkred"
        sva_object <- sm(sva::sva(log2_mtrx, conditional_model,
                                  null_model, n.sv = chosen_surrogates))
        surrogate_result <- sva::fsva(log2_mtrx, conditional_model,
                                      sva_object, newdat = as.matrix(log2_mtrx),
                                      method = "exact")
        model_adjust <- as.matrix(surrogate_result[["newsv"]])
        source_counts <- surrogate_result[["new"]]
      },
      "isva" = {
        mesg("Attempting isva surrogate estimation with ",
             chosen_surrogates, " ", sword, ".")
        warning("isva, in my estimation, performs incredibly poorly, or I misread the documentation.")
        type_color <- "darkgreen"
        condition_vector <- as.numeric(conditions)
        batch_vector <- as.numeric(batches)

        confounder_lst <- list()
        if (is.null(confounders)) {
          confounder_lst[["batch"]] <- as.numeric(batches)
        } else {
          for (c in 1:length(confounders)) {
            name <- confounders[c]
            confounder_lst[[name]] <- as.numeric(my_design[[name]])
          }
        }
        confounder_mtrx <- matrix(data = confounder_lst[[1]], ncol = 1)
        colnames(confounder_mtrx) <- names(confounder_lst)[1]
        if (length(confounder_lst) > 1) {
          for (i in 2:length(confounder_lst)) {
            confounder_mtrx <- cbind(confounder_mtrx, confounder_lst[[i]])
            names(confounder_mtrx)[i] <- names(confounder_lst)[i]
          }
        }
        mesg("Attempting isva surrogate estimation with ",
             chosen_surrogates, " ", sword, ".")
        surrogate_result <- my_isva(data.m = log2_mtrx, pheno.v = condition_vector,
                                    ncomp = chosen_surrogates,
                                    icamethod = "JADE")
        model_adjust <- as.matrix(surrogate_result[["isv"]])
        output_scale <- "log2"
      },
      "limma" = {
        if (is.null(batch2)) {
          mesg("batch_counts: Using limma's removeBatchEffect to remove batch effect.")
          new_counts <- limma::removeBatchEffect(log2_mtrx, batch = batches)
        } else {
          batches2 <- as.factor(design[[batch2]])
          new_counts <- limma::removeBatchEffect(log2_mtrx, batch = batches, batch2 = batches2)
        }
        message(strwrap(prefix = " ", initial = "", "If you receive a warning: 'NANs produced', one
 potential reason is that the data was quantile normalized."))
        if (expt_state[["transform"]] == "raw") {
          new_counts <- (2 ^ new_counts) - 1
        }
      },
 "limmaresid" = {
   ## ok a caveat:  voom really does require input on the base 10 scale and returns
   ## log2 scale data.  Therefore we need to make sure that the input is
   ## provided appropriately.
   mesg("batch_counts: Using residuals of limma's lmfit to remove batch effect.")
   batch_model <- model.matrix(~batches)
   batch_voom <- limma::voom(linear_mtrx, batch_model,
                             normalize.method = "quantile",
                             plot = FALSE)
   batch_fit <- limma::lmFit(batch_voom, design = batch_model)
   new_counts <- residuals(batch_fit, batch_voom[["E"]])
   ## This is still fubar!
   ##new_counts <- limma::residuals.MArrayLM(batch_fit, batch_voom)
   if (expt_state[["transform"]] == "raw") {
     new_counts <- (2 ^ new_counts) - 1
   }
 },
 "pca" = {
   mesg("Attempting pca surrogate estimation with ",
        chosen_surrogates, " ", sword, ".")
   type_color <- "green"
   data_vs_means <- as.matrix(log2_mtrx - rowMeans(log2_mtrx))
   surrogate_result <- corpcor::fast.svd(data_vs_means)
   model_adjust <- as.matrix(surrogate_result[["v"]][, 1:chosen_surrogates])
 },
 "ruvg" = {
   mesg("Using RUVSeq and edgeR for batch correction (similar to lmfit residuals.)")
   ## Adapted from: http://jtleek.com/svaseq/simulateData.html -- but not quite correct yet
   ruv_input <- edgeR::DGEList(counts = linear_mtrx, group = conditions)
   ruv_input_norm <- ruv_input
   if (expt_state[["normalization"]] == "raw") {
     ruv_input_norm <- edgeR::calcNormFactors(ruv_input, method = "upperquartile")
   }
   ruv_input_glm <- edgeR::estimateGLMCommonDisp(ruv_input_norm, conditional_model)
   ruv_input_tag <- edgeR::estimateGLMTagwiseDisp(ruv_input_glm, conditional_model)
   ruv_fit <- edgeR::glmFit(ruv_input_tag, conditional_model)
   ## Use RUVSeq with empirical controls
   ## The previous instance of ruv_input should work here, and the ruv_input_norm
   ## Ditto for _glm and _tag, and indeed ruv_fit
   ## Thus repeat the first 7 lines of the previous RUVSeq before anything changes.
   ruv_lrt <- edgeR::glmLRT(ruv_fit, coef = 2)
   ruv_control_table <- ruv_lrt[["table"]]
   ranked <- as.numeric(rank(ruv_control_table[["LR"]]))
   bottom_third <- (summary(ranked)[[2]] + summary(ranked)[[3]]) / 2
   ruv_controls <- ranked <= bottom_third  ## what is going on here?!
   ## ruv_controls = rank(ruv_control_table$LR) <= 400  ## some data sets
   ## fail with 400 hard-set
   surrogate_result <- RUVSeq::RUVg(linear_mtrx, ruv_controls, k = chosen_surrogates)
   model_adjust <- surrogate_result[["W"]]
   source_counts <- surrogate_result[["normalizedCounts"]]
 },
 "ruv_empirical" = {
   mesg("Attempting ruvseq empirical surrogate estimation with ",
        chosen_surrogates, " ", sword, ".")
   type_color <- "orange"
   ruv_input <- edgeR::DGEList(counts = linear_mtrx, group = conditions)
   ruv_input_norm <- edgeR::calcNormFactors(ruv_input, method = "upperquartile")
   ruv_input_glm <- edgeR::estimateGLMCommonDisp(ruv_input_norm, conditional_model)
   ruv_input_tag <- edgeR::estimateGLMTagwiseDisp(ruv_input_glm, conditional_model)
   ruv_fit <- edgeR::glmFit(ruv_input_tag, conditional_model)
   ## Use RUVSeq with empirical controls
   ## The previous instance of ruv_input should work here, and the ruv_input_norm
   ## Ditto for _glm and _tag, and indeed ruv_fit
   ## Thus repeat the first 7 lines of the previous RUVSeq before anything changes.
   ruv_lrt <- edgeR::glmLRT(ruv_fit, coef = 2)
   ruv_control_table <- ruv_lrt[["table"]]
   ranked <- as.numeric(rank(ruv_control_table[["LR"]]))
   bottom_third <- (summary(ranked)[[2]] + summary(ranked)[[3]]) / 2
   ruv_controls <- ranked <= bottom_third  ## what is going on here?!
   ## ruv_controls = rank(ruv_control_table$LR) <= 400  ## some data sets
   ## fail with 400 hard-set
   surrogate_result <- RUVSeq::RUVg(round(linear_mtrx), ruv_controls, k = chosen_surrogates)
   model_adjust <- as.matrix(surrogate_result[["W"]])
   source_counts <- surrogate_result[["normalizedCounts"]]
 },
 "ruv_residuals" = {
   mesg("Attempting ruvseq residual surrogate estimation with ",
        chosen_surrogates, " ", sword, ".")
   type_color <- "purple"
   ## Use RUVSeq and residuals
   ruv_input <- edgeR::DGEList(counts = linear_mtrx, group = conditions)
   norm <- edgeR::calcNormFactors(ruv_input)
   ruv_input <- try(edgeR::estimateDisp(norm, design = conditional_model, robust = TRUE))
   ruv_fit <- edgeR::glmFit(ruv_input, conditional_model)
   ruv_res <- residuals(ruv_fit, type = "deviance")
   ruv_normalized <- EDASeq::betweenLaneNormalization(linear_mtrx, which = "upper")
   ## This also gets mad if you pass it a df and not matrix
   controls <- rep(TRUE, dim(linear_mtrx)[1])
   surrogate_result <- RUVSeq::RUVr(ruv_normalized, controls, k = chosen_surrogates, ruv_res)
   model_adjust <- as.matrix(surrogate_result[["W"]])
   source_counts <- as.matrix(surrogate_result[["normalizedCounts"]])
 },
 "ruv_supervised" = {
   mesg("Attempting ruvseq supervised surrogate estimation with ",
        chosen_surrogates, " ", sword, ".")
   type_color <- "black"
   ## Re-calculating the numer of surrogates with this modified data.
   surrogate_estimate <- sm(sva::num.sv(dat = log2_mtrx, mod = conditional_model))
   if (min(rowSums(linear_mtrx)) == 0) {
     warning("empirical.controls will likely fail because some rows are all 0.")
   }
   control_likelihoods <- sm(sva::empirical.controls(
                                      dat = log2_mtrx,
                                      mod = conditional_model,
                                      mod0=null_model,
                                      n.sv = surrogate_estimate))
   surrogate_result <- RUVSeq::RUVg(round(linear_mtrx),
                                    k = surrogate_estimate,
                                    cIdx = as.logical(control_likelihoods))
   model_adjust <- as.matrix(surrogate_result[["W"]])
   source_counts <- surrogate_result[["normalizedCounts"]]
 },
 "smartsva" = {
   mesg("Attempting svaseq estimation with ",
        chosen_surrogates, " ", sword, ".")
   surrogate_result <- SmartSVA::smartsva.cpp(
                                     dat = linear_mtrx,
                                     mod = conditional_model,
                                     mod0=null_model,
                                     n.sv = chosen_surrogates)
   model_adjust <- as.matrix(surrogate_result[["sv"]])
 },
 "svaseq" = {
   mesg("Attempting svaseq estimation with ",
        chosen_surrogates, " ", sword, ".")
   surrogate_result <- sm(sva::svaseq(dat = linear_mtrx,
                                      n.sv = chosen_surrogates,
                                      mod = conditional_model,
                                      mod0 = null_model))
   model_adjust <- as.matrix(surrogate_result[["sv"]])
 },
 "sva_supervised" = {
   mesg("Attempting sva supervised surrogate estimation with ",
        chosen_surrogates, " ", sword, ".")
   type_color <- "red"
   surrogate_result <- sm(sva::ssva(dat = log2_mtrx,
                                    controls = control_likelihoods,
                                    n.sv = chosen_surrogates))
   model_adjust <- as.matrix(surrogate_result[["sv"]])
 },
 "sva_unsupervised" = {
   mesg("Attempting sva unsupervised surrogate estimation with ",
        chosen_surrogates, " ", sword, ".")
   type_color <- "blue"
   if (min(rowSums(linear_mtrx)) == 0) {
     warning("sva will likely fail because some rowSums are 0.")
   }
   surrogate_result <- sm(sva::sva(dat = log2_mtrx,
                                   mod = conditional_model,
                                   mod0 = null_model,
                                   n.sv = chosen_surrogates))
   model_adjust <- as.matrix(surrogate_result[["sv"]])
 },
 "varpart" = {
   mesg("Taking residuals from a linear mixed model as suggested by variancePartition.")
   cl <- parallel::makeCluster(cpus)
   doParallel::registerDoParallel(cl)
   batch_model <- as.formula("~ (1|batch)")
   mesg("The function fitvarPartModel may take excessive memory, you have been warned.")
   batch_fit <- variancePartition::fitVarPartModel(linear_mtrx, formula = batch_model, design)
   new_counts <- residuals(batch_fit)
   rm(batch_fit)
   parallel::stopCluster(cl)
 },
 {
   type_color <- "grey"
   ## If given nothing to work with, use supervised sva
   mesg("Did not understand ", estimate_type, ", assuming supervised sva.")
   surrogate_result <- sva::svaseq(dat = linear_mtrx,
                                   mod = conditional_model,
                                   mod0 = null_model,
                                   n.sv = chosen_surrogates,
                                   controls = control_likelihoods)
   model_adjust <- as.matrix(surrogate_result[["sv"]])
 }) ## End of the switch.

  surrogate_plots <- NULL
  if (!is.null(model_adjust)) {
    rownames(model_adjust) <- sample_names
    sv_names <- glue("SV{1:ncol(model_adjust)}")
    colnames(model_adjust) <- sv_names
    if ("expt" %in% class(input)) {
      surrogate_plots <- plot_batchsv(input, model_adjust)
    }
  }
  ## Only use counts_from_surrogates if the method does not provide counts on its own
  if (is.null(new_counts)) {
    mesg("Adjusting counts with method: ", adjust_method, ".")
    new_counts <- counts_from_surrogates(data = surrogate_input, adjust = model_adjust,
                                         method = adjust_method, design = my_design)
    if (output_scale == "log2") {
      new_counts <- (2 ^ new_counts) - 1
    }
  }

  ret <- list(
      "surrogate_result" = surrogate_result,
      "null_model" = null_model,
      "model_adjust" = model_adjust,
      "source_counts" = source_counts,
      "new_counts" = new_counts,
      "sample_factor" = surrogate_plots[["sample_factor"]],
      "factor_svs" = surrogate_plots[["factor_svs"]],
      "svs_sample" = surrogate_plots[["svs_sample"]])
  return(ret)
}

#' Perform different batch corrections using limma, sva, ruvg, and cbcbSEQ.
#'
#' I found this note which is the clearest explanation of what happens with
#' batch effect data:
#' https://support.bioconductor.org/p/76099/
#' Just to be clear, there's an important difference between removing a batch
#' effect and modelling a batch effect. Including the batch in your design
#' formula will model the batch effect in the regression step, which means that
#' the raw data are not modified (so the batch effect is not removed), but
#' instead the regression will estimate the size of the batch effect and
#' subtract it out when performing all other tests. In addition, the model's
#' residual degrees of freedom will be reduced appropriately to reflect the fact
#' that some degrees of freedom were "spent" modelling the batch effects. This
#' is the preferred approach for any method that is capable of using it (this
#' includes DESeq2). You would only remove the batch effect (e.g. using limma's
#' removeBatchEffect function) if you were going to do some kind of downstream
#' analysis that can't model the batch effects, such as training a classifier.
#' I don't have experience with ComBat, but I would expect that you run it on
#' log-transformed CPM values, while DESeq2 expects raw counts as input. I
#' couldn't tell you how to properly use the two methods together.
#'
#' @param count_table Matrix of (pseudo)counts.
#' @param method Choose the method for batch/surrogate estimation.
#' @param expt_design Model matrix defining the experimental conditions/batches/etc.
#' @param batch1 String describing the method to try to remove the batch effect
#'  (or FALSE to leave it alone, TRUE uses limma).
#' @param current_state Current state of the expt in an attempt to avoid
#'  double-normalization.
#' @param current_design Redundant with expt_design above, but
#'  provides another place for normalize_expt() to send data.
#' @param expt_state Current state of the data
#' @param surrogate_method Also redundant for normalize_expt()
#' @param surrogates Number of surrogates or method to estimate them.
#' @param low_to_zero Send <0 entries to 0 to avoid shenanigans.
#' @param cpus Parallelize intensive operations.
#' @param batch2 Column in the design table describing the second covariant to
#'  remove (only used by limma at the moment).
#' @param noscale Used for combatmod, when true it removes the scaling
#'  parameter from the invocation of the modified combat.
#' @param ... More options for you!
#' @return The 'batch corrected' count table and new library size.  Please
#'  remember that the library size which comes out of this may not be what you
#'  want for voom/limma and would therefore lead to spurious differential
#'  expression values.
#' @seealso [limma] [edgeR] [RUVSeq] [sva]
#' @examples
#' \dontrun{
#'  limma_batch <- batch_counts(table, design, batch1='batch', batch2='strain')
#'  sva_batch <- batch_counts(table, design, batch='sva')
#' }
#' @export
batch_counts <- function(count_table, method = TRUE, expt_design = NULL, batch1 = "batch",
                         current_state = NULL, current_design = NULL, expt_state = NULL,
                         surrogate_method = NULL, surrogates = NULL, low_to_zero = FALSE,
                         cpus = 4, batch2 = NULL, noscale = TRUE, ...) {
  arglist <- list(...)
  if (!is.null(arglist[["batch"]])) {
    method <- arglist[["batch"]]
  }
  chosen_surrogates <- NULL
  if (is.null(surrogates) & is.null(surrogate_method)) {
    chosen_surrogates <- "be"
  } else if (!is.null(surrogates)) {
    chosen_surrogates <- surrogates
  } else {
    chosen_surrogates <- surrogate_method
  }

  if (is.null(cpus)) {
    cpus <- parallel::detectCores() - 2
  }
  if (cpus < 0) {
    cpus <- 1
  }

  prior.plots <- FALSE
  if (!is.null(arglist[["prior.plots"]])) {
    message("When using ComBat, using prior.plots may result in an error due to infinite ylim.")
    prior.plots <- arglist[["prior.plots"]]
  }

  ## Lets use expt_state to make sure we know if the data is already log2/cpm/whatever.
  ## We want to use this to back-convert or reconvert data to the appropriate
  ## scale on return.  Note that there are two possible variables, expt_state and current_state,
  ## I am using current_state to keep track of changes made on scale/etc during each step
  ## of the normalization process.
  used_state <- list(
      "filter" = "raw",
      "normalization" = "raw",
      "conversion" = "raw",
      "batch" = "raw",
      "transform" = "raw")
  if (is.null(expt_state) & is.null(current_state)) {
    mesg("Assuming a completely raw expressionset.")
  } else if (!is.null(current_state)) {
    used_state <- current_state
    mesg("Using the current state of normalization.")
  } else {
    used_state <- expt_state
    mesg("Using the initial state of the expressionset.")
  }

  design <- NULL
  if (is.null(expt_design) & is.null(current_design)) {
    stop("I require an experimental design.")
  } else if (!is.null(current_design)) {
    design <- current_design
  } else {
    design <- expt_design
  }

  ## These droplevels calls are required to avoid errors like 'confounded by batch'
  batches <- droplevels(as.factor(design[[batch1]]))
  conditions <- droplevels(as.factor(design[["condition"]]))

  if (isTRUE(method)) {
    method <- "limma"
  }

  count_df <- data.frame(count_table)
  count_mtrx <- as.matrix(count_df)
  conditional_model <- model.matrix(~conditions, data = count_df)
  null_model <- conditional_model[, 1]
  ## Set the number of surrogates for sva/ruv based methods.
  mesg("Passing the data to all_adjusters using the ", method, " estimate type.")
  new_material <- all_adjusters(count_table, design = design, estimate_type = method,
                                cpus = cpus, batch1 = batch1, batch2 = batch2,
                                expt_state = used_state, noscale = noscale,
                                chosen_surrogates = chosen_surrogates,
                                surrogates = surrogates,
                                low_to_zero = low_to_zero)
  count_table <- new_material[["new_counts"]]

  na_idx <- is.na(count_table)
  num_na <- sum(na_idx)
  if (num_na > 0) {
    message("Found ", num_na, " na entries in the new table, setting them to 0.")
    count_table[na_idx] <- 0
  }

  num_low <- sum(count_table <= 0)
  if (is.null(num_low)) {
    num_low <- 0
  }
  if (num_low > 0) {
    elements <- nrow(count_table) * ncol(count_table)
    low_pct <- scales::percent(num_low / elements)
    mesg("There are ", num_low, " (", low_pct,
         ") elements which are < 0 after batch correction.")
    if (isTRUE(low_to_zero)) {
      message("Setting ", num_low, " low elements to zero.")
      count_table[count_table < 0] <- 0
    }
  }
  libsize <- colSums(count_table)
  counts <- list(count_table = count_table, libsize = libsize, result = new_material)
  return(counts)
}

#' A function suggested by Hector Corrada Bravo and Kwame Okrah for batch
#' removal.
#'
#' During a lab meeting, the following function was suggested as a quick and
#' dirty batch removal tool.  It takes data and a model including a 'batch'
#' factor, invokes limma on them, removes the batch factor, does a cross
#' product of the fitted data and modified model and uses that with residuals to
#' get a new data set.
#'
#' @param normalized_counts Data frame of log2cpm counts.
#' @param model Balanced experimental model containing condition and batch
#'  factors.
#' @param conditional_model Experimental model with the conditional factor.
#' @param batch_model Experimental model with the batch factor.
#' @param batch1 Column containing the first batch's metadata in the experimental design.
#' @param condition Column containing the condition information in the metadata.
#' @param matrix_scale Is the data on a linear or log scale?
#' @param return_scale Do you want the data returned on the linear or log scale?
#' @param method I found a couple ways to apply the surrogates to the data.  One
#'  method subtracts the residuals of a batch model, the other adds the
#'  conditional.
#' @return Dataframe of residuals after subtracting batch from the model.
#' @seealso [limma::voom()] [limma::lmFit()]
#' @examples
#' \dontrun{
#'  newdata <- cbcb_batch_effect(counts, expt_model)
#' }
#' @export
cbcb_batch <- function(normalized_counts, model,
                       conditional_model = NULL,
                       batch_model = NULL,
                       batch1="batch", condition = "condition",
                       matrix_scale = "linear", return_scale = "linear",
                       method = "subtract") {
  batch_idx <- grep(pattern = batch1, x = colnames(model))
  cond_idx <- grep(pattern = condition, x = colnames(model))
  if (is.null(batch_model)) {
    batch_model <- model[, batch_idx] <- 0
  }
  if (is.null(conditional_model)) {
    conditional_model <- model[, cond_idx] <- 0
  }
  ## Voom takes counts on the linear scale, so change them if they are log.
  ## It also does a log2(counts + 0.5), so take that into account as well.
  if (matrix_scale == "log2") {
    normalized_counts <- (2 ^ normalized_counts) - 0.5
  } else if (matrix_scale != "linear") {
    stop("I do not understand the scale: ", matrix_scale, ".")
  }
  normal_voom <- limma::voom(normalized_counts, design = model, plot = FALSE)
  cond_voom <- limma::voom(normalized_counts, design = conditional_model, plot = FALSE)
  batch_voom <- limma::voom(normalized_counts, design = batch_model, plot = FALSE)
  if (method == "subtract") {
    modified_fit <- limma::lmFit(batch_voom)
    new_data <- residuals(modified_fit, batch_voom)
  } else if (method == "add") {
    fit <- limma::lmFit(normal_voom)
    ## I got confusered here, this might be incorrect.
    ##new_data <- tcrossprod(normal_voom[["coefficient"]], cond_modified_model) +
    ##  residuals(normal_voom, normalized_counts)
    new_data <- tcrossprod(normal_voom[["coefficient"]], conditional_model) +
      residuals(normal_voom, normalized_counts)
  } else {
    stop("This currently only understands 'add' or 'subtract'.")
  }
  ## The new_data is on the log scale, so switch back assuming linear is required.
  if (return_scale == "linear") {
    ## 0.5 is used here due to limma::voom()
    new_data <- (2 ^ new_data) - 0.5
  } else if (return_scale != "log2") {
    stop("I do not understand the return scale: ", return_scale, ".")
  }
  return(new_data)
}


#' Attempt to compare the results from the various batch/sv methods.
#'
#' Given an expressionset and list of methods, try to find out how
#' well the various methods agree via correlation.
#'
#' @param expt Input expressionset
#' @param methods Set of methods to try out.
compare_batches <- function(expt = NULL, methods = NULL) {
  if (is.null(methods)) {
    ## writing this oddly until I work out which ones to include
    methods <- c(
        "combat",
        "combatmod",
        "combat_notnull",
        ## "combat_noprior",
        ## "combat_noprior_scale",
        "combat_scale",
        "fsva",
        "isva",
        "limma",
        ## "limmaresid",
        "pca",
        "ruv_empirical",
        "ruvg",
        "ruv_residuals",
        "ruv_supervised",
        "smartsva",
        "svaseq",
        "sva_supervised",
        "sva_unsupervised")
    ## "varpart")
  }
  if (is.null(expt)) {
    message("Going to make an spombe expressionSet from the fission data set.")
    expt <- make_pombe_expt()
  }
  combined <- data.frame()
  lst <- list()
  for (m in 1:length(methods)) {
    method <- methods[m]
    res <- exprs(normalize_expt(expt, filter = TRUE, batch = method))
    lst[[method]] <- res
    column <- c()
    names <- c()
    for (c in 1:length(colnames(res))) {
      names <- c(names, rownames(res))
      column <- c(column, res[, c])
    }
    column <- as.data.frame(column)
    rownames(column) <- make.unique(names)
    if (m == 1) {
      combined <- column
    } else {
      combined <- merge(combined, column, by = "row.names")
      rownames(combined) <- combined[["Row.names"]]
      combined <- combined[, -1]
      colnames(combined)[m] <- method
    }
  }
  ## fun <- corrplot::corrplot(cor(combined), method = "ellipse", type = "lower", tl.pos = "d")
  ## tt <- cor(combined)
  something <- plot_disheat(combined)
}

#' Perform a comparison of the surrogate estimators demonstrated by Jeff Leek.
#'
#' This is entirely derivative, but seeks to provide similar estimates for one's
#' own actual data and catch corner cases not taken into account in that
#' document (for example if the estimators don't converge on a surrogate
#' variable). This will attempt each of the surrogate estimators described by
#' Leek: pca, sva supervised, sva unsupervised, ruv supervised, ruv residuals,
#' ruv empirical. Upon completion it will perform the same limma expression
#' analysis and plot the ranked t statistics as well as a correlation plot
#' making use of the extracted estimators against condition/batch/whatever
#' else. Finally, it does the same ranking plot against a linear fitting Leek
#' performed and returns the whole pile of information as a list.
#'
#' @param expt Experiment containing a design and other information.
#' @param extra_factors Character list of extra factors which may be included in
#'  the final plot of the data.
#' @param filter_it Most of the time these surrogate methods get mad if there
#'  are 0s in the data.  Filter it?
#' @param filter_type Type of filter to use when filtering the input data.
#' @param do_catplots Include the catplots?  They don't make a lot of sense yet,
#'  so probably no.
#' @param surrogates Use 'be' or 'leek' surrogate estimates, or choose a
#'  number.
#' @param ... Extra arguments when filtering.
#' @return List of the results.
#' @seealso [normalize_expt()] [plot_pca()] [all_adjuster()] [corrplot] [ffpe]
#' @export
compare_surrogate_estimates <- function(expt, extra_factors = NULL,
                                        filter_it = TRUE, filter_type = TRUE,
                                        do_catplots = FALSE, surrogates = "be", ...) {
  arglist <- list(...)
  design <- pData(expt)
  do_batch <- TRUE
  if (length(levels(design[["batch"]])) == 1) {
    mesg("There is 1 batch in the data, fitting condition+batch will fail.")
    do_batch <- FALSE
  }

  if (isTRUE(filter_it) & expt[["state"]][["filter"]] == "raw") {
    mesg("The expt has not been filtered, ",
         "set filter_type/filter_it if you want other options.")
    expt <- sm(normalize_expt(expt, filter = filter_type,
                              ...))
  }
  pca_plots <- list()
  pca_plots[["null"]] <- plot_pca(expt)[["plot"]]

  pca_adjust <- all_adjusters(expt, estimate_type = "pca",
                              surrogates = surrogates)
  pca_plots[["pca"]] <- plot_pca(pca_adjust[["new_counts"]],
                                 design = design,
                                 plot_colors = expt[["colors"]])[["plot"]]

  sva_supervised <- sm(all_adjusters(expt, estimate_type = "sva_supervised",
                                     surrogates = surrogates))
  pca_plots[["svasup"]] <- plot_pca(sva_supervised[["new_counts"]],
                                    design = design,
                                    plot_colors = expt[["colors"]])[["plot"]]

  sva_unsupervised <- sm(all_adjusters(expt, estimate_type = "sva_unsupervised",
                                       surrogates = surrogates))
  pca_plots[["svaunsup"]] <- plot_pca(sva_unsupervised[["new_counts"]],
                                      design = design,
                                      plot_colors = expt[["colors"]])[["plot"]]

  ruv_supervised <- sm(all_adjusters(expt, estimate_type = "ruv_supervised",
                                     surrogates = surrogates))
  pca_plots[["ruvsup"]] <- plot_pca(ruv_supervised[["new_counts"]],
                                    design = design,
                                    plot_colors = expt[["colors"]])[["plot"]]

  ruv_residuals <- sm(all_adjusters(expt, estimate_type = "ruv_residuals",
                                    surrogates = surrogates))
  pca_plots[["ruvresid"]] <- plot_pca(ruv_residuals[["new_counts"]],
                                      design = design,
                                      plot_colors = expt[["colors"]])[["plot"]]

  ruv_empirical <- sm(all_adjusters(expt, estimate_type = "ruv_empirical",
                                    surrogates = surrogates))
  pca_plots[["ruvemp"]] <- plot_pca(ruv_empirical[["new_counts"]],
                                    design = design,
                                    plot_colors = expt[["colors"]])[["plot"]]

  first_svs <- data.frame(
      "condition" = as.numeric(as.factor(expt[["conditions"]])),
      "batch" = as.numeric(as.factor(expt[["batches"]])),
      "pca_adjust" = pca_adjust[["model_adjust"]][, 1],
      "sva_supervised" = sva_supervised[["model_adjust"]][, 1],
      "sva_unsupervised" = sva_unsupervised[["model_adjust"]][, 1],
      "ruv_supervised" = ruv_supervised[["model_adjust"]][, 1],
      "ruv_residuals" = ruv_residuals[["model_adjust"]][, 1],
      "ruv_empirical" = ruv_empirical[["model_adjust"]][, 1])
  batch_adjustments <- list(
      "condition" = as.factor(expt[["conditions"]]),
      "batch" = as.factor(expt[["batches"]]),
      "pca_adjust" = pca_adjust[["model_adjust"]],
      "sva_supervised" = sva_supervised[["model_adjust"]],
      "sva_unsupervised" = sva_unsupervised[["model_adjust"]],
      "ruv_supervised" = ruv_supervised[["model_adjust"]],
      "ruv_residuals" = ruv_residuals[["model_adjust"]],
      "ruv_empirical" = ruv_empirical[["model_adjust"]])
  batch_names <- c("condition", "batch", "pca", "sva_sup", "sva_unsup",
                   "ruv_sup", "ruv_resid", "ruv_emp")
  first_samples <- data.frame(
      "pca_adjust" = pca_adjust[["new_counts"]][, 1],
      "sva_supervised" = sva_supervised[["new_counts"]][, 1],
      "sva_unsupervised" = sva_unsupervised[["new_counts"]][, 1],
      "ruv_supervised" = ruv_supervised[["new_counts"]][, 1],
      "ruv_residuals" = ruv_residuals[["new_counts"]][, 1],
      "ruv_empirical" = ruv_empirical[["new_counts"]][, 1])

  if (!is.null(extra_factors)) {
    for (fact in extra_factors) {
      if (!is.null(design[, fact])) {
        batch_names <- append(x = batch_names, values = fact)
        first_svs[[fact]] <- as.numeric(as.factor(design[, fact]))
        batch_adjustments[[fact]] <- as.numeric(as.factor(design[, fact]))
      }
    }
  }
  correlations <- cor(first_svs)
  corrplot::corrplot(correlations, type = "lower", method = "ellipse", tl.pos = "d")
  ret_plot <- grDevices::recordPlot()
  sample_correlations <- cor(first_samples)
  corrplot::corrplot(sample_correlations, method = "ellipse", type = "upper", tl.pos = "d")
  sample_dist <- plot_disheat(first_samples)

  sample_corplot <- grDevices::recordPlot()
  adjustments <- c("+ batch_adjustments$batch", "+ batch_adjustments$pca",
                   "+ batch_adjustments$sva_sup", "+ batch_adjustments$sva_unsup",
                   "+ batch_adjustments$ruv_sup", "+ batch_adjustments$ruv_resid",
                   "+ batch_adjustments$ruv_emp")
  adjust_names <- gsub(
      pattern = "^.*adjustments\\$(.*)$", replacement = "\\1", x = adjustments)
  starter <- edgeR::DGEList(counts = exprs(expt))
  norm_start <- edgeR::calcNormFactors(starter)

  ## Create a baseline to compare against.
  null_formula <- as.formula("~ condition ")
  null_limma_design <- model.matrix(null_formula, data = design)
  null_voom_result <- limma::voom(norm_start, null_limma_design, plot = FALSE)
  null_limma_fit <- limma::lmFit(null_voom_result, null_limma_design)
  null_fit <- limma::eBayes(null_limma_fit)
  null_tstat <- null_fit[["t"]]
  null_catplot <- NULL
  if (isTRUE(do_catplots)) {
    if (!isTRUE("ffpe" %in% .packages(all.available = TRUE))) {
      ## ffpe has some requirements which do not install all the time.
      tt <- please_install("ffpe")
    }
    if (isTRUE("ffpe" %in% .packages(all.available = TRUE))) {
      null_catplot <- ffpe::CATplot(-rank(null_tstat), -rank(null_tstat),
                                    maxrank = 1000, make.plot = TRUE)
    } else {
      catplots[[adjust_name]] <- NULL
    }
  }

  catplots <- vector("list", length(adjustments))  ## add 1 for a null adjustment
  names(catplots) <- adjust_names
  tstats <- list()
  oldpar <- par(mar = c(5, 5, 5, 5))
  num_adjust <- length(adjust_names)
  ## Now perform other adjustments
  for (a in 1:num_adjust) {
    adjust_name <- adjust_names[a]
    adjust <- adjustments[a]
    if (adjust_name == "batch" & !isTRUE(do_batch)) {
      message("A friendly reminder that there is only 1 batch in the data.")
      tstats[[adjust_name]] <- null_tstat
      catplots[[adjust_name]] <- null_catplot
    } else {
      mesg(a, "/", num_adjust, ": Performing lmFit(data) etc. with ",
           adjust_name, " in the model.")
      modified_formula <- as.formula(glue("~ condition {adjust}"))
      limma_design <- model.matrix(modified_formula, data = design)
      voom_result <- limma::voom(norm_start, limma_design, plot = FALSE)
      limma_fit <- limma::lmFit(voom_result, limma_design)
      modified_fit <- limma::eBayes(limma_fit)
      tstats[[adjust_name]] <- modified_fit[["t"]]
      ##names(tstats[[counter]]) <- as.character(1:dim(data)[1])
      catplot_together <- NULL
      if (isTRUE(do_catplots)) {
        if (!isTRUE("ffpe" %in% .packages(all.available = TRUE))) {
          ## ffpe has some requirements which do not install all the time.
          tt <- please_install("ffpe")
        }
        if (isTRUE("ffpe" %in% .packages(all.available = TRUE))) {
          catplots[[adjust_name]] <- ffpe::CATplot(
                                               rank(tstats[[adjust_name]]), rank(null_tstat),
                                               maxrank = 1000, make.plot = TRUE)
        } else {
          catplots[[adjust_name]] <- NULL
        }
      }
    }
  } ## End for a in 2:length(adjustments)

  ## Final catplot plotting, if necessary.
  if (isTRUE(do_catplots)) {
    catplot_df <- as.data.frame(catplots[[1]][[2]])
    for (c in 2:length(catplots)) {
      cat <- catplots[[c]]
      catplot_df <- cbind(catplot_df, cat[["concordance"]])
    }
    colnames(catplot_df) <- names(catplots)
    catplot_df[["x"]] <- rownames(catplot_df)
    gg_catplot <- reshape2::melt(data = catplot_df, id.vars = "x")
    colnames(gg_catplot) <- c("x", "adjust", "y")
    gg_catplot[["x"]] <- as.numeric(gg_catplot[["x"]])
    gg_catplot[["y"]] <- as.numeric(gg_catplot[["y"]])

    cat_plot <- ggplot(data = gg_catplot, mapping = aes_string(x = "x", y = "y", color = "adjust")) +
      ggplot2::geom_point() +
      ggplot2::geom_jitter() +
      ggplot2::geom_line() +
      ggplot2::xlab("Rank") +
      ggplot2::ylab("Concordance") +
      ggplot2::theme_bw()
  } else {
    cat_plot <- NULL
  }

  ret <- list(
      "pca_adjust" = pca_adjust,
      "sva_supervised_adjust" = sva_supervised,
      "sva_unsupervised_adjust" = sva_unsupervised,
      "ruv_supervised_adjust" = ruv_supervised,
      "ruv_residual_adjust" = ruv_residuals,
      "ruv_empirical_adjust" = ruv_empirical,
      "adjustments" = batch_adjustments,
      "correlations" = correlations,
      "plot" = ret_plot,
      "sample_corplot" = sample_corplot,
      "pca_plots" = pca_plots,
      "catplots" = cat_plot)
  return(ret)
}

#' A single place to extract count tables from a set of surrogate variables.
#'
#' Given an initial set of counts and a series of surrogates, what would the
#' resulting count table look like? Hopefully this function answers that
#' question.
#'
#' @param data Original count table, may be an expt/expressionset or df/matrix.
#' @param adjust Surrogates with which to adjust the data.
#' @param design Experimental design if it is not included in the expressionset.
#' @param method Which methodology to follow, ideally these agree but that seems untrue.
#' @param cond_column design column containing the condition data.
#' @param batch_column design column with the batch data, used for
#'  subtractive methods.
#' @param matrix_scale Was the input for the surrogate estimator on a log or linear scale?
#' @param return_scale Does one want the output linear or log?
#' @param ... Arguments passed to downstream functions.
#' @return A data frame of adjusted counts.
#' @seealso [sva] [RUVSeq] [crossprod()] [tcrossprod()] [solve()]
#' @export
counts_from_surrogates <- function(data, adjust = NULL, design = NULL, method = "ruv",
                                   cond_column = "condition", batch_column = "batch",
                                   matrix_scale = "linear", return_scale = "linear", ...) {
  arglist <- list(...)
  data_mtrx <- NULL
  my_design <- NULL
  conditions <- NULL
  batches <- NULL
  if (class(data)[1] == "expt") {
    my_design <- pData(data)
    conditions <- droplevels(as.factor(pData(data)[[cond_column]]))
    batches <- droplevels(as.factor(pData(data)[[batch_column]]))
    data_mtrx <- exprs(data)
  } else if (class(data)[1] == "ExpressionSet") {
    my_design <- pData(data)
    conditions <- droplevels(as.factor(pData(data)[[cond_column]]))
    batches <- droplevels(as.factor(pData(data)[[batch_column]]))
    data_mtrx <- exprs(data)
  } else {
    my_design <- design
    conditions <- droplevels(as.factor(design[[cond_column]]))
    batches <- droplevels(as.factor(design[[batch_column]]))
    data_mtrx <- as.matrix(data)
  }
  conditional_model <- model.matrix(~ conditions, data = my_design)

  new_model <- conditional_model
  ## Explicitly append columns of the adjust matrix to the conditional model.
  ## In the previous code, this was: 'X <- cbind(conditional_model, sva$sv)'
  ## new_model <- cbind(conditional_model, adjust)
  new_colnames <- colnames(conditional_model)
  if (is.null(adjust)) {
    message("No adjust was provided, leaving the data alone.")
    adjust <- data.frame(row.names = rownames(my_design))
    adjust[["SV1"]] <- 1
  }
  adjust_mtrx <- as.matrix(adjust)
  ##for (col in 1:ncol(adjust_mtrx)) {
  ##  new_model <- cbind(new_model, adjust_mtrx[, col])
  ##  new_colname <- glue("sv{col}")
  ##  new_colnames <- append(new_colnames, new_colname)
  ##}
  ##colnames(new_model) <- new_colnames
  full_models <- sm(choose_model(data, conditions = conditions, batches = batches,
                                 model_batch = adjust_mtrx))
  full_model <- full_models[["noint_model"]]
  ## FIXME: Reverse this, it looks weird.
  batch_model <- full_models[["model_batch"]]
  cond_model <- sm(choose_model(data, conditions = conditions,
                                batches = batches)[["noint_model"]])

  switchret <- switch(
      method,
      "cbcb_add" = {
        ## we also have conditional_model which may make this easier
        new_counts <- cbcb_batch(data_mtrx, full_model,
                                 conditional_model = cond_model,
                                 batch_model = batch_model,
                                 method = "add",
                                 matrix_scale = matrix_scale,
                                 return_scale = return_scale,
                                 ...)
      },
      "cbcb_subtract" = {
        ##new_counts <- cbcb_batch(data_mtrx, my_design, method = "subtract",
        ##                         matrix_scale = matrix_scale, return_scale = return_scale,
        ##                         ...)
        new_counts <- cbcb_batch(data_mtrx, new_model,
                                 conditional_model = cond_model,
                                 batch_model = batch_model,
                                 method = "subtract",
                                 matrix_scale = matrix_scale, return_scale = return_scale,
                                 ...)
      },
      "ruv" = {
        ## Here is the original code, as a reminder: W is the matrix of surrogates.
        ## W <- svdWa$u[, (first:k), drop = FALSE]
        ## alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
        ## correctedY <- Y - W %*% alpha
        ## if(!isLog & all(.isWholeNumber(x))) {
        ##   if(round) {
        ##     correctedY <- round(exp(correctedY) - epsilon)
        ##     correctedY[correctedY<0] <- 0
        ##   } else {
        ##     correctedY <- exp(correctedY) - epsilon
        ##   }
        ## }
        ## colnames(W) <- paste("W", seq(1, ncol(W)), sep = "_")
        ## return(list(W = W, normalizedCounts = t(correctedY)))

        ##alpha <- try(solve(t(adjust_mtrx) %*% adjust_mtrx))
        ## Y <- t(data_mtrx)
        ## W <- ruv_model

        log_data_mtrx <- NULL
        if (matrix_scale != "log2") {
          ## Note that we are actually going on scale e.
          log_data_mtrx <- log(data_mtrx + 1)
        }

        ## original_alpha <- solve(t(adjust_mtrx) %*% adjust_mtrx) %*% t(adjust_mtrx) %*% t(data_mtrx)
        alpha <- try(solve(crossprod(adjust_mtrx)))
        if (class(alpha)[1] == "try-error") {
          warning("Data modification by the model failed. Leaving counts untouched.")
          return(data_mtrx)
        }
        beta <- tcrossprod(t(adjust_mtrx), log_data_mtrx)
        gamma <- alpha %*% beta
        delta <- t(log_data_mtrx) - (adjust_mtrx %*% gamma)
        new_counts <- t(delta)

        if (return_scale == "log2") {
          new_counts <- new_counts / log(2)
        } else if (return_scale == "linear") {
          new_counts <- exp(new_counts) - 1
        } else {
          stop("I do not understand the scale: ", return_scale, ".")
        }
      },
      "solve_crossproducts" = {
        ## I think that if I did the math out, this is the same as ruv above.

        ##data_modifier <- try(solve(t(new_model) %*% new_model) %*% t(new_model))
        ## In the previous code, this was: 'Hat <- solve(t(X) %*% X) %*% t(X)'
        ## Now it is in two separate lines, first the solve operation:

        if (matrix_scale != "log2") {
          ## Note that we are actually going on scale e.
          data_mtrx <- log2(data_mtrx + 1)
        }

        data_solve <- try(solve(t(new_model) %*% new_model), silent = TRUE)
        if (class(data_solve)[1] == "try-error") {
          warning("Data modification by the model failed. Leaving counts untouched.")
          return(data_mtrx)
        }
        ## If the solve operation passes, then the '%*% t(X)' is allowed to happen.
        data_modifier <- data_solve %*% t(new_model)
        transformation <- (data_modifier %*% t(data_mtrx))
        conds <- ncol(conditional_model)
        new_counts <- data_mtrx - t(as.matrix(new_model[, -c(1:conds)]) %*%
                                    transformation[-c(1:conds), ])
        if (return_scale == "linear") {
          new_counts <- (2 ^ new_counts) - 1
        } else if (return_scale != "log2") {
          stop("I do not understand the return scale: ", return_scale, ".")
        }
      },
      {
        stop("I do not understand method: ", method, ".")
      }
  ) ## End the switch

  ## If the matrix state and return state are not the same, fix it.
  ## It appears to me that the logic of this is wrong, but I am not yet certain why.
  return(new_counts)
}

#' A modified version of comBatMod.
#'
#' This is a hack of Kwame Okrah's combatMod to make it not fail on corner-cases.
#' This was mostly copy/pasted from
#' https://github.com/kokrah/cbcbSEQ/blob/master/R/transform.R
#'
#' @param dat Df to modify.
#' @param batch Factor of batches.
#' @param mod Factor of conditions.
#' @param noscale The normal 'scale' option squishes the data too much, so this
#'   defaults to TRUE.
#' @param prior.plots Print out prior plots?
#' @param ... Extra options are passed to arglist
#' @return Df of batch corrected data
#' @seealso [sva] [sva::ComBat()]
#' @examples
#' \dontrun{
#'  df_new = cbcb_combat(df, batches, model)
#' }
#' @export
cbcb_combat <- function(dat, batch, mod, noscale = TRUE, prior.plots = FALSE, ...) {
  arglist <- list(...)
  par.prior <- TRUE
  numCovs <- NULL
  mod <- cbind(mod, batch)
  check <- apply(mod, 2, function(x) all(x == 1))
  mod <- as.matrix(mod[, !check])
  colnames(mod)[ncol(mod)] <- "Batch"
  if (sum(check) > 0 & !is.null(numCovs)) {
    numCovs <- numCovs - 1
  }
  design <- survJamda::design.mat(mod)
  batches <- survJamda::list.batch(mod)
  n.batch <- length(batches)
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  NAs <- any(is.na(dat))
  B.hat <- NULL
  ## This is taken from sva's github repository in helper.R
  Beta.NA <- function(y, X) {
    des <- X[!is.na(y), ]
    y1 <- y[!is.na(y)]
    B <- solve(t(des)%*%des)%*%t(des)%*%y1
    B
  }
  var.pooled <- NULL
  mesg("Standardizing data across genes\n")
  if (NAs) {
    warning(glue("Found {sum(is.na(dat)} missing data values."))
    warning("The original combatMod uses an undefined variable Beta.NA here,
I set it to 1 not knowing what its purpose is.")
    B.hat <- apply(dat, 1, Beta.NA)
  } else {
    ## There are no NAs in the data, this is a good thing(Tm)!
    B.hat <- solve(t(design) %*% design) %*% t(design) %*% t(as.matrix(dat))
  }
  grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]

  if (NAs) {
    var.pooled <- apply(dat - t(design %*% B.hat), 1, var, na.rm = TRUE)
  } else {
    transposed <- t(design %*% B.hat)
    subtracted <- dat - transposed
    second_half <- rep(1 / n.array, n.array)
    var.pooled <- as.matrix(subtracted ^ 2) %*% as.matrix(second_half)
  }
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1, n.array)))
  if (isTRUE(noscale)) {
    m.data <- dat - stand.mean
    second_half <- as.matrix(rep(1 / (n.array - ncol(design)), n.array))
    mse <- as.matrix((dat - t(design %*% B.hat)) ^ 2) %*% second_half
    hld <- NULL
    bayesdata <- dat
    for (k in 1:n.batch) {
      mesg("Fitting 'shrunk' batch ", k, " effects.")
      sel <- batches[[k]]
      gammaMLE <- rowMeans(m.data[, sel])
      mprior <- mean(gammaMLE, na.rm = TRUE)
      vprior <- var(gammaMLE, na.rm = TRUE)
      prop <- vprior / (mse / (length(sel)) + vprior)
      gammaPost <- prop * gammaMLE + (1 - prop) * mprior
      for (i in sel) {
        bayesdata[, i] <- bayesdata[, i] - gammaPost
      }
      stats <- data.frame(gammaPost = gammaPost, gammaMLE = gammaMLE, prop = prop)
      hld[[paste("Batch", k, sep = ".")]] <- list(
          "stats" = stats,
          "indices" = sel,
          "mprior" = mprior,
          "vprior" = vprior)
    }
    mesg("Adjusting data for batch effects.")
    return(bayesdata)
  } else {
    mesg("Fitting L/S model and finding priors.")
    batch.design <- design[, 1:n.batch]
    if (NAs) {
      gamma.hat <- apply(s.data, 1, Beta.NA, batch.design)
    } else {
      gamma.hat <- solve(t(batch.design) %*% batch.design) %*%
        t(batch.design) %*% t(as.matrix(s.data))
    }
    delta.hat <- NULL
    for (i in batches) {
      delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var, na.rm = TRUE))
    }
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, aprior)
    b.prior <- apply(delta.hat, 1, bprior)
    if (prior.plots & par.prior) {
      oldpar <- par(mfrow = c(2, 2))
      tmp <- density(gamma.hat[1, ], na.rm = TRUE)
      plot(tmp, type = "l", main = "Density Plot")
      xx <- seq(min(tmp$x), max(tmp$x), length = 100)
      lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
      stats::qqnorm(gamma.hat[1, ])
      stats::qqline(gamma.hat[1, ], col = 2)
      tmp <- stats::density(delta.hat[1, ], na.rm = TRUE)
      invgam <- 1 / stats::rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
      tmp1 <- try(stats::density(invgam, na.rm = TRUE))
      plot(tmp, typ = "l", main = "Density Plot", ylim = c(0, max(tmp$y, tmp1$y)))
      if (class(tmp1)[1] != "try-error") {
        lines(tmp1, col = 2)
      }
      try(stats::qqplot(delta.hat[1, ], invgam,
                        xlab = "Sample Quantiles",
                        ylab = "Theoretical Quantiles"))
      lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
      title("Q-Q Plot")
      newpar <- par(oldpar)
    }
    gamma.star <- delta.star <- NULL
    if (par.prior) {
      mesg("Finding parametric adjustments.")
      for (i in 1:n.batch) {
        temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                       delta.hat[i, ], gamma.bar[i],
                       t2[i], a.prior[i], b.prior[i])
        gamma.star <- rbind(gamma.star, temp[1, ])
        delta.star <- rbind(delta.star, temp[2, ])
      }
    } else {
      mesg("Finding nonparametric adjustments.")
      for (i in 1:n.batch) {
        temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                           gamma.hat[i, ], delta.hat[i, ])
        gamma.star <- rbind(gamma.star, temp[1, ])
        delta.star <- rbind(delta.star, temp[2, ])
      }
    }
    mesg("Adjusting the Data.")
    bayesdata <- s.data
    j <- 1
    for (i in batches) {
      bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, ] %*% gamma.star)) /
        (sqrt(delta.star[j, ]) %*% t(rep(1, n.batches[j])))
      j <- j + 1
    }
    bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, n.array)))) + stand.mean
    return(bayesdata)
  }
}

#' There are some funky scoping problems in isva::DoISVA().
#'
#' Thus I copy/pasted the function and attempted to address them here.
#'
#' @param data.m Input matrix.
#' @param pheno.v Vector of conditions of interest in the data.
#' @param cf.m Matrix of confounded conditions in the data.
#' @param factor.log I forget.
#' @param pvthCF Minimal p-value for considering.
#' @param th threshold for inclusion.
#' @param ncomp Number of SVA components to estimate.
#' @param icamethod Which ICA implementation to use?
#' @seealso [isva]
my_isva <- function(data.m, pheno.v, cf.m = NULL, factor.log = FALSE, pvthCF = 0.01,
                    th = 0.05, ncomp = NULL, icamethod = "fastICA") {
  isva.o <- isva::isvaFn(data.m, pheno.v, ncomp, icamethod)
  ## The default values of selisv.idx and pv.m
  selisv.idx <- 1:ncol(isva.o[["isv"]])
  pv.m <- NULL
  ## These change if cf.m is set.
  if (!is.null(cf.m)) {
    tmp.m <- cbind(pheno.v, cf.m)
    treatfactor <- c(FALSE, factor.log)
    pv.m <- matrix(nrow = ncol(isva.o$isv), ncol = 1 + ncol(cf.m))
    colnames(pv.m) <- c("POI", colnames(cf.m))
    for (c in 1:ncol(tmp.m)) {
      if (treatfactor[c] == FALSE) {
        for (sv in 1:ncol(isva.o[["isv"]])) {
          lm.o <- lm(isva.o[["isv"]][, sv] ~ as.numeric(tmp.m[, c]))
          pv.m[sv, c] <- summary(lm.o)[["coefficients"]][2, 4]
        }
      } else {
        for (sv in 1:ncol(isva.o$isv)) {
          lm.o <- lm(
              isva.o$isv[, sv] ~ as.factor(tmp.m[, c]))
          pv.m[sv, c] <- stats::pf(summary(lm.o)$fstat[1], summary(lm.o)$fstat[2],
                                   summary(lm.o)$fstat[3], lower.tail = FALSE)
        }
      }
    }
    selisv.idx <- vector()
    for (sv in 1:nrow(pv.m)) {
      ncf <- length(which(pv.m[sv, 2:ncol(pv.m)] < pvthCF))
      minpv <- min(pv.m[sv, 2:ncol(pv.m)])
      phpv <- pv.m[sv, 1]
      if (ncf > 0) {
        if (minpv < phpv) {
          selisv.idx <- c(selisv.idx, sv)
        }
      }
    }
    if (length(selisv.idx) == 0) {
      mesg("No ISVs selected because none correlated with the given confounders.
Rerun ISVA with cf.m = NULL option")
      stop()
    }
  }

  mesg("Running final multivariate regressions with selected ISVs")
  selisv.m <- matrix(isva.o[["isv"]][, selisv.idx], ncol = length(selisv.idx))
  ## print(selisv.m)
  mod <- model.matrix(~ pheno.v + selisv.m)
  modNULL <- model.matrix(~ selisv.m)
  isv_fit <- limma::lmFit(data.m, design = modNULL)
  new_counts <- limma::residuals.MArrayLM(isv_fit, data.m)
  df1 <- dim(mod)[2]
  df0 <- dim(modNULL)[2]
  pv.v <- rep(0, nrow(data.m))
  Id <- diag(ncol(data.m))
  resid <- data.m %*%
    (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  rss1 <- rowSums(resid * resid)
  residNULL <- data.m %*%
    (Id - modNULL %*% solve(t(modNULL) %*% modNULL) %*% t(modNULL))
  rssNULL <- rowSums(residNULL * residNULL)
  fstats <- ((rssNULL - rss1) / (df1 - df0)) /
    (rss1 / (ncol(data.m) - df1))
  pv.v <- 1 - isva::pf(fstats,
                       df1=(df1 - df0),
                       df2=(ncol(data.m) - df1))
  pv.s <- sort(pv.v, decreasing = FALSE, index.return = TRUE)
  qv.v <- p.adjust(pv.s[["x"]])
  ##qv.v <- qvalue(pv.s[["x"]])[["qvalue"]]
  ntop <- length(which(qv.v < th))
  sig_mtrx <- as.matrix(data.m[pred.idx, ])
  mesg("Number of DEGs after ISV adjustment = ", ntop)
  if (ntop > 0) {
    pred.idx <- pv.s[["ix"]][1:ntop]
    lm.o <- lm(t(data.m) ~ pheno.v + selisv.m)
    lm_sum <- summary(lm.o)
    tstats.v <- unlist(lapply(lm_sum,
                              function(x) {
                                x[["coefficients"]][2, 3]
                              }))
    lm.m <- cbind(tstats.v, pv.s$x, qv.v)
    colnames(lm.m) <- c("t-stat", "P-value", "q-value")
  }
  retlist <- list(
      "spv" = pv.s[["x"]],
      "qv" = qv.v,
      "rk" = pv.s[["ix"]],
      "ndeg" = ntop,
      "deg" = pred.idx,
      "lm" = lm.m,
      "isv" = selisv.m,
      "nsv" = length(selisv.idx),
      "pvCF" = pv.m,
      "selisv" = selisv.idx)
  return(retlist)
  ##  return(list(spv = pv.s$x, qv = qv.v, rk = pv.s$ix, ndeg = ntop,
  ##              deg = pred.idx, lm = lm.m, isv = selisv.m, nsv = length(selisv.idx),
  ##              pvCF = pv.m, selisv = selisv.idx))
}

## EOF
