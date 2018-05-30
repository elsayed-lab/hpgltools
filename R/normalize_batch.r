#' A function suggested by Hector Corrada Bravo and Kwame Okrah for batch removal
#'
#' During a lab meeting, the following function was suggested as a quick and
#' dirty batch removal tool
#'
#' @param normalized_counts Data frame of log2cpm counts.
#' @param model Balanced experimental model containing condition and batch factors.
#' @return Dataframe of residuals after subtracting batch from the model.
#' @seealso \pkg{limma}
#'  \code{\link[limma]{voom}} \code{\link[limma]{lmFit}}
#' @examples
#' \dontrun{
#'  newdata <- cbcb_batch_effect(counts, expt_model)
#' }
#' @export
cbcb_batch_effect <- function(normalized_counts, model) {
  ## model = model.matrix(~ condition + batch)
  voomed <- hpgl_voom(normalized_counts, model)
  voomed_fit <- limma::lmFit(voomed)
  modified_model <- model
  ## Drop batch from the model
  modified_model <- modified_model[, grep("batch", colnames(modified_model))] <- 0
  new_data <- tcrossprod(voomed_fit$coefficient, modified_model) +
    residuals(voomed_fit, normalized_counts)
  return(new_data)
}

#' Perform different batch corrections using limma, sva, ruvg, and cbcbSEQ.
#'
#' I found this note which is the clearest explanation of what happens with batch effect data:
#' https://support.bioconductor.org/p/76099/
#' Just to be clear, there's an important difference between removing a batch effect and modelling a
#' batch effect. Including the batch in your design formula will model the batch effect in the
#' regression step, which means that the raw data are not modified (so the batch effect is not
#' removed), but instead the regression will estimate the size of the batch effect and subtract it
#' out when performing all other tests. In addition, the model's residual degrees of freedom will
#' be reduced appropriately to reflect the fact that some degrees of freedom were "spent"
#' modelling the batch effects. This is the preferred approach for any method that is capable of
#' using it (this includes DESeq2). You would only remove the batch effect (e.g. using limma's
#' removeBatchEffect function) if you were going to do some kind of downstream analysis that can't
#' model the batch effects, such as training a classifier.
#' I don't have experience with ComBat, but I would expect that you run it on log-transformed CPM
#' values, while DESeq2 expects raw counts as input. I couldn't tell you how to properly use the
#' two methods together.
#'
#' @param count_table  Matrix of (pseudo)counts.
#' @param design  Model matrix defining the experimental conditions/batches/etc.
#' @param batch  String describing the method to try to remove the batch effect
#'  (or FALSE to leave it alone, TRUE uses limma).
#' @param expt_state  Current state of the expt in an attempt to avoid double-normalization.
#' @param batch1  Column in the design table describing the presumed covariant to remove.
#' @param batch2  Column in the design table describing the second covariant to remove
#'  (only used by limma at the moment).
#' @param noscale  Used for combatmod, when true it removes the scaling parameter
#'  from the invocation of the modified combat.
#' @param ...  More options for you!
#' @return The 'batch corrected' count table and new library size.  Please remember that the
#'  library size which comes out of this may not be what you want for voom/limma and would
#'  therefore lead to spurious differential expression values.
#' @seealso \pkg{limma} \pkg{edgeR} \pkg{RUVSeq} \pkg{sva} \pkg{cbcbSEQ}
#' @examples
#' \dontrun{
#'  limma_batch <- batch_counts(table, design, batch1='batch', batch2='strain')
#'  sva_batch <- batch_counts(table, design, batch='sva')
#' }
#' @export
batch_counts <- function(count_table, design, batch=TRUE, batch1="batch", expt_state=NULL,
                         batch2=NULL, noscale=TRUE, ...) {
  arglist <- list(...)
  low_to_zero <- FALSE
  if (!is.null(arglist[["low_to_zero"]])) {
    low_to_zero <- arglist[["low_to_zero"]]
  }
  num_surrogates <- NULL
  surrogate_method <- NULL
  if (is.null(arglist[["num_surrogates"]]) & is.null(arglist[["surrogate_method"]])) {
    surrogate_method <- "be"
  } else if (!is.null(arglist[["num_surrogates"]])) {
    if (class(arglist[["num_surrogates"]]) == "character") {
      surrogate_method <- arglist[["num_surrogates"]]
    } else {
      num_surrogates <- arglist[["num_surrogates"]]
    }
  } else if (!is.null(arglist[["surrogate_method"]])) {
    if (class(arglist[["surrogate_method"]]) == "numeric") {
      num_surrogates <- arglist[["surrogate_method"]]
    } else {
      surrogate_method <- arglist[["surrogate_method"]]
    }
  } else {
    warning("Both num_surrogates and surrogate_method were defined.
This will choose the number of surrogates differently depending on method chosen.")
  }
  message("In norm_batch, after testing logic of surrogate method/number, the
number of surrogates is: ", num_surrogates, " and the method is: ", surrogate_method, ".")

  cpus <- 4
  if (!is.null(arglist[["cpus"]])) {
    cpus <- arglist[["cpus"]]
  }
  prior.plots <- FALSE
  if (!is.null(arglist[["prior.plots"]])) {
    message("When using ComBat, using prior.plots may result in an error due to infinite ylim.")
    prior.plots <- arglist[["prior.plots"]]
  }

  ## Lets use expt_state to make sure we know if the data is already log2/cpm/whatever.
  ## We want to use this to back-convert or reconvert data to the appropriate scale on return.
  if (is.null(expt_state)) {
    expt_state <- list(
      "filter" = "raw",
      "normalization" = "raw",
      "conversion" = "raw",
      "batch" = "raw",
      "transform" = "raw")
  }
  ## Use current_state to keep track of changes made on scale/etc during batch correction
  ## This is pointed directly at limmaresid for the moment, which converts to log2
  current_state <- expt_state
  ## These droplevels calls are required to avoid errors like 'confounded by batch'
  batches <- droplevels(as.factor(design[[batch1]]))
  conditions <- droplevels(as.factor(design[["condition"]]))

  message("Note to self:  If you get an error like 'x contains missing values'; I think this
 means that the data has too many 0's and needs to have a better low-count filter applied.")

  num_low <- sum(count_table < 1 & count_table > 0)
  if (is.null(num_low)) {
    num_low <- 0
  }
  if (num_low > 0) {
    message("batch_counts: Before batch correction, ", num_low, " entries 0<x<1.")
  }
  num_zero <- sum(count_table <= 0)
  if (num_zero > 0) {
    message("batch_counts: Before batch correction, ", num_zero, " entries are >= 0.")
  }
  if (isTRUE(batch)) {
    batch <- "limma"
  }

  count_df <- data.frame(count_table)
  count_mtrx <- as.matrix(count_df)
  conditional_model <- model.matrix(~conditions, data=count_df)
  null_model <- conditional_model[, 1]
  ## Set the number of surrogates for sva/ruv based methods.
  if (!is.null(surrogate_method)) {
    num_surrogates <- sm(sva::num.sv(count_mtrx, conditional_model, method=surrogate_method))
  }
  if (num_surrogates < 1) {
    message("0 surrogates were detected by the ", surrogate_method, " method.")
    message("This will end badly, so setting num_surrogates to 1.")
    num_surrogates <- 1
  }
  message("After checking/setting the number of surrogates, it is: ", num_surrogates, ".")

  switchret <- switch(
    batch,
    "limma" = {
      if (expt_state[["transform"]] == "raw") {
        count_table <- log2(count_table + 1)
      }
      if (is.null(batch2)) {
        ## A reminder of removeBatchEffect usage
        ## adjusted_batchdonor = removeBatchEffect(data,
        ##                             batch=as.factor(as.character(des$donor)),
        ##                             batch2=as.factor(as.character(des$batch)))
        message("batch_counts: Using limma's removeBatchEffect to remove batch effect.")
        ## limma's removeBatchEffect() function does in fact take and return log2.
        ## Thus we need to check the transform state and make sure the input and output
        ## is base 10.
        count_table <- limma::removeBatchEffect(count_table, batch=batches)
      } else {
        batches2 <- as.factor(design[[batch2]])
        count_table <- limma::removeBatchEffect(count_table, batch=batches, batch2=batches2)
      }
      message(strwrap(prefix=" ", initial="", "If you receive a warning: 'NANs produced', one
 potential reason is that the data was quantile normalized."))
      if (expt_state[["transform"]] == "raw") {
        count_table <- (2 ^ count_table) - 1
      }
    },
    "limmaresid" = {
      ## ok a caveat:  voom really does require input on the base 10 scale and returns
      ## log2 scale data.  Therefore we need to make sure that the input is provided appropriately.
      message("batch_counts: Using residuals of limma's lmfit to remove batch effect.")
      batch_model <- model.matrix(~batches)
      if (expt_state[["transform"]] == "log2") {
        ## For now just assume only log2/base10.
        count_table <- (2 ^ count_table) - 1
      }
      batch_voom <- NULL
      if (expt_state[["normalization"]] == "raw") {
        batch_voom <- limma::voom(count_df, batch_model,
                                  normalize.method="quantile",
                                  plot=FALSE)
      } else {
        batch_voom <- limma::voom(count_df, batch_model, plot=FALSE)
      }
      batch_fit <- limma::lmFit(batch_voom, design=batch_model)
      ## count_table <- residuals(batch_fit, batch_voom[["E"]])
      ## This is still fubar!
      count_table <- limma::residuals.MArrayLM(batch_fit, batch_voom)
      ## Make sure to change this soon to take into account whether we are working on the log
      ## or non-log scale. Perhaps switch out the call from limma::voom to my own voom -- though
      ## I think I would prefer to use their copy and have a check that way if they change
      ## something important I will pick up on it.
    },
    "combatmod" = {
      ## normalized_data = hpgl_combatMod(dat=data.frame(counts), batch=batches,
      ##                                  mod=conditions, noScale=noscale, ...)
      message("batch_counts: Using a modified cbcbSEQ combatMod for batch correction.")
      count_table <- hpgl_combatMod(dat=count_df, batch=batches,
                                    mod=conditions, noScale=noscale, ...)
    },
    "fsva" = {
      message("batch_counts: Using sva::fsva for batch correction.")
      sva_object <- sm(sva::sva(count_mtrx,
                                conditional_model,
                                null_model,
                                n.sv=num_surrogates))
      fsva_result <- sm(sva::fsva(count_mtrx,
                                  conditional_model,
                                  sva_object,
                                  newdat=as.matrix(count_mtrx),
                                  method="exact"))
      count_table <- fsva_result[["db"]]
    },
    "combat" = ,  ## This peculiar syntax should match combat and combat_noscale to the same result.
    "combat_noscale" = {
      message("batch_counts: Using sva::combat with a prior for batch correction and no scaling.")
      count_table <- sm(sva::ComBat(count_table, batches, mod=NULL,
                                    par.prior=TRUE, prior.plots=prior.plots, mean.only=TRUE))
    },
    "combat_noprior" = {
      message("batch_counts: Using sva::combat without a prior for batch correction and no scaling.")
      message("This takes a long time!")
      count_table <- sm(sva::ComBat(count_table, batches, mod=conditions,
                                    par.prior=FALSE, prior.plots=prior.plots, mean.only=TRUE))
    },
    "combat_scale" = {
      message("batch_counts: Using sva::combat with a prior for batch correction and with scaling.")
      count_table <- sm(sva::ComBat(count_table, batches, mod=conditions,
                                    par.prior=TRUE, prior.plots=prior.plots, mean.only=FALSE))
    },
    "combat_noprior_scale" = {
      message("batch_counts: Using sva::combat without a prior for batch correction and with scaling.")
      count_table <- sm(sva::ComBat(count_table, batches, mod=conditions,
                                    par.prior=FALSE, prior.plots=prior.plots, mean.only=FALSE))
    },
    "svaseq" = {
      message("batch_counts: Using sva::svaseq for batch correction.")
      message("Note to self:  If you feed svaseq a data frame you will get an error like:")
      message("data %*% (Id - mod %*% blah blah requires numeric/complex arguments.")
      svaseq_result <- sm(sva::svaseq(count_mtrx, conditional_model, null_model, n.sv=num_surrogates))
      count_table <- counts_from_surrogates(count_mtrx, svaseq_result[["sv"]], design=design)
    },
    "varpart" = {
      message("Taking residuals from a linear mixed model as suggested by the variancePartition package.")
      cl <- parallel::makeCluster(cpus)
      doParallel::registerDoParallel(cl)
      batch_model <- as.formula("~ (1|batch)")
      message("The function fitvarPartModel may take excessive memory, you have been warned.")
      batch_fit <- variancePartition::fitVarPartModel(count_df, batch_model, design)
      count_table <- residuals(batch_fit)
      rm(batch_fit)
      parallel::stopCluster(cl)
    },
    "ruvg" = {
      message("Using RUVSeq and edgeR for batch correction (similar to lmfit residuals.)")
      ## Adapted from: http://jtleek.com/svaseq/simulateData.html -- but not quite correct yet
      ruv_input <- edgeR::DGEList(counts=count_df, group=conditions)
      ruv_input_norm <- ruv_input
      if (expt_state[["normalization"]] == "raw") {
        ruv_input_norm <- edgeR::calcNormFactors(ruv_input, method="upperquartile")
      }
      ruv_input_glm <- edgeR::estimateGLMCommonDisp(ruv_input_norm, conditional_model)
      ruv_input_tag <- edgeR::estimateGLMTagwiseDisp(ruv_input_glm, conditional_model)
      ruv_fit <- edgeR::glmFit(ruv_input_tag, conditional_model)
      ## Use RUVSeq with empirical controls
      ## The previous instance of ruv_input should work here, and the ruv_input_norm
      ## Ditto for _glm and _tag, and indeed ruv_fit
      ## Thus repeat the first 7 lines of the previous RUVSeq before anything changes.
      ruv_lrt <- edgeR::glmLRT(ruv_fit, coef=2)
      ruv_control_table <- ruv_lrt[["table"]]
      ranked <- as.numeric(rank(ruv_control_table[["LR"]]))
      bottom_third <- (summary(ranked)[[2]] + summary(ranked)[[3]]) / 2
      ruv_controls <- ranked <= bottom_third  ## what is going on here?!
      ## ruv_controls = rank(ruv_control_table$LR) <= 400  ## some data sets fail with 400 hard-set
      ruv_result <- RUVSeq::RUVg(count_mtrx, ruv_controls, k=num_surrogates)
      count_table <- ruv_result[["normalizedCounts"]]
    },
    {
      message("Passing the batch method to get_model_adjust().")
      message("It understands a few additional batch methods.")
      surrogate_result <- try(get_model_adjust(
        count_table, design=design, estimate_type=batch, surrogates=num_surrogates,  ...))
      if (class(surrogate_result) != "try-error") {
        count_table <- surrogate_result[["new_counts"]]
      }
    }
  )  ## Ending the switch statement.
  count_table[is.na(count_table)] <- 0
  num_low <- sum(count_table <= 0)
  if (is.null(num_low)) {
    num_low <- 0
  }
  if (num_low > 0) {
    message("The number of elements which are < 0 after batch correction is: ", num_low)
    message("The variable low_to_zero sets whether to change <0 values to 0 and is: ", low_to_zero)
    if (isTRUE(low_to_zero)) {
      count_table[count_table < 0] <- 0
    }
  }
  libsize <- colSums(count_table)
  counts <- list(count_table=count_table, libsize=libsize)
  return(counts)
}

#' A single place to extract count tables from a set of surrogate variables.
#'
#' Given an initial set of counts and a series of surrogates, what would the resulting count table
#' look like? Hopefully this function answers that question.
#'
#' @param data  Original count table, may be an expt/expressionset or df/matrix.
#' @param adjust  Surrogates with which to adjust the data.
#' @param design  Experimental design if it is not included in the expressionset.
#' @return A data frame of adjusted counts.
#' @seealso \pkg{Biobase}
#' @export
counts_from_surrogates <- function(data, adjust, design=NULL) {
  data_mtrx <- NULL
  my_design <- NULL
  if (class(data) == "expt") {
    my_design <- pData(data)
    conditions <- droplevels(as.factor(pData(data)[["condition"]]))
    data_mtrx <- exprs(data)
  } else if (class(data) == "ExpressionSet") {
    my_design <- pData(data)
    conditions <- droplevels(as.factor(pData(data)[["condition"]]))
    data_mtrx <- exprs(data)
  } else {
    my_design <- design
    conditions <- droplevels(as.factor(design[["condition"]]))
    data_mtrx <- as.matrix(data)
  }
  conditional_model <- model.matrix(~ conditions, data=my_design)

  new_model <- conditional_model
  ## Explicitly append columns of the adjust matrix to the conditional model.
  ## In the previous code, this was: 'X <- cbind(conditional_model, sva$sv)'
  ## new_model <- cbind(conditional_model, adjust)
  new_colnames <- colnames(conditional_model)
  adjust_mtrx <- as.matrix(adjust)
  for (col in 1:ncol(adjust_mtrx)) {
    new_model <- cbind(new_model, adjust_mtrx[, col])
    new_colname <- paste0("sv", col)
    new_colnames <- append(new_colnames, new_colname)
  }
  colnames(new_model) <- new_colnames

  ##data_modifier <- try(solve(t(new_model) %*% new_model) %*% t(new_model))
  ## In the previous code, this was: 'Hat <- solve(t(X) %*% X) %*% t(X)'
  ## Now it is in two separate lines, first the solve operation:
  data_solve <- try(solve(t(new_model) %*% new_model), silent=TRUE)
  if (class(data_solve) == "try-error") {
    message("Data modification by the model failed.")
    message("Leaving counts untouched.")
    return(data_mtrx)
  }
  ## If the solve operation passes, then the '%*% t(X)' is allowed to happen.
  data_modifier <- data_solve %*% t(new_model)
  transformation <- (data_modifier %*% t(data_mtrx))
  conds <- ncol(conditional_model)
  new_counts <- data_mtrx - t(as.matrix(new_model[, -c(1:conds)]) %*%
                              transformation[-c(1:conds), ])
  return(new_counts)
}

#' A modified version of comBatMod.
#'
#' This is a hack of Kwame Okrah's combatMod to make it not fail on corner-cases.
#' This was mostly copy/pasted from https://github.com/kokrah/cbcbSEQ/blob/master/R/transform.R
#'
#' @param dat Df to modify.
#' @param batch Factor of batches.
#' @param mod Factor of conditions.
#' @param noScale The normal 'scale' option squishes the data too much, so this defaults to TRUE.
#' @param prior.plots Print out prior plots?
#' @param ... Extra options are passed to arglist
#' @return Df of batch corrected data
#' @seealso \pkg{sva}
#'  \code{\link[sva]{ComBat}}
#' @examples
#' \dontrun{
#'  df_new = hpgl_combatMod(df, batches, model)
#' }
#' @export
hpgl_combatMod <- function(dat, batch, mod, noScale=TRUE, prior.plots=FALSE, ...) {
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
  message("Standardizing data across genes\n")
  if (NAs) {
    warning(paste0("Found ", sum(is.na(dat)), " missing data values."))
    warning("The original combatMod uses an undefined variable Beta.NA here,
I set it to 1 not knowing what its purpose is.")
    B.hat <- apply(dat, 1, Beta.NA)
  } else {
    ## There are no NAs in the data, this is a good thing(Tm)!
    B.hat <- solve(t(design) %*% design) %*% t(design) %*% t(as.matrix(dat))
  }
  grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]

  if (NAs) {
    var.pooled <- apply(dat - t(design %*% B.hat), 1, var, na.rm=TRUE)
  }
  else {
    var.pooled <- ((dat - t(design %*% B.hat)) ^ 2) %*% rep(1 / n.array, n.array)
  }
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1, n.array)))
  if (noScale) {
    m.data <- dat - stand.mean
    mse <- ((dat - t(design %*% B.hat)) ^ 2) %*% rep(1 / (n.array - ncol(design)), n.array)
    hld <- NULL
    bayesdata <- dat
    for (k in 1:n.batch) {
      message("Fitting 'shrunk' batch ", k, " effects.")
      sel <- batches[[k]]
      gammaMLE <- rowMeans(m.data[, sel])
      mprior <- mean(gammaMLE, na.rm = TRUE)
      vprior <- var(gammaMLE, na.rm = TRUE)
      prop <- vprior / (mse / (length(sel)) + vprior)
      gammaPost <- prop * gammaMLE + (1 - prop) * mprior
      for (i in sel) {
        bayesdata[, i] <- bayesdata[, i] - gammaPost
      }
      stats <- data.frame(gammaPost=gammaPost, gammaMLE=gammaMLE, prop=prop)
      hld[[paste("Batch", k, sep=".")]] <- list(
        "stats" = stats,
        "indices" = sel,
        "mprior" = mprior,
        "vprior" = vprior)
    }
    message("Adjusting data for batch effects.")
    return(bayesdata)
  } else {
    message("Fitting L/S model and finding priors.")
    batch.design <- design[, 1:n.batch]
    if (NAs) {
      gamma.hat <- apply(s.data, 1, Beta.NA, batch.design)
    } else {
      gamma.hat <- solve(t(batch.design) %*% batch.design) %*% t(batch.design) %*% t(as.matrix(s.data))
    }
    delta.hat <- NULL
    for (i in batches) {
      delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var, na.rm=TRUE))
    }
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, sva:::aprior)
    b.prior <- apply(delta.hat, 1, sva:::bprior)
    if (prior.plots & par.prior) {
      oldpar <- par(mfrow = c(2, 2))
      tmp <- density(gamma.hat[1, ])
      plot(tmp, type="l", main="Density Plot")
      xx <- seq(min(tmp$x), max(tmp$x), length = 100)
      lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
      stats::qqnorm(gamma.hat[1, ])
      stats::qqline(gamma.hat[1, ], col = 2)
      tmp <- stats::density(delta.hat[1, ])
      invgam <- 1 / stats::rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
      tmp1 <- stats::density(invgam)
      plot(tmp, typ="l", main="Density Plot", ylim=c(0, max(tmp$y, tmp1$y)))
      lines(tmp1, col = 2)
      stats::qqplot(delta.hat[1, ], invgam, xlab="Sample Quantiles", ylab="Theoretical Quantiles")
      lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
      title("Q-Q Plot")
      newpar <- par(oldpar)
    }
    gamma.star <- delta.star <- NULL
    if (par.prior) {
      message("Finding parametric adjustments.")
      for (i in 1:n.batch) {
        temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                             delta.hat[i, ], gamma.bar[i],
                             t2[i], a.prior[i], b.prior[i])
        gamma.star <- rbind(gamma.star, temp[1, ])
        delta.star <- rbind(delta.star, temp[2, ])
      }
    } else {
      message("Finding nonparametric adjustments.")
      for (i in 1:n.batch) {
        temp <- sva:::int.eprior(as.matrix(s.data[, batches[[i]]]), gamma.hat[i, ], delta.hat[i, ])
        gamma.star <- rbind(gamma.star, temp[1, ])
        delta.star <- rbind(delta.star, temp[2, ])
      }
    }
    message("Adjusting the Data.")
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

## EOF
