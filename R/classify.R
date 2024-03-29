## Some functions to help with classification

#' Rerun a model generator and classifier on a training/testing set multiple times.
#'
#' @param full_df The matrix of preProcessed() data.
#' @param interesting_meta dataframe of metadata of potential interest.
#' @param outcome_column metadata column of interest.
#' @param p The proportion of training/testing
#' @param list How to return the partitions.
#' @param formula_string Optional formula string, otherwise genrated on thee fly.
#' @param run_times How many times to repeat this process
#' @param method Modelling method to employ with caret.
#' @param sampler Sampler to employ, bootstrap or cv right now.
#' @param sample_number How many times to use the sampler
#' @param tuner Tuning arguments for the method above.
#' @param state When provided, passes the state of the data to the return
#'  so it may be reported later.
#' @param ... Others, currently unused I think
#' @export
classify_n_times <- function(full_df, interesting_meta, outcome_column = "finaloutcome",
                             p = 0.4, list = FALSE, formula_string = NULL,
                             run_times = 10, method = "xgbTree",
                             sampler = "cv", sample_number = 10,
                             tuner = NULL, state = NULL, ...) {
  ## Note, the function declaration assumes outcome_factor is a single element character,
  ## but in my actual document it is a factor of 1 element per sample,
  ## e.g. 109 cures and 57 fails.  This should be addressed by create_partitions!()
  parameter_lst <- list(
    "proportion_trained" = p,
    "method_used" = method,
    "model_iterations" = run_times,
    "outcome_column" = outcome_column,
    "sampler_used" = sampler,
    "sampler_iterations" = sample_number)
  if (!is.null(state)) {
    parameter_lst <- append(state, parameter_lst)
  }

  if (! outcome_column %in% colnames(interesting_meta)) {
    stop("The outcome column: ", outcome_column, " is not in the metadata.")
  }
  outcome_factor <- as.factor(interesting_meta[[outcome_column]])
  partitions <- create_partitions(full_df, interesting_meta, outcome_factor = outcome_factor,
                                  p = p, list = list, times = run_times)

  if (is.null(formula_string)) {
    formula_string <- "outcome ~ ."
  }

  train_method <- NULL
  if (sampler == "cv") {
    sampling_method <- caret::trainControl(method = "cv", number = sample_number)
    parameter_lst[["return_resample"]] <- "unused"
  } else if (sampler == "bootstrap") {
    sampling_method <- caret::trainControl(method = "boot", number = sample_number,
                                    returnResamp = "all")
    parameter_lst[["return_resample"]] <- "all"
  } else {
    sampling_method <- sampler
  }

  trainer <- NULL
  train_args <- list()
  if (method == "xgbTree") {
    if (is.null(tuner)) {
      tuner <- data.frame(
        "nrounds" = 200,
        "eta" = c(0.05, 0.1, 0.3),
        "max_depth" = 4,
        "gamma" = 0,
        "colsample_bytree" = 1,
        "subsample" = 0.5,
        "min_child_weight" = 1)
    }
  } else if (method == "ranger") {
    if (is.null(tuner)) {
      tuner <- data.frame(
          mtry = 200,
          min.node.size = 1,
          splitrule = "gini")
      }
    train_args <- list(
      "importance" = "permutation")
    parameter_lst[["importance"]] <- "permutation"
  } else if (method == "glmnet") {
    if (is.null(tuner)) {
      tuner <- data.frame(
        alpha = 0.5,
        lambda = seq(0.1, 0.7, 0.05))
      train_args <- list(
        "family" = "binomial",
        "importance" = "permutation",
        "verbose" = TRUE)
      parameter_lst <- append(parameter_lst, train_args)
    }
  } else if (method == "knn") {
    if (is.null(tuner)) {
      tuner <- data.frame("k" = 1:10)
    }
  } else {
    stop("I do not yet know this method.")
  }

  for (key in names(tuner)) {
    new_key <- paste0("tuner_", key)
    if (length(tuner[[key]]) > 1) {
      parameter_lst[[new_key]] <- toString(tuner[[key]])
    } else {
      parameter_lst[[new_key]] <- tuner[[key]]
    }
  }

  ## I learned something important when setting this up,
  ## the trainer I am using has parallelism built in, so if I try and use doParallel,
  ## badness is likely to ensue.

  ## 20231129: FIXME
  ## !!!! Another _VERY_IMPORTANT_ note: do _NOT_ set verbose=TRUE!!!!!
  ## For reasons passing all understanding, it causes the trainer to fail.

  all_results <- list()
  test_eval_df <- data.frame()
  train_eval_df <- data.frame()
  for (d in seq_len(sample_number)) {
    message("Working on iteration: ", d, ".")
    ## This might require .packages = c("hpgltools" ...)
    train_all <- partitions[["trainers"]][[d]]
    train_df <- partitions[["trainers_stripped"]][[d]]
    train_idx <- partitions[["train_idx"]][[d]]
    train_outcomes <- partitions[["trainer_outcomes"]][[d]]
    test_df <- partitions[["testers"]][[d]]
    test_idx <- partitions[["test_idx"]][[d]]
    test_outcomes <- partitions[["tester_outcomes"]][[d]]

    all_train_args <- list(
      "form" = as.formula("outcome ~ ."),
      "data" = train_all,
      "method" = method,
      "trControl" = sampling_method,
      "tuneGrid" = tuner)
    for (a in seq_along(train_args)) {
      name <- names(train_args)[a]
      value <- train_args[[a]]
      all_train_args[[name]] <- value
    }

    trainer <- BiocGenerics::do.call(what = caret::train, args = all_train_args)
    trained <- stats::predict(trainer, train_df)
    message("Evaluating predictions.")
    trained_eval <- self_evaluate_model(trained, partitions,
                                        which_partition = d, type = "train")
    tested <- stats::predict(trainer, test_df)
    tested_eval <- self_evaluate_model(tested, partitions,
                                       which_partition = d, type = "test")
    this_result <- list(
      "trainer" = trainer,
      "trained" = trained,
      "trained_eval" = trained_eval,
      "tested" = tested,
      "tested_eval" = tested_eval)
    all_results[[d]] <- this_result
  }

  train_eval_df <- data.frame()
  test_eval_df <- data.frame()
  for (i in seq_along(all_results)) {
    train_eval <- all_results[[i]][["trained_eval"]]
    train_confused <- train_eval[["confusion_mtrx"]]
    ## Lets assume cure/fail for these categories.
    ## When the reference is cure and prediction is cure.
    true_positive <- train_confused[["table"]][1, 1]
    ## reference is fail and prediction is fail.
    true_negative <- train_confused[["table"]][2, 2]
    ## reference is fail and prediction is cure.
    false_positive <- train_confused[["table"]][1, 2]
    ## reference is cure and prediction is fail.
    false_negative <- train_confused[["table"]][2, 1]
    train_roc <- train_eval[["roc"]]
    truefalse <- c(true_positive, true_negative, false_positive, false_negative)
    names(truefalse) <- c("true_positive", "true_negative", "false_positive", "false_negative")
    train_numbers <- c(truefalse, train_confused[["overall"]], train_confused[["byClass"]])
    train_numbers <- c(train_numbers, train_eval[["auc"]])
    last <- length(train_numbers)
    names(train_numbers)[last] <- "auc"
    train_eval_df <- rbind(train_eval_df, train_numbers)
    colnames(train_eval_df) <- names(train_numbers)

    test_eval <- all_results[[i]][["tested_eval"]]
    test_confused <- test_eval[["confusion_mtrx"]]
    ## Lets assume cure/fail for these categories.
    ## When the reference is cure and prediction is cure.
    true_positive <- test_confused[["table"]][1, 1]
    ## reference is fail and prediction is fail.
    true_negative <- test_confused[["table"]][2, 2]
    ## reference is fail and prediction is cure.
    false_positive <- test_confused[["table"]][1, 2]
    ## reference is cure and prediction is fail.
    false_negative <- test_confused[["table"]][2, 1]
    test_roc <- test_eval[["roc"]]
    truefalse <- c(true_positive, true_negative, false_positive, false_negative)
    names(truefalse) <- c("true_positive", "true_negative", "false_positive", "false_negative")
    test_numbers <- c(truefalse, test_confused[["overall"]], test_confused[["byClass"]])
    test_numbers <- c(test_numbers, test_eval[["auc"]])
    last <- length(test_numbers)
    names(test_numbers)[last] <- "auc"
    test_eval_df <- rbind(test_eval_df, test_numbers)
    colnames(test_eval_df) <- names(test_numbers)
  }

  parameter_df <- t(as.data.frame(parameter_lst))
  colnames(parameter_df) <- c("Parameters Used")
  retlist <- list(
    "parameters" = parameter_df,
    "train_eval_summary" = train_eval_df,
    "test_eval_summary" = test_eval_df,
    "all_results" = all_results,
    "tuner" = tuner)
  class(retlist) <- "classified_n_times"
  return(retlist)
}

#' Use createDataPartition to create test/train sets and massage them a little.
#'
#' This will also do some massaging of the data to make it easier to work with for
#' downstream tasks.  Most notably, since I am mostly evaluating classifiers of
#' clinical data to see how well they agree with extant annotations, I want to make sure
#' the relevant columns are renamed in the testing sets.
#'
#' @param full_df Dataframe containing the measured data and relevant factors.
#' @param interesting_meta Other metadata (maybe not needed)
#' @param outcome_factor Name of the outcome column
#' @param p Ratio to split trainer and testers.
#' @param list Generate result as list or dataframe
#' @param times How many times to iterate
#' @seealso https://topepo.github.io/caret/data-splitting.html#simple-splitting-based-on-the-outcome
#'   and https://github.com/compgenomr/book/blob/master/05-supervisedLearning.Rmd
#' @export
create_partitions <- function(full_df, interesting_meta, outcome_factor = "condition",
                              p = 0.4, list = FALSE, times = 5) {
  if (length(outcome_factor) == 1) {
    outcome_fct <- as.factor(as.character(interesting_meta[[outcome_factor]]))
  } else {
    outcome_fct <- as.factor(outcome_factor)
  }
  training_mtrx <- caret::createDataPartition(outcome_fct, p = p,
                                              list = list, times = times)
  full_df <- as.data.frame(full_df)
  full_df[["outcome"]] <- outcome_fct
  cbind_df <- full_df %>%
    dplyr::select("outcome", tidyselect::everything())

  trainers <- list()
  trainers_stripped <- list()
  trainers_idx <- list()
  trainer_outcomes <- list()
  testers <- list()
  testers_idx <- list()
  tester_outcomes <- list()
  ## Each column of the training matrix is a set of training indices.
  ## I want extract from the cbind_df the relevant samples for them for
  ## testing/training and also move the outcome_factor to {outcome_factor}_bak
  ## for the testing data so that I can run predict but also find the data to
  ## create a ROC curve later.
  for (col in colnames(training_mtrx)) {
    train_idx <- training_mtrx[, col]
    train_rownames <- rownames(cbind_df)[train_idx]
    train_outcomes <- cbind_df[train_idx, "outcome"]
    names(train_outcomes) <- train_rownames

    test_idx <- ! rownames(cbind_df) %in% train_rownames
    test_rownames <- rownames(cbind_df)[test_idx]
    test_outcomes <- cbind_df[test_idx, "outcome"]
    names(test_outcomes) <- test_rownames

    train_df <- as.data.frame(cbind_df[train_rownames, ])
    test_df <- as.data.frame(cbind_df[test_rownames, ])

    trainers[[col]] <- train_df
    trainers_stripped[[col]] <- train_df
    trainers_stripped[[col]][["outcome"]] <- NULL
    trainers_idx[[col]] <- train_idx
    trainer_outcomes[[col]] <- train_outcomes
    ## Remove the outcome factor from the test data.
    test_df[["outcome"]] <- NULL
    testers[[col]] <- test_df
    testers_idx[[col]] <- test_idx
    tester_outcomes[[col]] <- test_outcomes


  }
  retlist <- list(
    "trainers" = trainers,
    "trainers_stripped" = trainers_stripped,
    "train_idx" = trainers_idx,
    "trainer_outcomes" = trainer_outcomes,
    "testers" = testers,
    "test_idx" = testers_idx,
    "tester_outcomes" = tester_outcomes,
    p = p,
    outcome_factor = outcome_factor,
    list = list,
    times = times)
  class(retlist) <- "partitioned_data"
  return(retlist)
}

#' Given an n-dimensional matrix, try some KNN-esque clustering on it.
#'
#' I want some functions to help me understand clustering.  This is a
#' first pass at that goal.
#'
#' @param mtrx Matrix to cluster, usually 2d from a point plot.
#' @param resolution Used after cluster generation for making neighbor
#'  groups.
#' @param k Used during cluster generation.
#' @param type Define the type of clustering to perform, currently
#'  only KNN/SNN
#' @param full Get the full set of metrics from bluster.
#' @param merge_to Use the neighborhood collapse function to set a
#'  hard ceiling on the number of clusters in the final result.
#' @param ... Extra args for bluster.
#' @return List containing the resulting groups and some information about them.
#' @export
generate_nn_groups <- function(mtrx, resolution = 1, k = 10, type = "snn",
                                full = TRUE, merge_to = NULL, ...) {
  params <- bluster::SNNGraphParam(k = k, ...)
  if (type == "knn") {
    params <- bluster::KNNGraphParam(k = k)
  }
  clusters <- bluster::clusterRows(mtrx, params, full = full)
  merged <- NULL
  groups <- as.factor(paste0("g", clusters$clusters))
  message("After clustering, there are: ", length(levels(groups)), " groups.")
  if (!is.null(merge_to)) {
    merged <- bluster::mergeCommunities(clusters[["objects"]][["graph"]],
                                        clusters[["clusters"]],
                                        steps = merge_to)
  }
  merged_groups <- as.factor(paste0("m", merged))
  message("After merging, there are: ", length(levels(merged_groups)), " groups.")

  retlist <- list(
    "clusters" = clusters,
    "merged" = merged,
    "start_groups" = groups,
    "merged_groups" = merged_groups)
  return(retlist)
}

#' Create a confusion matrix and ROC of a model against its training data. (and test data
#' if the annotations are known)
#'
#' This assumes a set of partitions from create_partitions() which
#' keeps the training metadata alongside the matrix of model
#' variables.  When available, that function also keeps the known
#' annotations of the testing data.  Given those annotations and the
#' model created/tested from them, this runs confusionMatrix and ROC,
#' collects the results, and provides them as a list.
#'
#' @param predictions Model created by train()
#' @param datasets Set of training/testing partitions along with
#'  associated metadata annotations.
#' @param which_partition Choose a paritiont to evaluate
#' @param type Use the training or testing data?
#' @export
self_evaluate_model <- function(predictions, datasets, which_partition = 1, type = "train") {
  stripped <- data.frame()
  idx <- numeric()
  outcomes <- factor()
  if (type == "train") {
    stripped <- datasets[["trainers_stripped"]][[which_partition]]
    idx <- datasets[["train_idx"]][[which_partition]]
    outcomes <- datasets[["trainer_outcomes"]][[which_partition]]
  } else {
    stripped <- datasets[["testers_stripped"]][[which_partition]]
    idx <- datasets[["test_idx"]][[which_partition]]
    outcomes <- datasets[["tester_outcomes"]][[which_partition]]
  }

  ## This assumes the input is a matrix of class probabilities.  First convert that to
  ## classes.
  predict_type <- "factor"
  predict_numeric <- NULL
  predict_df <- NULL
  if (class(predictions)[1] == "factor") {
    predict_class <- as.factor(predictions)
    names(predict_class) <- names(outcomes)
    predict_numeric <- as.numeric(predict_class)
  } else {
    predict_type <- "data.frame"
    predict_df <- as.data.frame(predictions)
    rownames(predict_df) <- names(outcomes)
    #predict_df <- predict_df %>%
    #  dplyr::mutate('class' = names(.)[apply(., 1, which.max)])
    ## I still do not fully understand the various implications of
    ## using . vs .data vs .env and when they are relevant for
    ## dplyr/magrittr.  As a result, this might be horribly wrong.
    predict_df <- predict_df %>%
      dplyr::mutate("class" = names(.data)[apply(.data, 1, which.max)])
    predict_class <- as.factor(predict_df[["class"]])
    names(predict_class) <- rownames(predict_df)
    predict_numeric <- predict_df[[1]]
  }

  confused <- caret::confusionMatrix(data = outcomes,
                                     reference = predict_class,
                                     mode = "everything")

  self_test <- predict_class == outcomes
  names(self_test) <- names(outcomes)
  wrong_sample_idx <- self_test == FALSE
  wrong_samples <- names(self_test)[wrong_sample_idx]
  self_summary <- summary(self_test)

  roc <- pROC::roc(response = outcomes, predictor = predict_numeric)
  roc_plot <- plot(roc)
  roc_record <- grDevices::recordPlot(roc_plot)
  auc <- pROC::auc(roc)
  retlist <- list(
    "self_test" = self_test,
    "self_summary" = self_summary,
    "wrong_samples" = wrong_samples,
    "confusion_mtrx" = confused,
    "reference" = predict_class,
    "roc" = roc,
    "roc_plot" = roc_record,
    "auc" = auc)
  class(retlist) <- "classifier_evaluation"
  return(retlist)
}

#' Write out the results of classify_n_times().
#'
#' @param result Ibid.
#' @param excel Output excel file
#' @param name Name of the sheet to write.
#' @export
write_classifier_summary <- function(result, excel = "ML_summary.xlsx", name = NULL) {
  ## FIXME: Make a generic and use a method here.
  if (is.null(name)) {
    name <- result[["parameters"]]["method_used", 1]
  }
  if (class(excel)[1] == "character" && !file.exists(excel)) {
    ## Then write a legend sheet
    legend <- data.frame(rbind(
      c("top_table", "Some information about the model generation."),
      c("bottom_table",
        "1 row for each iteration of the model followed by some meta summary statistics."),
      c("true_positive", "The raw number of true positives observed."),
      c("true_negative", "The raw number of true negatives observed."),
      c("false_positive", "The raw number of false positives observed."),
      c("false_negative", "The raw number of false negatives observed; the sum of these are the number of samples in the data."),
      c("accuracy", "Proportion of correct calls vs. all samples."),
      c("kappa", "Cohen's kappa: (accuracy - p(correct values vs. chance)) / (1 - p(correct values vs. chance))"),
      c("accuracy_lower", "Error rate if the model always chose the second class."),
      c("accuracy_upper", "Error rate if the model always chose the first class."),
      c("accuracy_null", "Error rate if the model always chose the majority class."),
      c("accuracy_pvalue", "Result of a binomial t-test to see if the accuracy is better than chance, taking into account imbalances in the samples/class."),
      c("mcnemar_pvalue", "McNemar's chi-squared test for row/column symmetry for 2 factors."),
      c("sensitivity", "(/ true_positive (+ true_positive false_positive))"),
      c("specificity", "(/ true_negative (+ true_negative false_negative))"),
      c("pos_pred_value", "(/ (* sensitivity prevelence) (+ (* sensitivity prevelance) (* (- specificity 1) (- prevalence 1))))"),
      c("neg_pred_value", "(/ (* specificity (- prevalence 1))  (+ (* (- specificity 1) * (- prevalence 1))  (* specificity (- prevalence 1))))"),
      c("precision", "(/ true_positive (+ true_positive true_negative))"),
      c("recall", "(/ true_positive (+ true_positive false_negative))"),
      c("f1", "2 * (/ 1 (+ (/ 1 precision) (/ 1 recall)))"),
      c("prevalence", "(/ (+ true_positive false_negative) (sum all))"),
      c("detection_rate", "(/ true_positive (sum all))"),
      c("detection_prevalence", "(/ (+ true_positive true_negative) (sum all))"),
      c("balanced_accuracy", "(/ (+ sensitivity specificity) 2)"),
      c("auc", "Area under the ROC curve.")))
    excel <- write_xlsx(data = legend, start_col = 1,
                          start_row = 1, sheet = "legend",
                          excel = excel)
  }
  summary_df <- rbind_summary_rows(result[["test_eval_summary"]])
  written <- write_xlsx(data = result[["parameters"]], sheet = name,
                        excel = excel, title = "parameters used")
  written <- write_xlsx(data = summary_df, start_col = 1,
                        start_row = written[["end_row"]] + 1,
                        excel = written)
  return(written)
}
