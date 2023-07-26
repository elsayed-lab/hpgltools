## Some functions to help with classification

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
                              p = 0.6, list = FALSE, times = 5) {
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
    "tester_outcomes" = tester_outcomes)
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
#' @param which Choose a paritiont to evaluate
#' @param type Use the training or testing data?
#' @export
self_evaluate_model <- function(predictions, datasets, which = 1, type = "train") {
  stripped <- data.frame()
  idx <- numeric()
  outcomes <- factor()
  if (type == "train") {
    stripped <- datasets[["trainers_stripped"]][[which]]
    idx <- datasets[["train_idx"]][[which]]
    outcomes <- datasets[["trainer_outcomes"]][[which]]
  } else {
    stripped <- datasets[["testers_stripped"]][[which]]
    idx <- datasets[["test_idx"]][[which]]
    outcomes <- datasets[["tester_outcomes"]][[which]]
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
    predict_df <- predict_df %>%
      dplyr::mutate('class' = names(.)[apply(., 1, which.max)])
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
    "roc" = roc,
    "roc_plot" = roc_record,
    "auc" = auc)
  class(retlist) <- "classifier_evaluation"
  return(retlist)
}
