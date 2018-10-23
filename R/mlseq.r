
simple_mlseq <- function(expt, comparison="condition", number_by_var=100,
                         ceiling_factor=1/3, training_number=2, training_repeats=2,
                         training_method="repeatedcv",
                         classify_method="svmRadial", classify_preprocess="deseq-rlog",
                         reference_factor=NULL, ...) {
  arglist <- list(...)
  tune_length <- 20
  if (!is.null(arglist[["tune_length"]])) {
    tune_length <- arglist[["tune_length"]]
  }

  expt <- make_pombe_expt()
  expt <- subset_expt(expt, subset="minute=='0'|minute=='120'")
  metadata <- pData(expt)
  expression <- exprs(expt)

  ## library(MLSeq)
  ## filepath <- system.file("extdata/cervical.txt", package = "MLSeq")
  ## expression <- read.table(filepath, header=TRUE)
  ## metadata <- S4Vectors::DataFrame(condition = factor(rep(c("N","T"), c(29, 29))))

  ## Simplify the experimental design to a single-column data frame
  ## with non-numeric factor elements.
  classification <- S4Vectors::DataFrame(
                                 "query" = as.factor(paste0("q", metadata[[comparison]])))

  ## Get the top-n genes by variance.
  variance_by_genes <- sort(apply(expression, 1, BiocGenerics::var, na.rm=TRUE), decreasing=TRUE)
  topn_variance_genes <- head(names(variance_by_genes), n=number_by_var)
  start_data <- expression[topn_variance_genes, ]
  num_to_test <- ceiling(ncol(start_data) * ceiling_factor)

  ## Extract a random slice of the columns with a caveat: there must be some columns
  ## from all factor elements in the original data.  This loop will resample either until
  ## that constraint is met, or after 20 tries it will give up and throw an error.
  try_count <- 0
  good_slice <- 0
  start_levels <- levels(classification[["query"]])
  if (length(start_levels) == 1) {
    stop("There is nothing to compare.")
  }
  while (good_slice == 0) {
    try_count <- try_count + 1
    test_slice <- sample(ncol(start_data), num_to_test, FALSE)
    training_class <- S4Vectors::DataFrame("query" = classification[-test_slice, ])
    testing_class <- S4Vectors::DataFrame("query" = classification[test_slice, ])
    training_levels <- levels(training_class[["query"]])
    testing_levels <- levels(testing_class[["query"]])
    if (identical(training_levels, start_levels) &
        identical(testing_levels, start_levels)) {
      message("Got suitable training/testing columns after ", try_count, " tries.")
      good_slice <- 1
    } else if (try_count > 20) {
      stop("Tried 20 times to get samples with all levels.")
    }
  }

  ## Now we have some slices for training/testing along with our set of high-variance genes.
  ## At this point we want to set up to perform classification, this + 1 still makes no sense
  ## to me though.
  training_data <- as.matrix(start_data[, -test_slice] + 1)
  testing_data <- as.matrix(start_data[, test_slice] + 1)
  chosen_design <- formula(~query)
  ## Create DESeq data sets from the matrix subsets.
  training_deseq <- DESeq2::DESeqDataSetFromMatrix(
                              countData=training_data, colData=training_class,
                              design=chosen_design)
  testing_deseq <- DESeq2::DESeqDataSetFromMatrix(
                             countData=testing_data, colData=testing_class,
                             design=chosen_design)


  ## If not provided, set up a reference factor from the first element in the metadata.
  if (is.null(reference_factor)) {
    reference_factor <- as.character(training_levels[1])
    message("Choosing the first observed factor in the metadata as the reference: ",
            reference_factor, ".")
  } else {
    if (! reference_factor %in% training_levels) {
      stop("Could not find the reference factor: ", reference_factor,
           " in the set of available factors in the design: ", training_levels, ".")
    }
  }

  controllers <- list(
    "chosen" = caret::trainControl(method=training_method, number=training_number,
                                   repeats=training_repeats, classProbs=TRUE),
    "continuous" = caret::trainControl(method="repeatedcv", number=training_number,
                                       repeats=training_repeats),
    "discrete" = MLSeq::discreteControl(method="repeatedcv", number=training_number,
                                        repeats=training_repeats, tuneLength=tune_length),
    "voom_nsc" = MLSeq::voomControl(tuneLength=tune_length)
    )


  ## Set up a trainer.
  \training_control <- NULL
  if (grepl(pattern="voom", x=classify_method ) {
    training_control <-
  } else {
    training_control <- caret::trainControl(method=training_method, number=training_number,
                                            repeats=training_repeats, classProbs=TRUE)
  }


  ## I think I should make a list of the trainers/controllers/classifiers.
  deseq_fit <- MLSeq::classify(data=training_deseq, method=classify_method,
                               preProcessing=classify_preprocess,
                               ref=reference_factor, control=training_control)
  training_summary <- MLSeq::trained(deseq_fit)
  deseq_plot <- MLSeq::plot(deseq_fit)

  svm_control <- caret::trainControl(method="repeatedcv", number=5,
                                     repeats=10, classProbs=TRUE)
  svm_fit <- MLSeq::classify(data=training_deseq, method="svmRadial",
                             preProcessing="deseq-vst", ref=reference_factor,
                             tuneLength=10, control=svm_control)
  svm_plot <- MLSeq::plot(svm_fit)

  plda_control <- MLSeq::discreteControl(method="repeatedcv", number=5, repeats=1,
                                         tuneLength=tune_length)
  plda_fit <- MLSeq::classify(data=training_deseq, method="PLDA", normalize="deseq",
                              ref=reference_factor, control=plda_control)
  plda_plot <- MLSeq::plot(plda_fit)

  voom_control <- MLSeq::voomControl(method="repeatedcv", number=5, repeats=1,
                                     tuneLength=tune_length)
  voom_fit <- MLSeq::classify(data=training_deseq, method="voomDLDA",
                              normalize="deseq", ref=reference_factor, control=voom_control)

  deseq_predict <- relevel(MLSeq::predict(deseq_fit, testing_deseq), ref=reference_factor)
  svm_predict <- relevel(MLSeq::predict(svm_fit, testing_deseq), ref=reference_factor)
  plda_predict <- relevel(MLSeq::predict(plda_fit, testing_deseq), ref=reference_factor)
  voom_predict <- relevel(MLSeq::predict(voom_fit, testing_deseq), ref=reference_factor)

  actual_values <- relevel(testing_class[["query"]], ref=reference_factor)
  deseq_table <- table(Predicted=deseq_predict, Actual=actual_values)
  svm_table <- table(Predicted=svm_predict, Actual=actual_values)
  plda_table <- table(Predicted=plda_predict, Actual=actual_values)
  voom_table <- table(Predicted=voom_predict, Actual=actual_values)

  deseq_confused <- caret::confusionMatrix(deseq_table, positive=reference_factor)
  svm_confused <- caret::confusionMatrix(svm_table, positive=reference_factor)
  plda_confused <- caret::confusionMatrix(plda_table, positive=reference_factor)
  voom_confused <- caret::confusionMatrix(voom_table, positive=reference_factor)

  deseq_genes <- MLSeq::selectedGenes(deseq_fit)
  svm_genes <- MLSeq::selectedGenes(svm_fit)
  plda_genes <- MLSeq::selectedGenes(plda_fit)
  voom_genes <- MLSeq::selectedGenes(voom_fit)



  svm_fit <- MLSeq::classify(data=training_data, method="svmRadial",
                             preProcessing="deseq-vst", ref="T",
                             tuneLength=10, control=svm_control)
  svm_trained <- MLSeq::trained(svm_fit)
  show(svm_fit)
  plot(svm_fit)

  ctrl_svm <- caret::trainControl(method="repeatedcv", number=5, repeats=1)
  ctrl_plda <- MLSeq::discreteControl(method="repeatedcv", number=5,
                                      repeats=1, tuneLength=10)
  ctrl_voomDLDA <- MLSeq::voomControl(method = "repeatedcv", number = 5, repeats = 1,
                                      tuneLength = 10)
  ## Support vector machines with radial basis function kernel
  fit_svm <- MLSeq::classify(data=data_trainS4, method="svmRadial",
                             preProcessing="deseq-vst", ref="T", tuneLength=10,
                             control=ctrl_svm)
  ## Poisson linear discriminant analysis
  fit_plda <- MLSeq::classify(data=data_trainS4, method="PLDA", normalize="deseq",
                              ref="T", control=ctrl_plda)
  ## Voom-based diagonal linear discriminant analysis
  fit_voomDLDA <- MLSeq::classify(data=data_trainS4, method="voomDLDA",
                                  normalize="deseq", ref="T", control=ctrl_voomDLDA)
  voomdlda_training <- MLSeq::trained(fit_voomDLDA)

  pred_svm <- MLSeq::predict(fit_svm, data_testS4)
  pred_svm

  pred_svm <- relevel(pred_svm, ref="T")
  actual <- relevel(classts$condition, ref="T")
  tbl <- table(Predicted = pred_svm, Actual=actual)
  caret::confusionMatrix(tbl, positive="T")

  ## Define control lists.
  ctrl_continuous <- caret::trainControl(method="repeatedcv", number=5, repeats=10)
  ctrl_discrete <- MLSeq::discreteControl(method="repeatedcv", number=5, repeats=10,
                                   tuneLength=10)
  ctrl_voom <- MLSeq::voomControl(method="repeatedcv", number=5, repeats=10,
                                  tuneLength=10)
  ## 1. aContinuous classifiers, SVM and NSC
  fit_svm <- MLSeq::classify(data=data_trainS4, method="svmRadial",
                             preProcessing="deseq-vst", ref="T", tuneLength=10,
                             control=ctrl_continuous)
  fit_NSC <- MLSeq::classify(data=data_trainS4, method="pam",
                             preProcessing="deseq-vst", ref="T", tuneLength=10,
                             control=ctrl_continuous)
  ## 2. Discrete classifiers
  fit_plda <- MLSeq::classify(data=data_trainS4, method="PLDA", normalize="deseq",
                              ref="T", control=ctrl_discrete)

  fit_plda2 <- MLSeq::classify(data=data_trainS4, method="PLDA2", normalize="deseq",
                               ref="T", control=ctrl_discrete)
  fit_nblda <- MLSeq::classify(data=data_trainS4, method="NBLDA", normalize="deseq",
                               ref="T", control=ctrl_discrete)
  ## 3. voom-based classifiers
  fit_voomDLDA <- MLSeq::classify(data=data_trainS4, method="voomDLDA",
                                  normalize="deseq", ref="T", control=ctrl_voom)
  fit_voomNSC <- MLSeq::classify(data=data_trainS4, method="voomNSC",
                                 normalize="deseq", ref="T", control=ctrl_voom)
  ## 4. Predictions
  pred_svm <- MLSeq::predict(fit_svm, data_testS4)
  pred_NSC <- MLSeq::predict(fit_NSC, data_testS4)

  nblda_selected <- MLSeq::selectedGenes(fit_nblda)
  voom_selected <- MLSeq::selectedGenes(fit_voomNSC)
}



mlseq_example <- function() {
  library(S4Vectors)
  library(DESeq2)
  filepath <- system.file("extdata/cervical.txt", package = "MLSeq")
  cervical <- read.table(filepath, header=TRUE)
  class <- DataFrame(condition = factor(rep(c("N","T"), c(29, 29))))
  vars <- sort(apply(cervical, 1, var, na.rm = TRUE), decreasing = TRUE)
  data <- cervical[names(vars)[1:100], ]
  nTest <- ceiling(ncol(data) * 0.3)
  ind <- sample(ncol(data), nTest, FALSE)
                                        # Minimum count is set to 1 in order to prevent 0 division problem within
                                        # classification models.
  data.train <- as.matrix(data[ ,-ind] + 1)
  data.test <- as.matrix(data[ ,ind] + 1)
  classtr <- DataFrame(condition = class[-ind, ])
  classts <- DataFrame(condition = class[ind, ])
  data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,
                                        design = formula(~condition))
  data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,
                                       design = formula(~condition))
  fit <- classify(data = data.trainS4, method = "svmRadial",
                  preProcessing = "deseq-rlog", ref = "T",
                  control = trainControl(method = "repeatedcv", number = 2,
                                         repeats = 2, classProbs = TRUE))
  show(fit)

}
