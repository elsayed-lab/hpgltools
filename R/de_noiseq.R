## NOISeq looks to me like it does not use a statistical model in the same way as
## limma/deseq/edger, but instead implements its own (un)guided batch correction method
## and expects the user to use that data as input for its set of comparisons.

## Thus, I think setting model_batch to TRUE is likely the place to use their default function
## 'ARSyNseq()' and if we use svaseq, simply send the modified counts to the various noiseq()
## functions which follow.

noiseq_pairwise <- function(input = NULL, conditions = NULL,
                            batches = NULL, model_cond = TRUE,
                            model_batch = TRUE, annot_df = NULL,
                            ...) {
  arglist <- list(...)

  message("Starting noiseq pairwise comparisons.")
  input <- sanitize_expt(input)
  input_data <- choose_binom_dataset(input, force = force)
  design <- pData(input)
  conditions <- design[["condition"]]
  conditions_table <- table(conditions)
  batches <- design[["batch"]]
  batches_table <- table(batches)
  data <- input_data[["data"]]
  conditions <- as.factor(conditions)
  batches <- as.factor(batches)

  noiseq_input <- NOISeq::readData(input_data[["data"]], factors = pData(input))
  norm <- NOISeq::ARSyNseq(noiseq_input, factor = "condition", batch = FALSE,
                           norm = "rpkm", logtransf = FALSE)


  ## Yes I know NOISeq doesn't use models in the same way as other
  ## methods I have applied, but this will make it easier to set up
  ## the contrasts.
  model_choice <- choose_model(input, conditions = conditions,
                               batches = batches,
                               model_batch = model_batch,
                               model_cond = model_cond)
  model_including <- model_choice[["including"]]
  if (class(model_choice[["model_batch"]])[1] == "matrix") {
    model_batch <- model_choice[["model_batch"]]
  }
  model_data <- model_choice[["chosen_model"]]
  model_string <- model_choice[["chosen_string"]]
  apc <- make_pairwise_contrasts(model_data, conditions)

  contrast_list <- list()
  result_list <- list()
  lrt_list <- list()
  sc <- vector("list", length(apc[["names"]]))
  end <- length(apc[["names"]])
  if (isTRUE(verbose)) {
    bar <- utils::txtProgressBar(style = 3)
  }
  for (con in seq_along(apc[["names"]])) {
    name <- apc[["names"]][[con]]
    numerator <- apc[["numerators"]][[name]]
    denominator <- apc[["denominators"]][[name]]
    if (isTRUE(verbose)) {
      pct_done <- con / length(apc[["names"]])
      utils::setTxtProgressBar(bar, pct_done)
    }
    noiseq_table <- NOISeq::noiseqbio(
      norm, k = 0.5, norm = "rpkm", factor = "condition", lc = 1,
      r = 20, adj = 1.5, plot = TRUE, a0per = 0.9, filter = 1,
      conditions = c(numerator, denominator))
    result_list[[name]] <- noiseq_table@results[[1]]
  }
  if (isTRUE(verbose)) {
    close(bar)
  }

  retlist <- list(
      "all_tables" = result_list,
      "batches" = batches,
      "batches_table" = batches_table,
      "conditions" = conditions,
      "conditions_table" = conditions_table,
      "contrast_list" = contrast_list,
      "contrasts" = apc,
      "contrast_string" = contrast_string,
      "contrasts_performed" = apc[["names"]],
      "input_data" = input,
      "method" = "noiseq",
      "model" = model_data,
      "model_string" = model_string)
  class(retlist) <- c("noiseq_result", "list")
  return(retlist)
}

## EOF
