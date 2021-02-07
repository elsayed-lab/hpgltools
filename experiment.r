Experiment <- R6Class(
    "Experiment",
    public = list(
        "expressionset" = NA,
        "metadata" = NA,
        "gene_info" = NA,
        "count_data" = NA,
        "sample_colors" = NULL,
        "title" = NULL,
        "notes" = NULL,
        "file_column" = NULL,
        "id_column" = NULL,
        "savefile" = NULL,
        "conditions" = NULL,
        "batches" = NULL,
        "sanitize_rownames" = FALSE,
        "state" = list("norm" = "raw", "transform" = "raw", "convert" = "raw",
                       "batch" = "raw", "filter" = "raw"),

        initialize = function(expressionset = NA, metadata = NA,
                              gene_info = NA, count_data = NA) {
          self$metadata <- metadata
          self$gene_info <- gene_info
          self$count_data <- count_data
          if (is.na(expressionset) & is.na(metadata)) {
            stop("This requires either an expressionset or some metadata at minimum.")
          } else if (is.na(metadata)) {
            expt <- coerce_expressionset(expressionset)
          }

          invisible(expt)
        },
        set_colors = function(conditions = NA) {
        },
        print_experiment = function() {
        },

        coerce_expressionset(expressionset) <- function(expressionset, condition_column = "condition",
                                                        batch_column = "batch", id_column = "row.names") {
          metadata <- pData(expressionset)
          self$conditions <- as.factor(metadata[[condition_column]])
          self$batches <- as.factor(metadata[[batch_column]])
          self$colors <- colors_from_factor(self$conditions)
        }
    ))


colors_from_factor <- function(fact, palette="Dark2") {
  num_colors <- length(levels(fact))
  sample_names <- names(fact)
  sample_colors <- sm(
      grDevices::colorRampPalette(
                     RColorBrewer::brewer.pal(num_colors, palette))(num_colors))
  mapping <- setNames(sample_colors, unique(chosen_colors))
  chosen_colors <- mapping[chosen_colors]
}
