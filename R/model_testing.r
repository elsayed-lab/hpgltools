#' Make sure a given experimental factor and design will play together.
#'
#' Have you ever wanted to set up a differential expression analysis and after minutes of the
#' computer churning away it errors out with some weird error about rank?  Then this is the function
#' for you!
#'
#' @param design Dataframe describing the design of the experiment.
#' @param goal Experimental factor you actually want to learn about.
#' @param factors Experimental factors you rather wish would just go away.
#' @param ... I might decide to add more options from other functions.
#' @return List of booleans telling if the factors + goal will work.
#' @seealso \code{\link{model.matrix}} \code{\link{qr}}
#' @export
model_test <- function(design, goal="condition", factors=NULL, ...) {
  arglist <- list(...)
  ## For testing, use some existing matrices/data
  message("There are ", length(levels(as.factor(design[, goal]))),
          " levels in the goal: ", goal, ".")
  ret_list <- list()
  if (is.null(factors)) {
    message("Testing an experimental design with only ", goal, ".")
    matrix_all_formula <- as.formula(glue("~ 0 + {goal}"))
    matrix_test <- model.matrix(matrix_all_formula, data=design)
    num_columns <- ncol(matrix_test)
    matrix_decomp <- qr(matrix_test)
    message("The model of ", goal, "has ", num_columns,
            " levels and rank ", matrix_decomp[["rank"]], ".")
    if (matrix_decomp[["rank"]] < num_columns) {
      message("This will not work, a different factor should be used.")
      ret_list[[goal]] <- 0
    } else {
      ret_list[[goal]] <- 1
    }
    for (factor in colnames(design)) {
      if (factor == goal) {
        next
      }
      message("Testing an experimental design with ", goal, " and ", factor, ".")
      matrix_goal <- design[, goal]
      matrix_factor <- design[, factor]
      if (length(levels(as.factor(matrix_factor))) == 1) {
        message("Factor ", factor, " has only 1 level, skipping it.")
        next
      }
      matrix_all_formula <- as.formula(glue("~ 0 + {goal} + {factor}"))
      matrix_test <- model.matrix(matrix_all_formula, data=design)
      num_columns <- ncol(matrix_test)
      matrix_decomp <- qr(matrix_test)
      message("The model of ", goal, " and ", factor, " has ", num_columns,
              " and rank ", matrix_decomp[["rank"]])
      if (matrix_decomp[["rank"]] < num_columns) {
        message("This will not work, a different factor should be used.")
        ret_list[[factor]] <- 0
      } else {
        ret_list[[factor]] <- 1
      }
    } ## End for loop
  } else {
    for (factor in factors) {
      matrix_goal <- design[, goal]
      matrix_factor <- design[, factor]
      matrix_all_formula <- as.formula(glue("~ 0 + {goal} + {factor}"))
      matrix_test <- model.matrix(matrix_all_formula, data=design)
      num_columns <- ncol(matrix_test)
      matrix_decomp <- qr(matrix_test)
      message("The model of ", goal, " and ", factor, " has ", num_columns,
              " and rank ", matrix_decomp[["rank"]])
      if (matrix_decomp[["rank"]] < num_columns) {
        message("This will not work, a different factor should be used.")
        ret_list[[factor]] <- 0
      } else {
        ret_list[[factor]] <- 1
      }
    }
  }
  return(ret_list)
}

## EOF
