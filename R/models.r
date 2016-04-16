## Time-stamp: <Sat Apr 16 00:12:42 2016 Ashton Trey Belew (abelew@gmail.com)>

model_test <- function(design, goal="condition", factors=NULL, ...) {
    arglist <- list(...)
    ## For testing, use some existing matrices/data
    message(paste0("There are ", length(levels(as.factor(design[, goal]))), " levels in the goal."))
    ret_list <- list()
    if (is.null(factors)) {
        for (factor in colnames(design)) {
            matrix_goal <- design[, goal]
            matrix_factor <- design[, factor]
            matrix_all_formula <- as.formula(paste0("~ 0 + ", goal, " + ", factor))
            matrix_test <- model.matrix(matrix_all_formula, data=design)
            num_columns <- ncol(matrix_test)
            matrix_decomp <- qr(matrix_test)
            message(paste0("The model of ", goal, " and ", factor, " has ", num_columns, " and rank ", matrix_decomp[["rank"]]))
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
            matrix_all_formula <- as.formula(paste0("~ 0 + ", goal, " + ", factor))
            matrix_test <- model.matrix(matrix_all_formula, data=design)
            num_columns <- ncol(matrix_test)
            matrix_decomp <- qr(matrix_test)
            message(paste0("The model of ", goal, " and ", factor, " has ", num_columns, " and rank ", matrix_decomp[["rank"]]))
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

