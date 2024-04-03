start <- as.POSIXlt(Sys.time())
context("105helpers_misc.R")

## This function actually makes untenable assumptions about where the repository lives.
testing <- get_git_commit(gitdir = "")
expected <- c("glue", "character")
actual <- class(testing)
test_that("get_git_commit() gave me a commit id?", {
  expect_equal(expected, actual)
})

sp <- plot_hypotrochoid()
test_that("We can make fun spirograph plots?", {
  expect_equal(class(sp)[1], "gg")
})
print_file <- "spirograph.png"
printed <- pp(file = print_file, image = sp)
test_that("We can print them easily to disk?", {
  expect_equal(class(printed)[1], "gg")
  expect_true(file.exists(print_file))
  expect_true(file.remove(print_file))
})

## So, a bunch of functions exported in helpers are difficult to test because
## they are pretty specific to their little domains.  I probably therefore will
## skip testing some (many) of them.
end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 105helpers_misc.R in ", elapsed,  " seconds.")
