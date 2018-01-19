start <- as.POSIXlt(Sys.time())
context("105helpers_misc.R:\n")
## 2017-12, exported functions in helpers_misc:
## bioc_all() clear_session() cordist() get_git_commit()
## hpgl_arescore() hpgl_cor() install_packrat_globally() loadme()
## make_report() pp() require.auto() rex() saveme() sillydist()
## sm() unAsIs() ymxb_print()

## bioc_all()

## clear_session()

## cordist()

## get_git_commit()
testing <- sm(get_git_commit())
expected <- "character"
actual <- class(testing)
test_that("get_git_commit() gave me a commit id?", {
  expect_equal(expected, actual)
})

## hpgl_arescore()

## hpgl_cor()

## install_packrat_globally()

## loadme()

## make_report()

## pp()

## require.auto()

## rex()

## saveme()

## sillydist()

## sm()

## unAsIs()

## ymxb_print()

## So, a bunch of functions exported in helpers are difficult to test because
## they are pretty specific to their little domains.  I probably therefore will
## skip testing some (many) of them.
end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 105helpers_misc.R in ", elapsed,  " seconds."))
