start <- as.POSIXlt(Sys.time())
context("030annotation_microbesonline.R\n")
## 2017-12, exported functions in annotation_microbesonline:
## get_microbesonline_ids(), get_microbesonline_name(), load_microbesonline_annotations()
## load_microbesonline_go(), load_microbesonline_kegg().

## get_microbesonline_ids()
## ecoli_ids <- get_microbesonline_ids()
## actual <- head(ecoli_ids[[1]])  ## Taxonomy IDs
## expected <- c(155864, 199310, 316407, 364106, 362663, 386585)
## test_that("We get expected microbesonline IDs?", {
##   expect_equal(expected, actual)
## })

## actual <- head(ecoli_ids[[2]])  ## short species names
## expected <- c("Escherichia coli O157:H7 EDL933", "Escherichia coli CFT073",
##               "Escherichia coli W3110", "Escherichia coli UTI89",
##               "Escherichia coli 536", "Escherichia coli O157:H7 str. Sakai")
## test_that("We get expected microbesonline names?", {
##   expect_equal(expected, actual)
## })

## ## get_microbesonline_name()
## expected <- expected[1]
## id <- get_microbesonline_ids(expected)[1,1]
## name <- get_microbesonline_name(id)
## actual <- name[1, 1]
## test_that("Can we convert a name to id back to a name?", {
##   expect_equal(expected, actual)
## })

## ## load_microbesonline_annotations()
## spy_annot <- sm(load_microbesonline_annotations())
## actual <- dim(spy_annot[[1]])
## expected <- c(1871, 18)
## test_that("Do we get the expected amount of pyogenes data?", {
##   expect_equal(expected[1], actual[1])
##   expect_equal(expected[2], actual[2])
## })

## ## load_microbesonline_go()
## ecoli_go <- sm(load_microbesonline_go())
## actual <- dim(ecoli_go)
## expected <- c(3709, 5)
## test_that("Do we get the expected amount of ecoli GO data?", {
##   expect_equal(expected[1], actual[1])
##   expect_equal(expected[2], actual[2])
## })

## ## load_microbesonline_kegg().
## spy_kegg <- sm(load_microbesonline_kegg())
## actual <- dim(spy_kegg)
## expected <- c(1801, 6)
## test_that("Do we get the expected amount of spyogenes KEGG data?", {
##   expect_equal(expected[1], actual[1])
##   expect_equal(expected[2], actual[2])
## })

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 030nnotation_microbesonline.R in ", elapsed,  " seconds."))
