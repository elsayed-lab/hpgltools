start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
##context("42ann_microbes.R: May I used microbesonline?\n")
##
##ids <- sm(get_microbesonline_ids())
##actual <- head(sort(ids[["shortName"]]))
##expected <- c("Escherichia albertii TW07627", "Escherichia coli",
##              "Escherichia coli 'BL21-Gold(DE3)pLysS AG'", "Escherichia coli 101-1",
##              "Escherichia coli 1520", "Escherichia coli 536")
## test_that("Do we get an expected set of species names?", {
##     expect_equal(expected, actual)
## })

## actual <- head(sort(ids[["taxonomyId"]]))
## expected <- c(562, 37762, 155864, 199310, 299586, 316385)
## test_that("Do we get an expected set of microbesonline IDs?", {
##     expect_equal(expected, actual)
## })

## shortname <- sm(get_microbesonline_name())
## actual <- shortname[["shortName"]]
## expected <- "Escherichia coli str. K-12 substr. DH10B"
## test_that("Do we get the correct E. coli strain name?", {
##     expect_equal(expected, actual)
## })

## id_df <- sm(get_microbesonline_ids(name="pyogenes"))
## my_id <- id_df[["shortName"]] == "Streptococcus pyogenes MGAS5005"
## my_id <- id_df[my_id, 1]

## annotations <- sm(load_microbesonline_annotations(id=my_id))
## annotations <- annotations[[1]]
## expected <- c("dnaA", "dnaN", "M5005_Spy_0003",
##               "M5005_Spy_0004", "pth", "trcF")
## actual <- as.character(head(annotations[["name"]]))
## test_that("Do we get the correct gene names for Streptococcus pyogenes?", {
##     expect_equal(expected, actual)
## })

## ## The default is 160490, an E. coli strain.
## go_ids <- sm(load_microbesonline_go(name="pyogenes%5005"))
## actual <- head(sort(unique(go_ids[["acc"]])))
## expected <- c("GO:0000015", "GO:0000062", "GO:0000074",
##               "GO:0000105", "GO:0000150", "GO:0000154")
## test_that("Can we download GO accessions from microbesonline?", {
##     expect_equal(expected, actual)
## })
## end <- as.POSIXlt(Sys.time())
## elapsed <- round(x=as.numeric(end) - as.numeric(start))
## message(paste0("\nFinished 42ann_microbes.R in ", elapsed,  " seconds."))
## tt <- try(clear_session())
