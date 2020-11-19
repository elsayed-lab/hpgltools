start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("45ann_organdb.R: Is it possible to manipulate OrgDb/OrganismDbi instances?\n")

tt <- sm(please_install("org.Dm.eg.db"))
## I have a separate implementation of this function in EuPathDB which I may rename
fly_annotations <- sm(hpgltools::load_orgdb_annotations(
  "org.Dm.eg.db",
  keytype="ensembl",
  fields=c("symbol", "ensembl", "accnum", "entrezid")))
actual <- head(sort(fly_annotations[["genes"]][["accnum"]]))
##expected <- c("AAF45486", "AAL28774", "AAN09009",
##              "ADR83721", "AHN59203", "AY061226")
expected <- c("0908303C", "1205216A", "1304299A",
              "1410208A", "1505369A", "1505374A")
test_that("Can we extract annotation information from an orgdb instance with ease?", {
    expect_equal(expected, actual)
})

## New org.Dm.eg.db seems to have a different order of GO IDs.
gene_ids <- actual
orgdb_go <- sm(hpgltools::load_orgdb_go(
  org.Dm.eg.db,
  gene_ids=gene_ids,
  keytype="accnum",
  columns=c("go", "goall", "goid")))
actual <- sort(head(unique(orgdb_go[["GO"]])))
expected <- c("GO:0000122", "GO:0000413", "GO:0000976",
              "GO:0000978", "GO:0000981", "GO:0001222")
test_that("Can we use GO.db with Orgdb to extract gene ontology information?", {
    expect_equal(expected, actual)
})

orgdb_kegg <- sm(hpgltools::load_orgdb_go(org.Dm.eg.db, columns="path", keytype="entrezid"))
actual <- head(orgdb_kegg[["PATH"]])
expected <- c("00310", "00903", "01100", "00280", "00562", "00640")
test_that("Can we extract KEGG data from Orgdb?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 45ann_organdb.R in ", elapsed,  " seconds."))
