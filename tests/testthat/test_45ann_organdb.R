start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("45ann_organdb.R: Is it possible to manipulate OrgDb/OrganismDbi instances?\n")

tt <- sm(require.auto("org.Dm.eg.db"))
tt <- sm(library(org.Dm.eg.db))
fly_annotations <- sm(load_orgdb_annotations(org.Dm.eg.db, keytype="ENSEMBL", fields=c("SYMBOL", "ENSEMBL", "ACCNUM", "ENTREZID")))
actual <- head(fly_annotations[["genes"]][["accnum"]])
expected <- c("AAF45486", "AAL28774", "AAN09009",
              "ADR83721", "AHN59203", "AY061226")
test_that("Can we extract annotation information from an orgdb instance with ease?", {
    expect_equal(expected, actual)
})

orgdb_go <- sm(load_orgdb_go(org.Dm.eg.db, gene_ids=expected, keytype="ACCNUM"))
actual <- head(orgdb_go[["GO"]])
expected <- c("GO:0006486", "GO:0006486", "GO:0006486",
              "GO:0006486", "GO:0006486", "GO:0006486")
test_that("Can we use GO.db with Orgdb to extract gene ontology information?", {
    expect_equal(expected, actual)
})

orgdb_kegg <- sm(load_orgdb_kegg(org.Dm.eg.db, columns="PATH", keytype="ENTREZID"))
actual <- head(orgdb_kegg[["PATH"]])
expected <- c("00310", "00903", "01100", "00280", "00562", "00640")
test_that("Can we extract KEGG data from Orgdb?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 45ann_organdb.R in ", elapsed,  " seconds."))
tt <- clear_session()
