start <- as.POSIXlt(Sys.time())
context("200ontology_kegg.R:\n")
## 2017-12, exported functions in annotation_gff:

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 200ontology_kegg.R in ", elapsed,  " seconds."))
