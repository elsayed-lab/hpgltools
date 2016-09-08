## ----circos--------------------------------------------------------------
microbe_ids <- as.character(get_microbesonline_ids("pyogenes MGAS5005"))
mgas_df <- sm(get_microbesonline_annotation(microbe_ids[[1]])[[1]])
mgas_df$sysName <- gsub(pattern="Spy_", replacement="Spy", x=mgas_df$sysName)
rownames(mgas_df) <- make.names(mgas_df$sysName, unique=TRUE)


## First make a template configuration
circos_test <- circos_prefix()
## Fill it in with the data for s.pyogenes
circos_kary <- circos_karyotype("mgas", length=1895017)
## Fill in the gene category annotations by gene-strand
circos_plus <- sm(circos_plus_minus(mgas_df, circos_test))

circos_limma_hist <- sm(circos_hist(spyogenes_de$limma$all_tables[[1]], mgas_df, circos_test, outer=circos_plus))
circos_deseq_hist <- sm(circos_hist(spyogenes_de$deseq$all_tables[[1]], mgas_df, circos_test, outer=circos_limma_hist))
circos_edger_hist <- sm(circos_hist(spyogenes_de$edger$all_tables[[1]], mgas_df, circos_test, outer=circos_deseq_hist))
circos_suffix(cfgout=circos_test)
circos_made <- sm(circos_make(target="mgas"))

## ----sysinfo, results='asis'---------------------------------------------
library('pander')
pander(sessionInfo())

