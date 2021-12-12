## ----de_test------------------------------------------------------------------
spyogenes_de <- sm(all_pairwise(expt))
## Even the lowest correlations are quite high.

## ----keeper_example-----------------------------------------------------------
my_keepers <- list(
  ## name    =   numerator / denominator
  "wt_media" = c("wt_ll_cf", "wt_ll_cg"),
  "mga_media" = c("mga_ll_cf", "mga_ll_cg"))

## ----combine_test-------------------------------------------------------------
spyogenes_tables <- sm(combine_de_tables(spyogenes_de, excel = FALSE))
summary(spyogenes_tables)
## Try changing the p-adjustment
spyogenes_tables <- sm(combine_de_tables(spyogenes_de, excel = FALSE, padj_type = "BH"))
knitr::kable(head(spyogenes_tables$data[[1]]))

## ----sig_genes_test, fig.show = "hide"----------------------------------------
spyogenes_sig <- sm(extract_significant_genes(spyogenes_tables, excel = FALSE))
knitr::kable(head(spyogenes_sig$limma$ups[[1]]))

## ----circos-------------------------------------------------------------------
##microbe_ids <- as.character(sm(get_microbesonline_ids("pyogenes MGAS5005")))
## A caveat!  The new version of microbesonline changed the IDs so that they no longer
## match my old rnaseq analysis!!  Thus I put my old gff file used for mapping into inst/
## and will load the annotation data from that; but I will use this data to gather
## the COG information.
mgas_gff_df <- sm(load_gff_annotations(gff = system.file("share/gas.gff", package = "hpgltools")))
mgas_gff_df <- mgas_gff_df[-1, ]
mgas_gff_df[["sysName"]] <- gsub(pattern = "Spy_", replacement = "Spy",
                                 x = mgas_gff_df[["locus_tag"]])
rownames(mgas_gff_df) <- make.names(mgas_gff_df[["sysName"]], unique = TRUE)

mgas_microbes_df <- sm(load_microbesonline_annotations(id = 293653))
mgas_microbes_df$sysName <- gsub(pattern = "Spy_", replacement = "Spy",
                                 x = mgas_microbes_df$sysName)
rownames(mgas_microbes_df) <- make.names(mgas_microbes_df$sysName, unique = TRUE)

mgas_df <- merge(x = mgas_gff_df, y = mgas_microbes_df, by = "row.names")
rownames(mgas_df) <- mgas_df[["Row.names"]]
mgas_df <- mgas_df[, -1]
colnames(mgas_df) <- c("seqnames", "start", "end", "width", "strand", "source", "type",
                       "score", "phase", "ID", "Dbxref", "Is_circular", "gbkey", "genome",
                       "mol_type", "strain", "Name", "Note", "gene", "locus_tag",
                       "Parent", "product", "protein_id", "transl_table", "gene_synonym",
                       "sysName_again", "locusId", "accession", "GI", "scaffoldId",
                       "start_again", "stop", "strand_again", "sysName_again", "name",
                       "desc", "COG", "COGFun", "COGDesc", "TIGRFam", "TIGRRoles",
                       "GO", "EC", "ECDesc")

## First make a template configuration
circos_test <- circos_prefix(annotation = mgas_df)
## Fill it in with the data for s.pyogenes
lengths <- 1838600
names(lengths) <- "NC_007297"

circos_kary <- circos_karyotype(cfg = circos_test, lengths = lengths)
## Fill in the gene category annotations by gene-strand
circos_plus <- circos_plus_minus(cfg = circos_test)
circos_limma_hist <- circos_hist(cfg = circos_test,
                                 df = spyogenes_de$limma$all_tables[[1]],
                                 basename = "limma",
                                 colname = "logFC",
                                 outer = circos_plus)
circos_deseq_hist <- circos_hist(cfg = circos_test,
                                 df = spyogenes_de$deseq$all_tables[[1]],
                                 basename = "deseq",
                                 colname = "logFC",
                                 outer = circos_limma_hist)
circos_edger_hist <- circos_hist(cfg = circos_test,
                                 df = spyogenes_de$edger$all_tables[[1]],
                                 basename = "edger",
                                 colname = "logFC",
                                 outer = circos_deseq_hist)
circos_suffix(cfg = circos_test)
circos_made <- sm(circos_make(cfg = circos_test, target = "mgas"))
getwd()

## ----genoplot-----------------------------------------------------------------
genoplot_chromosome()

## ----wt_mga, fig.show = "hide"------------------------------------------------
wt_mga_expt <- set_expt_conditions(expt = expt, fact = "type")
wt_mga_plots <- sm(graph_metrics(wt_mga_expt))
wt_mga_norm <- sm(normalize_expt(wt_mga_expt, transform = "log2", convert = "raw", filter = TRUE, norm = "quant"))
wt_mga_nplots <- sm(graph_metrics(wt_mga_norm))
wt_mga_de <- sm(all_pairwise(input = wt_mga_expt,
                          combined_excel = "wt_mga.xlsx",
                          sig_excel = "wt_mga_sig.xlsx",
                          abundant_excel = "wt_mga_abundant.xlsx"))

## ----wt_mga_plots-------------------------------------------------------------
wt_mga_de$combined$comp_plot
## How well do the various DE tools agree on this data?

wt_mga_plots$tsne_plot
wt_mga_nplots$pc_plot
wt_mga_de$combined$limma_plots$WT_vs_mga$scatter
wt_mga_de$combined$limma_ma_plots$WT_vs_mga$plot
wt_mga_de$combined$limma_vol_plots$WT_vs_mga$plot

## ----sysinfo, results='asis'--------------------------------------------------
pander::pander(sessionInfo())

