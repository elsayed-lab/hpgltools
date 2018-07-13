## ----options, include=FALSE----------------------------------------------
## These are the options I tend to favor
library("hpgltools")
## tt <- devtools::load_all("~/hpgltools")
knitr::opts_knit$set(progress=TRUE,
                     verbose=TRUE,
                     width=90,
                     echo=TRUE)
knitr::opts_chunk$set(error=TRUE,
                      fig.width=8,
                      fig.height=8,
                      dpi=96)
old_options <- options(digits=4,
                       stringsAsFactors=FALSE,
                       knitr.duplicate.label="allow")
ggplot2::theme_set(ggplot2::theme_bw(base_size=10))
set.seed(1)
rmd_file <- "c-03_fission_differential_expression.Rmd"

## ----setup---------------------------------------------------------------
library(hpgltools)
tt <- sm(library(fission))
tt <- data(fission)

## ----data_import---------------------------------------------------------
## Extract the meta data from the fission dataset
meta <- as.data.frame(fission@colData)
## Make conditions and batches
meta$condition <- paste(meta$strain, meta$minute, sep=".")
meta$batch <- meta$replicate
meta$sample.id <- rownames(meta)
## Grab the count data
fission_data <- fission@assays$data$counts
## This will make an experiment superclass called 'expt' and it contains
## an ExpressionSet along with any arbitrary additional information one might want to include.
## Along the way it writes a Rdata file which is by default called 'expt.Rdata'
fission_expt <- create_expt(metadata=meta, count_dataframe=fission_data)

## ----simple_subset-------------------------------------------------------
fun_data <- subset_expt(fission_expt,
                        subset="condition=='wt.120'|condition=='wt.30'")
fun_norm <- sm(normalize_expt(fun_data, batch="limma", norm="quant",
                              transform="log2", convert="cpm"))

## ----simple_limma--------------------------------------------------------
limma_comparison <- sm(limma_pairwise(fun_data))
names(limma_comparison$all_tables)
summary(limma_comparison$all_tables$wt30_vs_wt120)
scatter_wt_mut <- extract_coefficient_scatter(limma_comparison, type="limma",
                                              x="wt30", y="wt120")
scatter_wt_mut$scatter
scatter_wt_mut$both_histogram$plot + ggplot2::scale_y_continuous(limits=c(0, 0.20))
ma_wt_mut <- extract_de_plots(limma_comparison, type="limma")
ma_wt_mut$ma$plot
ma_wt_mut$volcano$plot

## ----simple_deseq2-------------------------------------------------------
deseq_comparison <- sm(deseq2_pairwise(fun_data))
summary(deseq_comparison$all_tables$wt30_vs_wt120)
scatter_wt_mut <- extract_coefficient_scatter(deseq_comparison, type="deseq",
                                              x="wt30", y="wt120", gvis_filename=NULL)
scatter_wt_mut$scatter
plots_wt_mut <- extract_de_plots(deseq_comparison, type="deseq")
plots_wt_mut$ma$plot
plots_wt_mut$volcano$plot

## ----simple_edger--------------------------------------------------------
edger_comparison <- sm(edger_pairwise(fun_data, model_batch=TRUE))
plots_wt_mut <- extract_de_plots(edger_comparison, type="edger")
scatter_wt_mut <- extract_coefficient_scatter(edger_comparison, type="edger",
                                              x="wt30", y="wt120", gvis_filename=NULL)
scatter_wt_mut$scatter
plots_wt_mut$ma$plot
plots_wt_mut$volcano$plot

## ----simple_basic--------------------------------------------------------
basic_comparison <- sm(basic_pairwise(fun_data))
summary(basic_comparison$all_tables$wt30_vs_wt120)
scatter_wt_mut <- extract_coefficient_scatter(basic_comparison, type="basic",
                                              x="wt30", y="wt120", gvis_filename=NULL)
scatter_wt_mut$scatter
plots_wt_mut <- extract_de_plots(basic_comparison, type="basic")
plots_wt_mut$ma$plot
plots_wt_mut$volcano$plot

## ----simple_all----------------------------------------------------------
all_comparisons <- sm(all_pairwise(fun_data, model_batch=TRUE))
all_combined <- sm(combine_de_tables(all_comparisons, excel=FALSE))
head(all_combined$data[[1]])
sig_genes <- sm(extract_significant_genes(all_combined, excel=FALSE))
head(sig_genes$limma$ups[[1]])

## Here we see that edger and deseq agree the least:
all_comparisons$comparison$comp

## And here we can look at the set of 'significant' genes according to various tools:
yeast_sig <- extract_significant_genes(all_combined, excel=FALSE)
yeast_barplots <- sm(significant_barplots(combined=all_combined))
yeast_barplots$limma
yeast_barplots$edger
yeast_barplots$deseq

## ----ontology_setup------------------------------------------------------
limma_results <- limma_comparison$all_tables
## The set of comparisons performed
names(limma_results)
table <- limma_results$wt30_vs_wt120
dim(table)
gene_names <- rownames(table)

updown_genes <- get_sig_genes(table, p=0.05, lfc=0.4, p_column="P.Value")
tt <- please_install("GenomicFeatures")
tt <- please_install("biomaRt")
available_marts <- biomaRt::listMarts(host="fungi.ensembl.org")
available_marts
ensembl_mart <- biomaRt::useMart("fungi_mart", host="fungi.ensembl.org")
available_datasets <- biomaRt::listDatasets(ensembl_mart)
pombe_hit <- grep(pattern="pombe", x=available_datasets[["description"]])
pombe_name <- available_datasets[pombe_hit, "dataset"]
pombe_mart <- biomaRt::useDataset(pombe_name, mart=ensembl_mart)

pombe_goids <- biomaRt::getBM(attributes=c("pombase_transcript", "go_id"),
                              values=gene_names, mart=pombe_mart)
colnames(pombe_goids) <- c("ID", "GO")

## ----ontology_setup_hpgltools--------------------------------------------
## In theory, the above should work with a single function call:
pombe_goids_simple <- load_biomart_go(species="spombe", overwrite=TRUE,
                                      dl_rows=c("pombase_transcript", "go_id"),
                                      host="fungi.ensembl.org")
head(pombe_goids_simple)
head(pombe_goids)

## This used to work, but does so no longer and I do not know why.
## pombe <- sm(GenomicFeatures::makeTxDbFromBiomart(biomart="fungal_mart",
##                                                  dataset="spombe_eg_gene",
##                                                  host="fungi.ensembl.org"))

## I bet I can get all this information from ensembl now.
## This was found at the bottom of: https://www.biostars.org/p/232005/
link <- "ftp://ftp.ensemblgenomes.org/pub/release-34/fungi/gff3/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.34.gff3.gz"
pombe <- GenomicFeatures::makeTxDbFromGFF(link, format="gff3", organism="Schizosaccharomyces pombe",
                                          taxonomyId=4896)

pombe_transcripts <- as.data.frame(GenomicFeatures::transcriptsBy(pombe))
lengths <- pombe_transcripts[, c("group_name","width")]
colnames(lengths) <- c("ID","width")
## Something useful I didn't notice before:
## makeTranscriptDbFromGFF()  ## From GenomicFeatures, much like my own gff2df()
gff_from_txdb <- GenomicFeatures::asGFF(pombe)
## why is GeneID: getting prefixed to the IDs!?
gff_from_txdb$ID <- gsub(x=gff_from_txdb$ID, pattern="GeneID:", replacement="")
written_gff <- rtracklayer::export.gff3(gff_from_txdb, con="pombe.gff")

## ----test_goseq----------------------------------------------------------
summary(updown_genes)
test_genes <- updown_genes$down_genes
rownames(test_genes) <- paste0(rownames(test_genes), ".1")
lengths$ID <- paste0(lengths$ID, ".1")
goseq_result <- sm(simple_goseq(sig_genes=test_genes, go_db=pombe_goids, length_db=lengths))
head(goseq_result$alldata)
goseq_result$pvalue_plots$mfp_plot

test_genes <- updown_genes$up_genes
rownames(test_genes) <- paste0(rownames(test_genes), ".1")
goseq_result <- sm(simple_goseq(sig_genes=test_genes, go_db=pombe_goids, length_db=lengths))
head(goseq_result$alldata)
goseq_result$pvalue_plots$bpp_plot

## ----test_cp, eval=FALSE-------------------------------------------------
#  ## holy crap makeOrgPackageFromNCBI is slow, no slower than some of mine, so who am I to complain.
#  orgdb <- AnnotationForge::makeOrgPackageFromNCBI(version="0.1", author="atb <abelew@gmail.com>",
#                                                   maintainer="atb <abelew@gmail.com>", tax_id="4896",
#                                                   genus="Schizosaccharomyces", species="pombe")
#  ## This created the directory 'org.spombe.eg.db'
#  devtools::install_local("org.Spombe.eg.db")
#  library(org.Spombe.eg.db)
#  ## Don't forget to remove the terminal .1 from the gene names...
#  ## If you do forget this, it will fail for no easily visible reason until you remember
#  ## this and get really mad at yourself.
#  rownames(test_genes) <- gsub(pattern=".1$", replacement="", x=rownames(test_genes))
#  pombe_goids[["ID"]] <- gsub(pattern=".1$", replacement="", x=pombe_goids[["ID"]])
#  cp_result <- simple_clusterprofiler(sig_genes=test_genes, do_david=FALSE, do_gsea=FALSE,
#                                      de_table=all_combined$data[[1]],
#                                      orgdb=org.Spombe.eg.db, orgdb_to="ALIAS")
#  cp_result[["pvalue_plots"]][["ego_all_mf"]]
#  ## Yay bar plots!

## ----test_tp-------------------------------------------------------------
## Get rid of those stupid terminal .1s.
rownames(test_genes) <- gsub(pattern=".1$", replacement="", x=rownames(test_genes))
pombe_goids[["ID"]] <- gsub(pattern=".1$", replacement="", x=pombe_goids[["ID"]])
tp_result <- sm(simple_topgo(sig_genes=test_genes, go_db=pombe_goids, pval_column="limma_adjp"))

tp_result[["pvalue_plots"]][["mfp_plot_over"]]
tp_result[["pvalue_plots"]][["bpp_plot_over"]]

## ----gst_test, eval=FALSE------------------------------------------------
#  ## Get rid of those stupid terminal .1s.
#  rownames(test_genes) <- gsub(pattern=".1$", replacement="", x=rownames(test_genes))
#  pombe_goids[["ID"]] <- gsub(pattern=".1$", replacement="", x=pombe_goids[["ID"]])
#  ## universe_merge is the column in the final data frame when.
#  ## gff_type is the field in the gff file providing the id, this may be redundant with
#  ## universe merge, that is something to check on...
#  gst_result <- sm(simple_gostats(sig_genes=test_genes, go_db=pombe_goids, universe_merge="id",
#                                  gff_type="gene",
#                                  gff="pombe.gff", pval_column="limma_adjp"))

## ----sysinfo, results="asis"---------------------------------------------
pander::pander(sessionInfo())

