## ----options, include=FALSE---------------------------------------------------
library("hpgltools")
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
rmd_file <- "a-01_bacterial_example.Rmd"

## ----loading_data-------------------------------------------------------------
library(hpgltools)
data_file <- system.file("cdm_expt.rda", package="hpgltools")
cdm <- new.env()
load(data_file, envir=cdm)
rm(data_file)

ls()

## ----create_expt--------------------------------------------------------------
expt <- create_expt(count_dataframe=cdm$cdm_counts,
                    metadata=cdm$cdm_metadata,
                    gene_info=cdm$gene_info)

knitr::kable(head(expt$design))
summary(expt)

## ----test_se------------------------------------------------------------------
tt <- sm(library(SummarizedExperiment))
test_expr <- expt[["expressionset"]]

## The following in theory converts to a SE, but is a bit unhelpful in its result.
## Convert expressionSet to summarizedExperiment
test_se <- SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(test_expr)
## Holy crap long function name!
test_assay <- assay(test_se)
head(test_assay)
test_meta <- colData(test_se)
head(test_meta)
test_gene_info <- rowData(test_se)
head(test_gene_info)

## Here is a stupider version of the above, but one which will not make my fingers sad.
exprs_to_se <- function(exset) {
  mtrx <- Biobase::exprs(exset)
  annot <- Biobase::fData(exset)
  meta <- Biobase::pData(exset)
  se <- SummarizedExperiment(assays=mtrx, colData=meta, rowData=annot)
  return(se)
}

## and back?
test_ex <- as(test_se, "ExpressionSet")
head(exprs(test_ex))

## YAY!

## ----graph_original, fig.show="hide"------------------------------------------
raw_metrics <- sm(graph_metrics(expt, qq=TRUE, cis=NULL))

## ----show_original_plots------------------------------------------------------
## View a raw library size plot
raw_metrics$libsize
## Or boxplot to see the data distribution
raw_metrics$boxplot
## The warning is because it automatically uses a log scale and there are some 0 count genes.
## Perhaps you prefer density plots
raw_metrics$density
## quantile/quantile plots compared to the median of all samples
raw_metrics$qqrat
raw_metrics$tsne_plot
## Here we can see some samples are differently 'shaped' compared to the median than others
## There are other plots one may view, but this data set is a bit too crowded as is.
## The following summary shows the other available plots:
summary(raw_metrics)

## ----subset_data, fig.show='hide'---------------------------------------------
head(expt$design)
## elt stands for: "early/late in thy"
batch_a <- subset_expt(expt, subset="batch=='a'")
batch_b <- subset_expt(expt, subset="batch=='b'")

a_metrics <- sm(graph_metrics(batch_a, cis=NULL))
b_metrics <- sm(graph_metrics(batch_b, cis=NULL))

## ----subset_show_plots--------------------------------------------------------
a_metrics$pc_plot
b_metrics$pc_plot
a_metrics$tsne_plot
b_metrics$tsne_plot

## ----normalize_subset, fig.show="hide"----------------------------------------
## doing nothing to the data except log2 transforming it has a surprisingly large effect
norm_test <- normalize_expt(expt, transform="log2")
l2_metrics <- sm(graph_metrics(norm_test, cis=NULL))
## a quantile normalization alone affect some, but not all of the data
norm_test <- sm(normalize_expt(expt, norm="quant"))
q_metrics <- sm(graph_metrics(norm_test, cis=NULL))  ## q for quant, who quaffed nightshade.
## cpm alone brings out some samples, too
norm_test <- sm(normalize_expt(expt, convert="cpm"))
c_metrics <- sm(graph_metrics(norm_test, cis=NULL))  ## c for cpm, who could not see the train.
## low count filtering has some effect, too
norm_test <- sm(normalize_expt(expt, filter="pofa"))
f_metrics <- sm(graph_metrics(norm_test, cis=NULL))  ## f for filter, who was hit with a spade.
## how about if we mix and match methods?
norm_test <- sm(normalize_expt(expt, transform="log2", convert="cpm",
                               norm="quant", batch="combat_scale", filter=TRUE,
                               batch_step=4, low_to_zero=TRUE))
## Some metrics are not very useful on (especially quantile) normalized data
norm_graphs <- sm(graph_metrics(norm_test, cis=NULL))

## ----view_metrics-------------------------------------------------------------
l2_metrics$pc_plot
## Also viewable with plot_pca()$plot
## PCA plots seem (to me) to prefer log2 scale data.
q_metrics$pc_plot
## only normalizing on the quantiles leaves the data open to scale effects.
c_metrics$pc_plot
## but cpm alone is insufficient
f_metrics$pc_plot
## only filtering out low-count genes is helpful as well
norm_graphs$pc_plot
## The different batch effect testing methods have a pretty widely ranging effect on the clustering
## play with them by changing the batch= parameter to:
## "limma", "sva", "svaseq", "limmaresid", "ruvg", "combat", combatmod"
knitr::kable(norm_graphs$pc_summary)
## Thus we see a dramatic decrease in variance accounted for
## by batch after applying limma's 'removebatcheffect'
## (see batch.R2 here vs. above)
norm_graphs$smc
norm_graphs$disheat  ## svaseq's batch correction seems to draw out the signal quite nicely.
## It is worth noting that the wt, early log, thy, replicate c samples are still a bit weird.
norm_graphs$tsne_plot

## ----de_test------------------------------------------------------------------
spyogenes_de <- sm(all_pairwise(expt))
## Even the lowest correlations are quite high.

## ----keeper_example-----------------------------------------------------------
my_keepers <- list(
  ## name    =   numerator / denominator
  "wt_media" = c("wt_ll_cf", "wt_ll_cg"),
  "mga_media" = c("mga_ll_cf", "mga_ll_cg"))

## ----combine_test-------------------------------------------------------------
spyogenes_tables <- sm(combine_de_tables(spyogenes_de, excel=FALSE))
summary(spyogenes_tables)
## Try changing the p-adjustment
spyogenes_tables <- sm(combine_de_tables(spyogenes_de, excel=FALSE, padj_type="BH"))
knitr::kable(head(spyogenes_tables$data[[1]]))

## ----sig_genes_test, fig.show="hide"------------------------------------------
spyogenes_sig <- sm(extract_significant_genes(spyogenes_tables, excel=FALSE))
knitr::kable(head(spyogenes_sig$limma$ups[[1]]))

## ----circos-------------------------------------------------------------------
##microbe_ids <- as.character(sm(get_microbesonline_ids("pyogenes MGAS5005")))
## A caveat!  The new version of microbesonline changed the IDs so that they no longer
## match my old rnaseq analysis!!  Thus I put my old gff file used for mapping into inst/
## and will load the annotation data from that; but I will use this data to gather
## the COG information.
mgas_gff_df <- sm(load_gff_annotations(gff=system.file("gas.gff", package="hpgltools")))
mgas_gff_df <- mgas_gff_df[-1, ]
mgas_gff_df[["sysName"]] <- gsub(pattern="Spy_", replacement="Spy", x=mgas_gff_df[["locus_tag"]])
rownames(mgas_gff_df) <- make.names(mgas_gff_df[["sysName"]], unique=TRUE)

mgas_microbes_df <- sm(load_microbesonline_annotations(id=293653))
mgas_microbes_df$sysName <- gsub(pattern="Spy_", replacement="Spy", x=mgas_microbes_df$sysName)
rownames(mgas_microbes_df) <- make.names(mgas_microbes_df$sysName, unique=TRUE)

mgas_df <- merge(x=mgas_gff_df, y=mgas_microbes_df, by="row.names")
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
circos_test <- circos_prefix(annotation=mgas_df)
## Fill it in with the data for s.pyogenes
lengths <- 1838600
names(lengths) <- "NC_007297"
## I need to manually set the ID now, as my SQL-based selections are failing...
## I need to manually set the ID now, as my SQL-based selections are failing...

circos_kary <- circos_karyotype(cfg=circos_test, lengths=lengths)
## Fill in the gene category annotations by gene-strand
circos_plus <- circos_plus_minus(cfg=circos_test)
circos_limma_hist <- circos_hist(cfg=circos_test,
                                 df=spyogenes_de$limma$all_tables[[1]],
                                 basename="limma",
                                 colname="logFC",
                                 outer=circos_plus)
circos_deseq_hist <- circos_hist(cfg=circos_test,
                                 df=spyogenes_de$deseq$all_tables[[1]],
                                 basename="deseq",
                                 colname="logFC",
                                 outer=circos_limma_hist)
circos_edger_hist <- circos_hist(cfg=circos_test,
                                 df=spyogenes_de$edger$all_tables[[1]],
                                 basename="edger",
                                 colname="logFC",
                                 outer=circos_deseq_hist)
circos_suffix(cfg=circos_test)
circos_made <- sm(circos_make(cfg=circos_test, target="mgas"))
getwd()

## ----genoplot-----------------------------------------------------------------
genoplot_chromosome()

## ----wt_mga, fig.show="hide"--------------------------------------------------
wt_mga_expt <- set_expt_conditions(expt=expt, fact="type")
wt_mga_plots <- sm(graph_metrics(wt_mga_expt))
wt_mga_norm <- sm(normalize_expt(wt_mga_expt, transform="log2", convert="raw", filter=TRUE, norm="quant"))
wt_mga_nplots <- sm(graph_metrics(wt_mga_norm))
wt_mga_de <- sm(all_pairwise(input=wt_mga_expt,
                          combined_excel="wt_mga.xlsx",
                          sig_excel="wt_mga_sig.xlsx",
                          abundant_excel="wt_mga_abundant.xlsx"))

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

