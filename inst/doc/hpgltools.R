## ----startup-------------------------------------------------------------

## This block serves to load requisite libraries and set some options.
library("hpgltools")
## To set up an initial vignette, use the following line:
## devtools::use_vignette("hpgltools")
autoloads_all()
opts_knit$set(progress=TRUE, verbose=TRUE, purl=FALSE, error=TRUE, stop_on_error=FALSE, fig.width=7, fig.height=7)
options(java.parameters="-Xmx8g")  ## used for xlconnect -- damn 4g wasn't enough
theme_set(theme_bw(base_size=10))
set.seed(1)


## ----render_vignette, eval=FALSE-----------------------------------------
#  
#  load("RData")
#  rm(list=ls())
#  save(list=ls(all=TRUE), file="RData")
#  render("hpgltools.Rmd", output_format="pdf_document")
#  render("hpgltools.Rmd", output_format="html_document")
#  

## ----annotation_information, eval=FALSE----------------------------------
#  
#  tcruzi_annotations = import.gff3("reference/gff/clbrener_8.1_complete.gff.gz")
#  annotation_info = as.data.frame(tcruzi_annotations)
#  
#  genes = annotation_info[annotation_info$type=="gene",]
#  gene_annotations = genes
#  rownames(genes) = genes$Name
#  tooltip_data = genes
#  tooltip_data = tooltip_data[,c(11,12)]
#  tooltip_data$tooltip = paste(tooltip_data$Name, tooltip_data$description, sep=": ")
#  tooltip_data$tooltip = gsub("\\+", " ", tooltip_data$tooltip)
#  rownames(tooltip_data) = tooltip_data$Name
#  tooltip_data = tooltip_data[-1]
#  tooltip_data = tooltip_data[-1]
#  colnames(tooltip_data) = c("name.tooltip")
#  head(tooltip_data)
#  

## ----sample_sheet, results='asis'----------------------------------------

samples = read.csv("data/all_samples.csv")
knitr::kable(head(samples))


## ----creating_experiment-------------------------------------------------

example_data = counts(make_exampledata(ngenes=10000, columns=24))
## create_expt() usually expects that there are a bunch of count tables
## from htseq in the directory: processed_data/count_tables/
## These may be organised in separate directories by condition(type)
## in one directory each by sample.  By default, this assumes they will be
## named sample_id.count.gz, but this may be changed with the suffix argument.
all_expt = create_expt("data/all_samples.csv", count_dataframe=example_data)


## ----examining_data------------------------------------------------------

## graph_metrics() performs the following:
## runs a libsize plot, non-zero genes plot, boxplot, correlation/distance heatmaps, and pca plots
## It performs a normalization of the data (log2(quantile(cpm)) by default), and does it again
## It then uses limma's removeBatchEffect() to make a stab at removing batch effect, and does it again.

## An important thing to remember: the data from makeExampleData() is not very interesting, so the resulting
## plots are also not interesting...
fun = graph_metrics(expt=all_expt)
fun
## The following are some examples of other ways to make use of these plots:

##fun_boxplot = hpgl_boxplot(df=fun)
##print(fun_boxplot)
##log_boxplot = hpgl_boxplot(df=fun, scale="log")
##print(log_boxplot)
##hpgl_corheat(df=fun, colors=hpgl_colors)
##hpgl_disheat(df=fun, colors=hpgl_colors)
##hpgl_smc(df=fun, colors=hpgl_colors)
##hpgl_libsize(df=fun)
##hpgl_qq_all(df=fun)


## ----normalize_data------------------------------------------------------

## normalize_expt will do this on the expt class, replace the expressionset therein, and
## make a backup of the data inside the expt class.
norm_expt = normalize_expt(all_expt)
head(exprs(norm_expt$expressionset))
## size factor, tmm, rle, upperQuartile all require a design matrix.
norm_boxplot = hpgl_boxplot(expt=norm_expt)
print(norm_boxplot)
norm_disheat = hpgl_disheat(expt=norm_expt)
print(norm_disheat)


## ----etc-----------------------------------------------------------------

## el_subset means to pull out only those samples which represent 'Early Log' growth.
el_subset = expt_subset(norm_expt, "stage=='EL'")
## Conversely, one may pull samples which are early log and also wild type
elwt_subset = expt_subset(norm_expt, "stage=='EL'&type=='WT'")
## These subsets may be characterized with the plots as above
## Here is a qq plot as an example.
elwt_qqs = hpgl_qq_all(expt=elwt_subset)
## Simple comparison will take the first condition as control and the second
## as experimental, if we look at el_subset, we will see that means conditions
## 'a' and 'b'.  Thus performing simple_comparison will look for differentially
## expressed genes between them.
head(el_subset$design)
ab_comparison = simple_comparison(el_subset)
## A summary of the data will show the data provided:
## The following plots and pieces of data show the output provided by simple_comparison()
## This function isn't really intended to be used, but provides a reference point for performing other analyses.
summary(ab_comparison)
print(ab_comparison$amean_histogram)  ## A histogram of the per-gene mean values
print(ab_comparison$coef_amean_cor)   ## The correlation of the means (should not be significant)
print(ab_comparison$coefficient_scatter) ## A scatter plot of condition b with respect to a
print(ab_comparison$coefficient_x) ## A histogram of the gene abundances of a
print(ab_comparison$coefficient_y) ## A histogram of the gene abundances of b
print(ab_comparison$coefficient_both) ## A histogram of the gene abundances of a and b
## Note to self, I keep meaning to change the colors of that to match the others
print(ab_comparison$coefficient_lm) ## The description of the line which describes the relationship
## of all of the genes in a to those in b
print(ab_comparison$coefficient_lmsummary) ## A summary of the robust linear model in coefficient_lm
## This has some neat things like the R-squared value and the parameters used to arrive at the linear model.
## ab_comparison$coefficient_weights ## a list of weights by gene, bigger weights mean closer to the linear model.
## ab_comparison$comparisons ## the raw output from limma
print(ab_comparison$contrasts)  ## The output from limma's makeContrasts()
print(ab_comparison$contrast_histogram)  ## A histogram of the values of b-a for each gene
head(ab_comparison$downsignificant)  ## The list of genes which are significantly down in b vs a
dim(ab_comparison$downsignificant)
## ab_comparison$fit ## the result from lmFit()
print(ab_comparison$ma_plot)  ## An ma plot of b vs a
print(ab_comparison$pvalue_histogram) ## A histogram of the p-values, one would hope to see a spike in the low numbers
head(ab_comparison$table) ## The full contrast table
head(ab_comparison$upsignificant)  ## The list of genes which are significantly up in b vs a
dim(ab_comparison$upsignificant)
print(ab_comparison$volcano_plot) ## A Volcano plot of b vs a
## ab_comparison$voom_data  ## The output from voom()
print(ab_comparison$voom_plot) ## A ggplot2 version of the mean/variance trend provided by voom()

## The data structure ab_comparison$comparisons contains the output from eBayes() which comprises the last
## limma step...
funkytown = write_limma(data=ab_comparison$comparisons, excel=FALSE, csv=FALSE)
## Lets make up some gene lengths
gene_lengths = funkytown[[1]]
gene_lengths$width = sample(nrow(gene_lengths))
gene_lengths$ID = rownames(gene_lengths)
gene_lengths = gene_lengths[,c("ID","width")]

## And some GO categories
goids=funkytown[[1]]
all_go_categories = AnnotationDbi::keys(GO.db)
goids$GO = sample(all_go_categories, nrow(gene_lengths))
goids$ID = rownames(goids)
goids = goids[,c("ID","GO")]

ontology_fun = limma_ontology(funkytown, gene_lengths=gene_lengths, goids=goids, n=100, overwrite=TRUE)

testme = head(funkytown[[1]], n=40)
tt = simple_clusterprofiler(testme, goids=goids, gff=goids)
ttt = cluster_trees(testme, tt)
tttt = simple_topgo(testme)


## ----example_acb, eval=FALSE---------------------------------------------
#  
#  
#  ## acb stands for "kept_conditions_batches"  which takes too long to
#  ## type when setting up the contrasts.
#  acb = paste0(kept_qcpml2$conditions, kept_qcpml2$batches)
#  kept_data = exprs(kept_qcpml2$expressionset)
#  table(acb)
#  ## The invocation of table() keptows me to count up the contribution of
#  ## each condition/batch combination to the whole data set.
#  
#  ## Doing this (as I understand it) means I do nothave to worry about
#  ## balanced samples so much, but must be more careful to understand
#  ## the relative contribution of each sample type to the entire data
#  ## set.
#  
#  complete_model = model.matrix(~0 + acb)
#  complete_fit = lmFit(kept_data, complete_model)
#  complete_voom = hpgl_voom(kept_data, complete_model)
#  complete_voom$plot
#  complete_model
#  ## This is an example of what happens when I have heterogenous numbers of samples
#  ## on each side of a contrast, so that a normal design matrix of conditions + batches
#  ## would not work, so instead I add up the contributions of each batch (capital letters)
#  ## and average them out, then use the resulting terms in the various contrasts below.
#  epi_cl14 = "acbcl14_epiF"
#  epi_clbr = "acbclbr_epiE"
#  tryp_cl14 = "(acbcl14_trypB + acbcl14_trypD + acbcl14_trypG) / 3"
#  tryp_clbr = "acbclbr_trypG"
#  a60_cl14 =  "(acbcl14_a60A * 2/3) + (acbcl14_a60B * 1/3)"
#  a60_clbr = "acbclbr_a60A"
#  a96_cl14 = "acbcl14_a96C"
#  a96_clbr = "acbclbr_a96C"
#  epi_cl14clbr = paste0("(",epi_cl14,")", "  -  ", "(",epi_clbr,")")
#  tryp_cl14clbr = paste0("(",tryp_cl14,")", "  -  ", "(",tryp_clbr,")")
#  a60_cl14clbr = paste0("(",a60_cl14,")", "  -  ", "(",a60_clbr,")")
#  a96_cl14clbr = paste0("(",a96_cl14,")", "  -  ", "(",a96_clbr,")")
#  epitryp_cl14 = paste0("(",tryp_cl14,")", "  -  ", "(",epi_cl14,")")
#  epitryp_clbr = paste0("(",tryp_clbr,")", "  -  ", "(",epi_clbr,")")
#  epia60_cl14 = paste0("(",a60_cl14,")", "  -  ", "(",epi_cl14,")")
#  epia60_clbr = paste0("(",a60_clbr,")", "  -  ", "(",epi_clbr,")")
#  a60a96_cl14 = paste0("(",a96_cl14,")", "  -  ", "(",a60_cl14,")")
#  a60a96_clbr = paste0("(",a96_clbr,")", "  -  ", "(",a60_clbr,")")
#  a60tryp_cl14 = paste0("(",tryp_cl14,")", "  -  ", "(",a60_cl14,")")
#  a60tryp_clbr = paste0("(",tryp_clbr,")", "  -  ", "(",a60_clbr,")")
#  ## The following contrast is messed up in some as of yet unknown way.
#  epitryp_cl14clbr = paste0("(",epitryp_cl14,")", "  -  ", "(",epitryp_clbr,")")
#  ## So I will add some more contrasts using data which doesn't get screwed up
#  epia60_cl14clbr = paste0("(",epia60_cl14,")", "  -  ", "(",epia60_clbr,")")
#  a60tryp_cl14clbr = paste0("(",a60tryp_cl14,")", "  -  ", "(",a60tryp_clbr,")")
#  a60a96_cl14clbr = paste0("(",a60a96_cl14,")", "  -  ", "(",a60a96_clbr,")")
#  
#  complete_contrasts_v2 = makeContrasts(
#      epi_cl14=epi_cl14,
#      epi_clbr=epi_clbr,
#      tryp_cl14=tryp_cl14,
#      tryp_clbr=tryp_clbr,
#      a60_cl14=a60_cl14,
#      a60_clbr=a60_clbr,
#      a96_cl14=a96_cl14,
#      a96_clbr=a96_clbr,
#      epi_cl14clbr=epi_cl14clbr,
#      tryp_cl14clbr=tryp_cl14clbr,
#      a60_cl14clbr=a60_cl14clbr,
#      a96_cl14clbr=a96_cl14clbr,
#      epitryp_cl14=epitryp_cl14,
#      epitryp_clbr=epitryp_clbr,
#      epia60_cl14=epia60_cl14,
#      epia60_clbr=epia60_clbr,
#      a60a96_cl14=a60a96_cl14,
#      a60a96_clbr=a60a96_clbr,
#      a60tryp_cl14=a60tryp_cl14,
#      a60tryp_clbr=a60tryp_clbr,
#      epitryp_cl14clbr=epitryp_cl14clbr,
#      epia60_cl14clbr=epia60_cl14clbr,
#      a60tryp_cl14clbr=a60tryp_cl14clbr,
#      a60a96_cl14clbr=a60a96_cl14clbr,
#      levels=complete_voom$design)
#  ## This colnames() is annoyingly necessary to avoid really obnoxious contrast names.
#  colnames(complete_contrasts_v2) = c("epi_cl14","epi_clbr","tryp_cl14","tryp_clbr","a60_cl14","a60_clbr","a96_cl14","a96_clbr","epi_cl14clbr","tryp_cl14clbr","a60_cl14clbr","a96_cl14clbr","epitryp_cl14","epitryp_clbr","epia60_cl14","epia60_clbr","a60tryp_cl14","a60tryp_clbr","a60a96_cl14","a60a96_clbr","epitryp_cl14clbr","epia60_cl14clbr","a60tryp_cl14clbr","a60a96_cl14clbr")
#  kept_fits = contrasts.fit(complete_fit, complete_contrasts_v2)
#  kept_comparisons = eBayes(kept_fits)
#  

## ----acb_balanced, eval=FALSE--------------------------------------------
#  
#  all_data = exprs(norm_expt$expressionset)
#  complete_model = model.matrix(~0 + all_human_expt$conditions + all_human_expt$batches)
#  ## Shorten the column names of the model so I don't have to type so much later...
#  tmpnames = colnames(complete_model)
#  tmpnames = gsub("all_human_expt[[:punct:]]","", tmpnames)
#  tmpnames = gsub("conditions","", tmpnames)
#  colnames(complete_model) = tmpnames
#  rm(tmpnames)
#  
#  complete_voom = hpgl_voom(all_data, complete_model)
#  complete_voom$plot
#  complete_fit = lmFit(complete_voom, complete_model)
#  
#  all_contrasts = makeContrasts(
#      ## Start with the simple coefficient groupings for each condition
#      none4=none4,
#      none24=none24,
#      none48=none48,
#      none72=none72,
#      bead4=bead4,
#      bead24=bead24,
#      bead48=bead48,
#      bead72=bead72,
#      maj4=maj4,
#      maj24=maj24,
#      maj48=maj48,
#      maj72=maj72,
#      ama4=ama4,
#      ama24=ama24,
#      ama48=ama48,
#      ama72=ama72,
#      ## Now do a few simple comparisons
#      ## compare beads to uninfected
#      beadnone_4=bead4-none4,
#      beadnone_24=bead24-none24,
#      beadnone_48=bead48-none48,
#      beadnone_72=bead72-none72,
#      majnone_4=maj4-none4,
#      majnone_24=maj24-none24,
#      majnone_48=maj48-none48,
#      majnone_72=maj72-none72,
#      amanone_4=ama4-none4,
#      amanone_24=ama24-none24,
#      amanone_48=ama48-none48,
#      amanone_72=ama72-none72,
#      ## compare samples to beads
#      majbead_4=maj4-bead4,
#      majbead_24=maj24-bead24,
#      majbead_48=maj48-bead48,
#      majbead_72=maj72-bead72,
#      amabead_4=ama4-bead4,
#      amabead_24=ama24-bead24,
#      amabead_48=ama48-bead48,
#      amabead_72=ama72-bead72,
#      ## (x-z)-(a-b)
#      ## Use this to compare major and amazonensis
#      amamaj_bead_4=(ama4-bead4)-(maj4-bead4),
#      amamaj_bead_24=(ama24-bead24)-(maj24-bead24),
#      amamaj_bead_48=(ama48-bead48)-(maj48-bead48),
#      amamaj_bead_72=(ama72-bead72)-(maj72-bead72),
#      ## (c-d)-(e-f) where c/d are: (amazon|major/none)/(beads/none)
#      majbead_none_4=(maj4-none4)-(bead4-none4),
#      majbead_none_24=(maj24-none24)-(bead24-none24),
#      majbead_none_48=(maj48-none48)-(bead48-none48),
#      majbead_none_72=(maj72-none72)-(bead72-none72),
#      amabead_none_4=(ama4-none4)-(bead4-none4),
#      amabead_none_24=(ama24-none24)-(bead24-none24),
#      amabead_none_48=(ama48-none48)-(bead48-none48),
#      amabead_none_72=(ama72-none72)-(bead72-none72),
#      levels=complete_voom$design)
#  all_fits = contrasts.fit(complete_fit, all_contrasts)
#  all_comparisons = eBayes(all_fits)
#  limma_list = write_limma(data=all_comparisons)
#  
#  all_table = topTable(all_comparisons, adjust="fdr", n=nrow(all_data))
#  write.csv(all_comparisons, file="excel/all_tables.csv")
#  ## write_limma() is a shortcut for writing out all the data structures
#  all_comparison_tables = write_limma(all_comparisons, excel=FALSE)
#  

## ----ontology_searches---------------------------------------------------

ontology_info = read.csv(file="data/trinotate_go_trimmed.csv.gz", header=FALSE, sep="\t")
##ontology_info = read.csv(file="data/transcript_go.csv.gz", header=FALSE, sep="\t")
colnames(ontology_info) = c("gene_id","transcript_id","group","startend","blast_go","pfam_go")
## Drop any entries which don't have a putative length
ontology_info = subset(ontology_info, startend != 0)
## Split the column 'startend' into two columns by the '-' sign
ontology_info = as.data.frame(transform(ontology_info, startend=reshape::colsplit(startend, split="\\-", names=c("start","end"))))
## Make the resulting pieces into two separate columns, start and end.
ontology_info$start = ontology_info$startend$start
ontology_info$end = ontology_info$startend$end
## Use start and end to make length
ontology_info$length = abs(ontology_info$start - ontology_info$end)
## Drop the unneeded columns
ontology_info = ontology_info[,c("gene_id","transcript_id","group","start","end","length","blast_go","pfam_go")]
head(ontology_info)

## goseq() requires mappings between ID/length and ID/GO category
## Currently I have my toy set to assume column names, which is admittedly stupid.
gene_lengths = ontology_info[,c("transcript_id","length")]
colnames(gene_lengths) = c("ID","width")
split_go = ontology_info[,c("transcript_id","blast_go")]
split_go$blast_go = as.character(split_go$blast_go)

## The following few lines were pulled from the internet
## they serve to generate a data structure in the format expected by goseq()
## It simply splits all space separated GO categories into separate rows
## with the same ID
require.auto("splitstackshape")
id_go = concat.split.multiple(split_go, "blast_go", seps=" ", "long")
id_go = as.data.frame(id_go)
colnames(id_go) = c("ID","GO")
go_ids = subset(id_go, GO != 0)

## Pull out all entries from group 1
group_one = subset(ontology_info, group == "1")
group_one = group_one[,c("transcript_id","start","end")]
colnames(group_one) = c("ID","start","end")

## Perform the goseq() analysis
group_one_go = simple_goseq(group_one, lengths=gene_lengths, goids=go_ids)
group_one_go$pvalue_histogram
head(group_one_go$godata_interesting)
head(group_one_go$mf_subset)
group_one_go$mfp_plot
group_one_go$bpp_plot
group_one_go$ccp_plot

## Print trees of the goseq() data
initial_trees = goseq_trees(group_one, group_one_go, goids_df=go_ids)
initial_trees$MF
initial_trees$BP
initial_trees$CC


## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

