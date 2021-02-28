#' Pipe operator
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' data.table's funky column assignment operator
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name :=
#' @rdname column_assignment
#' @keywords internal
#' @export
#' @importFrom data.table :=
NULL

#' dopar
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name %dopar%
#' @rdname dopar
#' @keywords internal
#' @export
#' @importFrom foreach %dopar%
NULL

#' hpgltools: a suite of tools to make our analyses easier
#'
#' This provides a series of helpers for working with sequencing data
#'
#' It falls under a few main topics
#'
#' \itemize{
#' \item Data exploration, look for trends in sequencing data and identify batch
#'       effects or skewed distributions.
#' \item Differential expression analyses, use DESeq2/limma/EdgeR in a hopefully
#'       robust and flexible fashion.
#' \item Ontology analyses, use goseq/clusterProfiler/topGO/GOStats/gProfiler in
#'       hopefully robust ways.
#' \item Perform some simple TnSeq analyses.
#' }
#'
#' To see examples of this in action, check out the vignettes:
#' \code{browseVignettes(package = 'hpgltools')}
#'
#' @docType package
#' @name hpgltools
#' @importFrom Biobase exprs pData fData notes sampleNames
#' @importFrom data.table data.table
#' @importFrom dplyr filter group_by n summarise
#' @importFrom foreach foreach
#' @importFrom ggplot2 aes aes_string ggplot
#' @importFrom glue glue glue_data
#' @importFrom grDevices recordPlot
#' @importFrom rlang abort sym
#' @importFrom stats
#'  aggregate as.dendrogram as.formula ave biplot coef coefficients complete.cases
#'  cor cor.test density dist dnorm formula glm hclust lm lowess median
#'  model.matrix na.omit order.dendrogram p.adjust p.adjust.methods pnorm
#'  princomp quantile relevel reorder resid residuals rnbinom sd setNames
#'  t.test var
#' @import graphics
#' @import grDevices
#' @import methods
#' @import utils
NULL

#' The following sets the ggplot2 default text size.
base_size <- 16

#' Set the xlsx table style
table_style <- "TableStyleMedium9"

#' R CMD check is super annoying about :::.
#'
#' In a fit of pique, I did a google search to see if anyone else has been
#' annoyed in the same way as was I.  Yihui Xie was, and in his email to r-devel
#' in 2013 he proposed a game of hide-and-seek; which I am repeating here.
#'
#' This just implements ::: as an infix operator that will not trip check.
#'
#' @param pkg on the left hand side
#' @param fun on the right hand side
`%:::%` <- function(pkg, fun) {
  get(fun, envir = asNamespace(pkg), inherits = FALSE)
}

getMaintainer <- "GenomicFeatures" %:::% ".getMaintainer"
getMetaDataValue <- "GenomicFeatures" %:::% ".getMetaDataValue"
getTxDbVersion <- "GenomicFeatures" %:::% ".getTxDbVersion"
normAuthor <- "GenomicFeatures" %:::% ".normAuthor"
sortCols <- "variancePartition" %:::% ".sortCols"
aprior <- "sva" %:::% "aprior"
bprior <- "sva" %:::% "bprior"
it.sol <- "sva" %:::% "it.sol"
int.eprior <- "sva" %:::% "int.eprior"
getGOLevel <- "clusterProfiler" %:::% "getGOLevel"

## Set up some datasets/variables for working with the examples/tests/vignettes.

#' Group A Streptococcus pyogenes strain 5005 GFF file.
#'
#' This gff contains the annotations which correspond to NC_007297.
#' I am not sure if I want to include the fasta file.
example_gff <- system.file("share", "gas.gff", package = "hpgltools")

#' Group A Streptococcus pyogenes strain 5005 genome.
#'
#' Yeah, what the heck, I am not likely to try to include this in CRAN, I can
#' waste some space.
example_fasta <- system.file("share", "spyogenes_5005.fasta", package = "hpgltools")

#' Data file containing an expressionset of some GAS 5005 RNASeq data.
#'
#' This rda file contains an expt which comprises a subset of a RNASeq
#' experiment performed with the McIver lab.  The goal was to examine changes in
#' transcription across two media types (the very permissive THY and restrictive
#' CDM).  In addition, this used two strains, one which is more and less
#' pathogenic.
cdm_expt_rda <- system.file("share", "cdm_expt.rda", package = "hpgltools")

#' Pseudomonas areuginosa strain PA14 sample sheet.
#'
#' The sample sheet from an experiment with another ESKAPE pathogen,
#' Pesudomonas!  This is emblematic of how I like to organize samples.  The most
#' relevant columns for creating an expressionset with create_expt() include:
#' 'Sample ID', 'Condition', 'Batch', and 'file'.  This actually provides a
#' subset of an experiment in which we were looking simultaneously at the
#' 'large' and 'small' RNA populations in two PA strains, one of which is
#' deficient in an oligonucleotide degradation enzyme 'orn'.  We were also
#' seeking to find changes from exponential growth to stationary.  The portions
#' of the experiment included in this sample sheet are only 3 replicates of the
#' large RNA samples.
pa_sample_sheet <- system.file("share", "pa_samples.xlsx", package = "hpgltools")

#' PA14 annotations in a GFF file
#'
#' Along with the sample sheet, this provides everything required to make an
#' expressionsheet via create_expt() for our PA 14 experiment.
pa_gff <- system.file("share", "paeruginosa_pa14.gff", package = "hpgltools")

#' PA14 genome
#'
#' This is the genome which corresponds to pa_gff.
pa_fasta <- system.file("share", "paeruginosa_pa14.fasta", package = "hpgltools")

#' PA14 genbank flat file
#'
#' This may be used instead to extract genomic information and annotations via
#' load_genbank_annotations().
pa_gb <- system.file("share", "paeruginosa.pa14.gb", package = "hpgltools")

#' The count directory from the Pseudomonas experiment subset.
#'
#' In other examples I include an archive file of the counts, here are the raw
#' tables.  These were created via bowtie2 -> samtools -> htseq-count.
pa_counts <- system.file("share", "counts", package = "hpgltools")

#' Subset of a human RNASeq expressionset.
#'
#' This is a portion of an expressionset used to examine changes caused by
#' infection with Leishmania panamensis.
hs_expt_rda <- system.file("share", "hs_expt.rda", package = "hpgltools")

#' Extra fonts and useful bits and bobs for working with circos.
#'
#' This is a tarball of the circos etc/ directory from my Debian linux
#' installation.  I discovered to my annoyance that other systems were missing
#' the fonts required to make circos plots work properly along with something
#' else.  In a fit of pique I tarred them up and left them here.
circos_etc <- system.file("share", "circos", "circos_etc.tar.xz", package = "hpgltools")

#' Sample sheet from a portion of an RNASeq experiment of Trypanosoma cruzi CL-Brener.
#'
#' This contains a portion of an experiment performed with Santuza in which we
#' were comparing two closely related T.cruzi strains: CL-14 and CL-Brener,
#' which are super-similar but vastly different in terms of their
#' pathogenicity.  (I would much rather be infected with CL-14!).
clbr_sample_sheet <- system.file("share", "clbr", "clbr_samples_combined.xlsx",
                                 package = "hpgltools")

#' Archive of the count tables from the CL-Brener experiment.
#'
#' These tables were created via TopHat2 mapping of the CL-Br and CL-14 samples
#' using the CL-Brener genome.
clbr_count_tables <- system.file("share", "clbr", "clbr_counts.tar", package = "hpgltools")

#' CL-Brener GFF file containing all haplotypes.
#'
#' One of the interesting and bizarre things about CL-Brener: it is a
#' multi-haplotype strain, containing bits and pieces from two lineages, named
#' 'Esmeraldo' and 'Non-Esmeraldo'.  In addition, there is a large portion which
#' has not been characterized and is therefore called 'Unassigned'.  Thus, when
#' mapping the data, we create a concatenated genome with all three haplotypes.
clbr_gff <- system.file("share", "clbr", "clbrener_8.1_complete_genes.gff.gz",
                        package = "hpgltools")

#' CL-Brener/CL-14 output from vcfutils.
#'
#' One thing we did not include in the paper (because I didn't think of it until
#' later), was an analysis of the single nucleotide variants between the two
#' strains.  RNASeq data is of course not an ideal format for performing these
#' analyses, but I think I figured out a reasonable method to extract mostly
#' robust differences (in snp.r).
clbr_vcf_output <- system.file("share", "clbr", "vcfutils_output.tar",
                               package = "hpgltools")

#' Sample sheet used for a portion of a Group B Streptococcal TNSeq experiment.
#'
#' TNSeq is sort of the inverse of RNASeq, one is instead looking for the genes
#' _not_ represented in the dataset.  This sample sheet lays out the
#' experimental design for an in vitro TNSeq experiment from Streptococcus
#' agalactiae strain CJB111.
gbs_sample_sheet <- system.file("share", "gbs_tnseq", "sagalactiae_samples.xlsx",
                                package = "hpgltools")

#' GFF annotations used for the GBS TNSeq experiment.
#'
#' At the time of the experiment, there was not a very good genome for this
#' strain.  We therefore chose to use strain A909 as the reference.  Since then,
#' Lindsey's lab made a complete genome, though the annotations remain a bit
#' sparse.
gbs_gff <- system.file("share", "gbs_tnseq", "sagalactiae_a909.gff",
                       package = "hpgltools")

#' Genome for the GBS TNSeq experiment.
#'
#' If you read this far, you know what this is.
gbs_fasta <- system.file("share", "gbs_tnseq", "sagalactiae_a909.fasta",
                         package = "hpgltools")

#' Count tables of the mapped GBS TNSeq DNA fragments.
#'
#' One thing I like to do with TNSeq data is to treat it similarly to RNASeq
#' data in order to get a sense of the changing 'fitness' of each gene.  The
#' data has all the same distribution attributes of a RNASeq dataset, after all;
#' so why not use the same plots and tests to see if it is valid?
gbs_counts <- system.file("share", "gbs_tnseq", "gbs_essentiality_counts.tar",
                          package = "hpgltools")

#' Outputs from the DeJesus Essentiality package from the GBS TNSeq experiment.
#'
#' The Essentiality package from the DeJesus lab uses a Bayesian framework to
#' look for genes which are essential in a TNSeq experiment.  In my pipeline, I
#' invoke this tool with multiple parameters in an attempt to find the
#' parameters which provide the most likely 'true' result.  This archive
#' contains those results for our GBS TNSeq experiment.
gbs_essentiality <- system.file("share", "gbs_tnseq", "gbs_essentiality.tar",
                                package = "hpgltools")

#' Wig files used as input for the DeJesus Essentiality package.
#'
#' My preprocessing pipeline uses the bam alignments from bowtie to extract all
#' reads which start/end on a mariner insertion site (TA) and count how many
#' occured at every position of the genome.  These files are the result of that
#' process, thus each line is the position of a 'T' in the 'TA' followed by the
#' number of reads which start/end with it.
gbs_wig <- system.file("share", "gbs_tnseq", "gbs_essentiality_wig.tar",
                       package = "hpgltools")

#' The UTR regions of every highly translated L.major gene.
#'
#' Once upon a time I performed a ribosome profiling experiment in Leishmania
#' major.  Sadly, we still have not published it.  This file contains the UTRs
#' of every highly translated gene in procyclic promastigotes.  I used it in
#' some motif analyses for fun.
lm_high_utr <- system.file("share", "motif", "pro_high.fasta",
                           package = "hpgltools")

#' Portion of a sample sheet used in a DIA-SWATH proteomics experiment.
#'
#' In this experiment, Dr. Briken sought to learn about proteins which are
#' exported by Mycobacterium tuberculosis.  He therefore performed a DIA SWATH
#' experiment using two strains and collected the supernatant fraction (to get
#' the exported proteins) and the intracellular fraction.  This file contains
#' the metadata for a portion of that experiment.
dia_samples <- system.file("share", "mtb_prot", "dia_samples.ods", package = "hpgltools")

#' A few scored tsv files from Dr. Briken's DIA-SWATH experiment.
#'
#' I used a fairly exhaustive set of open source tools to interpret Dr. Briken's
#' data.  These files comprise the endpoint of the preprocessing and the inputs
#' for the R package 'SWATH2stats'.  This file is excessively large, but the
#' smallest by far of the various inputs I wanted to include to test my various
#' proteomics functions.
dia_counts <- system.file("share", "mtb_prot", "sb_prot.tar", package = "hpgltools")

#' The Rmats alternative splicing results for a Mycobacterium tuberculosis experiment.
#'
#' This is the comparison of the splicing isoforms observed in the host
#' following infection with a few Mtb strains.  Our goal was to learn how
#' reliable the various alternative splicing quantification tools are. (spoiler:
#' not very).
mtb_rmats <- system.file("share", "mtb_rmats.tar.xz", package = "hpgltools")

#' The suppa alternative splicing results for a Mycobacterium tuberculosis experiment.
#'
#' This is the comparison of the splicing isoforms observed in the host
#' following infection with a few Mtb strains.  Our goal was to learn how
#' reliable the various alternative splicing quantification tools are. (spoiler:
#' not very).
mtb_suppa <- system.file("share", "mtb_suppa.tar.xz", package = "hpgltools")

#' Portion of the RNASeq results from Solanum betaceum.
#'
#' I had the opportunity to work the Sandra Correia, she was awesome.  She was
#' seeking to learn about differences among embryogenic cells in the Tree
#' tomato.  I therefore got to learn first-hand a tiny portion of what is meant
#' when one says 'plant genetics are hard.'  I had it far easier than Sandra.  I
#' just used Trinity to make some de-novo transcriptomes and attempt to provide
#' some metrics about which ones are real and really different across conditions
#' in her experiment.  Her work was many thousands of times more difficult.
sb_data <- system.file("share", "sb", "preprocessing.tar", package = "hpgltools")

#' Portion of the raw trinotate 'annotation' output for Solanum betaceum.
#'
#' The full S.betaceum trinotate annotation is quite large, so I just pulled a
#' portion as an example for this package.
sb_annot <- system.file("share", "sb", "trinotate_head.csv.xz", package = "hpgltools")

#' Name of the Drosphila melanogaster orgdb.
#'
#' Some of my examples use this and I forgot the make it explicit...
dm_orgdb <- "org.Dm.eg.db"
