## The karyotype file is circos/data/5005_5448_karyotype.txt
## The 5005 genome is 1838562 nt. long (looking at reference/genbank/mgas_5005.gb)
## The 5448 genome is 1829516 nt. long
## (+ 1838562 1829516) 3668078
## I have been writing these text files manually, maybe I should write a function to do it...
## Except there are so many too many ways of expressing them
## So I will write a function to write individual chromosomes

#' Create the description of (a)chromosome(s) for circos.
#'
#' This function tries to save me from having to get the lengths of arcs for bacterial chromosomes
#' manually correct, and writes them as a circos compatible karyotype file.
#' The outfile parameter was chosen to match the configuration directive
#' outlined in circos_prefix(), however that will need to be changed in order for this to work
#' in variable conditions.  Next time I make one of these graphs I will do that I suspect.
#' In addition, this currently only understands how to write bacterial chromosomes, that will
#' likely be fixed when I am asked to write out a L.major karyotype.
#' These defaults were chosen because I have a chromosome of this length that is correct.
#'
#' @param name Name of the chromosome (This currently assumes a bacterial chromosome).
#' @param conf_dir Where to put the circos configuration file(s).
#' @param length Length of the chromosome (the default is mgas5005).
#' @param chr_name Short name of the chromosome.
#' @param segments How many segments to cut the chromosome into?
#' @param color Color segments of the chromosomal arc?
#' @param chr_num Number to record for each chromosome.
#' @param fasta Fasta file to use to create the karyotype.
#' @return The output filename.
#' @export
circos_karyotype <- function(name="default", conf_dir="circos/conf", length=NULL,
                             chr_name="chr1", segments=6, color="white",
                             chr_num=1, fasta=NULL) {
  genome_length <- 0
  if (is.null(length) & is.null(fasta)) {
    stop("circos_karyotype() requires a chromosome length or fasta file to gather sequence data from.")
  } else if (!is.null(length)) {
    genome_length <- length
  } else {
    raw_seq <- Rsamtools::FaFile(fasta)
    all_seq <- Biostrings::getSeq(raw_seq)
    genome_length <- sum(as.data.frame(all_seq@ranges)[["width"]])
  }

  outfile <- paste0(conf_dir, "/karyotypes/", name, ".conf")
  out <- file(outfile, open="w+")
  ## First write the summary line
  start_string <- sprintf("chr - %s %d 0 %d %s", chr_name, chr_num, genome_length, color)
  cat(start_string, file=out, sep="\n")
  for (segment in 1:segments) {
    current <- segment - 1
    begin <- floor(current * (genome_length / segments))
    end <- floor(segment * (genome_length / segments))
    string <- sprintf("band %s %d.1 %d.1 %d %d %s", chr_name, chr_num, chr_num, begin, end, color)
    cat(string, file=out, sep="\n")
  }
  close(out)
  message("Wrote karyotype to ", outfile)
  message("This should match the karyotype= line in ", name, ".conf")
  return(outfile)
}

#' Create the description of chromosome markings.
#'
#' This function writes ideogram files for circos.
#'
#' @param name Name of the configuration file to which to add the ideogram.
#' @param conf_dir Where does the configuration live?
#' @param band_url Provide a url for making these imagemaps?
#' @return The file to which the ideogram configuration was written.
#' @export
circos_ideogram <- function(name="default", conf_dir="circos/conf", band_url=NULL,
                            fill="yes", stroke_color="black",
                            thickness="20", stroke_thickness="2", fill_color="black", radius="0.85",
                            label_size="36", band_stroke_thickness="2"
                            ) {
  ideogram_outfile <- paste0(conf_dir, "/ideograms/", name, ".conf")
  out <- file(ideogram_outfile, open="w+")
  show_label <- "no"
  ideogram_string <- sprintf("## The following plot stanza describes the ideograms
<ideogram>
 <spacing>
  default = 0u
  break = 0u
 </spacing>
 # thickness (px) of chromosome ideogram
 thickness = %sp
 stroke_thickness = %s
 # ideogram border color
 stroke_color = %s
 fill = %s
 # the default chromosome color is set here and any value
 # defined in the karyotype file overrides it
 fill_color = %s
 # fractional radius position of chromosome ideogram within image
 radius = %sr
 show_label = %s
 label_font = condensedbold
 label_radius = dims(ideogram,radius) + 0.05r
 label_size = %s
 # cytogenetic bands
 band_stroke_thickness = %s
 # show_bands determines whether the outline of cytogenetic bands
 # will be seen
 show_bands = yes
 # in order to fill the bands with the color defined in the karyotype
 # file you must set fill_bands
 fill_bands = yes
", thickness, stroke_thickness, stroke_color, fill, fill_color, radius,
show_label, label_size, band_stroke_thickness)
  cat(ideogram_string, file=out, sep="")
  ideogram_band <- ""
  if (!is.null(band_url)) {
    ideogram_band <- sprintf("
 band_url = %s
 # image_map_missing_parameter=removeparam
 ideogram_url = %s
", band_url)
    cat(ideogram_band, file=out, sep="")
  }
  end_string <- "\n</ideogram>\n"
  cat(end_string, file=out, sep="")
  close(out)
  message("Wrote karyotype to ", ideogram_outfile)
  message("This should match the karyotype= line in ", name, ".conf")
  return(ideogram_outfile)
}

#' Create the ticks for a circos plot.
#'
#' This function writes ticks for circos.  This has lots of options, the defaults are all taken
#' from the circos example documentation for a bacterial genome.
#'
#' @param name Name of the configuration file to which to add the ideogram.
#' @param conf_dir Where does the configuration live?
#' @param tick_separation Top-level separation between tick marks.
#' @param min_label_distance distance to the edge of the plot for labels.
#' @param label_separations radial distance between labels.
#' @param label_offset  The offset for the labels.
#' @param label_size  Top-level label size.
#' @param multiplier  When writing the position, by what factor to lower the numbers?
#' @param main_color  Color for top-level labels?
#' @param main_thickness Top-level thickness of lines etc.
#' @param main_size Top-level size of text.
#' @param first_size Second level size of text.
#' @param first_spacing Second level spacing of ticks.
#' @param first_color Second-level text color.
#' @param first_show_label  Show a label for the second level ticks?
#' @param first_label_size Text size for second level labels?
#' @param second_size  Size of ticks for the third level.
#' @param second_spacing third-level spacing
#' @param second_color  Text color for the third level.
#' @param second_show_label  Give them a label?
#' @param second_label_size  And a size.
#' @param third_size  Now for the size of the almost-largest ticks
#' @param third_spacing  How far apart?
#' @param third_color  and their color
#' @param third_show_label  give a label?
#' @param third_label_size  and a size.
#' @param fourth_spacing The largest ticks!
#' @param fourth_color  The largest color.
#' @param fourth_show_label  Provide a label?
#' @param fourth_label_size  They are big!
#' @return The file to which the ideogram configuration was written.
#' @export
circos_ticks <- function(name="default", conf_dir="circos/conf",
                         tick_separation=2, min_label_distance=0, label_separation=5, label_offset=5,
                         label_size=8, multiplier=0.001, main_color="black", main_thickness=3,
                         main_size=20, first_size=10, first_spacing=1, first_color="black",
                         first_show_label="no", first_label_size=12, second_size=15,
                         second_spacing=5, second_color="black", second_show_label="yes",
                         second_label_size=16, third_size=18, third_spacing=10, third_color="black",
                         third_show_label="yes", third_label_size=16, fourth_spacing=100,
                         fourth_color="black", fourth_show_label="yes", fourth_label_size=36) {
  tick_outfile <- paste0(conf_dir, "/ticks_", name, ".conf")
  out <- file(tick_outfile, open="w")
  show_label <- "no"
  tick_string <- sprintf("## The following plot stanza describes the ticks
show_ticks = yes
show_tick_labels = yes
show_grid = no
grid_start = dims(ideogram,radius_inner)-0.5r
grid_end = dims(ideogram,radius_inner)
<ticks>
  skip_first_label = yes
  skip_last_label = no
  radius = dims(ideogram,radius_outer)
  tick_separation = %sp
  min_label_distance_to_edge = %sp
  label_separation = %sp
  label_offset = %sp
  label_size = %sp
  multiplier = %s
  color = %s
  thickness = %sp
  size = %sp
  <tick>
    size = %sp
    spacing = %su
    color = %s
    show_label = %s
    label_size = %sp
    format = %%.2f
    grid = no
    grid_color = lblue
    grid_thickness = 1p
  </tick>
  <tick>
    size = %sp
    spacing = %su
    color = %s
    show_label = %s
    label_size = %sp
    format = %%s
    grid = yes
    grid_color = lgrey
    grid_thickness = 1p
  </tick>
  <tick>
    size = %sp
    spacing = %su
    color = %s
    show_label = %s
    label_size = %sp
    format = %%s
    grid = yes
    grid_color = grey
    grid_thickness = 1p
  </tick>
  <tick>
    spacing = %su
    color = %s
    show_label = %s
    suffix = \" kb\"
    label_size = %sp
    format = %%s
    grid = yes
    grid_color = dgrey
    grid_thickness = 1p
  </tick>
</ticks>
",
tick_separation, min_label_distance, label_separation, label_offset, label_size, multiplier,
main_color, main_thickness, main_size, first_size, first_spacing, first_color, first_show_label,
first_label_size, second_size, second_spacing, second_color, second_show_label, second_label_size,
third_size, third_spacing, third_color, third_show_label, third_label_size, fourth_spacing,
fourth_color, fourth_show_label, fourth_label_size)
  cat(tick_string, file=out, sep="")
  close(out)
  message("Wrote ticks to ", tick_outfile)
  return(tick_outfile)
}

#' Write tiles of bacterial ontology groups using the categories from microbesonline.org.
#'
#' This function tries to save me from writing out ontology
#' definitions and likely making mistakes.  It uses the start/ends
#' from the gff annotation along with the 1 letter GO-like categories
#' from microbesonline.org.  It then writes two data files
#' circos/data/bob_plus_go.txt, circos/data/bob_minus_go.txt along
#' with two configuration files circos/conf/bob_minus_go.conf and
#' circos/conf/bob_plus_go.conf and finally adds an include to
#' circos/bob.conf
#'
#' @param table Dataframe with starts/ends and categories.
#' @param cfgout Master configuration file to write.
#' @param chr Name of the chromosome.
#' @param outer Floating point radius of the circle into which to place the plus-strand data.
#' @param width Radial width of each tile.
#' @param spacing Radial distance between outer, inner, and inner to whatever follows.
#' @param acol A color: RNA processing and modification.
#' @param bcol B color: Chromatin structure and dynamics.
#' @param ccol C color: Energy production conversion.
#' @param dcol D color: Cell cycle control, mitosis and meiosis.
#' @param ecol E color: Amino acid transport metabolism.
#' @param fcol F color: Nucleotide transport and metabolism.
#' @param gcol G color: Carbohydrate transport and metabolism.
#' @param hcol H color: Coenzyme transport and metabolism.
#' @param icol I color: Lipid transport and metabolism.
#' @param jcol J color: Translation, ribosome structure and biogenesis.
#' @param kcol K color: Transcription.
#' @param lcol L color: Replication, recombination, and repair.
#' @param mcol M color: Cell wall/membrane biogenesis.
#' @param ncol N color: Cell motility
#' @param ocol O color: Posttranslational modification, protein turnover, chaperones.
#' @param pcol P color: Inorganic ion transport and metabolism.
#' @param qcol Q color: Secondary metabolite biosynthesis, transport, and catabolism.
#' @param rcol R color: General function prediction only.
#' @param scol S color: Function unknown.
#' @param tcol T color: Signal transduction mechanisms.
#' @param ucol U color: Intracellular trafficking(sp?) and secretion.
#' @param vcol V color: Defense mechanisms.
#' @param wcol W color: Extracellular structures.
#' @param xcol X color: Not in COG.
#' @param ycol Y color: Nuclear structure.
#' @param zcol Z color: Cytoskeleton.
#' @return Radius after adding the plus/minus information and the spacing between them.
#' @export
circos_plus_minus <- function(table, cfgout="circos/conf/default.conf", chr="chr1",
                              outer=1.0, width=0.08, spacing=0.0,
                              acol="orange", bcol="reds-9-seq", ccol="yellow", dcol="vlpurple",
                              ecol="vlgreen", fcol="dpblue", gcol="vlgreen", hcol="vlpblue",
                              icol="vvdpgreen", jcol="dpred", kcol="orange", lcol="vvlorange",
                              mcol="dpgreen", ncol="vvlpblue", ocol="vvlgreen", pcol="vvdpred",
                              qcol="ylgn-3-seq", rcol="vlgrey", scol="grey", tcol="vlpurple",
                              ucol="greens-3-seq", vcol="vlred", wcol="vvdppurple", xcol="black",
                              ycol="lred", zcol="vlpblue") {
  end_col <- "end"
  if (!is.null(table[["stop"]])) {
    end_col <- "stop"
  }
  number_pluses <- sum(table[["strand"]] == "+")
  number_ones <- sum(table[["strand"]] == 1)
  plus_string <- "+"
  minus_string <- "-"
  if (number_pluses + number_ones == 0) {
    stop("This function requires some way of understanding plus/minus strand.")
  } else if (number_ones > 0) {
    plus_string <- 1
    minus_string <- -1
  }
  if (is.null(table[["start"]]) | is.null(table[[end_col]]) |
      is.null(table[["strand"]]) | is.null(table[["COGFun"]])) {
    stop("This function assumes an input table including the columns: 'start', 'stop', 'strand', and 'COGFun'")
  }
  plus_cfg_file <- cfgout
  minus_cfg_file <- cfgout
  plus_cfg_file <- gsub(".conf$", "_plus_go.conf", plus_cfg_file)
  minus_cfg_file <- gsub(".conf$", "_minus_go.conf", minus_cfg_file)
  ## What I should do is spend some time thinking and reformat this to handle
  ## and arbitrary number of arbitrary columns so that I have some flexibility later.
  if (is.null(table[["id"]])) {
    table <- table[, c("start", "stop", "strand", "COGFun")]
  } else {
    table <- table[, c("start", "stop", "strand", "COGFun", "id")]
  }
  go_plus <- as.data.frame(table[table[["strand"]] == plus_string, ])
  go_minus <- as.data.frame(table[table[["strand"]] == minus_string, ])
  go_plus[["chr"]] <- chr
  go_minus[["chr"]] <- chr
  if (is.null(table[["id"]])) {
    go_plus <- go_plus[, c("chr", "start", "stop", "COGFun")]
    go_minus <- go_minus[, c("chr", "start", "stop", "COGFun")]
    go_plus[["go"]] <- paste0("value=", go_plus[["COGFun"]], "0")
    go_minus[["go"]] <- paste0("value=", go_minus[["COGFun"]], "0")
  } else {
    go_plus <- go_plus[, c("chr", "start", "stop", "COGFun", "id")]
    go_minus <- go_minus[, c("chr", "start", "stop", "COGFun", "id")]
    go_plus[["id"]] <- paste0("id=", go_plus[["id"]])
    go_minus[["id"]] <- paste0("id=", go_minus[["id"]])
  }
  data_prefix <- cfgout
  data_prefix <- gsub("/conf/", "/data/", data_prefix)
  data_prefix <- gsub(".conf$", "", data_prefix)

  plus_file <- paste0(data_prefix, "_plus_go.txt")
  message("Writing data file: ", plus_file, " with the + strand GO data.")
  write.table(go_plus, file=plus_file, quote=FALSE, row.names=FALSE, col.names=FALSE, na="no_go")

  minus_file <- paste0(data_prefix, "_minus_go.txt")
  message("Writing data file: ", minus_file, " with the - strand GO data.")
  write.table(go_minus, file=minus_file, quote=FALSE, row.names=FALSE, col.names=FALSE, na="no_go")

  first_outer <- outer
  first_inner <- first_outer - width
  plus_cfg_out <- file(plus_cfg_file, open="w+")
  plus_cfg_filename <- gsub("^circos/", "", plus_file)
  plus_cfg_string <- sprintf("## The following plot stanza describes the broad ontological categories
## Of genes on the plus strand.
 <plot>
  type = tile
  file = %s
  layers = 1
  layers_overflow = hide
  margin = 0.00u
  thickness = 95
  padding = 0
  orientation = out
  stroke_thickness = 0
  stroke_color = black
  color = green
  r1 = %sr
  r0 = %sr
  url = script?type=label&value=[id]&color=[color]
  <rules>
   <rule>
    condition = var(value) =~ \"^A\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^B\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^C\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^D\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^E\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^F\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^G\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^H\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^I\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^J\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^K\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^L\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^M\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^N\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^O\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^P\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^Q\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^R\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^S\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^T\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^U\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^V\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^W\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^X\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^Y\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^Z\"
    fill_color = %s
    color = %s
   </rule>
  </rules>
 </plot>
", plus_cfg_filename, first_outer, first_inner,
acol, acol, bcol, bcol, ccol, ccol, dcol, dcol,
ecol, ecol, fcol, fcol, gcol, gcol, hcol, hcol,
icol, icol, jcol, jcol, kcol, kcol, lcol, lcol,
mcol, mcol, ncol, ncol, ocol, ocol, pcol, pcol,
qcol, qcol, rcol, rcol, scol, scol, tcol, tcol,
ucol, ucol, vcol, vcol, wcol, wcol, xcol, xcol,
ycol, ycol, zcol, zcol)
  cat(plus_cfg_string, file=plus_cfg_out, sep="")
  close(plus_cfg_out)

  ## Now move the ring in one width and print the minus strand.
  second_outer <- first_inner - spacing
  second_inner <- second_outer - width
  minus_cfg_out <- file(minus_cfg_file, open="w+")
  minus_cfg_filename <- gsub("^circos/", "", minus_file)
  minus_cfg_string <- sprintf("## The following plot stanza describes the broad ontological categories
## Of genes on the minus strand.
 <plot>
  type = tile
  file = %s
  layers = 1
  layers_overflow = hide
  margin = 0.00u
  thickness = 95
  padding = 0
  orientation = out
  stroke_thickness = 0
  stroke_color = black
  color = green
  r1 = %sr
  r0 = %sr
  url = script?type=label&value=[id]&color=[color]
  <rules>
   <rule>
    condition = var(value) =~ \"^A\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^B\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^C\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^D\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^E\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^F\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^G\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^H\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^I\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^J\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^K\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^L\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^M\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^N\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^O\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^P\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^Q\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^R\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^S\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^T\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^U\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^V\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^W\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^X\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^Y\"
    fill_color = %s
    color = %s
   </rule>
   <rule>
    condition = var(value) =~ \"^Z\"
    fill_color = %s
    color = %s
   </rule>
  </rules>
 </plot>
", minus_cfg_filename, second_outer, second_inner,
acol, acol, bcol, bcol, ccol, ccol, dcol, dcol,
ecol, ecol, fcol, fcol, gcol, gcol, hcol, hcol,
icol, icol, jcol, jcol, kcol, kcol, lcol, lcol,
mcol, mcol, ncol, ncol, ocol, ocol, pcol, pcol,
qcol, qcol, rcol, rcol, scol, scol, tcol, tcol,
ucol, ucol, vcol, vcol, wcol, wcol, xcol, xcol,
ycol, ycol, zcol, zcol)
  cat(minus_cfg_string, file=minus_cfg_out, sep="")
  close(minus_cfg_out)
  message("Wrote the +/- config files.  Appending their inclusion to the master file.")

  master_cfg_out <- file(cfgout, open="a+")
  plus_cfg_include <- plus_cfg_filename
  plus_cfg_include <- gsub("^circos/", "", plus_cfg_file)
  minus_cfg_include <- minus_cfg_filename
  minus_cfg_include <- gsub("^circos/", "", minus_cfg_file)
  master_cfg_string <- sprintf("## The +/- ontology rings.
## The plus strand ring extends from %sr to %sr
<<include %s>>
## The minus strand ring extends from %sr to %sr
<<include %s>>

", first_outer, first_inner, plus_cfg_include, second_outer, second_inner, minus_cfg_include)
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  message("Returning the inner width: ", second_inner, ".  Use it as the outer for the next ring.")
  new_outer <- second_inner - spacing
  return(new_outer)
}

#' Write tiles of arbitrary categorical point data in circos.
#'
#' This function tries to make the writing circos tiles easier.  Like
#' circos_plus_minus() and circos_hist() it works in 3 stages,
#' It writes out a data file using cfgout as a basename and the data
#' from df in the circos histogram format into
#' circos/data/bob_tile.txt
#' It then writes out a configuration plot stanza in
#' circos/conf/bob_tile.conf
#' and finally adds an include to circos/bob.conf
#'
#' @param df Dataframe with starts/ends and the floating point information.
#' @param annot_df Annotation data frame defining starts/stops.
#' @param cfgout Master configuration file to write.
#' @param colname Name of the column with the data of interest.
#' @param chr Name of the chromosome (This currently assumes a bacterial chromosome)
#' @param colors Colors of the data.
#' @param outer Floating point radius of the circle into which to place the categorical data.
#' @param width Width of each tile.
#' @param spacing Radial distance between outer, inner, and inner to whatever follows.
#' @return Radius after adding the histogram and the spacing.
#' @export
circos_tile <- function(df, annot_df=NULL, cfgout="circos/conf/default.conf", colname="logFC",
                        chr="chr1", basename="", colors=NULL,
                        outer=0.9, width=0.08, spacing=0.0) {
  ## I am going to have this take as input a data frame with genes as rownames
  ## starts, ends, and functional calls
  ## I will tell R to print out a suitable stanza for circos while I am at it
  ## because I am tired of mistyping something stupid.
  if (is.null(annot_df)) {
    full_table <- df
  } else {
    full_table <- merge(df, annot_df, by.x="row.names", by.y="row.names")
  }

  start_name <- "start"
  stop_name <- "stop"
  if (is.null(full_table[[stop_name]])) {
    stop_name <- "end"
  }

  if (! start_name %in% colnames(full_table)) {
    stop("This requires a column named start.")
  }
  if (! stop_name %in% colnames(full_table)) {
    stop("This requires a column named ", stop_name, ".")
  }
  if (! colname %in% colnames(full_table)) {
    stop("This requires a column named ", colname, ".")
  }
  if (is.null(rownames(full_table))) {
    stop("This requires rownames.")
  }
  full_table <- full_table[, c(start_name, stop_name, colname)]
  non_empty <- full_table[[colname]] != ""
  full_table <- full_table[non_empty, ]
  if (is.null(full_table[[start_name]]) | is.null(full_table[[stop_name]]) |
      is.null(rownames(full_table)) | is.null(full_table[[colname]])) {
    stop("This requires columns: start, stop, rownames, and datum")
  }

  datum_cfg_file <- cfgout
  datum_cfg_file <- gsub(".conf$", "", datum_cfg_file)
  datum_cfg_file <- paste0(datum_cfg_file, "_", colname, "_tile.conf")
  full_table[["chr"]] <- chr
  full_table <- full_table[, c("chr", start_name, stop_name, colname)]
  data_prefix <- cfgout
  data_prefix <- gsub("/conf/", "/data/", data_prefix)
  data_prefix <- gsub(".conf$", "", data_prefix)
  data_filename <- paste0(data_prefix, "_", basename, colname, "_tile.txt")
  message("Writing data file: ", data_filename, " with the ", colname, " column.")
  write.table(full_table, file=data_filename, quote=FALSE, row.names=FALSE, col.names=FALSE)

  num_colors <- 1
  if (is.null(colors)) {
    conditions <- levels(as.factor(full_table[[colname]]))
    num_colors <- length(conditions)
    colors <- suppressWarnings(grDevices::colorRampPalette(
                                            RColorBrewer::brewer.pal(num_colors, "Dark2"))(num_colors))
    names(colors) <- conditions
  } else {
    num_colors <- length(colors)
  }

  ## Now write the config stanza
  inner <- outer - width
  data_cfg_out <- file(datum_cfg_file, open="w+")
  data_cfg_filename <- gsub("^circos/", "", datum_cfg_file)
  data_relative_filename <- gsub(pattern="circos/", replacement="", x=data_filename)
  data_cfg_string <- sprintf("## This is a circos tile configuration stanza.
 <plot>
  type = tile
  file = %s
  layers = 1
  layers_overflow = hide
  margin = 0.00u
  thickness = 90
  padding = 1
  orientation = out
  stroke_thickness = 0
  stroke_color = black
  color = black
  fill_color = black
  r1 = %fr
  r0 = %fr
  <rules>
", data_relative_filename, outer, inner)
  cat(data_cfg_string, file=data_cfg_out, sep="")
  for (c in 1:num_colors) {
    red_component <- paste0("0x", substr(colors[[c]], 2, 3))
    green_component <- paste0("0x", substr(colors[[c]], 4, 5))
    blue_component <- paste0("0x", substr(colors[[c]], 5, 6))
    red_component <- strtoi(red_component)
    green_component <- strtoi(green_component)
    blue_component <- strtoi(blue_component)
    color_string <- paste0(red_component, ",", blue_component, ",", green_component)
    new_string <- sprintf("   <rule>
    condition = var(value) =~ \"^%s\"
    fill_color = %s
    color = %s
   </rule>
", names(colors)[[c]], color_string, color_string)
    cat(new_string, file=data_cfg_out, sep="")
  }
  end_string <- sprintf("  </rules>\n </plot>", sep="")
  cat(end_string, file=data_cfg_out, sep="")
  close(data_cfg_out)

  ## Now add to the master configuration file.
  master_cfg_out <- file(cfgout, open="a+")
  data_cfg_include <- data_cfg_filename
  data_cfg_include <- gsub("^circos/", "", data_cfg_include)
  master_cfg_string <- sprintf("
## The tile ring for %s
<<include %s>>
", colname, data_cfg_include)
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  new_outer <- inner - spacing
  return(new_outer)
}

#' Write tiles of arbitrary heat-mappable data in circos.
#'
#' This function tries to make the writing circos heatmaps easier.  Like
#' circos_plus_minus() and circos_hist() it works in 3 stages,
#' It writes out a data file using cfgout as a basename and the data
#' from df in the circos histogram format into
#' circos/data/bob_heatmap.txt
#' It then writes out a configuration plot stanza in
#' circos/conf/bob_heatmap.conf
#' and finally adds an include to circos/bob.conf
#'
#' @param df Dataframe with starts/ends and the floating point information.
#' @param annot_df Annotation data frame with starts/ends.
#' @param cfgout Master configuration file to write.
#' @param colname Name of the column with the data of interest.
#' @param chr Name of the chromosome (This currently assumes a bacterial chromosome).
#' @param basename Make sure the written configuration files get different names with this.
#' @param colors Colors of the heat map.
#' @param outer Floating point radius of the circle into which to place the heatmap.
#' @param width Width of each tile in the heatmap.
#' @param spacing Radial distance between outer, inner, and inner to whatever follows.
#' @return Radius after adding the histogram and the spacing.
#' @export
circos_heatmap <- function(df, annot_df, cfgout="circos/conf/default.conf", colname="logFC",
                           chr="chr1", basename="", colors=NULL, color_choice="spectral-9-div",
                           scale_log_base=1, outer=0.9, width=0.08, spacing=0.02) {
  ## I am going to have this take as input a data frame with genes as rownames
  ## starts, ends, and functional calls
  ## I will tell R to print out a suitable stanza for circos while I am at it
  ## because I am tired of mistyping something stupid.
  full_table <- merge(df, annot_df, by="row.names")
  if (nrow(full_table) == 0) {
    stop("Merging the annotations and data failed.")
  }
  start_colnames <- colnames(full_table)
  new_colnames <- gsub(x=start_colnames, pattern="\\.x$", replacement="")
  colnames(full_table) <- new_colnames

  start_name <- "start"
  stop_name <- "stop"
  if (is.null(full_table[[stop_name]])) {
    stop_name <- "end"
  }

  if (! start_name %in% colnames(full_table)) {
    stop("This requires a column named start.")
  }
  if (! stop_name %in% colnames(full_table)) {
    stop("This requires a column named ", stop_name, ".")
  }
  if (! colname %in% colnames(full_table)) {
    stop("This requires a column named ", colname, ".")
  }
  if (is.null(rownames(full_table))) {
    stop("This requires rownames.")
  }

  if (is.null(full_table[["id"]])) {
    full_table <- full_table[, c(start_name, stop_name, colname)]
  } else {
    full_table <- full_table[, c(start_name, stop_name, colname, "id")]
    full_table[["id"]] <- paste0("id=", full_table[["id"]])
  }
  full_table[[start_name]] <- as.numeric(full_table[[start_name]])
  na_drop <- ! is.na(full_table[[start_name]])
  full_table <- full_table[na_drop, ]

  if (is.null(full_table[[start_name]]) | is.null(full_table[[stop_name]]) |
      is.null(rownames(full_table)) | is.null(full_table[[colname]])) {
    stop("This requires columns: start, stop, rownames, and datum")
  }
  datum_cfg_file <- cfgout
  datum_cfg_file <- gsub(".conf$", "", datum_cfg_file)
  datum_cfg_file <- paste0(datum_cfg_file, "_", basename, colname, "_heatmap.conf")
  full_table[["chr"]] <- chr
  if (is.null(full_table[["id"]])) {
    full_table <- full_table[, c("chr", start_name, stop_name, colname)]
  } else {
    full_table <- full_table[, c("chr", start_name, stop_name, colname, "id")]
  }
  data_prefix <- cfgout
  data_prefix <- gsub("/conf/", "/data/", data_prefix)
  data_prefix <- gsub(".conf$", "", data_prefix)
  data_filename <- paste0(data_prefix, "_", basename, colname, "_heatmap.txt")
  message("Writing data file: ", data_filename, " with the ", basename, colname, " column.")
  write.table(full_table, file=data_filename, quote=FALSE, row.names=FALSE, col.names=FALSE)

  num_colors <- 1
  if (is.null(colors)) {
    conditions <- levels(as.factor(full_table[[colname]]))
    num_colors <- length(conditions)
    colors <- suppressWarnings(grDevices::colorRampPalette(
                                            RColorBrewer::brewer.pal(num_colors, "Dark2"))(num_colors))
    names(colors) <- conditions
  } else {
    num_colors <- length(conditions)
  }

  ## Now write the config stanza
  minval <- min(full_table[[colname]])
  maxval <- max(full_table[[colname]])
  inner <- outer - width
  data_cfg_out <- file(datum_cfg_file, open="w+")
  data_cfg_filename <- gsub("^circos/", "", datum_cfg_file)
  data_relative_filename <- gsub(pattern="circos/", replacement="", x=data_filename)
  data_cfg_string <- sprintf("## This is a circos heatmap configuration stanza.
 <plot>
  type = heatmap
  file = %s
##  pattern = hline,vline  # Also solid
  color_mapping = 0 ## also 1, 2
  min = %s
  max = %s
  stroke_thickness = 0
  r1 = %fr
  r0 = %fr
  color = %s
  scale_log_base = %s
  url = script?type=label&value=[id]
 </plot>", data_relative_filename, minval, maxval, outer, inner, color_choice, scale_log_base)
  cat(data_cfg_string, file=data_cfg_out, sep="")
  close(data_cfg_out)

  ## Now add to the master configuration file.
  master_cfg_out <- file(cfgout, open="a+")
  data_cfg_include <- data_cfg_filename
  data_cfg_include <- gsub("^circos/", "", data_cfg_include)
  master_cfg_string <- sprintf("
## The tile ring for %s
<<include %s>>
", colname, data_cfg_include)
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  new_outer <- inner - spacing
  return(new_outer)
}

#' Write histograms of arbitrary floating point data in circos.
#'
#' This function tries to make the writing of histogram data in circos
#' easier.  Like circos_plus_minus() it works in 3 stages,
#' It writes out a data file using cfgout as a basename and the data
#' from df in the circos histogram format into
#' circos/data/bob_hist.txt
#' It then writes out a configuration plot stanza in
#' circos/conf/bob_hist.conf
#' and finally adds an include to circos/bob.conf
#'
#' @param df Dataframe with starts/ends and the floating point information.
#' @param annot_df Annotation data frame containing starts/ends.
#' @param cfgout Master configuration file to write.
#' @param colname Name of the column with the data of interest.
#' @param chr Name of the chromosome (This currently assumes a bacterial chromosome).
#' @param basename Location to write the circos data (usually cwd).
#' @param color Color of the plotted data.
#' @param fill_color Guess!
#' @param outer Floating point radius of the circle into which to place the data.
#' @param width Radial width of each tile.
#' @param spacing Distance between outer, inner, and inner to whatever follows.
#' @return Radius after adding the histogram and the spacing.
#' @export
circos_hist <- function(df, annot_df, cfgout="circos/conf/default.conf", colname="logFC",
                        chr="chr1", basename="", color="blue", fill_color="blue",
                        outer=0.9, width=0.08, spacing=0.0) {
  ## I am going to have this take as input a data frame with genes as rownames
  ## starts, ends, and functional calls
  ## I will tell R to print out a suitable stanza for circos while I am at it
  ## because I am tired of mistyping something stupid.
  full_table <- merge(df, annot_df, by="row.names")
  if (nrow(full_table) == 0) {
    stop("Merging the annotations and data failed.")
  }
  start_colnames <- colnames(full_table)
  new_colnames <- gsub(x=start_colnames, pattern="\\.x$", replacement="")
  colnames(full_table) <- new_colnames

  start_name <- "start"
  stop_name <- "stop"
  if (is.null(full_table[[stop_name]])) {
    stop_name <- "end"
  }

  if (! start_name %in% colnames(full_table)) {
    stop("This requires a column named start.")
  }
  if (! stop_name %in% colnames(full_table)) {
    stop("This requires a column named ", stop_name, ".")
  }
  if (! colname %in% colnames(full_table)) {
    stop("This requires a column named ", colname, ".")
  }
  if (is.null(rownames(full_table))) {
    stop("This requires rownames.")
  }

  full_table <- full_table[, c(start_name, stop_name, colname)]
  full_table[[start_name]] <- as.integer(full_table[[start_name]])
  full_table[[stop_name]] <- as.integer(full_table[[stop_name]])
  keep_idx <- !is.na(full_table[[start_name]])
  full_table <- full_table[keep_idx, ]

  datum_cfg_file <- cfgout
  datum_cfg_file <- gsub(".conf$", "", datum_cfg_file)
  datum_cfg_file <- paste0(datum_cfg_file, "_", basename, colname, "_hist.conf")
  full_table[["chr"]] <- chr
  full_table <- full_table[, c("chr", start_name, stop_name, colname)]
  data_prefix <- cfgout
  data_prefix <- gsub("/conf/", "/data/", data_prefix)
  data_prefix <- gsub(".conf$", "", data_prefix)
  data_filename <- paste0(data_prefix, "_", basename, colname, "_hist.txt")
  message("Writing data file: ", data_filename, " with the ", basename, colname, " column.")
  write.table(full_table, file=data_filename, quote=FALSE, row.names=FALSE, col.names=FALSE)

  num_colors <- 1
  ## if (is.null(colors)) {
  ##     conditions <- levels(as.factor(df[["call"]]))
  ##     num_colors <- length(conditions)
  ##     colors <- suppressWarnings(grDevices::colorRampPalette(
  ##                                RColorBrewer::brewer.pal(num_colors, "Dark2"))(num_colors))
  ##     names(colors) <- conditions
  ## }

  ## Now write the config stanza
  inner <- outer - width
  data_cfg_out <- file(datum_cfg_file, open="w+")
  data_cfg_filename <- gsub("^circos/", "", datum_cfg_file)
  data_relative_filename <- gsub(pattern="circos/", replacement="", x=data_filename)
  data_cfg_string <- sprintf("
 <plot>
  type = histogram
  file = %s
  extend_bin = no
  fill_under = yes
  fill_color = %s
  color = %s
  thickness = 0
  r1 = %fr
  r0 = %fr
  orientation = out
 </plot>
", data_relative_filename, fill_color, color, outer, inner)
  cat(data_cfg_string, file=data_cfg_out, sep="")
  close(data_cfg_out)

  ## Now add to the master configuration file.
  master_cfg_out <- file(cfgout, open="a+")
  data_cfg_include <- data_cfg_filename
  data_cfg_include <- gsub("^circos/", "", data_cfg_include)
  master_cfg_string <- sprintf("
## The histogram ring for %s
<<include %s>>

", colname, data_cfg_include)
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  new_outer <- inner - spacing
  return(new_outer)
}

#' Write a simple makefile for circos.
#'
#' I regenerate all my circos pictures with make(1).  This is my makefile.
#'
#' @param target Default make target.
#' @param output Makefile to write.
#' @param circos Location of circos.  I have a copy in home/bin/circos and use that sometimes.
#' @return a kitten
#' @export
circos_make <- function(target="", output="circos/Makefile", circos="circos") {
  circos_dir <- dirname(output)
  if (!file.exists(circos_dir)) {
    message("The circos directory does not exist, creating: ", circos_dir)
    dir.create(circos_dir, recursive=TRUE)
  }
  if (!file.exists("circos/etc")) {
    system("ln -s /etc/circos circos/etc")
  }
  out <- file(output, open="w+")
  makefile_string <- sprintf("
.PHONY:\tclean
CIRCOS=\"%s\"

clean:
\trm -rf conf data *.conf *.png *.svg *.html

%%.png:\t%%.conf
\t$(CIRCOS) -conf $< -outputfile $*.png

%%.svg:\t%%.conf
\t$(CIRCOS) -conf $< -outputfile $*.svg

%%:\t%%.conf
\t$(CIRCOS) -conf $< -outputfile $*.png
\t$(CIRCOS) -conf $< -outputfile $*.svg
\techo '<img src=\"$*.svg\" usemap=\"#$*\">' > map.html
\tcat $*.html >> map.html
\tmv map.html $*.html

", circos)
  cat(makefile_string, file=output, sep="")
  ## close(output)

  make_target <- gsub(pattern="circos/conf/", replacement="", x=target)
  make_target <- gsub(pattern="\\.conf", replacement="", x=make_target)
  make_target_svg <- paste0(make_target, ".svg")
  make_target_png <- paste0(make_target, ".png")

  make_command <- paste0("cd circos && touch Makefile && make ", make_target_svg, " 2>>make.out 1>&2 && make ",
                         make_target, "")
  result <- system(make_command) ##, show.output.on.console=FALSE)

  return(result)
}

#' Write arcs between chromosomes in circos.
#'
#' Ok, so when I said I only do 1 chromosome images, I lied.
#' This function tries to make writing arcs between chromosomes easier.
#' It too works in 3 stages,
#' It writes out a data file using cfgout as a basename and the data
#' from df in the circos arc format into
#' circos/data/bob_arc.txt
#' It then writes out a configuration plot stanza in
#' circos/conf/bob_arc.conf
#' and finally adds an include to circos/bob.conf
#'
#' In its current implementation, this only understands two
#' chromosomes.  A minimal amount of logic and data organization will
#' address this weakness.
#'
#' @param df Dataframe with starts/ends and the floating point information.
#' @param cfgout Master configuration file to write.
#' @param first_col Name of the first chromosome.
#' @param second_col Name of the second chromosome.
#' @param color Color of the chromosomes.
#' @param radius Outer radius at which to add the arcs.
#' @param thickness Integer thickness of the arcs.
#' @return The file to which the arc configuration information was written.
#' @export
circos_arc <- function(df, cfgout="circos/conf/default.conf", first_col="chr1", second_col="chr2",
                       color="blue", radius=0.75, thickness=3) {
  if (is.null(df$start) | is.null(df$end) | is.null(rownames(df)) |
      is.null(df[[first_col]]) | is.null(df[[second_col]])) {
    stop("This requires columns: start, end, rownames, and datum")
  }
  datum_cfg_file <- cfgout
  datum_cfg_file <- gsub(".conf$", "", datum_cfg_file)
  datum_cfg_file <- paste0(datum_cfg_file, "_arc.conf")
  first_name <- paste0(first_col, "_name")
  second_name <- paste0(second_col, "_name")
  first_start_name <- paste0(first_col, "_start")
  first_end_name <- paste0(first_col, "_end")
  second_start_name <- paste0(second_col, "_start")
  second_end_name <- paste0(second_col, "_end")
  message("This function assumes an input table including columns: ", first_start_name,
                 ",", first_end_name, ",", second_start_name, ",", second_end_name, ",",
                 first_name, ",", second_name, ".")
  df <- df[, c(first_name, second_name, first_start_name,
               first_end_name, second_start_name, second_end_name)]
  data_prefix <- cfgout
  data_prefix <- gsub("/conf/", "/data/", data_prefix)
  data_prefix <- gsub(".conf$", "", data_prefix)
  data_filename <- paste0(data_prefix, "_", first_col, "_arc.txt")
  message("Writing data file: ", data_filename, " with the ", first_col, " column.")
  print_arc <- function(x) {
    cat(x[5], " chr5005 ", x[1], " ", x[2], "\n", x[5], " chr5448 ", x[3], " ", x[4], "\n\n",
        file="circos/data/crossref_5005_5448.txt", append=TRUE, sep="")
  }
  file.remove(data_filename, showWarnings=FALSE) ## To avoid appending forever.
  apply(df, 1, print_arc)

  ## Now write the config stanza
  ## I just realized that there are the possibility of multiple
  ## link stanzas just like multiple plot stanzas...
  ## well, deal with that later.
  data_cfg_out <- file(datum_cfg_file, open="w+")
  data_cfg_filename <- gsub("^circos/", "", data_filename)
  data_cfg_string <- sprintf("
<links>
 z = 0
 radius = %sr
 <link>
  ribbon = yes
  show = yes
  file = %s
  color = %s
  thickness = %s
  </link>
 </links>
", radius, data_cfg_string, color, thickness)
  cat(data_cfg_string, file=data_cfg_out, sep="")
  close(data_cfg_out)

  ## Now add to the master configuration file.
  master_cfg_out <- file(cfgout, open="a+")
  data_cfg_include <- data_cfg_filename
  data_cfg_include <- gsub("^circos/", "", data_cfg_include)
  master_cfg_string <- sprintf("
## The histogram ring for %s
<<include %s>>

", first_col, data_cfg_include)
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  return(radius)
}

#' Write the beginning of a circos configuration file.
#'
#' A few parameters need to be set when starting circos.  This sets
#' some of them and gets ready for plot stanzas.
#'
#' In its current implementation, this really assumes that there will
#' be no highlight stanzas and at most 1 link stanza.
#' chromosomes.  A minimal amount of logic and data organization will
#' address these weaknesses.
#'
#' @param name Name of the map, called with 'make name'.
#' @param conf_dir Directory containing the circos configuration data.
#' @param radius Size of the image.
#' @param band_url Place to imagemap link.
#' @return The master configuration file name.
#' @export
circos_prefix <- function(name="mgas", conf_dir="circos/conf", radius=1800, band_url=NULL) {
  message("This assumes you have a colors.conf in circos/colors/ and fonts.conf in circos/fonts/")
  message("It also assumes you have conf/ideogram.conf, conf/ticks.conf, and conf/housekeeping.conf")
  cfgout <- paste0(conf_dir, "/", name, ".conf")
  message("It will write ", cfgout, " with a reasonable first approximation config file.")

  data_dir <- gsub(pattern="conf", replacement="data", x=conf_dir)
  if (!file.exists(data_dir)) {
    message("Creating the data directory: ", data_dir)
    dir.create(data_dir, recursive=TRUE)
  }
  if (!file.exists(conf_dir)) {
    message("The circos directory does not exist, creating: ", conf_dir)
    dir.create(conf_dir, recursive=TRUE)
  }
  karyotype_dir <- paste0(conf_dir, "/karyotypes")
  ideogram_dir <- paste0(conf_dir, "/ideograms")
  if (!file.exists(karyotype_dir)) {
    message("The karyotype directory does not exist, creating: ", karyotype_dir)
    dir.create(karyotype_dir, recursive=TRUE)
  }
  if (!file.exists(paste0(conf_dir, "/ideograms"))) {
    message("The ideogram directory does not exist, creating: ", ideogram_dir)
    dir.create(ideogram_dir, recursive=TRUE)
  }
  karyotype_file <- gsub("circos/conf", "conf/karyotypes", cfgout)
  ideogram_file <- gsub("circos/conf", "conf/ideograms", cfgout)
  etc_file <- paste0(path.package("hpgltools"), "/circos/circos_etc.tar.xz")
  etc_cmd <- paste0("tar -C ", dirname(conf_dir), " -xavf ", etc_file, " 2>/dev/null 1>&2")
  system(command=etc_cmd)

  ## If you want clickable ideograms, add band_url='script?start=[start]&end=[end]&label=[label]
  ideogram_filename <- circos_ideogram(name=name, conf_dir=conf_dir, band_url=band_url)
  tick_filename <- circos_ticks(name=name, conf_dir=conf_dir)
  tick_file <- gsub(x=tick_filename, pattern="^circos/", replacement="")

  out <- file(cfgout, open="w+")
  prefix_string <- sprintf("## This is the prefix of a circos configuration file written by hpgltools.
<colors>
 <<include colors.conf>>
</colors>
<fonts>
 <<include fonts/fonts.conf>>
</fonts>

## The ideograms are generated by circos_ideogram() in R.
<<include %s>>
## The ticks are generated by circos_ticks() in R.
<<include %s>>
<<include housekeeping.conf>>
karyotype = %s

<image>
 image_map_use = yes
 image_map_missing_parameter = removeurl
 dir = .
 radius = %sp
 background = white
 angle_offset = -90
</image>

chromosomes_units = 1000
chromosomes_display_default = yes

<highlights>
</highlights>

<plots>
", ideogram_file, tick_file, karyotype_file, radius)
  cat(prefix_string, file=out, sep="")
  close(out)
  to_path <- paste0(name, ".conf")
  if (!file.exists(to_path)) {
    wd <- getwd()
    tmpwd <- paste0(wd, "/circos")
    setwd(tmpwd)
    from <- gsub("circos/", "", cfgout)
    to <- "."
    file.symlink(from, to)
    setwd(wd)
  }
  return(cfgout)
}

#' Write the end of a circos master configuration.
#'
#' circos configuration files need an ending.  This writes it.
#'
#' @param cfgout Master configuration file to write.
#' @return The filename of the configuration.
#' @export
circos_suffix <- function(cfgout="circos/conf/default.conf") {
  out <- file(cfgout, open="a+")
  suffix_string <- "</plots>"
  cat(suffix_string, file=out, sep="\n")
  close(out)
}

## EOF
