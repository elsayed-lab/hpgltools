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
#' @param cfg Result of circos_prefix(), contains a bunch of useful material.
#' @param df Dataframe with starts/ends and the floating point information.
#' @param first_col Name of the first chromosome.
#' @param second_col Name of the second chromosome.
#' @param color Color of the chromosomes.
#' @param radius Outer radius at which to add the arcs.
#' @param thickness Integer thickness of the arcs.
#' @param ribbon Print as a ribbon?
#' @param show Show these arcs?
#' @param z Correction parameter.
#' @return The file to which the arc configuration information was written.
#' @export
circos_arc <- function(cfg, df, first_col="seqnames", second_col="seqnames.2",
                       color="blue", radius=0.75, thickness=3, ribbon="yes", show="yes", z="0") {
  annot <- cfg[["annot"]]

  arc_cfg_file <- cfg[["cfg_file"]]
  arc_cfg_file <- gsub(pattern=".conf$", replacement="", x=arc_cfg_file)
  arc_cfg_file <- paste0(arc_cfg_file, first_col, "_arc.conf")
  arc_data_file <- file.path(cfg[["data_dir"]], basename(arc_cfg_file))
  arc_data_file <- gsub(pattern=".conf$", replacement=".txt", x=arc_data_file)

  first_name <- glue::glue("{first_col}_name")
  second_name <- glue::glue("{second_col}_name")
  first_start_name <- glue::glue("{first_col}_start")
  first_end_name <- glue::glue("{first_col}_end")
  second_start_name <- glue::glue("{second_col}_start")
  second_end_name <- glue::glue("{second_col}_end")
  message("This function assumes an input table including columns: ", first_start_name,
          ",", first_end_name, ",", second_start_name, ",", second_end_name, ",",
          first_name, ",", second_name, ".")
  df <- df[, c(first_name, second_name, first_start_name,
               first_end_name, second_start_name, second_end_name)]

  message("Writing data file: ", arc_data_file, " with the ", first_col, " column.")
  print_arc <- function(x) {
    cat(x[5], " chr5005 ", x[1], " ", x[2], "\n", x[5], " chr5448 ", x[3], " ", x[4], "\n\n",
        file="circos/data/crossref_5005_5448.txt", append=TRUE, sep="")
  }
  file.remove(arc_data_file, showWarnings=FALSE) ## To avoid appending forever.
  apply(df, 1, print_arc)

  ## Now write the config stanza
  ## I just realized that there are the possibility of multiple
  ## link stanzas just like multiple plot stanzas...
  ## well, deal with that later.
  data_cfg_string <- glue::glue("
<links>
 z = {z}
 radius = {radius}r
 <link>
  ribbon = {ribbon}
  show = {show}
  file = {data_cfg_string}
  color = {color}
  thickness = {thickness}
  </link>
 </links>

")
  data_cfg_out <- file(arc_cfg_file, open="w+")
  cat(data_cfg_string, file=data_cfg_out, sep="")
  close(data_cfg_out)

  rel_cfg_file <- file.path("conf", basename(arc_cfg_file))
  rel_data_file <- file.path("data", basename(arc_data_file))
  master_cfg_string <- glue::glue("
  ## The histogram ring for {first_col}
  <<include {rel_cfg_file}>>

")
  master_cfg_out <- file(cfg[["cfg_file"]], open="a+")
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  return(radius)
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
#' @param cfg Result of circos_prefix(), contains a bunch of useful material.
#' @param df Dataframe with starts/ends and the floating point information.
#' @param colname Name of the column with the data of interest.
#' @param color_mapping 0 means no overflows for min/max, 1 means overflows
#'        of min get a chosen color, 2 means overflows of both min/max get
#'        chosen colors.
#' @param min_value Minimum value for the data.
#' @param max_value Maximum value for the data.
#' @param basename Make sure the written configuration files get different names with this.
#' @param colors Colors of the heat map.
#' @param color_choice  Name of the heatmap to use, I forget how this interacts with color...
#' @param scale_log_base  Defines how the range of colors will be ranged with
#'   respect to the values in the data.
#' @param outer Floating point radius of the circle into which to place the heatmap.
#' @param rules some extra rules?
#' @param width Width of each tile in the heatmap.
#' @param spacing Radial distance between outer, inner, and inner to whatever follows.
#' @return Radius after adding the histogram and the spacing.
#' @export
circos_heatmap <- function(cfg, df, colname="logFC",
                           color_mapping=0, min_value=NULL, max_value=NULL,
                           basename="", colors=NULL,
                           color_choice="spectral-9-div", scale_log_base=1, outer=0.9, rules=NULL,
                           width=0.08, spacing=0.02) {
  annot <- cfg[["annot"]]
  full_table <- merge(df, annot, by="row.names")
  if (nrow(full_table) == 0) {
    stop("Merging the annotations and data failed.")
  }
  start_colnames <- colnames(full_table)
  new_colnames <- gsub(x=start_colnames, pattern="\\.x$", replacement="")
  colnames(full_table) <- new_colnames
  rownames(full_table) <- full_table[["Row.names"]]
  full_table[["Row.names"]] <- NULL

  full_table <- full_table[, c("chr", "start", "stop", colname)]
  start_undefined_idx <- full_table[["start"]] == "undefined"
  full_table <- full_table[!start_undefined_idx, ]
  stop_undefined_idx <- full_table[["stop"]] == "undefined"
  full_table <- full_table[!stop_undefined_idx, ]
  full_table[["start"]] <- as.numeric(full_table[["start"]])
  full_table[["stop"]] <- as.numeric(full_table[["stop"]])
  keep_idx <- !is.na(full_table[["start"]])
  full_table <- full_table[keep_idx, ]
  keep_idx <- !is.na(full_table[["stop"]])
  full_table <- full_table[keep_idx, ]

  heat_cfg_file <- cfg[["cfg_file"]]
  heat_cfg_file <- gsub(pattern=".conf$", replacement="", x=heat_cfg_file)
  heat_cfg_file <- paste0(heat_cfg_file, colname, "_heatmap.conf")
  heat_data_file <- file.path(cfg[["data_dir"]], basename(heat_cfg_file))
  heat_data_file <- gsub(pattern=".conf$", replacement=".txt", x=heat_data_file)
  message("Writing data file: ", heat_data_file, " with the ", basename, colname, " column.")
  write.table(full_table, file=heat_data_file, quote=FALSE, row.names=FALSE, col.names=FALSE)

  num_colors <- 1
  if (is.null(colors)) {
    conditions <- levels(as.factor(full_table[[colname]]))
    num_colors <- length(conditions)
    colors <- sm(grDevices::colorRampPalette(
                              RColorBrewer::brewer.pal(num_colors, "Dark2"))(num_colors))
    names(colors) <- conditions
  } else {
    num_colors <- length(conditions)
  }

  ## Now write the config stanza
  minval <- -1
  maxval <- 1
  if (is.null(min_value)) {
    minval <- min(full_table[[colname]])
  } else if (is.numeric(min_value)) {
    minval <- min_value
  }
  if (is.null(max_value)) {
    maxval <- max(full_table[[colname]])
  } else if (is.numeric(max_value)) {
    maxval <- max_value
  }
  extra_rules <- ""
  if (!is.null(rules)) {
    extra_rules <- "<rules>"
    for (rule in rules) {
      extra_rules <- glue::glue("
{extra_rules}
   <rule>
    {rule}
   </rule>

")
    }
    extra_rules <- glue::glue("
{extra_rules}
  </rules>

")
  }
  inner <- outer - width
  rel_cfg_file <- file.path("conf", basename(heat_cfg_file))
  rel_data_file <- file.path("data", basename(heat_data_file))
  heat_cfg_string <- glue::glue("
## This is a circos heatmap configuration stanza.
 <plot>
  type = heatmap
  file = {rel_data_file}
  ## pattern = hline,vline  ## Also solid
  color_mapping = {color_mapping}
  min = {minval}
  max = {maxval}
  stroke_thickness = 0
  r1 = {outer}r
  r0 = {inner}r
  color = {color_choice}
  scale_log_base = {scale_log_base}
  url = script?type=label&value=[id]
  {extra_rules}
 </plot>

")
  heat_cfg_out <- file(heat_cfg_file, open="w+")
  cat(heat_cfg_string, file=heat_cfg_out, sep="")
  close(heat_cfg_out)

  ## Now add to the master configuration file.
  master_cfg_string <- glue::glue("
## The tile ring for {colname}
<<include {rel_cfg_file}>>

")
  master_cfg_out <- file(cfg[["cfg_file"]], open="a+")
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  new_outer <- inner - spacing
  message("Returning the inner width: ", new_outer,
          ".  Use it as the outer for the next ring.")
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
#' @param cfg Result of circos_prefix(), contains a bunch of useful material.
#' @param df Dataframe with starts/ends and the floating point information.
#' @param colname Name of the column with the data of interest.
#' @param basename Location to write the circos data (usually cwd).
#' @param color Color of the plotted data.
#' @param fill_color Guess
#' @param fill_under The circos histogram fill under parameter
#' @param extend_bin Extend bins?
                              #' @param thickness histogram thickness.
#' @param orientation facing in or out?
                                     #' @param outer Floating point radius of the circle into which to place the data.
#' @param width Radial width of each tile.
#' @param spacing Distance between outer, inner, and inner to whatever follows.
#' @return Radius after adding the histogram and the spacing.
#' @export
circos_hist <- function(cfg, df, colname="logFC", basename="", color="blue", fill_color="blue",
                        fill_under="yes", extend_bin="no", thickness="0", orientation="out",
                        outer=0.9, width=0.08, spacing=0.0) {
  ## I am going to have this take as input a data frame with genes as rownames
  ## starts, ends, and functional calls
  ## I will tell R to print out a suitable stanza for circos while I am at it
  ## because I am tired of mistyping something stupid.
  annot <- cfg[["annot"]]
  full_table <- merge(df, annot, by="row.names")
  if (nrow(full_table) == 0) {
    stop("Merging the annotations and data failed.")
  }
  start_colnames <- colnames(full_table)
  new_colnames <- gsub(x=start_colnames, pattern="\\.x$", replacement="")
  colnames(full_table) <- new_colnames
  rownames(full_table) <- full_table[["Row.names"]]
  full_table[["Row.names"]] <- NULL

  full_table <- full_table[, c("chr", "start", "stop", colname)]
  start_undefined_idx <- full_table[["start"]] == "undefined"
  full_table <- full_table[!start_undefined_idx, ]
  stop_undefined_idx <- full_table[["stop"]] == "undefined"
  full_table <- full_table[!stop_undefined_idx, ]
  full_table[["start"]] <- as.numeric(full_table[["start"]])
  full_table[["stop"]] <- as.numeric(full_table[["stop"]])
  keep_idx <- !is.na(full_table[["start"]])
  full_table <- full_table[keep_idx, ]
  keep_idx <- !is.na(full_table[["stop"]])
  full_table <- full_table[keep_idx, ]

  ## This can be simplified using cfg[["data_dir"]] and cfg[["conf_dir"]]
  hist_cfg_file <- cfg[["cfg_file"]]
  hist_cfg_file <- gsub(pattern=".conf$", replacement="", x=hist_cfg_file)
  hist_cfg_file <- paste0(hist_cfg_file, "_", basename, colname, "_hist.conf")
  hist_data_file <- file.path(cfg[["data_dir"]], basename(hist_cfg_file))
  hist_data_file <- gsub(pattern=".conf$", replacement=".txt", x=hist_data_file)
  message("Writing data file: ", hist_data_file, " with the ", basename, colname, " column.")
  write.table(full_table, file=hist_data_file, quote=FALSE, row.names=FALSE, col.names=FALSE)

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
  rel_cfg_file <- file.path("conf", basename(hist_cfg_file))
  rel_data_file <- file.path("data", basename(hist_data_file))
  hist_cfg_string <- glue::glue("
 <plot>
  type = histogram
  file = {rel_data_file}
  extend_bin = {extend_bin}
  fill_under = {fill_under}
  fill_color = {fill_color}
  color = {color}
  thickness = {thickness}
  r1 = {outer}r
  r0 = {inner}r
  orientation = {orientation}
 </plot>

")
  hist_cfg_out <- file(hist_cfg_file, open="w+")
  cat(hist_cfg_string, file=hist_cfg_out, sep="")
  close(hist_cfg_out)

  ## Now add to the master configuration file.
  master_cfg_string <- glue::glue("
  ## The histogram ring for {colname}
  <<include {rel_cfg_file}>>

")
  master_cfg_out <- file(cfg[["cfg_file"]], open="a+")
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  new_outer <- inner - spacing
  message("Returning the inner width: ", new_outer,
          ".  Use it as the outer for the next ring.")
  return(new_outer)
}

#' Create the description of chromosome markings.
#'
#' This function writes ideogram files for circos.
#'
#' @param name Name of the configuration file to which to add the ideogram.
#' @param conf_dir Where does the configuration live?
#' @param band_url Provide a url for making these imagemaps?
#' @param fill Fill in the strokes?
#' @param stroke_color What color?
#' @param show_bands Show the bands for the ideogram?
#' @param fill_bands and fill them in?
#' @param thickness How thick to color the lines
#' @param stroke_thickness How much of them to fill in
#' @param label_font What font to use.
#' @param spacing_default How much space between elements.
#' @param spacing_break Space between breaks.
#' @param fill_color What color to fill
#' @param radius Where on the circle to put them
#' @param radius_padding How much to pad between radii.
#' @param label_size How large to make the labels in px.
#' @param band_stroke_thickness How big to make the strokes!
                                                      #' @return The file to which the ideogram configuration was written.
#' @export
circos_ideogram <- function(name="default", conf_dir="circos/conf", band_url=NULL,
                            fill="yes", stroke_color="black",
                            show_bands="yes", fill_bands="yes",
                            thickness="20", stroke_thickness="2",
                            label_font="condensedbold",
                            spacing_default="0", spacing_break="0",
                            fill_color="black", radius="0.85", radius_padding="0.05",
                            label_size="36", band_stroke_thickness="2") {
  ideogram_outfile <- glue::glue("{conf_dir}/ideograms/{name}.conf")
  out <- file(ideogram_outfile, open="w+")
  show_label <- "no"
  ideogram_string <- glue::glue("
## The following plot stanza describes the ideograms
<ideogram>
  <spacing>
    default = {spacing_default}u
    break = {spacing_break}u
  </spacing>
  thickness = {thickness}p
  stroke_thickness = {stroke_thickness}
  stroke_color = {stroke_color}
  fill = {fill}
  fill_color = {fill_color}
  radius = {radius}r
  show_label = {show_label}
  label_font = {label_font}
  label_radius = dims(ideogram,radius) + {radius_padding}r
  label_size = {label_size}
  band_stroke_thickness = {band_stroke_thickness}
  show_bands = {show_bands}
  fill_bands = {fill_bands}

")
  cat(ideogram_string, file=out, sep="")
  ideogram_band <- ""
  if (!is.null(band_url)) {
    ideogram_band <- glue::glue("
 band_url = {band_url}
 ## image_map_missing_parameter=removeparam
 ideogram_url = {band_url}

")
    cat(ideogram_band, file=out, sep="")
  }
  end_string <- "</ideogram>\n"
  cat(end_string, file=out, sep="")
  close(out)
  message("Wrote karyotype to ", ideogram_outfile)
  message("This should match the karyotype= line in ", name, ".conf")
  return(ideogram_outfile)
}

#' Create the description of (a)chromosome(s) for circos.
#'
#' This function tries to save me from having to get the lengths of arcs for
#' bacterial chromosomes manually correct, and writes them as a circos
#' compatible karyotype file. The outfile parameter was chosen to match the
#' configuration directive outlined in circos_prefix(), however that will need
#' to be changed in order for this to work in variable conditions.  Next time I
#' make one of these graphs I will do that I suspect. In addition, this
#' currently only understands how to write bacterial chromosomes, that will
#' likely be fixed when I am asked to write out a L.major karyotype. These
#' defaults were chosen because I have a chromosome of this length that is
#' correct.
#'
#' @param cfg Result from circos_prefix(), contains a bunch of useful things.
#' @param segments How many segments to cut the chromosome into?
                                                             #' @param color Color segments of the chromosomal arc?
                                                                                                               #' @param fasta Fasta file to use to create the karyotype.
#' @param lengths If no sequence file is provided, use a named numeric vector to provide them.
#' @return The output filename.
#' @export
circos_karyotype <- function(cfg, segments=6, color="white", fasta=NULL, lengths=NULL) {
  conf_dir <- cfg[["conf_dir"]]
  name <- cfg[["name"]]
  ## genome_length <- 0
  chr_df <- data.frame()
  if (is.null(length) & is.null(fasta)) {
    stop("circos_karyotype() requires chromosome length or fasta file.")
  } else if (!is.null(lengths)) {
    ## genome_length <- sum(lengths)
    chr_df <- data.frame("width" = lengths, "names" = names(lengths))
  } else {
    raw_seq <- Rsamtools::FaFile(fasta)
    all_seq <- Biostrings::getSeq(raw_seq)
    ## genome_length <- sum(as.data.frame(all_seq@ranges)[["width"]])
    chr_df <- data.frame("width" = BiocGenerics::width(all_seq), "names" = names(all_seq))
    chr_df[["names"]] <- gsub(x=chr_df[["names"]], pattern="^(\\w+) .*", replacement="\\1")
  }

  chr_num <- nrow(chr_df)
  outfile <- glue::glue("{conf_dir}/karyotypes/{name}.conf")
  out <- file(outfile, open="w+")
  ## First write the summary line
  for (ch in 1:chr_num) {
    chr_name <- chr_df[ch, "names"]
    chr_width <- chr_df[ch, "width"]
    start_string <- glue::glue("chr - {chr_name} {chr_num} 0 {chr_width} {color}")
    cat(start_string, file=out, sep="\n")

    individual_segments <- segments
    if (chr_width < 100000) {
      individual_segments <- 1
    }
    for (segment in 1:individual_segments) {
      current <- segment - 1
      begin <- floor(current * (chr_width / individual_segments))
      end <- floor(segment * (chr_width / individual_segments))
      string <- glue::glue("band {chr_name} {chr_num}.1 {chr_num}.1 {begin} {end} {color}")
      cat(string, file=out, sep="\n")
    }
  }
  close(out)
  message("Wrote karyotype to ", outfile)
  message("This should match the karyotype= line in ", name, ".conf")
  return(outfile)
}

#' Write a simple makefile for circos.
#'
#' I regenerate all my circos pictures with make(1).  This is my makefile.
#'
#' @param cfg Configuration from circos_prefix().
#' @param target Default make target.
#' @param circos Location of circos.  I have a copy in home/bin/circos and use that sometimes.
#' @return a kitten
#' @export
circos_make <- function(cfg, target="", circos="circos") {
  circos_dir <- cfg[["basedir"]]
  output <- file.path(circos_dir, "Makefile")
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

%%.png:\t%%.conf
\t$(CIRCOS) -conf $< -outputfile $*.png

clean:
\trm -rf conf data *.conf *.png *.svg *.html

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
  close(out)

  make_target <- gsub(pattern="circos/conf/", replacement="", x=target)
  make_target <- gsub(pattern="\\.conf", replacement="", x=make_target)
  make_target_svg <- glue::glue("{make_target}.svg")
  make_target_png <- glue::glue("{make_target}.png")

  make_command <- glue::glue("cd circos && touch Makefile && make {make_target} 2>&1 | grep -v Redundant")
  result <- system(make_command) ##, show.output.on.console=FALSE)

  return(result)
}

#' Write tiles of bacterial ontology groups using the categories from
#' microbesonline.org.
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
#' @param cfg Result from circos_prefix().
#' @param outer Floating point radius of the circle into which to place the
#'   plus-strand data.
#' @param width Radial width of each tile.
#' @param thickness How wide to make the bars.
#' @param spacing Radial distance between outer, inner, and inner to whatever follows.
#' @param padding How much space between them.
#' @param margin Margin between elements.
#' @param plus_orientation Orientation of the plus pieces.
#' @param minus_orientation Orientation of the minus pieces.
#' @param layers How many layers to use
#' @param layers_overflow How to handle too many layers.
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
circos_plus_minus <- function(cfg, outer=1.0, width=0.08, thickness=95,
                              spacing=0.0, padding=1, margin=0.00,
                              plus_orientation="out", minus_orientation="in",
                              layers=1, layers_overflow="hide",
                              acol="orange", bcol="reds-9-seq", ccol="yellow",
                              dcol="vlpurple", ecol="vlgreen", fcol="dpblue",
                              gcol="vlgreen", hcol="vlpblue", icol="vvdpgreen",
                              jcol="dpred", kcol="orange", lcol="vvlorange",
                              mcol="dpgreen", ncol="vvlpblue", ocol="vvlgreen",
                              pcol="vvdpred", qcol="ylgn-3-seq", rcol="vlgrey",
                              scol="grey", tcol="vlpurple", ucol="greens-3-seq",
                              vcol="vlred", wcol="vvdppurple", xcol="black",
                              ycol="lred", zcol="vlpblue") {

  message("Writing data file: ", cfg[["plus_data_file"]], " with the + strand GO data.")
  write.table(cfg[["plus_df"]], file=cfg[["plus_data_file"]], quote=FALSE,
              row.names=FALSE, col.names=FALSE, na="no_go")
  message("Writing data file: ", cfg[["minus_data_file"]], " with the - strand GO data.")
  write.table(cfg[["minus_df"]], file=cfg[["minus_data_file"]], quote=FALSE,
              row.names=FALSE, col.names=FALSE, na="no_go")

  first_outer <- outer
  first_inner <- first_outer - width
  plus_cfg_out <- file(cfg[["plus_cfg_file"]], open="w+")
  rel_plus_file <- file.path("data", basename(cfg[["plus_data_file"]]))
  rules_string <- glue::glue("
  <rules>
   <rule>
    condition = var(value) =~ \"^A\"
    fill_color = {acol}
    color = {acol}
   </rule>
   <rule>
    condition = var(value) =~ \"^B\"
    fill_color = {bcol}
    color = {bcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^C\"
    fill_color = {ccol}
    color = {ccol}
   </rule>
   <rule>
    condition = var(value) =~ \"^D\"
    fill_color = {dcol}
    color = {dcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^E\"
    fill_color = {ecol}
    color = {ecol}
   </rule>
   <rule>
    condition = var(value) =~ \"^F\"
    fill_color = {fcol}
    color = {fcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^G\"
    fill_color = {gcol}
    color = {gcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^H\"
    fill_color = {hcol}
    color = {hcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^I\"
    fill_color = {icol}
    color = {icol}
   </rule>
   <rule>
    condition = var(value) =~ \"^J\"
    fill_color = {jcol}
    color = {jcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^K\"
    fill_color = {kcol}
    color = {kcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^L\"
    fill_color = {lcol}
    color = {lcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^M\"
    fill_color = {mcol}
    color = {mcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^N\"
    fill_color = {ncol}
    color = {ncol}
   </rule>
   <rule>
    condition = var(value) =~ \"^O\"
    fill_color = {ocol}
    color = {ocol}
   </rule>
   <rule>
    condition = var(value) =~ \"^P\"
    fill_color = {pcol}
    color = {pcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^Q\"
    fill_color = {qcol}
    color = {qcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^R\"
    fill_color = {rcol}
    color = {rcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^S\"
    fill_color = {scol}
    color = {scol}
   </rule>
   <rule>
    condition = var(value) =~ \"^T\"
    fill_color = {tcol}
    color = {tcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^U\"
    fill_color = {ucol}
    color = {ucol}
   </rule>
   <rule>
    condition = var(value) =~ \"^V\"
    fill_color = {vcol}
    color = {vcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^W\"
    fill_color = {wcol}
    color = {wcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^X\"
    fill_color = {xcol}
    color = {xcol}
   </rule>
   <rule>
    condition = var(value) =~ \"^Y\"
    fill_color = {ycol}
    color = {ycol}
   </rule>
   <rule>
    condition = var(value) =~ \"^Z\"
    fill_color = {zcol}
    color = {zcol}
   </rule>
  </rules>
")

  plus_cfg_string <- glue::glue("

## The following plot stanza describes the broad ontological categories
## Of genes on the plus strand.
 <plot>
  type = tile
  file = {rel_plus_file}
  layers = {layers}
  layers_overflow = {layers_overflow}
  margin = {margin}u
  thickness = {thickness}
  padding = {padding}
  orientation = {plus_orientation}
  stroke_thickness = 0
  stroke_color = black
  color = green
  r1 = {first_outer}r
  r0 = {first_inner}r
  url = script?type=label&value=[id]&color=[color]
{rules_string}
 </plot>

")
  cat(plus_cfg_string, file=plus_cfg_out, sep="")
  close(plus_cfg_out)

  ## Now move the ring in one width and print the minus strand.
  second_outer <- first_inner - spacing
  second_inner <- second_outer - width
  minus_cfg_out <- file(cfg[["minus_cfg_file"]], open="w+")
  rel_minus_file <- file.path("data", basename(cfg[["minus_data_file"]]))
  minus_cfg_string <- glue::glue("

## The following plot stanza describes the broad ontological categories
## Of genes on the minus strand.
 <plot>
  type = tile
  file = {rel_minus_file}
  layers = {layers}
  layers_overflow = {layers_overflow}
  margin = {margin}u
  thickness = {thickness}
  padding = {padding}
  orientation = {minus_orientation}
  stroke_thickness = 0
  stroke_color = black
  color = green
  r1 = {second_outer}r
  r0 = {second_inner}r
  url = script?type=label&value=[id]&color=[color]
{rules_string}
 </plot>

")
  cat(minus_cfg_string, file=minus_cfg_out, sep="")
  close(minus_cfg_out)
  message("Wrote the +/- config files.  Appending their inclusion to the master file.")

  rel_plus_cfg <- file.path("conf", basename(cfg[["plus_cfg_file"]]))
  rel_minus_cfg <- file.path("conf", basename(cfg[["minus_cfg_file"]]))

  master_cfg_out <- file(cfg[["cfg_file"]], open="a+")
  master_cfg_string <- glue::glue("

  ## The +/- ontology rings.
  ## The plus strand ring extends from {first_outer}r to {first_inner}r
  <<include {rel_plus_cfg}>>
  ## The minus strand ring extends from {second_outer}r to {second_inner}r
  <<include {rel_minus_cfg}>>

")
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  message("Returning the inner width: ", second_inner,
          ".  Use it as the outer for the next ring.")
  new_outer <- second_inner - spacing
  return(new_outer)
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
#' @param annotation Annotation data frame.
#' @param name Name of the map, called with 'make name'.
#' @param basedir Base directory for writing the data.
#' @param chr_column Name of the column containing the chromosome names in the annotations.
#' @param cog_column Name of the column containing the COG groups in the annotations.
#' @param start_column Name of the column containing the starts in the annotations.
#' @param stop_column Name of the column containing the stops in the annotations.
#' @param strand_column Name of the column containing the strand information.
#' @param id_column Where do the gene IDs live? NULL means rownames.
#' @param cog_map Not yet used, but used to provide an alternate map of groups/colors.
#' @param radius Size of the image.
#' @param chr_units How often to print chromosome in 'prefix' units.
#' @param band_url Place to imagemap link.
#' @param ... Extra arguments passed to the tick/karyotype makers.
#' @return The master configuration file name.
#' @export
circos_prefix <- function(annotation, name="mgas", basedir="circos",
                          chr_column="seqnames", cog_column="COGFun",
                          start_column="start", stop_column="end",
                          strand_column="strand", id_column=NULL,
                          cog_map=NULL,
                          radius=1800, chr_units=1000, band_url=NULL, ...) {
  message("This assumes you have a colors.conf in circos/colors/ ",
          "and fonts.conf in circos/fonts/")
  message("It also assumes you have conf/ideogram.conf, conf/ticks.conf, ",
          "and conf/housekeeping.conf")
  conf_dir <- file.path(basedir, "conf")
  data_dir <- file.path(basedir, "data")
  cfgout <- paste0(file.path(conf_dir, name), ".conf")
  message("It will write ", cfgout, " with a reasonable first approximation config file.")

  if (!file.exists(data_dir)) {
    message("Creating the data directory: ", data_dir)
    dir.create(data_dir, recursive=TRUE)
  }
  if (!file.exists(conf_dir)) {
    message("The circos directory does not exist, creating: ", conf_dir)
    dir.create(conf_dir, recursive=TRUE)
  }

  ## Set up some data which will be shared by all the other functions.
  number_pluses <- sum(annotation[[strand_column]] == "+")
  number_ones <- sum(annotation[[strand_column]] == 1)
  plus_string <- "+"
  minus_string <- "-"
  if (number_pluses + number_ones == 0) {
    stop("This function requires some way of understanding plus/minus strand.")
  } else if (number_ones > 0) {
    plus_string <- 1
    minus_string <- -1
  }
  if (is.null(annotation[[start_column]]) | is.null(annotation[[stop_column]]) |
      is.null(annotation[[strand_column]]) | is.null(annotation[[cog_column]])) {
    stop("This function assumes columns for start, stop, strand, chromosome names, and cog.")
  }
  plus_cfg_file <- cfgout
  minus_cfg_file <- cfgout
  plus_cfg_file <- gsub(pattern=".conf$", replacement="_plus_go.conf", x=plus_cfg_file)
  plus_data_file <- file.path(data_dir, basename(plus_cfg_file))
  plus_data_file <- gsub(pattern=".conf$", replacement=".txt", x=plus_data_file)
  minus_cfg_file <- gsub(pattern=".conf$", replacement="_minus_go.conf", x=minus_cfg_file)
  minus_data_file <- file.path(data_dir, basename(minus_cfg_file))
  minus_data_file <- gsub(pattern=".conf$", replacement=".txt", x=minus_data_file)
  ## What I should do is spend some time thinking and reformat this to handle
  ## and arbitrary number of arbitrary columns so that I have some flexibility later.
  gids <- NULL
  if (is.null(id_column)) {
    gids <- rownames(annotation)
  } else {
    gids <- annotation[[id_column]]
  }

  plus_idx <- annotation[[strand_column]] == plus_string
  plus_gids <- gids[plus_idx]
  minus_idx <- annotation[[strand_column]] == minus_string
  minus_gids <- gids[minus_idx]

  needed_columns <- c(chr_column, start_column, stop_column, cog_column)
  plus_df <- as.data.frame(annotation[plus_idx, needed_columns])
  rownames(plus_df) <- plus_gids
  colnames(plus_df) <- c("chr", "start", "stop", "cog")
  minus_df <- as.data.frame(annotation[minus_idx, needed_columns])
  rownames(minus_df) <- minus_gids
  colnames(minus_df) <- c("chr", "start", "stop", "cog")
  plus_df[["value"]] <- glue::glue("value={plus_df[['cog']]}0")
  minus_df[["value"]] <- glue::glue("value={minus_df[['cog']]}0")

  needed_columns <- c(chr_column, start_column, stop_column, strand_column)
  annot <- annotation[, needed_columns]
  rownames(annot) <- gids
  colnames(annot) <- c("chr", "start", "stop", "strand")

  karyotype_dir <- file.path(conf_dir, "karyotypes")
  if (!file.exists(karyotype_dir)) {
    message("The karyotype directory does not exist, creating: ", karyotype_dir)
    dir.create(karyotype_dir, recursive=TRUE)
  }
  ideogram_dir <- file.path(conf_dir, "ideograms")
  if (!file.exists(ideogram_dir)) {
    message("The ideogram directory does not exist, creating: ", ideogram_dir)
    dir.create(ideogram_dir, recursive=TRUE)
  }

  etc_file <- file.path(path.package('hpgltools'), "circos", "circos_etc.tar.xz")
  etc_cmd <- glue::glue("tar -C {dirname(conf_dir)} -xavf {etc_file} 2>/dev/null 1>&2")
  system(command=etc_cmd)

  karyotype_cfg_file <- paste0(file.path(basedir, "karyotypes", name), ".conf")
  rel_karyotype_file <- paste0(file.path("conf", "karyotypes", name), ".conf")
  ## If you want clickable ideograms, add
  ## band_url='script?start=[start]&end=[end]&label=[label]
  ideogram_cfg_file <- circos_ideogram(name=name, conf_dir=conf_dir, band_url=band_url)
  rel_ideogram_file <- paste0(file.path("conf", "ideograms", name), ".conf")
  tick_cfg_file <- circos_ticks(name=name, conf_dir=conf_dir,
                                ...)
  rel_tick_file <- file.path("conf", basename(tick_cfg_file))

  out <- file(cfgout, open="w+")
  prefix_string <- glue::glue("
## This is the prefix of a circos configuration file written by hpgltools.
<colors>
 <<include colors.conf>>
</colors>
<fonts>
 <<include fonts.conf>>
</fonts>

<<include {rel_ideogram_file}>>
<<include {rel_tick_file}>>
<<include housekeeping.conf>>
karyotype = {rel_karyotype_file}

<image>
 image_map_use = yes
 image_map_missing_parameter = removeurl
 dir = .
 radius = {radius}p
 background = white
 angle_offset = -90
</image>

chromosomes_units = {chr_units}
chromosomes_display_default = yes

<highlights>
</highlights>

<plots>

")
  cat(prefix_string, file=out, sep="")
  close(out)
  to_path <- glue::glue("{name}.conf")
  wd <- getwd()
  final_cfg <- file.path(wd, basedir, to_path)
  if (!file.exists(final_cfg)) {
    tmpwd <- glue::glue("{wd}/circos")
    setwd(file.path(wd, basedir))
    from <- gsub(pattern="circos/", replacement="", x=cfgout)
    file.symlink(from, to_path)
    setwd(wd)
  }

  retlist <- list(
    "name" = name,
    "basedir" = basedir,
    "cfg_file" = cfgout,
    "conf_dir" = conf_dir,
    "data_dir" = data_dir,
    "annotation" = annotation,
    "annot" = annot,
    "karyotype_cfg_file" = karyotype_cfg_file,
    "ideogram_cfg_file" = ideogram_cfg_file,
    "tick_cfg_file" = tick_cfg_file,
    "plus_df" = plus_df,
    "plus_cfg_file" = plus_cfg_file,
    "plus_data_file" = plus_data_file,
    "minus_df" = minus_df,
    "minus_cfg_file" = minus_cfg_file,
    "minus_data_file" = minus_data_file)
  return(retlist)
}

#' Write the end of a circos master configuration.
#'
#' circos configuration files need an ending.  This writes it.
#'
#' @param cfg Result from circos_prefix()
#' @return The filename of the configuration.
#' @export
circos_suffix <- function(cfg) {
  cfgout <- cfg[["cfg_file"]]
  out <- file(cfgout, open="a+")
  suffix_string <- "
</plots>"
  cat(suffix_string, file=out, sep="\n")
  close(out)
}

#' Create the ticks for a circos plot.
#'
#' This function writes ticks for circos.  This has lots of options, the
#' defaults are all taken from the circos example documentation for a bacterial
#' genome.
#'
#' @param name Name of the configuration file to which to add the ideogram.
#' @param conf_dir Where does the configuration live.
#' @param show_ticks Show them or not.
#' @param show_tick_labels Show the tick labels, or do not.
#' @param show_grid Print a grid behind.
#' @param skip_first_label Like a clock.
#' @param skip_last_label Ditto.
#' @param tick_separation Top-level separation between tick marks.
#' @param min_label_distance distance to the edge of the plot for labels.
#' @param label_separation radial distance between labels.
#' @param label_offset The offset for the labels.
#' @param label_size Top-level label size.
#' @param multiplier When writing the position, by what factor to lower the numbers?
                                                                              #' @param main_color Color for top-level labels?
                                                                                                                       #' @param main_thickness Top-level thickness of lines etc.
#' @param main_size Top-level size of text.
#' @param first_size Second level size of text.
#' @param first_spacing Second level spacing of ticks.
#' @param first_color Second-level text color.
#' @param first_show_label Show a label for the second level ticks?
                                                               #' @param first_label_size Text size for second level labels?
                                                                                                                     #' @param second_size Size of ticks for the third level.
#' @param second_spacing third-level spacing
#' @param second_color Text color for the third level.
#' @param second_show_label Give them a label?
                                          #' @param second_label_size And a size.
#' @param third_size Now for the size of the almost-largest ticks
#' @param third_spacing How far apart?
                                  #' @param third_color and their color
#' @param third_show_label give a label?
                                    #' @param third_label_size and a size.
#' @param fourth_spacing The largest ticks!
                                       #' @param fourth_color The largest color.
#' @param fourth_show_label Provide a label?
                                        #' @param suffix String for printing chromosome distances.
#' @param fourth_label_size They are big!
                                       #' @param include_first_label Provide the smallest labels?
                                                                                          #' @param include_second_label Second smallest labels?
                                                                                                                                         #' @param include_third_label Second biggest labels?
                                                                                                                                                                                      #' @param include_fourth_label Largest labels?
                                                                                                                                                                                                                             #' @param ... Extra arguments from circos_prefix().
#' @return The file to which the ideogram configuration was written.
#' @export
circos_ticks <- function(name="default", conf_dir="circos/conf",
                         show_ticks="yes", show_tick_labels="yes",
                         show_grid="no", skip_first_label="yes",
                         skip_last_label="no",
                         tick_separation=2, min_label_distance=0,
                         label_separation=5, label_offset=5,
                         label_size=8, multiplier=0.001, main_color="black",
                         main_thickness=3, main_size=20, first_size=10,
                         first_spacing=1, first_color="black",
                         first_show_label="no", first_label_size=12, second_size=15,
                         second_spacing=5, second_color="black",
                         second_show_label="yes", second_label_size=16,
                         third_size=18, third_spacing=10, third_color="black",
                         third_show_label="yes", third_label_size=16, fourth_spacing=100,
                         fourth_color="black", fourth_show_label="yes",
                         suffix=" kb", fourth_label_size=36,
                         include_first_label=TRUE, include_second_label=TRUE,
                         include_third_label=TRUE, include_fourth_label=TRUE,
                         ...) {

  tick_outfile <- file.path(conf_dir, paste0("ticks_", name, ".conf"))
  out <- file(tick_outfile, open="w")
  show_label <- "no"
  tick_string <- glue::glue("
## The following plot stanza describes the ticks
show_ticks = {show_ticks}
show_tick_labels = {show_tick_labels}
show_grid = {show_grid}
grid_start = dims(ideogram,radius_inner) - 0.5r
grid_end = dims(ideogram,radius_inner)
<ticks>
  skip_first_label = {skip_first_label}
  skip_last_label = {skip_last_label}
  radius = dims(ideogram,radius_outer)
  tick_separation = {tick_separation}p
  min_label_distance_to_edge = {min_label_distance}p
  label_separation = {label_separation}p
  label_offset = {label_offset}p
  label_size = {label_size}p
  multiplier = {multiplier}
  color = {main_color}
  thickness = {main_thickness}p
  size = {main_size}p

")
  if (isTRUE(include_first_label)) {
    tick_string <- glue::glue("
{tick_string}
  <tick>
    size = {first_size}p
    spacing = {first_spacing}u
    color = {first_color}
    show_label = {first_show_label}
    label_size = {first_label_size}p
    format = %.2f
    grid = no
    grid_color = lblue
    grid_thickness = 1p
  </tick>

")
  }
  if (isTRUE(include_second_label)) {
    tick_string <- glue::glue("
{tick_string}
  <tick>
    size = {second_size}p
    spacing = {second_spacing}u
    color = {second_color}
    show_label = {second_show_label}
    label_size = {second_label_size}p
    format = %s
    grid = yes
    grid_color = lgrey
    grid_thickness = 1p
  </tick>

")
  }
  if (isTRUE(include_third_label)) {
    tick_string <- glue::glue("
{tick_string}
  <tick>
    size = {third_size}p
    spacing = {third_spacing}u
    color = {third_color}
    show_label = {third_show_label}
    label_size = {third_label_size}p
    format = %s
    grid = yes
    grid_color = grey
    grid_thickness = 1p
  </tick>

")
  }
  if (isTRUE(include_fourth_label)) {
    tick_string <- glue::glue("
{tick_string}
  <tick>
    spacing = {fourth_spacing}u
    color = {fourth_color}
    show_label = {fourth_show_label}
    suffix = \"{suffix}\"
    label_size = {fourth_label_size}p
    format = %s
    grid = yes
    grid_color = dgrey
    grid_thickness = 1p
  </tick>

")
  }
  tick_string <- glue::glue("
{tick_string}
</ticks>

")
  cat(tick_string, file=out, sep="")
  close(out)
  message("Wrote ticks to ", tick_outfile)
  return(tick_outfile)
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
#' @param cfg Result from circos_prefix().
#' @param df Dataframe with starts/ends and the floating point information.
#' @param colname Name of the column with the data of interest.
#'   chromosome)
#' @param basename  Used to make unique filenames for the data/conf files.
#' @param colors Colors of the data.
#' @param thickness How thick to make the tiles in radial units.
#' @param margin How much space between other rings and the tiles?
                                                              #' @param stroke_thickness Size of the tile outlines.
#' @param padding Space between tiles.
#' @param orientation Facing in or out.
#' @param outer Floating point radius of the circle into which to place the
#'   categorical data.
#' @param width Width of each tile.
#' @param spacing Radial distance between outer, inner, and inner to whatever
#'   follows.
#' @return Radius after adding the histogram and the spacing.
#' @export
circos_tile <- function(cfg, df, colname="logFC", basename="", colors=NULL,
                        thickness=80, padding=1, margin=0.00, stroke_thickness=0.00,
                        orientation="out",
                        outer=0.9, width=0.08, spacing=0.0) {
  annot <- cfg[["annot"]]
  full_table <- merge(df, annot, by="row.names")
  if (nrow(full_table) == 0) {
    stop("Merging the annotations and data failed.")
  }
  start_colnames <- colnames(full_table)
  new_colnames <- gsub(x=start_colnames, pattern="\\.x$", replacement="")
  colnames(full_table) <- new_colnames
  rownames(full_table) <- full_table[["Row.names"]]
  full_table[["Row.names"]] <- NULL

  full_table <- full_table[, c("chr", "start", "stop", colname)]
  start_undefined_idx <- full_table[["start"]] == "undefined"
  full_table <- full_table[!start_undefined_idx, ]
  stop_undefined_idx <- full_table[["stop"]] == "undefined"
  full_table <- full_table[!stop_undefined_idx, ]
  full_table[["start"]] <- as.numeric(full_table[["start"]])
  full_table[["stop"]] <- as.numeric(full_table[["stop"]])
  keep_idx <- !is.na(full_table[["start"]])
  full_table <- full_table[keep_idx, ]
  keep_idx <- !is.na(full_table[["stop"]])
  full_table <- full_table[keep_idx, ]

  tile_cfg_file <- cfg[["cfg_file"]]
  tile_cfg_file <- gsub(pattern=".conf$", replacement="", x=tile_cfg_file)
  tile_cfg_file <- paste0(tile_cfg_file, colname, "_tile.conf")
  tile_data_file <- file.path(cfg[["data_dir"]], basename(tile_cfg_file))
  tile_data_file <- gsub(pattern=".conf$", replacement=".txt", x=tile_data_file)
  message("Writing data file: ", tile_data_file, " with the ", basename, colname, " column.")
  write.table(full_table, file=tile_data_file, quote=FALSE, row.names=FALSE, col.names=FALSE)

  num_colors <- 1
  if (is.null(colors)) {
    conditions <- levels(as.factor(full_table[[colname]]))
    num_colors <- length(conditions)
    colors <- sm(grDevices::colorRampPalette(
                              RColorBrewer::brewer.pal(num_colors, "Dark2"))(num_colors))
    names(colors) <- conditions
  } else {
    num_colors <- length(colors)
  }

  ## Now write the config stanza
  inner <- outer - width
  rel_cfg_file <- file.path("conf", basename(tile_cfg_file))
  rel_data_file <- file.path("data", basename(tile_data_file))
  tile_cfg_string <- glue::glue("
 ## This is a circos tile configuration stanza.
 <plot>
  type = tile
  file = {rel_data_file}
  layers = 1
  layers_overflow = hide
  margin = {margin}u
  thickness = {thickness}
  padding = {padding}
  orientation = {orientation}
  stroke_thickness = {stroke_thickness}
  stroke_color = black
  color = black
  fill_color = black
  r1 = {outer}r
  r0 = {inner}r
  <rules>

")
  tile_cfg_out <- file(tile_cfg_file, open="w+")
  cat(tile_cfg_string, file=tile_cfg_out, sep="")
  for (c in 1:num_colors) {
    red_component <- "0x00"
    green_component <- "0x00"
    blue_compnent <- "0x00"
    this_color <- gsub(pattern="^#", replacement="", x=colors[[c]])
    red_component <- strtoi(glue::glue("0x{substr(this_color, 1, 2)}"))
    green_component <- strtoi(glue::glue("0x{substr(this_color, 3, 4)}"))
    blue_component <- strtoi(glue::glue("0x{substr(colors[[c]], 5, 6)}"))
    color_string <- glue::glue("{red_component},{green_component},{blue_component}")
    color_name <- names(colors)[[c]]
    new_string <- glue::glue("
   <rule>
    condition = var(value) =~ \"^{color_name}\"
    fill_color = {color_string}
    color = {color_string}
   </rule>

")
    cat(new_string, file=tile_cfg_out, sep="")
  }
  end_string <- glue::glue("
  </rules>
 </plot>

")
  cat(end_string, file=tile_cfg_out, sep="")
  close(tile_cfg_out)

  ## Now add to the master configuration file.
  master_cfg_string <- glue::glue("
  ## The tile ring for {colname}
  <<include {rel_cfg_file}>>

")
  master_cfg_out <- file(cfg[["cfg_file"]], open="a+")
  cat(master_cfg_string, file=master_cfg_out, sep="")
  close(master_cfg_out)

  new_outer <- inner - spacing
  message("Returning the inner width: ", new_outer,
          ".  Use it as the outer for the next ring.")
  return(new_outer)
}

## EOF
