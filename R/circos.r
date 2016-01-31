## Time-stamp: <Fri Jan 29 12:35:31 2016 Ashton Trey Belew (abelew@gmail.com)>

## The karyotype file is circos/data/5005_5448_karyotype.txt
## The 5005 genome is 1838562 nt. long (looking at reference/genbank/mgas_5005.gb)
## The 5448 genome is 1829516 nt. long
## (+ 1838562 1829516) 3668078
## I have been writing these text files manually, maybe I should write a function to do it...
## Except there are so many too many ways of expressing them
## So I will write a function to write individual chromosomes

#' circos_karyotype()  Create the description of (a)chromosome(s) for circos.
#'
#' This function tries to save me from having to get the lengths of arcs for bacterial chromosomes
#' manually correct, and writes them as a circos compatible karyotype file.
#' The outfile parameter was chosen to match the configuration directive
#' outlined in circos_prefix(), however that will need to be changed in order for this to work
#' in variable conditions.  Next time I make one of these graphs I will do that I suspect.
#' In addition, this currently only understands how to write bacterial chromosomes, that will
#' likely be fixed when I am asked to write out a L.major karyotype.
#'
#' These defaults were chosen because I have a chromosome of this length that is correct.
#' @param outfile default='circos/conf/karyotypes/default.txt'  a file to write
#' @param length default=1838554  the default length of the chromosome (That is mgas5005)
#' @param name default='chr1'  the name of the chromosome (This currently assumes a bacterial chromosome)
#' @param segments default=6  how many segments to cut it into
#' @param color default='white'  how to colors the chromosomal arc. (circos images are cluttered enough)
#' @param chr_num default=1  the number to record (This and name above should change for multi-chromosomal species)
#'
#' @export
#' @return undef
circos_karyotype <- function(name='default', conf_dir='circos/conf', length=NULL,
                            chr_name='chr1', segments=6, color='white',
                            chr_num=1, fasta=NULL) {
    genome_length <- 0
    if (is.null(length) & is.null(fasta)) {
        stop("circos_karyotype() requires a chromosome length or fasta file to gather sequence data from.")
    } else if (!is.null(length)) {
        genome_length <- length
    } else {
        raw_seq <- Rsamtools::FaFile(fasta)
        all_seq <- Biostrings::getSeq(raw_seq)
        genome_length <- sum(as.data.frame(all_seq@ranges)$width)
    }

    outfile <- paste0(conf_dir, '/karyotypes/', name, ".conf")
    out <- file(outfile, open='w+')
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
    message(paste0("Wrote karyotype to ", outfile))
    message(paste0("This should match the karyotype= line in ", name, ".conf"))
    return(outfile)
}

#' circos_ideogram()  Create the description of chromosome markings
#'
#' This function writes ideogram files for circos.
#' Currently it only has a single format.
#'
#' @param name default='default'  the name of the configuration
#' @param conf_dir default='circos/conf'  where does the configuration live?
#' @param band_url default=NULL  provide a url for making these imagemaps.
#'
#' @export
#' @return undef
circos_ideogram <- function(name='default', conf_dir='circos/conf', band_url=NULL) {
    ideogram_outfile <- paste0(conf_dir, '/ideograms/', name, ".conf")
    out <- file(ideogram_outfile, open='w+')
    show_label <- "no"
    ideogram_string <- sprintf("## The following plot stanza describes the ideograms
<ideogram>
 <spacing>
 default = 0u
 break = 0u
 </spacing>
 # thickness (px) of chromosome ideogram
 thickness = 20p
 stroke_thickness = 2
 # ideogram border color
 stroke_color = black
 fill = yes
 # the default chromosome color is set here and any value
 # defined in the karyotype file overrides it
 fill_color = black
 # fractional radius position of chromosome ideogram within image
 radius = 0.85r
 show_label = %s
 label_font = condensedbold
 label_radius = dims(ideogram,radius) + 0.05r
 label_size = 36
 # cytogenetic bands
 band_stroke_thickness = 2
 # show_bands determines whether the outline of cytogenetic bands
 # will be seen
 show_bands = yes
 # in order to fill the bands with the color defined in the karyotype
 # file you must set fill_bands
 fill_bands = yes
", show_label)
    cat(ideogram_string, file=out, sep="")
    ideogram_band = ""
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
    message(paste0("Wrote karyotype to ", out))
    message(paste0("This should match the karyotype= line in ", name, ".conf"))
    return(out)
}

#' circos_plus_minus()  Write tiles of bacterial ontology groups using the categories from microbesonline.org
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
#' @param go_table  a dataframe with starts/ends and categories
#' @param cfgout default='circos/conf/default.conf'  The master
#'     configuration file to write.
#' @param chr default='chr1'  the name of the chromosome (This currently assumes a bacterial chromosome)
#' @param outer default=1.0  the floating point radius of the circle
#'     into which to place the plus-strand data
#' @param width default=0.08  the radial width of each tile
#' @param spacing default=0.0  the radial distance between outer,inner
#'     and inner,whatever follows.
#'
#' @export
#' @return the radius after adding the plus/minus information and the spacing between them.
circos_plus_minus <- function(go_table, cfgout="circos/conf/default.conf", chr='chr1',
                             outer=1.0, width=0.08, spacing=0.0) {
    plus_cfg_file <- cfgout
    minus_cfg_file <- cfgout
    plus_cfg_file <- gsub(".conf$", "_plus_go.conf", plus_cfg_file)
    minus_cfg_file <- gsub(".conf$", "_minus_go.conf", minus_cfg_file)
    go_plus <- as.data.frame(go_table[go_table$strand == "+",])
    go_minus <- as.data.frame(go_table[go_table$strand == "-",])
    go_plus$chr <- chr
    go_minus$chr <- chr
    ##    go_plus = go_plus[,c(5,1,2,4)]
    message("This function assumes an input go table including the columns: 'start', 'end', 'strand', and 'go'")
    go_plus <- go_plus[,c("chr","start","end","go")]
    ##go_minus = go_minus[,c(5,1,2,4)]
    go_minus <- go_minus[,c("chr","start","end","go")]
    go_plus$go <- paste0("value=", go_plus$go, "0")
    go_minus$go <- paste0("value=", go_minus$go, "0")
    data_prefix <- cfgout
    data_prefix <- gsub("/conf/", "/data/", data_prefix)
    data_prefix <- gsub(".conf$", "", data_prefix)

    plus_file <- paste0(data_prefix, "_plus_go.txt")
    message(paste0("Writing data file: ", plus_file, " with the + strand GO data."))
    write.table(go_plus, file=plus_file, quote=FALSE, row.names=FALSE, col.names=FALSE, na="no_go")

    minus_file <- paste0(data_prefix, "_minus_go.txt")
    message(paste0("Writing data file: ", minus_file, " with the - strand GO data."))
    write.table(go_minus, file=minus_file, quote=FALSE, row.names=FALSE, col.names=FALSE, na="no_go")

    first_outer <- outer
    first_inner <- first_outer - width
    plus_cfg_out <- file(plus_cfg_file, open='w+')
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
  <rules>
   <rule>
    condition = var(value) =~ \"^A\"
    fill_color = orange
    color = orange
   </rule>
   <rule>
    condition = var(value) =~ \"^B\"
    fill_color = reds-9-seq
    color = reds-9-seq
   </rule>
   <rule>
    condition = var(value) =~ \"^C\"
    fill_color = yellow
    color = yellow
   </rule>
   <rule>
    condition = var(value) =~ \"^D\"
    fill_color = vlpurple
    color = vlpurple
   </rule>
   <rule>
    condition = var(value) =~ \"^E\"
    fill_color = vlgreen
    color = vlgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^F\"
    fill_color = dpblue
    color = dpblue
   </rule>
   <rule>
    condition = var(value) =~ \"^G\"
    fill_color = vlgreen
    color = vlgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^H\"
    fill_color = vlpblue
    color = vlpblue
   </rule>
   <rule>
    condition = var(value) =~ \"^I\"
    fill_color = vvdpgreen
    color = vvdpgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^J\"
    fill_color = dpred
    color = dpred
   </rule>
   <rule>
    condition = var(value) =~ \"^K\"
    fill_color = orange
    color = orange
   </rule>
   <rule>
    condition = var(value) =~ \"^L\"
    fill_color = vvlorange
    color = vvlorange
   </rule>
   <rule>
    condition = var(value) =~ \"^M\"
    fill_color = dpgreen
    color = dpgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^N\"
    fill_color = vvlpblue
    color = vvlpblue
   </rule>
   <rule>
    condition = var(value) =~ \"^O\"
    fill_color = vvlgreen
    color = vvlgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^P\"
    fill_color = vvdpred
    color = vvdpred
   </rule>
   <rule>
    condition = var(value) =~ \"^Q\"
    fill_color = ylgn-3-seq
    color = ylgn-3-seq
   </rule>
   <rule>
    condition = var(value) =~ \"^R\"
    fill_color = vlgrey
    color = vlgrey
   </rule>
   <rule>
    condition = var(value) =~ \"^S\"
    fill_color = grey
    color = grey
   </rule>
   <rule>
    condition = var(value) =~ \"^T\"
    fill_color = vlpurple
    color = vlpurple
   </rule>
   <rule>
    condition = var(value) =~ \"^U\"
    fill_color = greens-3-seq
    color = greens-3-seq
   </rule>
   <rule>
    condition = var(value) =~ \"^V\"
    fill_color = vlred
    color = vlred
   </rule>
   <rule>
    condition = var(value) =~ \"^W\"
    fill_color = vvdppurple
    color = vvdppurple
   </rule>
   <rule>
    condition = var(value) =~ \"^X\"
    fill_color = black
    color = black
   </rule>
   <rule>
    condition = var(value) =~ \"^Y\"
    fill_color = lred
    color = lred
   </rule>
   <rule>
    condition = var(value) =~ \"^Z\"
    fill_color = vlpblue
    color = vlpblue
   </rule>
  </rules>
 </plot>
", plus_cfg_filename, first_outer, first_inner)
    cat(plus_cfg_string, file=plus_cfg_out, sep="")
    close(plus_cfg_out)

    ## Now move the ring in one width and print the minus strand.
    second_outer <- first_inner - spacing
    second_inner <- second_outer - width
    minus_cfg_out <- file(minus_cfg_file, open='w+')
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
  <rules>
   <rule>
    condition = var(value) =~ \"^A\"
    fill_color = orange
    color = orange
   </rule>
   <rule>
    condition = var(value) =~ \"^B\"
    fill_color = reds-9-seq
    color = reds-9-seq
   </rule>
   <rule>
    condition = var(value) =~ \"^C\"
    fill_color = yellow
    color = yellow
   </rule>
   <rule>
    condition = var(value) =~ \"^D\"
    fill_color = vlpurple
    color = vlpurple
   </rule>
   <rule>
    condition = var(value) =~ \"^E\"
    fill_color = vlgreen
    color = vlgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^F\"
    fill_color = dpblue
    color = dpblue
   </rule>
   <rule>
    condition = var(value) =~ \"^G\"
    fill_color = vlgreen
    color = vlgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^H\"
    fill_color = vlpblue
    color = vlpblue
   </rule>
   <rule>
    condition = var(value) =~ \"^I\"
    fill_color = vvdpgreen
    color = vvdpgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^J\"
    fill_color = dpred
    color = dpred
   </rule>
   <rule>
    condition = var(value) =~ \"^K\"
    fill_color = orange
    color = orange
   </rule>
   <rule>
    condition = var(value) =~ \"^L\"
    fill_color = vvlorange
    color = vvlorange
   </rule>
   <rule>
    condition = var(value) =~ \"^M\"
    fill_color = dpgreen
    color = dpgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^N\"
    fill_color = vvlpblue
    color = vvlpblue
   </rule>
   <rule>
    condition = var(value) =~ \"^O\"
    fill_color = vvlgreen
    color = vvlgreen
   </rule>
   <rule>
    condition = var(value) =~ \"^P\"
    fill_color = vvdpred
    color = vvdpred
   </rule>
   <rule>
    condition = var(value) =~ \"^Q\"
    fill_color = ylgn-3-seq
    color = ylgn-3-seq
   </rule>
   <rule>
    condition = var(value) =~ \"^R\"
    fill_color = vlgrey
    color = vlgrey
   </rule>
   <rule>
    condition = var(value) =~ \"^S\"
    fill_color = grey
    color = grey
   </rule>
   <rule>
    condition = var(value) =~ \"^T\"
    fill_color = vlpurple
    color = vlpurple
   </rule>
   <rule>
    condition = var(value) =~ \"^U\"
    fill_color = greens-3-seq
    color = greens-3-seq
   </rule>
   <rule>
    condition = var(value) =~ \"^V\"
    fill_color = vlred
    color = vlred
   </rule>
   <rule>
    condition = var(value) =~ \"^W\"
    fill_color = vvdppurple
    color = vvdppurple
   </rule>
   <rule>
    condition = var(value) =~ \"^X\"
    fill_color = black
    color = black
   </rule>
   <rule>
    condition = var(value) =~ \"^Y\"
    fill_color = lred
    color = lred
   </rule>
   <rule>
    condition = var(value) =~ \"^Z\"
    fill_color = vlpblue
    color = vlpblue
   </rule>
  </rules>
 </plot>
", minus_cfg_filename, second_outer, second_inner)
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

    message(paste0("Returning the inner width: ", second_inner, ".  Use it as the outer for the next ring."))
    new_outer = second_inner - spacing
    return(new_outer)
}

#' circos_tile()  Write tiles of arbitrary categorical point data in
#' circos.
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
#' @param df  a dataframe with starts/ends and the floating point information
#' @param cfgout default='circos/conf/default.conf'  The master
#'     configuration file to write.
#' @param colname default='datum'  The name of the column with the
#'     data of interest.
#' @param chr default='chr1'  the name of the chromosome (This
#'     currently assumes a bacterial chromosome)
#' @param color default='blue'  the color of the histogram
#' @param fill_color default='blue'  guess
#' @param outer default=1.0  the floating point radius of the circle
#'     into which to place the plus-strand data
#' @param width default=0.08  the radial width of each tile
#' @param spacing default=0.0  the radial distance between outer,inner
#'     and inner,whatever follows.
#'
#' @export
#' @return the radius after adding the histogram and the spacing.
circos_tile <- function(df, cfgout="circos/conf/default.conf", colname="datum",
                       chr='chr1', colors=NULL, outer=0.9, width=0.08, spacing=0.0) {
    ## I am going to have this take as input a data frame with genes as rownames
    ## starts, ends, and functional calls
    ## I will tell R to print out a suitable stanza for circos while I am at it
    ## because I am tired of mistyping something stupid.
    if (is.null(df$start) | is.null(df$end) | is.null(rownames(df)) | is.null(df[[colname]])) {
        stop("This requires columns: start, end, rownames, and datum")
    }
    datum_cfg_file <- cfgout
    datum_cfg_file <- gsub(".conf$", "", datum_cfg_file)
    datum_cfg_file <- paste0(datum_cfg_file, "_", colname, "_tile.conf")
    df$chr <- chr
    df <- df[,c("chr","start","end", colname)]
    data_prefix <- cfgout
    data_prefix <- gsub("/conf/", "/data/", data_prefix)
    data_prefix <- gsub(".conf$", "", data_prefix)
    data_filename <- paste0(data_prefix, "_", colname, "_tile.txt")
    message(paste0("Writing data file: ", data_filename, " with the ", colname, " column."))
    write.table(df, file=data_filename, quote=FALSE, row.names=FALSE, col.names=FALSE)

    num_colors <- 1
    if (is.null(colors)) {
        conditions <- levels(as.factor(df[[colname]]))
        num_colors <- length(conditions)
        colors <- suppressWarnings(grDevices::colorRampPalette(RColorBrewer::brewer.pal(num_colors, "Dark2"))(num_colors))
        names(colors) <- conditions
    }

    ## Now write the config stanza
    inner <- outer - width
    data_cfg_out <- file(datum_cfg_file, open="w+")
    data_cfg_filename <- gsub("^circos/", "", datum_cfg_file)
    data_file_string <- gsub("^circos/", "", data_filename)

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
", data_file_string, outer, inner)
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
    new_outer <- inner - spacing
    master_cfg_string <- sprintf("## The tile ring for column:%s, extending from %sr to %sr
<<include %s>>
", colname, outer, new_outer, data_cfg_include)
    cat(master_cfg_string, file=master_cfg_out, sep="")
    close(master_cfg_out)

    return(new_outer)
}

#' circos_heatmap()  Write tiles of arbitrary heat-mappable data in
#' circos.
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
#' @param df  a dataframe with starts/ends and the floating point information
#' @param cfgout default='circos/conf/default.conf'  The master
#'     configuration file to write.
#' @param colname default='datum'  The name of the column with the
#'     data of interest.
#' @param chr default='chr1'  the name of the chromosome (This
#'     currently assumes a bacterial chromosome)
#' @param color default='blue'  the color of the histogram
#' @param fill_color default='blue'  guess
#' @param outer default=1.0  the floating point radius of the circle
#'     into which to place the plus-strand data
#' @param width default=0.08  the radial width of each tile
#' @param spacing default=0.0  the radial distance between outer,inner
#'     and inner,whatever follows.
#'
#' @export
#' @return the radius after adding the histogram and the spacing.
circos_heatmap <- function(df, cfgout="circos/conf/default.conf", colname="datum",
                          chr='chr1', colors=NULL, outer=0.9, width=0.08, spacing=0.0) {
    ## I am going to have this take as input a data frame with genes as rownames
    ## starts, ends, and functional calls
    ## I will tell R to print out a suitable stanza for circos while I am at it
    ## because I am tired of mistyping something stupid.
    if (is.null(df$start) | is.null(df$end) | is.null(rownames(df)) | is.null(df[[colname]])) {
        stop("This requires columns: start, end, rownames, and datum")
    }
    datum_cfg_file <- cfgout
    datum_cfg_file <- gsub(".conf$", "", datum_cfg_file)
    datum_cfg_file <- paste0(datum_cfg_file, "_", colname, "_heatmap.conf")
    df$chr <- chr
    df <- df[,c("chr","start","end", colname)]
    data_prefix <- cfgout
    data_prefix <- gsub("/conf/", "/data/", data_prefix)
    data_prefix <- gsub(".conf$", "", data_prefix)
    data_filename <- paste0(data_prefix, "_", colname, "_heatmap.txt")
    message(paste0("Writing data file: ", data_filename, " with the ", colname, " column."))
    write.table(df, file=data_filename, quote=FALSE, row.names=FALSE, col.names=FALSE)

    num_colors <- 1
    if (is.null(colors)) {
        conditions <- levels(as.factor(df[[colname]]))
        num_colors <- length(conditions)
        colors <- suppressWarnings(RColorBrewer::colorRampPalette(grDevices::brewer.pal(num_colors, "Dark2"))(num_colors))
        names(colors) <- conditions
    }

    ## Now write the config stanza
    inner <- outer - width
    data_cfg_out <- file(datum_cfg_file, open="w+")
    data_cfg_filename <- gsub("^circos/", "", datum_cfg_file)
    data_cfg_string <- sprintf("## This is a circos heatmap configuration stanza.
 <plot>
  type = heatmap
  r1 = %fr
  r0 = %fr
  color = conf(plots, color_alt)
  file = %s
  pattern = hline,vline  # Also solid
  color_mapping = 0 ## also 1, 2
  min = 0
  max = 10
  stroke_thickness = 0
##  layers = 1
##  layers_overflow = hide
##  margin = 0.00u
##  thickness = 90
##  padding = 1
##  orientation = out
##  stroke_thickness = 0
##  stroke_color = black
##  color = black
##  fill_color = black
 </plot>", sep="")
    cat(data_cfg_string, file=data_cfg_out, sep="")
    close(data_cfg_out)

    ## Now add to the master configuration file.
    master_cfg_out <- file(cfgout, open="a+")
    data_cfg_include <- data_cfg_filename
    data_cfg_include <- gsub("^circos/", "", data_cfg_include)
    new_outer <- inner - spacing
    master_cfg_string <- sprintf("## The tile ring for column:%s, extending from %sr to %sr
<<include %s>>
", colname, outer, new_outer, data_cfg_include)
    cat(master_cfg_string, file=master_cfg_out, sep="")
    close(master_cfg_out)
    return(new_outer)
}

#' circos_hist()  Write histograms of arbitrary floating point data in
#' circos.
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
#' @param df  a dataframe with starts/ends and the floating point information
#' @param cfgout default='circos/conf/default.conf'  The master
#'     configuration file to write.
#' @param colname default='datum'  The name of the column with the
#'     data of interest.
#' @param chr default='chr1'  the name of the chromosome (This
#'     currently assumes a bacterial chromosome)
#' @param color default='blue'  the color of the histogram
#' @param fill_color default='blue'  guess
#' @param outer default=1.0  the floating point radius of the circle
#'     into which to place the plus-strand data
#' @param width default=0.08  the radial width of each tile
#' @param spacing default=0.0  the radial distance between outer,inner
#'     and inner,whatever follows.
#'
#' @export
#' @return the radius after adding the histogram and the spacing.
circos_hist <- function(df, cfgout="circos/conf/default.conf", colname="datum", chr='chr1',
                       color="blue", fill_color="blue", outer=0.9, width=0.08, spacing=0.0) {
    ## I am going to have this take as input a data frame with genes as rownames
    ## starts, ends, and functional calls
    ## I will tell R to print out a suitable stanza for circos while I am at it
    ## because I am tired of mistyping something stupid.
    if (is.null(df$start) | is.null(df$end) | is.null(rownames(df)) | is.null(df[[colname]])) {
        stop("This requires columns: start, end, rownames, and datum")
    }
    datum_cfg_file <- cfgout
    datum_cfg_file <- gsub(".conf$", "", datum_cfg_file)
    datum_cfg_file <- paste0(datum_cfg_file, "_", colname, "_hist.conf")
    df$chr <- chr
    df <- df[,c("chr","start","end", colname)]
    data_prefix <- cfgout
    data_prefix <- gsub("/conf/", "/data/", data_prefix)
    data_prefix <- gsub(".conf$", "", data_prefix)
    data_filename <- paste0(data_prefix, "_", colname, "_hist.txt")
    message(paste0("Writing data file: ", data_filename, " with the ", colname, " column."))
    write.table(df, file=data_filename, quote=FALSE, row.names=FALSE, col.names=FALSE)

    num_colors <- 1
    if (is.null(colors)) {
        conditions <- levels(as.factor(df$call))
        num_colors <- length(conditions)
        colors <- suppressWarnings(grDevices::colorRampPalette(RColorBrewer::brewer.pal(num_colors, "Dark2"))(num_colors))
        names(colors) <- conditions
    }

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

    new_outer = inner - spacing
    return(new_outer)
}

#' circos_make()  Write a simple makefile for circos.
#'
#' I regenerate all my circos pictures with make(1).  This is my
#' makefile.
#'
#' @param output default='circos/Makefile' the makefile
#' @param circos default='/usr/bin/circos'  the location of circos. (I
#'     have a copy in home/bin/circos and use that sometimes.
#'
#' @export
#' @return a kitten
circos_make <- function(target="", output="circos/Makefile", circos="/usr/bin/circos") {
    circos_dir <- dirname(output)
    if (!file.exists(circos_dir)) {
        message(paste0("The circos directory does not exist, creating: ", circos_dir))
        dir.create(circos_dir, recursive=TRUE)
    }
    out <- file(output, open='w+')
    makefile_string <- sprintf("
.PHONY:\tclean
CIRCOS=\"%s\"

%%.png:\t%%.conf
\t$(CIRCOS) -conf $< -outputfile $*.png

%%.svg:\t%%.conf
\t$(CIRCOS) -conf $< -outputfile $*.svg

%%:\t%%.conf
\t$(CIRCOS) -conf $< -outputfile $*.png && \t$(CIRCOS) -conf $< -outputfile $*.svg
", circos)
    cat(makefile_string, file=out, sep="")
    close(out)

    make_target <- gsub(pattern="circos/conf/", replacement="", x=target)
    make_target_svg <- gsub(pattern="\\.conf", replacement="\\.svg", x=make_target)
    make_target_png <- gsub(pattern="\\.conf", replacement="\\.png", x=make_target)
    make_command <- paste0("cd circos && make ", make_target_svg, " && make ", make_target_png)
    message(paste0("Running: ", make_command))
    system(make_command)
}

#' circos_arc()  Write arcs between chromosomes in circos.
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
#' @param df  a dataframe with starts/ends and the floating point information
#' @param cfgout default='circos/conf/default.conf'  The master
#'     configuration file to write.
#' @param first_col default='chr1'  The name of the first chromosome
#' @param second_col default='chr2'  The name of the second chromosome
#' @param color default='blue'  the color of the histogram
#' @param radius default=0.75  the radius at which to add the arcs
#' @param thickness default=3  integer thickness of the arcs
#'
#' @export
#' @return undef
circos_arc <- function(df, cfgout="circos/conf/default.conf", first_col='chr1', second_col='chr2',
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
    message(paste0("This function assumes an input table including columns: ",first_start_name,",",first_end_name,",", second_start_name,",",second_end_name,",",first_name,",",second_name,"."))
    df <- df[,c(first_name, second_name, first_start_name, first_end_name, second_start_name, second_end_name)]
    data_prefix <- cfgout
    data_prefix <- gsub("/conf/", "/data/", data_prefix)
    data_prefix <- gsub(".conf$", "", data_prefix)
    data_filename <- paste0(data_prefix, "_", first_col, "_arc.txt")
    message(paste0("Writing data file: ", data_filename, " with the ", first_col, " column."))
    print_arc <- function(x) {
        cat(x[5], " chr5005 ", x[1], " ", x[2], "\n", x[5], " chr5448 ", x[3], " ", x[4], "\n\n", file="circos/data/crossref_5005_5448.txt", append=TRUE, sep="")
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

#' circos_prefix()  Write the beginning of a circos configuration file.
#'
#' A few parameters need to be set when starting circos.  This sets
#' some of them and gets ready for plot stanzas.
#'
#' In its current implementation, this really assumes that there will
#' be no highlight stanzas and at most 1 link stanza.
#' chromosomes.  A minimal amount of logic and data organization will
#' address these weaknesses.
#'
#' @param name default='default'  The name of the map, called with 'make name'
#' @param circos_dir default='circos/conf'  The directory containing the circos configuration data.
#' @param radius default=1800  The size of the image.
#'
#' @export
#' @return undef
circos_prefix <- function(name="default", conf_dir="circos/conf", radius=1800, band_url=NULL) {
    message("This assumes you have a colors.conf in circos/colors/ and fonts.conf in circos/fonts/")
    message("It also assumes you have conf/ideogram.conf, conf/ticks.conf, and conf/housekeeping.conf")
    cfgout <- paste0(conf_dir, "/", name, ".conf")
    message(paste0("It will write ", cfgout, " with a reasonable first approximation config file."))

    data_dir <- gsub(pattern="conf", replacement="data", x=conf_dir)
    if (!file.exists(data_dir)) {
        message(paste0("Creating the data directory: ", data_dir))
        dir.create(data_dir, recursive=TRUE)
    }
    if (!file.exists(conf_dir)) {
        message(paste0("The circos directory does not exist, creating: ", conf_dir))
        dir.create(conf_dir, recursive=TRUE)
    }
    karyotype_dir <- paste0(conf_dir, "/karyotypes")
    ideogram_dir <- paste0(conf_dir, "/ideograms")
    if (!file.exists(karyotype_dir)) {
        message(paste0("The karyotype directory does not exist, creating: ", karyotype_dir))
        dir.create(karyotype_dir, recursive=TRUE)
    }
    if (!file.exists(paste0(conf_dir, '/ideograms'))) {
        message(paste0("The ideogram directory does not exist, creating: ", ideogram_dir))
        dir.create(ideogram_dir, recursive=TRUE)
    }
    karyotype_file <- gsub("circos/conf", "conf/karyotypes", cfgout)
    ideogram_file <- gsub("circos/conf", "conf/ideograms", cfgout)
    etc_file <- paste0(path.package("hpgltools"), "/inst/circos/circos_etc.tar.xz")
    etc_cmd <- paste0("tar -C ", dirname(conf_dir), " -xavf ", etc_file)
    system(command=etc_cmd)

    ## If you want clickable ideograms, add band_url='script?start=[start]&end=[end]&label=[label]
    ## for example, thus one could bring up
    circos_ideogram(name=name, conf_dir=conf_dir, band_url=band_url)

    out <- file(cfgout, open='w+')
    prefix_string <- sprintf("## This is the prefix of a circos configuration file written by hpgltools.
<colors>
 <<include etc/colors.conf>>
</colors>
<fonts>
 <<include etc/fonts/fonts.conf>>
</fonts>

## The ideograms are generated by circos_ideogram() in R.
<<include %s>>

<<include etc/ticks.conf>>
<<include etc/housekeeping.conf>>

karyotype = %s

<image>
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
", ideogram_file, karyotype_file, radius)
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

#' circos_suffix()  Write the end of a circos master configuration.
#'
#' circos configuration files need an ending.  This writes it.
#'
#' @param cfgout default='circos/conf/default.conf'  The master
#'     configuration file to write.
#'
#' @export
#' @return undef
circos_suffix <- function(cfgout="circos/conf/default.conf") {
    out <- file(cfgout, open='a+')
    suffix_string <- "</plots>"
    cat(suffix_string, file=out, sep="\n")
    close(out)
}

## EOF
