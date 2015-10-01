## The karyotype file is circos/data/5005_5448_karyotype.txt
## The 5005 genome is 1838562 nt. long (looking at reference/genbank/mgas_5005.gb)
## The 5448 genome is 1829516 nt. long
## (+ 1838562 1829516) 3668078
## I have been writing these text files manually, maybe I should write a function to do it...
## Except there are so many too many ways of expressing them
## So I will write a function to write individual chromosomes

## These defaults were chosen because I have a chromosome of this length that is correct.
write_karyotype = function(outfile='circos/conf/karyotypes/default.txt', length=1838554, name='chr1', segments=6, color='white', chr_num=1) {
    out = file(outfile, open='w+')
    ## First write the summary line
    start_string = sprintf("chr - %s %d 0 %d %s", name, chr_num, length, color)
    cat(start_string, file=out, sep="\n")
    for (segment in 1:segments) {
        current = segment - 1
        begin = floor(current * (length / segments))
        end = floor(segment * (length / segments))
        string = sprintf("band %s %d.1 %d.1 %d %d %s", name, chr_num, chr_num, begin, end, color)
        cat(string, file=out, sep="\n")
    }
    close(out)
}

write_plus_minus = function(go_table, cfgout="circos/conf/default.conf", chr='chr1', outer=1.0, width=0.08) {
    plus_cfg_file = cfgout
    minus_cfg_file = cfgout
    plus_cfg_file = gsub(".conf$", "_plus_go.conf", plus_cfg_file)
    minus_cfg_file = gsub(".conf$", "_minus_go.conf", minus_cfg_file)
    go_plus = as.data.frame(go_table[go_table$strand == "+",])
    go_minus = as.data.frame(go_table[go_table$strand == "-",])
    go_plus$chr = chr
    go_minus$chr = chr
    go_plus = go_plus[,c(5,1,2,4)]
    go_minus = go_minus[,c(5,1,2,4)]
    go_plus$COGFun = paste0("value=", go_plus$COGFun, "0")
    go_minus$COGFun = paste0("value=", go_minus$COGFun, "0")
    data_prefix = cfgout
    data_prefix = gsub("/conf/", "/data/", data_prefix)
    data_prefix = gsub(".conf$", "", data_prefix)

    plus_file = paste0(data_prefix, "_plus_go.txt")        
    message(paste0("Writing data file: ", plus_file, " with the + strand GO data."))
    write.table(go_plus, file=plus_file, quote=FALSE, row.names=FALSE, col.names=FALSE, na="no_go")

    minus_file = paste0(data_prefix, "_minus_go.txt")
    message(paste0("Writing data file: ", minus_file, " with the - strand GO data."))
    write.table(go_minus, file=minus_file, quote=FALSE, row.names=FALSE, col.names=FALSE, na="no_go")

    inner = outer - width
    plus_cfg_out = file(plus_cfg_file, open='w+')
    plus_cfg_filename = gsub("^circos/", "", plus_file)
    plus_cfg_string = sprintf("
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
", plus_cfg_filename, outer, inner)
    cat(plus_cfg_string, file=plus_cfg_out, sep="")
    close(plus_cfg_out)

    ## Now move the ring in one width and print the minus strand.
    outer = inner
    inner = outer - width
    minus_cfg_out = file(minus_cfg_file, open='w+')
    minus_cfg_filename = gsub("^circos/", "", minus_file)    
    minus_cfg_string = sprintf("
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
", minus_cfg_filename, outer, inner)
    cat(minus_cfg_string, file=minus_cfg_out, sep="")
    close(minus_cfg_out)
    message("Wrote the +/- config files.  Appending their inclusion to the master file.")
    master_cfg_out = file(cfgout, open="a+")
    plus_cfg_include = plus_cfg_filename
    plus_cfg_include = gsub("^circos/", "", plus_cfg_file)
    minus_cfg_include = minus_cfg_filename
    minus_cfg_include = gsub("^circos/", "", minus_cfg_file)
    master_cfg_string = sprintf("

## The +/- ontology rings.
<<include %s>>
<<include %s>>

", plus_cfg_include, minus_cfg_include)
    cat(master_cfg_string, file=master_cfg_out, sep="")
    close(master_cfg_out)

    message(paste0("Returning the inner width: ", inner, ".  Use it as the outer for the next ring."))
    return(inner)
}

circos_calls = function(df, outfile="calls.txt", name='chr1', stanza_file="circos_stanza.txt", colors=NULL, outer=0.9, width=0.08) {
    ## I am going to have this take as input a data frame with genes as rownames
    ## starts, ends, and functional calls
    ## I will tell R to print out a suitable stanza for circos while I am at it
    ## because I am tired of mistyping something stupid.
    if (is.null(df$start) | is.null(df$end) | is.null(rownames(df)) | is.null(df$call)) {
        stop("This requires columns: start, end, rownames, and call")
    }
    df$chr = name
    df = df[,c("chr","start","end","call")]
    write.table(df, file=outfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
    message(paste0("Writing ", stanza_file, " with a suitable circos configuration stanza for this data."))
    ## If colors is not provided, use colorbrewer to make up colors for every condition found in the factor of
    ## df$call.
    ## If colors is provided, it needs to have 1 element for every found condition.
    num_colors = 1
    if (is.null(colors)) {
        conditions = levels(as.factor(df$call))
        num_colors = length(conditions)
        colors = suppressWarnings(colorRampPalette(brewer.pal(num_colors, "Dark2"))(num_colors))
        names(colors) = conditions
    }
    out = file(stanza_file, open='w+')
    inner = outer - width
    start_string = sprintf("
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
", outfile, outer, inner)
    cat(start_string, file=out, sep="")
    for (c in 1:num_colors) {
        red_component = paste0("0x", substr(colors[[c]], 2, 3))
        green_component = paste0("0x", substr(colors[[c]], 4, 5))
        blue_component = paste0("0x", substr(colors[[c]], 5, 6))
        red_component = strtoi(red_component)
        green_component = strtoi(green_component)
        blue_component = strtoi(blue_component)
        color_string = paste0(red_component, ",", blue_component, ",", green_component)
        new_string = sprintf("   <rule>
    condition = var(value) =~ \"^%s\"
    fill_color = %s
    color = %s
   </rule>
", names(colors)[[c]], color_string, color_string)
        cat(new_string, file=out, sep="")
    }
    end_string = sprintf("  </rules>\n </plot>", sep="")
    cat(end_string, file=out, sep="")
    close(out)
    return(inner)
}

circos_makefile = function(output="circos/Makefile", circos="/usr/bin/circos") {
    circos_dir = dirname(output)
    if (!file.exists(circos_dir)) {
        message(paste0("The circos directory does not exist, creating: ", circos_dir))
        dir.create(circos_dir, recursive=TRUE)
    }
    out = file(output, open='w+')    
    makefile_string = sprintf("
.PHONY:\tclean
CIRCOS=\"%s\"
%%.png:\t%%.conf
\t$(CIRCOS) -conf $< -outputfile $*.png

%%.svg:\t%%.conf
\t$(CIRCOS) -conf $< -outputfile $*.svg
", circos)
    cat(makefile_string, file=out, sep="")
    close(out)
}

circos_prefix = function(cfgout="circos/conf/default.conf", radius=1800) {
    message("This assumes you have a colors.conf in circos/colors/ and fonts.conf in circos/fonts/")
    message("It also assumes you have conf/ideogram.conf, conf/ticks.conf, and conf/housekeeping.conf")
    message(paste0("It will write ", cfgout, " with a reasonable first approximation config file."))
    circos_dir = dirname(cfgout)
    if (!file.exists(circos_dir)) {
        message(paste0("The circos directory does not exist, creating: ", circos_dir))
        dir.create(circos_dir, recursive=TRUE)
    }
    karyotype_dir = paste0(circos_dir, "/karyotypes")
    if (!file.exists(karyotype_dir)) {
        message(paste0("The karyotype directory does not exist, creating: ", karyotype_dir))
        dir.create(karyotype_dir, recursive=TRUE)
    }
    karyotype_file = cfgout
    karyotype_file = gsub("circos/conf", "conf/karyotypes", karyotype_file)
    karyotype_file = gsub("\\.conf$", "\\.txt", karyotype_file)
    
    out = file(cfgout, open='w+')        
    prefix_string = sprintf("## This is the prefix of a circos configuration file written by hpgltools.
<colors>
 <<include colors/colors.conf>>
</colors>
<fonts>
 <<include fonts/fonts.conf>>
</fonts>

<<include conf/ideogram.conf>>
<<include conf/ticks.conf>>
<<include conf/housekeeping.conf>>

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
", karyotype_file, radius)
    cat(prefix_string, file=out, sep="")
    close(out)
    wd = getwd()
    tmpwd = paste0(wd, "/circos")
    setwd(tmpwd)
    from = gsub("circos/", "", cfgout)
    to = "."
    file.symlink(from, to)
    setwd(wd)
}


circos_suffix = function(cfgout="circos/conf/default.conf") {
    out = file(cfgout, open='a+')        
    suffix_string = "</plots>"
    cat(suffix_string, file=out, sep="\n")
    close(out)
}

