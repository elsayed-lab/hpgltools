## The karyotype file is circos/data/5005_5448_karyotype.txt
## The 5005 genome is 1838562 nt. long (looking at reference/genbank/mgas_5005.gb)
## The 5448 genome is 1829516 nt. long
## (+ 1838562 1829516) 3668078
## I have been writing these text files manually, maybe I should write a function to do it...
## Except there are so many too many ways of expressing them
## So I will write a function to write individual chromosomes

## These defaults were chosen because I have a chromosome of this length that is correct.
circos_chr = function(outfile='chr1_karyotype.txt', length=1838554, name='chr1', segments=6, color='white', chr_num=1) {

    out = file(outfile, open='a+')
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


write_plus_minus = function(go_table, prefix, chr='chr1') {
    go_plus = as.data.frame(go_table[go_table$strand == "+",])
    go_minus = as.data.frame(go_table[go_table$strand == "-",])
    go_plus$chr = chr
    go_minus$chr = chr
    go_plus = go_plus[,c(5,1,2,4)]
    go_minus = go_minus[,c(5,1,2,4)]
    go_plus$COGFun = paste0("value=", go_plus$COGFun, "0")
    go_minus$COGFun = paste0("value=", go_minus$COGFun, "0")
    plus_file = paste0(prefix, "_plus.txt")
    minus_file = paste0(prefix, "_minus.txt")
    write.table(go_plus, file=plus_file, quote=FALSE, row.names=FALSE, col.names=FALSE, na="no_go")
    write.table(go_minus, file=minus_file, quote=FALSE, row.names=FALSE, col.names=FALSE, na="no_go")
}

circos_calls = function(df, outfile="calls.txt", name='chr1', stanza_file="circos_stanza.txt", colors=NULL, outer_radius=0.9, inner_radius=0.81) {
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
    start_string = sprintf(" <plot>\n  type = tile\n  file = %s\n  layers = 1\n  layers_overflow = hide\n  margin = 0.00u\n  thickness = 90\n  padding = 1\n  orientation = out\n  stroke_thickness = 0\n  stroke_color = black\n  color = black\n  fill_color = black\n  r1 = %fr\n  r0 = %fr\n  <rules>\n", outfile, outer_radius, inner_radius)
    cat(start_string, file=out, sep="")
    for (c in 1:num_colors) {
        red_component = paste0("0x", substr(colors[[c]], 2, 3))
        green_component = paste0("0x", substr(colors[[c]], 4, 5))
        blue_component = paste0("0x", substr(colors[[c]], 5, 6))
        red_component = strtoi(red_component)
        green_component = strtoi(green_component)
        blue_component = strtoi(blue_component)
        color_string = paste0(red_component, ",", blue_component, ",", green_component)
        new_string = sprintf("   <rule>\n    condition = var(value) =~ \"^%s\"\n    fill_color = %s\n    color = %s\n   </rule>\n", names(colors)[[c]], color_string, color_string)
        cat(new_string, file=out, sep="")
    }
    end_string = sprintf("  </rules>\n </plot>", sep="")
    cat(end_string, file=out, sep="")
    close(out)
}
