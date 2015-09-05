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
