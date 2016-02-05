
#' run the rGADEM suite
#'
#' This is another function I started but never had cause to finish
#' for the test sequences it works though
#'
#' @export
simple_gadem <- function() {
    ##    require.auto('rGADEM')
    Hsapiens <- NULL
    require.auto('BSgenome.Hsapiens.UCSC.hg19')
##    require.auto("Biostrings")
    ## The following is for using bed files.
    pwd <- "" #INPUT FILES- BedFiles, FASTA, etc.
    path <- system.file("extdata/Test_100.bed",package="rGADEM")
    BedFile <- paste(pwd,path,sep="")
    BED <- read.table(BedFile,header=FALSE,sep="\t")
    BED <- data.frame("chr"=as.factor(BED[,1]),
                      "start"=as.numeric(BED[,2]),
                      "end"=as.numeric(BED[,3]))
    rgBED <- IRanges::IRanges(start=BED[,2], end=BED[,3])
    Sequences <- IRanges::RangedData(rgBED, space=BED[,1])
    ## The following is for using fasta files
    pwd <- "" #INPUT FILES- BedFiles, FASTA, etc.
    path <- system.file("extdata/Test_100.fasta",package="rGADEM")
    FastaFile <- paste(pwd,path,sep="")
    Sequences <- Biostrings::readDNAStringSet(FastaFile, "fasta")
    ## The actual gadem analysis follows
    gadem <- rGADEM::GADEM(Sequences,verbose=1,genome=Hsapiens)
    ## Or use a seeded analysis to go faster...
    path <- system.file("extdata/jaspar2009.txt",package="rGADEM")
    seededPwm <- rGADEM::readPWMfile(path)
    grep("STAT1", names(seededPwm))
    STAT1.PWM <- seededPwm[103]
    ##    gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens,Spwm=STAT1.PWM, fixSeeded=TRUE)
    gadem <- rGADEM::GADEM(Sequences, verbose=1, genome=Hsapiens, Spwm=STAT1.PWM)
    ## gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens,Spwm=STAT1.PWM)
    ## viewing results
    rGADEM::nOccurrences(gadem)
    rGADEM::consensus(gadem)

    ## require.auto('MotIV')
    motifs <- rGADEM::getPWM(gadem)
    ## plot(foxa1.filter.combine ,ncol=2,top=5, rev=FALSE, main="Logo", bysim=TRUE)
    ## foxa1.alignment <- viewAlignments(foxa1.filter.combine )
    ## print(foxa1.alignment[[1]] )
    ## plot(foxa1.filter.combine ,gadem,ncol=2, type="distribution", correction=TRUE, group=FALSE, bysim=TRUE, strand=FALSE, sort=TRUE, main="Distribution of FOXA")
    ## plot(foxa1.filter.combine ,gadem,type="distance", correction=TRUE, group=TRUE, bysim=TRUE, main="Distance between FOXA and AP-1", strand=FALSE, xlim=c(-100,100), bw=8)
    return(motifs)
}

simple_motifRG <- function() {
    ## require.auto('motifRG')
    MD.motifs <- motifRG::findMotifFasta(system.file("extdata", "MD.peak.fa", package="motifRG"),
                                         system.file("extdata", "MD.control.fa", package="motifRG"),
                                         max.motif=3,enriched=T)
    motifRG::motifLatexTable(main="MyoD motifs", MD.motifs, prefix="myoD")
    YY1.peak <- NULL
    YY1.control <- NULL
    data(YY1.peak)
    data(YY1.control)
    Hsapiens <- NULL
    require.auto("BSgenome.Hsapiens.UCSC.hg19")
    YY1.peak.seq <- seqinr::getSequence(YY1.peak, genome=Hsapiens)
    YY1.control.seq <- seqinr::getSequence(YY1.control, genome=Hsapiens)
    YY1.motif.1 <- motifRG::findMotifFgBg(YY1.peak.seq, YY1.control.seq, enriched=T)
    motifRG::motifLatexTable(main="YY1 motifs", YY1.motif.1, prefix="YY1-1")
    summary(Biostrings::letterFrequency(YY1.peak.seq, "CG", as.prob=T))
    summary(Biostrings::letterFrequency(YY1.control.seq, "CG", as.prob=T))
    summary(BiocGenerics::width(YY1.peak.seq))
    YY1.narrow.seq <- XVector::subseq(YY1.peak.seq,
                                      pmax(round((BiocGenerics::width(YY1.peak.seq) - 200)/2), 1),
                                      width=pmin(200, BiocGenerics::width(YY1.peak.seq)))
    YY1.control.narrow.seq <- XVector::subseq(YY1.control.seq,
                                              pmax(round((BiocGenerics::width(YY1.control.seq) - 200)/2),1),
                                              width=pmin(200, BiocGenerics::width(YY1.control.seq)))
    category <- c(rep(1, length(YY1.narrow.seq)), rep(0, length(YY1.control.narrow.seq)))
    all.seq <- append(YY1.narrow.seq, YY1.control.narrow.seq)
    gc <- as.integer(cut(Biostrings::letterFrequency(all.seq, "CG", as.prob=T), c(-1, 0.4, 0.45, 0.5, 0.55, 0.6, 2)))
    all.weights <- c(YY1.peak$weight, rep(1, length(YY1.control.seq)))
    YY1.motif.2 <- motifRG::findMotif(all.seq,category, other.data=gc, max.motif=5,enriched=T, weights=all.weights)
    motifRG::motifLatexTable(main="Refined YY1 motifs", YY1.motif.2,prefix="YY1-2")
    ctcf.motifs <- NULL
    data(ctcf.motifs)
    ctcf.seq <- Biostrings::readDNAStringSet(system.file("extdata", "ctcf.fa",package="motifRG"))
    pwm.match <- motifRG::refinePWMMotif(ctcf.motifs$motifs[[1]]@match$pattern, ctcf.seq)
    ## require.auto('seqLogo')
    ## library(seqLogo)
    seqLogo::seqLogo(pwm.match$model$prob)
    pwm.match.extend <- motifRG::refinePWMMotifExtend(ctcf.motifs$motifs[[1]]@match$pattern, ctcf.seq)
    seqLogo::seqLogo(pwm.match.extend$model$prob)
    motifRG::plotMotif(pwm.match.extend$match$pattern)
}

simple_motifstack <- function() {
    ##require.auto('motifStack')
    pcm <- read.table(file.path(find.package("motifStack"), "extdata", "bin_SOLEXA.pcm"))
    pcm <- pcm[,3:ncol(pcm)]
    rownames(pcm) <- c("A","C","G","T")
    motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")
    ##pfm object
    ## motif <- pcm2pfm(pcm)
    ## motif <- new("pfm", mat=motif, name="bin_SOLEXA")
    opar <- par(mfrow=c(4,1))
    plot(motif)
    ## plot the logo with same height
    plot(motif, ic.scale=FALSE, ylab="probability")
    ## try a different font
    plot(motif, font="mono,Courier")
    ## try a different font and a different color group
    motif@color <- motifStack::colorset(colorScheme='basepairing')
    plot(motif,font="Times")
    par(opar)

    ##library(motifStack)
    protein <- read.table(file.path(find.package("motifStack"),"extdata","cap.txt"))
    protein <- t(protein[,1:20])
    motif <- motifStack::pcm2pfm(protein)
    motif <- new("pfm", mat=motif, name="CAP", color=motifStack::colorset(alphabet="AA",colorScheme="chemistry"))
    plot(motif)
}


flanking_sequence <- function(bsgenome, annotation, distance=200, type='gene', prefix='') {
    if (class(annotation) == 'character') { ## Assume it is a filename to a gff file
        annotations <- gff2df(annotation, type=type)
        name_key <- 'gene_id'
    } else if (class(annotation) == 'data.frame') {
        annotations <- annotation
        name_key <- 'tx_name'
    } else {
        ## Then assume it is a GenomicRanges from a TxDb or somesuch
        annotations <- as.data.frame(annotation)
        name_key <- 'tx_name'
    }
    seqinfo <- as.data.frame(bsgenome@seqinfo)
    annotations <- merge(annotations, seqinfo, by.x='seqnames', by.y='row.names', all.x=TRUE)

    before <- GenomicRanges::GRanges(seqnames=S4Vectors::Rle(annotations[,'seqnames']),
                                     ranges=IRanges::IRanges(ifelse(annotations[,'start'] <= distance,
                                                                    1,
                                                                    annotations[,'start'] - distance),
                                                             end=(annotations[,'start'] + 2)),
                      strand=S4Vectors::Rle(annotations[,'strand']),
                      name=S4Vectors::Rle(annotations[,name_key]))

    after <- GenomicRanges::GRanges(seqnames=S4Vectors::Rle(annotations[,'seqnames']),
                                    ranges=IRanges::IRanges(annotations[,'end'],
                                                            end=ifelse(annotations[,'seqlengths'] <= (annotations[,'end'] + distance),
                                                                       annotations[,'seqlengths'], annotations[,'end'] + distance)),
                                    strand=S4Vectors::Rle(annotations[,'strand']),
                                    name=S4Vectors::Rle(annotations[,name_key]))
    before_seq <- Biostrings::getSeq(bsgenome, before)
    after_seq <- Biostrings::getSeq(bsgenome, after)
    retlist <- list(before=before_seq, after=after_seq)
    return(retlist)
}

## EOF
