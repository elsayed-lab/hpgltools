gostats_kegg = function() {
    frame = toTable(org.Hs.egPATH)
    keggframedata = data.frame(frame$path_id, frame$gene_id)
    keggFrame=KEGGFrame(keggframeData,organism="Homo sapiens")
    gsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())
    universe = Lkeys(org.Hs.egGO)
    genes = universe[1:500]
    kparams <- GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params",
                                    geneSetCollection=gsc,
                                    geneIds = genes,
                                    universeGeneIds = universe,
                                    pvalueCutoff = 0.05,
                                    testDirection = "over")
    kOver <- hyperGTest(params)
    head(summary(kOver))
}

## Please note that the KGML parser fails if other XML parsers are loaded into R
#' Print some data onto KEGG pathways
#'
#' @param de_genes some differentially expressed genes
#' @param indir A directory into which the unmodified kegg images will be downloaded (or already exist).  Defaults to 'pathview_in'
#' @param outdir A directory which will contain the colored images.
#' @param pathway Perform the coloring for a specific pathway?  Defaults to 'all'
#' @param species The kegg identifier for the species of interest.
#' @param filenames Whether to give filenames as the kegg ID or pathway name.  Defaults to 'id'.
#' 
#' @return A list of some information for every KEGG pathway downloaded/examined.  This information includes:
#'   a. The filename of the final image for each pathway.
#'   b. The number of genes which were found in each pathway image.
#'   c. The number of genes in the 'up' category
#'   d. The number of genes in the 'down' category
#' @seealso \code{\link{Ramigo}}
#' @export
#' @examples
#' thy_el_comp2_path = hpgl_pathview(thy_el_comp2_kegg, species="spz", indir="pathview_in", outdir="kegg_thy_el_comp2", string_from="_Spy", string_to="_Spy_", filenames="pathname")
hpgl_pathview = function(path_data, indir="pathview_in", outdir="pathview", pathway="all", species="lma", string_from="LmjF", string_to="LMJF", suffix="_colored", second_from=NULL, second_to=NULL, verbose=FALSE, filenames="id") {
    ##eh = new.env(hash=TRUE, size=NA)
    ## There is a weird namespace conflict when using pathview, so I will reload it here
    try(detach("package:Rgraphviz", unload=TRUE))
    try(detach("package:topGO", unload=TRUE))
    try(detach("package:pathview", unload=TRUE))
    try(detach("package:KEGGgraph", unload=TRUE))
    try(detach("package:RamiGO", unload=TRUE))
    try(detach("package:graph", unload=TRUE))
    require.auto("pathview")
    tmp_names = names(path_data)
    tmp_names = gsub(string_from, string_to, tmp_names)
    if (!is.null(second_from)) {
        tmp_names = gsub(second_from, second_to, tmp_names)
    }
##    tmp_names = gsub("\\.","_", tmp_names)
    names(path_data) = tmp_names
    rm(tmp_names)
    
    ## First check that the input pathview directory exists
    if (!file.exists(indir)) {
        dir.create(indir)
    }
    if (!file.exists(outdir)){
        dir.create(outdir)
    }
    paths = list()
    if (pathway == "all") {
        all_pathways = unique(KEGGREST::keggLink("pathway", species))
        paths = all_pathways
        paths = gsub("path:", "", paths)
        all_modules = unique(KEGGREST::keggLink("module", species))
    } else if (class(pathway) == "list") {
        paths = pathway
    } else {
        paths[1] = pathway
    }
    return_list = list()
    for (count in 1:length(paths)) {
        path = paths[count]
        path_name = keggGet(path)
        path_name = path_name[[1]]$NAME
        path_name = gsub("(.*) - .*", "\\1", path_name)
        path_name = tolower(path_name)
        path_name = gsub(" ", "_", path_name)
        gene_examples = try(keggLink(paste("path", path, sep=":"))[,2])  ## RCurl is crap and fails sometimes for no apparent reason.
        limits=c(min(path_data, na.rm=TRUE), max(path_data, na.rm=TRUE))
        if (isTRUE(verbose)) {
            print(paste("Here are some path gene examples: ", gene_examples, sep=""))
            print(paste("Here are your genes: ", head(names(path_data))), sep="")
            pv = try(pathview::pathview(gene.data=path_data, kegg.dir=indir, pathway.id=path, species=species, limit=list(gene=limits, cpd=limits), map.null=TRUE, gene.idtype="KEGG", out.suffix=suffix, split.group=TRUE, expand.node=TRUE, kegg.native=TRUE, map.symbol=TRUE, same.layer=FALSE, res=1200, new.signature=FALSE, cex=0.05, key.pos="topright"))
        } else {
            pv = suppressMessages(try(pathview::pathview(gene.data=path_data, kegg.dir=indir, pathway.id=path, species=species, limit=list(gene=limits, cpd=limits), map.null=TRUE, gene.idtype="KEGG", out.suffix=suffix, split.group=TRUE, expand.node=TRUE, kegg.native=TRUE, map.symbol=TRUE, same.layer=FALSE, res=1200, new.signature=FALSE, cex=0.05, key.pos="topright")))
        }
        if (class(pv) == "numeric") {
            colored_genes = NULL
            newfile = NULL
            up = NULL
            down = NULL
        } else {
            colored_genes = dim(pv$plot.data.gene)[1]
            ##        "lma04070._proeff.png"
            oldfile = paste(path, ".", suffix, ".png", sep="")
            ## An if-statement to see if the user prefers pathnames by kegg ID or pathway name
            ## Dr. McIver wants path names...
            if (filenames == "id") {
                newfile = paste(outdir,"/", path, suffix, ".png", sep="")
            } else {
                newfile = paste(outdir, "/", path_name, suffix, ".png", sep="")
            }
            if (isTRUE(verbose)) {
                message(paste("Moving file to: ", newfile, sep=""))
            }
            file.rename(from=oldfile, to=newfile)
            data_low = summary(path_data)[2]
            data_high = summary(path_data)[3]
            numbers_in_plot = as.numeric(pv$plot.data.gene$mol.data)
            up = sum(numbers_in_plot > data_high, na.rm=TRUE)
            down = sum(numbers_in_plot < data_low, na.rm=TRUE)
        }
        return_list[[path]]$file = newfile
        return_list[[path]]$genes = colored_genes
        return_list[[path]]$up = up
        return_list[[path]]$down = down
    }

    retdf = data.frame(rep(NA, length(names(return_list))))
    rownames(retdf) = names(return_list)
    retdf$genes = NA
    retdf$up = NA
    retdf$down = NA
    colnames(retdf) = c("file","genes","up","down")
    for (path in names(return_list)) {
        if (is.null(return_list[[path]]$genes)) {
            retdf[path,]$genes = 0
        } else {
            retdf[path,]$genes = as.numeric(return_list[[path]]$genes)
        }
        retdf[path,]$file = as.character(return_list[[path]]$file)
        retdf[path,]$up = as.numeric(return_list[[path]]$up)
        retdf[path,]$down = as.numeric(return_list[[path]]$down)        
    }
    retdf = retdf[with(retdf, order(up, down)), ]
    return(retdf)
}

#' Search the kegg identifier for a given species
#'
#' @param species A search string (Something like 'Homo sapiens')
#' 
#' @return a data frame of possible KEGG identifier codes, genome ID numbers, species, and phylogenetic classifications.
#' @seealso \code{\link{RCurl}}
#' @export
#' @examples
#' ## fun = kegg_get_orgn('Canis')
#' ## >     Tid     orgid      species                   phylogeny
#' ## >  17 T01007   cfa Canis familiaris (dog) Eukaryotes;Animals;Vertebrates;Mammals
kegg_get_orgn = function(species="Leishmania") {
    all_organisms = getURL("http://rest.kegg.jp/list/organism")
    org_tsv = textConnection(all_organisms)
    all = read.table(org_tsv, sep="\t", quote="", fill=TRUE)
    close(org_tsv)
    colnames(all) = c("Tid","orgid","species","phylogeny")
    candidates = all[grepl(species, all$species),]  ## Look for the search string in the species column
    return(candidates)
}


parseKGML2Graph2 <-function (file, ...) {
    pathway <- parseKGML2(file)
    gR <- KEGGpathway2Graph2(pathway, ...)
    return(gR)
}

hpgl_base_pathview = function (gene.data = NULL, cpd.data = NULL, xml.file = NULL, 
    pathway.id, species = "hsa", kegg.dir = ".", cpd.idtype = "kegg", 
    gene.idtype = "entrez", gene.annotpkg = NULL, min.nnodes = 3, 
    kegg.native = TRUE, map.null = TRUE, expand.node = FALSE, 
    split.group = FALSE, map.symbol = TRUE, map.cpdname = TRUE, 
    node.sum = "sum", discrete = list(gene = FALSE, cpd = FALSE), 
    limit = list(gene = 1, cpd = 1), bins = list(gene = 10, cpd = 10), 
    both.dirs = list(gene = T, cpd = T), trans.fun = list(gene = NULL, 
        cpd = NULL), low = list(gene = "green", cpd = "blue"), 
    mid = list(gene = "gray", cpd = "gray"), high = list(gene = "red", 
        cpd = "yellow"), na.col = "transparent", ...) {
    if (is.character(gene.data)) {
        gd.names = gene.data
        gene.data = rep(1, length(gene.data))
        names(gene.data) = gd.names
        both.dirs$gene = FALSE
        ng = length(gene.data)
        nsamp.g = 1
    }
    else if (!is.null(gene.data)) {
        if (length(dim(gene.data)) == 2) {
            gd.names = rownames(gene.data)
            ng = nrow(gene.data)
            nsamp.g = 2
        }
        else if (is.numeric(gene.data) & is.null(dim(gene.data))) {
            gd.names = names(gene.data)
            ng = length(gene.data)
            nsamp.g = 1
        }
        else stop("wrong gene.data format!")
    }
    else if (is.null(cpd.data)) {
        stop("gene.data and cpd.data are both NULL!")
    }
    gene.idtype = toupper(gene.idtype)
    data(bods)
    data(gene.idtype.list)
    if (species != "ko") {
        species.data = pathview::kegg.species.code(species, na.rm = T, 
            code.only = FALSE)
    }
    else {
        species.data = c(kegg.code = "ko", entrez.gnodes = "0", 
            kegg.geneid = "K01488", ncbi.geneid = "")
        gene.idtype = "KEGG"
        msg.fmt = "Only KEGG ortholog gene ID is supported, make sure it looks like \"%s\"!"
        msg = sprintf(msg.fmt, species.data["kegg.geneid"])
        message(msg)
    }
    if (length(dim(species.data)) == 2) {
        message("More than two valide species!")
        species.data = species.data[1, ]
    }
    species = species.data["kegg.code"]
    entrez.gnodes = species.data["entrez.gnodes"] == 1
    if (is.na(species.data["ncbi.geneid"])) {
        if (!is.na(species.data["kegg.geneid"])) {
            msg.fmt = "Only native KEGG gene ID is supported for this species,\nmake sure it looks like \"%s\"!"
            msg = sprintf(msg.fmt, species.data["kegg.geneid"])
            message(msg)
        }
        else {
            stop("This species is not annotated in KEGG!")
        }
    }
    if (is.null(gene.annotpkg)) 
        gene.annotpkg = bods[match(species, bods[, 3]), 1]
    if (length(grep("ENTREZ|KEGG", gene.idtype)) < 1 & !is.null(gene.data)) {
        if (is.na(gene.annotpkg)) 
            stop("No proper gene annotation package available!")
        if (!gene.idtype %in% gene.idtype.list) 
            stop("Wrong input gene ID type!")
        gene.idmap = id2eg(gd.names, category = gene.idtype, 
            pkg.name = gene.annotpkg)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "ENTREZ"
    }
    if (gene.idtype == "ENTREZ" & !entrez.gnodes & !is.null(gene.data)) {
        message("Getting gene ID data from KEGG...")
        gene.idmap = keggConv("ncbi-geneid", species)
        message("Done with data retrieval!")
        kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
        ncbi.ids = gsub("ncbi-geneid:", "", gene.idmap)
        gene.idmap = cbind(ncbi.ids, kegg.ids)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "KEGG"
    }
    if (is.character(cpd.data)) {
        cpdd.names = cpd.data
        cpd.data = rep(1, length(cpd.data))
        names(cpd.data) = cpdd.names
        both.dirs$cpd = FALSE
        ncpd = length(cpd.data)
    }
    else if (!is.null(cpd.data)) {
        if (length(dim(cpd.data)) == 2) {
            cpdd.names = rownames(cpd.data)
            ncpd = nrow(cpd.data)
        }
        else if (is.numeric(cpd.data) & is.null(dim(cpd.data))) {
            cpdd.names = names(cpd.data)
            ncpd = length(cpd.data)
        }
        else stop("wrong cpd.data format!")
    }
    if (length(grep("kegg", cpd.idtype)) < 1 & !is.null(cpd.data)) {
        data(rn.list)
        cpd.types = c(names(rn.list), "name")
        cpd.types = tolower(cpd.types)
        cpd.types = cpd.types[-grep("kegg", cpd.types)]
        if (!tolower(cpd.idtype) %in% cpd.types) 
            stop("Wrong input cpd ID type!")
        cpd.idmap = cpd2kegg(cpdd.names, in.type = cpd.idtype)
        cpd.data = mol.sum(cpd.data, cpd.idmap)
    }
    warn.fmt = "Parsing %s file failed, please check the file!"
    if (length(grep(species, pathway.id)) > 0) {
        pathway.name = pathway.id
        pathway.id = gsub(species, "", pathway.id)
    }
    else pathway.name = paste(species, pathway.id, sep = "")
    kfiles = list.files(path = kegg.dir, pattern = "[.]xml|[.]png")
    tfiles = paste(pathway.name, c("xml", "png"), sep = ".")
    if (!all(tfiles %in% kfiles)) {
        dstatus = download.kegg(pathway.id = pathway.id, species = species, 
            kegg.dir = kegg.dir)
        if (dstatus == "failed") {
            warn.fmt = "Failed to download KEGG xml/png files, %s skipped!"
            warn.msg = sprintf(warn.fmt, pathway.name)
            message(warn.msg)
            return(invisible(0))
        }
    }
    if (missing(xml.file)) 
        xml.file <- paste(kegg.dir, "/", pathway.name, ".xml", 
            sep = "")
    if (kegg.native) {
        node.data = try(node.info(xml.file), silent = T)
        if (class(node.data) == "try-error") {
            warn.msg = sprintf(warn.fmt, xml.file)
            message(warn.msg)
            return(invisible(0))
        }
        node.type = c("gene", "enzyme", "compound", "ortholog")
        sel.idx = node.data$type %in% node.type
        nna.idx = !is.na(node.data$x + node.data$y + node.data$width + 
            node.data$height)
        sel.idx = sel.idx & nna.idx
        if (sum(sel.idx) < min.nnodes) {
            warn.fmt = "Number of mappable nodes is below %d, %s skipped!"
            warn.msg = sprintf(warn.fmt, min.nnodes, pathway.name)
            message(warn.msg)
            return(invisible(0))
        }
        node.data = lapply(node.data, "[", sel.idx)
    }
    else {
        gR1 = try(parseKGML2Graph2(xml.file, genes = F, expand = expand.node, 
            split.group = split.group), silent = T)
        node.data = try(node.info(gR1), silent = T)
        if (class(node.data) == "try-error") {
            warn.msg = sprintf(warn.fmt, xml.file)
            message(warn.msg)
            return(invisible(0))
        }
    }
    if (species == "ko") 
        gene.node.type = "ortholog"
    else gene.node.type = "gene"
    if ((!is.null(gene.data) | map.null) & sum(node.data$type == 
        gene.node.type) > 1) {
        plot.data.gene = node.map(gene.data, node.data, node.types = gene.node.type, 
            node.sum = node.sum, entrez.gnodes = entrez.gnodes)
        kng = plot.data.gene$kegg.names
        kng.char = gsub("[0-9]", "", unlist(kng))
        if (any(kng.char > "")) 
            entrez.gnodes = FALSE
        if (map.symbol & species != "ko" & entrez.gnodes) {
            if (is.na(gene.annotpkg)) {
                warm.fmt = "No annotation package for the species %s, gene symbols not mapped!"
                warm.msg = sprintf(warm.fmt, species)
                message(warm.msg)
            }
            else {
                plot.data.gene$labels = eg2id(as.character(plot.data.gene$kegg.names), 
                  category = "SYMBOL", pkg.name = gene.annotpkg)[, 
                  2]
                mapped.gnodes = rownames(plot.data.gene)
                node.data$labels[mapped.gnodes] = plot.data.gene$labels
            }
        }
        cols.ts.gene = node.color(plot.data.gene, limit$gene, 
            bins$gene, both.dirs = both.dirs$gene, trans.fun = trans.fun$gene, 
            discrete = discrete$gene, low = low$gene, mid = mid$gene, 
            high = high$gene, na.col = na.col)
    }
    else plot.data.gene = cols.ts.gene = NULL
    if ((!is.null(cpd.data) | map.null) & sum(node.data$type == 
        "compound") > 1) {
        plot.data.cpd = node.map(cpd.data, node.data, node.types = "compound", 
            node.sum = node.sum)
        if (map.cpdname & !kegg.native) {
            plot.data.cpd$labels = cpdkegg2name(plot.data.cpd$labels)[, 
                2]
            mapped.cnodes = rownames(plot.data.cpd)
            node.data$labels[mapped.cnodes] = plot.data.cpd$labels
        }
        cols.ts.cpd = node.color(plot.data.cpd, limit$cpd, bins$cpd, 
            both.dirs = both.dirs$cpd, trans.fun = trans.fun$cpd, 
            discrete = discrete$cpd, low = low$cpd, mid = mid$cpd, 
            high = high$cpd, na.col = na.col)
    }
    else plot.data.cpd = cols.ts.cpd = NULL
    if (kegg.native) {
        pv.pars = keggview.native(plot.data.gene = plot.data.gene, 
            cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
            cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
            pathway.name = pathway.name, kegg.dir = kegg.dir, 
            limit = limit, bins = bins, both.dirs = both.dirs, 
            discrete = discrete, low = low, mid = mid, high = high, 
            na.col = na.col, ...)
    }
    else {
        pv.pars = keggview.graph(plot.data.gene = plot.data.gene, 
            cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
            cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
            path.graph = gR1, pathway.name = pathway.name, map.cpdname = map.cpdname, 
            split.group = split.group, limit = limit, bins = bins, 
            both.dirs = both.dirs, discrete = discrete, low = low, 
            mid = mid, high = high, na.col = na.col, ...)
    }
    plot.data.gene = cbind(plot.data.gene, cols.ts.gene)
    if (!is.null(plot.data.gene)) {
        cnames = colnames(plot.data.gene)[-(1:7)]
        nsamp = length(cnames)/2
        if (nsamp > 1) {
            cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                1):(2 * nsamp)], "col", sep = ".")
        }
        else cnames[2] = "mol.col"
        colnames(plot.data.gene)[-(1:7)] = cnames
    }
    plot.data.cpd = cbind(plot.data.cpd, cols.ts.cpd)
    if (!is.null(plot.data.cpd)) {
        cnames = colnames(plot.data.cpd)[-(1:7)]
        nsamp = length(cnames)/2
        if (nsamp > 1) {
            cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                1):(2 * nsamp)], "col", sep = ".")
        }
        else cnames[2] = "mol.col"
        colnames(plot.data.cpd)[-(1:7)] = cnames
    }
    return(invisible(list(plot.data.gene = plot.data.gene, plot.data.cpd = plot.data.cpd)))
}
