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
