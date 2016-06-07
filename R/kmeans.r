## This is for the moment just a code dump of some arbitrarily chosen kmeans clustering stuff
## Fill this in asap with real code for Ginger's search of gene sets which have similar profiles over time.
kmeans_testing <- function(gene_ids=get0("gene_ids")) {
    procyclic_bins  <- data.frame()
    id_count <- 0
    procyclic <- data.frame()
    metacyclic <- data.frame()
    bins <- 2
    for (id in gene_ids) {
        id_count <- id_count + 1
        h <- hist(procyclic[procyclic$name == id,]$position, breaks=bins, include.lowest=TRUE, plot=FALSE, na.rm=TRUE)
        ## QUESTION: Is it better to use h$density or h$counts?
        ## Using h$density would be less affected by scale of the profile..
        procyclic_bins <- rbind(procyclic_bins, h$density)
    }
    rownames(procyclic_bins) <- gene_ids
    set.seed(1)
    ## k = 1
    wss <- (nrow(procyclic_bins) - 1) * sum(apply(procyclic_bins, 2, var))

    ## k between 2 and 30
    for (i in 2:30) {
        wss[i] <- sum(kmeans(procyclic_bins, centers=i, iter.max=100, nstart=50)$withinss)
    }
    plot(1:30, wss, type="b", xlab="Number of Clusters",
         ylab="Within groups sum of squares")
    k <- 6
    clusters <- kmeans(procyclic_bins, centers=k, iter.max=100, nstart=50)

    merged_data <- data.frame()
    merged_clusters <- merged_data
    merged_clusters$pro_rank <- nrow(merged_clusters) - rank(merged_data$pro_eff)
    merged_clusters$meta_rank <- nrow(merged_clusters) - rank(merged_data$meta_eff)
    merged_clusters$name <- rownames(merged_clusters)

    df <- rbind(cbind(procyclic, stage='Procyclic'),
                cbind(metacyclic, stage='Metacyclic'))

    cluster_ids <- gene_ids[clusters$cluster == 4]

    cluster_list <- list()
    for (i in seq(1, k)) {
        ## Plot cluster histogram
        cluster_ids <- gene_ids[clusters$cluster == i]
        plt_title <- sprintf("Histogram for cluster %d/%d (%d genes total)",
                             i, k, length(cluster_ids))
        plt <- ggplot2::ggplot(df[df$name %in% cluster_ids,], ggplot2::aes_string(x="position", fill="stage")) +
            ggplot2::geom_histogram(breaks=bins, colour='#333333') +
            ggplot2::facet_wrap(~ stage) +
            ggplot2::ggtitle(plt_title)
        print(plt)
        ## hist(procyclic[procyclic$gene %in% cluster_ids,]$offset, breaks=bins, main=plt_title)
        ## Determine median translational efficiency rates
        eff_subset <- merged_clusters[rownames(merged_clusters) %in% cluster_ids,]
        ## Median translational efficiency rates and rankings
        avg_pro_rate <- median(eff_subset$pro_eff)
        avg_pro_rank <- round(median(eff_subset$pro_rank))
        avg_meta_rate <- median(eff_subset$meta_eff)
        avg_meta_rank <- round(median(eff_subset$meta_rank))
        averages <- data.frame("average_efficiency"=c(avg_pro_rate, avg_meta_rate),
                               "average_ranking"=c(avg_pro_rank, avg_meta_rank))
        cluster_list[i] <- list("averages"=averages, "subset"=eff_subset)
        rownames(averages) <- c("Procyclic", "Metacyclic")
        cat(sprintf("### Median translational efficiencies (cluster %d)\n", i))
        knitr::kable(head(averages))
        cat(sprintf("\n### Genes in cluster %d\n", i))
        knitr::kable(head(eff_subset))
        eff_subset_ids <- eff_subset[,c("name","pro_rank")]
        colnames(eff_subset_ids) <- c("ID","rank")
    }
    dev.off()
}

## EOF
