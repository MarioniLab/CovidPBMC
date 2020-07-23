#! /usr/bin/env Rscript

### quick and dirty visualisation
library(igraph)
library(SingleCellExperiment)
library(scran)
library(igraph)
library(irlba)
library(umap)
library(reshape2)
library(cowplot)
library(ggthemes)
library(ggplot2)

all.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE.RDS")


sigaa1.donor.id <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/demultiplexed/cellSNP_vireo/SIGAA1/donor_ids.tsv",
		              sep="\t", header=TRUE, stringsAsFactors=FALSE)
sigaa1.donor.id <- sigaa1.donor.id[sigaa1.donor.id$donor_id!="doublet",]
sigaa1.donor.id$donor_id <- paste0("SIGAA1_", sigaa1.donor.id$donor_id)
sigaa1.donor.id$CellID <- paste0("SIGAA1_", sigaa1.donor.id$cell)

donor.id <- sigaa1.donor.id
donor.id <- donor.id[donor.id$CellID %in% colnames(all.sce),]
all.sce <- all.sce[,donor.id$CellID]

message("Find HVGs")
hvg.stats <- modelGeneVar(all.sce)
table(hvg.stats$FDR <= 0.01)
# set NAs to 1
hvg.stats$FDR[is.na(hvg.stats$FDR)] <- 1

message("Do PCA")
all.pca <- prcomp_irlba(t(logcounts(all.sce[hvg.stats$FDR <= 0.01 ,])), n=50)

# message("Do graph building with k=21")
# set.seed(42)
# all.snn <- buildSNNGraph(all.pca$x[, c(1:30)], k=21,
#                          d=NA, transposed=TRUE)
# all.walktrap <- cluster_walktrap(all.snn, steps=4, membership=TRUE)
# all.cluster <- data.frame("CellID"=colnames(all.sce), "Graph.Cluster"=as.character(all.walktrap$membership))
# n.clusters <- length(unique(as.character(all.walktrap$membership)))
# all.cluster$Graph.Cluster <- factor(all.cluster$Graph.Cluster,
#                                     levels=c(1:n.clusters),
#                                     labels=c(1:n.clusters))
# all.meta.cluster <- merge(donor.id, all.cluster, by='CellID')
# all.meta.cluster <- all.cluster

# message("Making force-directed layout")
# all.fr.layout <- layout_with_fr(all.snn)
# 
# all.fr.df <- as.data.frame(all.fr.layout)
# all.fr.df$CellID <- colnames(all.sce)
# all.fr.merge <- as.data.frame(merge(all.fr.df, all.meta.cluster, by='CellID'))
# 
# get the edges of the layout
# all.edges <- get.data.frame(all.snn)
# all.edges$from.Sample <- all.fr.df$CellID[all.edges$from]
# all.edges$to.Sample <- all.fr.df$CellID[all.edges$to]
# 
# match edges to vertices and graph weights
# all.edges$from.x <- all.fr.df$V1[match(all.edges$from.Sample, all.fr.df$Sample)]
# all.edges$from.y <- all.fr.df$V2[match(all.edges$from.Sample, all.fr.df$Sample)]
# all.edges$to.x <- all.fr.df$V1[match(all.edges$to.Sample, all.fr.df$Sample)]
# all.edges$to.y <- all.fr.df$V2[match(all.edges$to.Sample, all.fr.df$Sample)]


message("Making UMAP")
set.seed(42)
# can I plug in the graph co-ordinates as an initial layout?
all.umap.map <- as.data.frame(umap(as.matrix(all.pca$x[, c(1:30)]),
                                   n_neighbors=21,
                                   init='spectral',
                                   metric='cosine',
                                   n_components=2,
                                   min_dist=0.2,
                                   method='naive')$layout)

colnames(all.umap.map) <- paste0("UMAP", c(1:2))
all.umap.map$CellID <- colnames(all.sce)

all.meta.umap <- merge(donor.id, all.umap.map, by='CellID')

write.table(all.meta.umap,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/SIGAA1_UMAP_woDoublets.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")
