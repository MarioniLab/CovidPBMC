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
sigaa1.abs.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE-SIGAA1_Ab.RDS")
sigaa2.abs.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE-SIGAA2_Ab.RDS")
sigaa3.abs.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE-SIGAA3_Ab.RDS")

abs.sce <- do.call(cbind, list("SIGAA1"=sigaa1.abs.sce, "SIGAA2"=sigaa2.abs.sce, "SIGAA3"=sigaa3.abs.sce))
saveRDS(abs.sce, "/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE_ABS.RDS")

print(head(rowData(all.sce)))
print(head(colData(all.sce)))

sigaa1.donor.id <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/demultiplexed/cellSNP_vireo/SIGAA1/donor_ids.tsv",
		              sep="\t", header=TRUE, stringsAsFactors=FALSE)
sigaa1.donor.id$donor_id <- paste0("SIGAA1_", sigaa1.donor.id$donor_id)
sigaa1.donor.id$CellID <- paste0("SIGAA1_", sigaa1.donor.id$cell)


sigaa2.donor.id <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/demultiplexed/cellSNP_vireo/SIGAA2/donor_ids.tsv",
		              sep="\t", header=TRUE, stringsAsFactors=FALSE)
sigaa2.donor.id$donor_id <- paste0("SIGAA2_", sigaa2.donor.id$donor_id)
sigaa2.donor.id$CellID <- paste0("SIGAA2_", sigaa2.donor.id$cell)


sigaa3.donor.id <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/demultiplexed/cellSNP_vireo/SIGAA3/donor_ids.tsv",
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)
sigaa3.donor.id$donor_id <- paste0("SIGAA3_", sigaa3.donor.id$donor_id)
sigaa3.donor.id$CellID <- paste0("SIGAA3_", sigaa3.donor.id$cell)

donor.id <- do.call(rbind.data.frame, list("SIGAA1"=sigaa1.donor.id "SIGAA2"=sigaa2.donor.id, "SIGAA3"=sigaa3.donor.id))

message("Find HVGs")
hvg.stats <- modelGeneVar(all.sce)
table(hvg.stats$FDR <= 0.01)
# set NAs to 1
hvg.stats$FDR[is.na(hvg.stats$FDR)] <- 1

message("Do PCA")
all.pca <- prcomp_irlba(t(logcounts(all.sce[hvg.stats$FDR <= 0.01 ,])), n=50)

message("Do graph building with k=21")
set.seed(42)
all.snn <- buildSNNGraph(all.pca$x[, c(1:30)], k=21,
                         d=NA, transposed=TRUE)
all.walktrap <- cluster_walktrap(all.snn, steps=4, membership=TRUE)
all.cluster <- data.frame("CellID"=colnames(all.sce), "Graph.Cluster"=as.character(all.walktrap$membership))
n.clusters <- length(unique(as.character(all.walktrap$membership)))
all.cluster$Graph.Cluster <- factor(all.cluster$Graph.Cluster,
                                    levels=c(1:n.clusters),
                                    labels=c(1:n.clusters))
all.meta.cluster <- merge(donor.id, all.cluster, by='CellID')
# all.meta.cluster <- all.cluster

message("Making force-directed layout")
all.fr.layout <- layout_with_fr(all.snn)

all.fr.df <- as.data.frame(all.fr.layout)
all.fr.df$CellID <- colnames(all.sce)
all.fr.merge <- as.data.frame(merge(all.fr.df, all.meta.cluster, by='CellID'))

# get the edges of the layout
all.edges <- get.data.frame(all.snn)
all.edges$from.Sample <- all.fr.df$CellID[all.edges$from]
all.edges$to.Sample <- all.fr.df$CellID[all.edges$to]

# match edges to vertices and graph weights
all.edges$from.x <- all.fr.df$V1[match(all.edges$from.Sample, all.fr.df$Sample)]
all.edges$from.y <- all.fr.df$V2[match(all.edges$from.Sample, all.fr.df$Sample)]
all.edges$to.x <- all.fr.df$V1[match(all.edges$to.Sample, all.fr.df$Sample)]
all.edges$to.y <- all.fr.df$V2[match(all.edges$to.Sample, all.fr.df$Sample)]


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

all.meta.umap <- merge(all.fr.merge, all.umap.map, by='CellID')

write.table(all.meta.umap,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/Covid-UMAP.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

## just use the Ab information - is it concordant?
message("Making projections from protein-only")
message("Do ADT PCA")
abs.pca <- prcomp_irlba(t(logcounts(abs.sce)), n=50)

message("Do ADT graph building with k=21")
set.seed(42)
abs.snn <- buildSNNGraph(abs.pca$x[, c(1:30)], k=21,
                         d=NA, transposed=TRUE)
abs.walktrap <- cluster_walktrap(abs.snn, steps=4, membership=TRUE)
abs.cluster <- data.frame("CellID"=colnames(abs.sce), "ADT.Cluster"=as.character(abs.walktrap$membership))
n.clusters <- length(unique(as.character(abs.walktrap$membership)))
abs.cluster$ADT.Cluster <- factor(abs.cluster$ADT.Cluster,
                                    levels=c(1:n.clusters),
                                    labels=c(1:n.clusters))
# abs.meta.cluster <- abs.cluster
abs.meta.cluster <- merge(donor.id, abs.cluster, by='CellID')

message("Making ADT force-directed layout")
abs.fr.layout <- layout_with_fr(abs.snn)

abs.fr.df <- as.data.frame(abs.fr.layout)
abs.fr.df$CellID <- colnames(abs.sce)
abs.fr.merge <- as.data.frame(merge(abs.fr.df, abs.meta.cluster, by='CellID'))

# get the edges of the layout
abs.edges <- get.data.frame(abs.snn)
abs.edges$from.Sample <- abs.fr.df$CellID[abs.edges$from]
abs.edges$to.Sample <- abs.fr.df$CellID[abs.edges$to]

# match edges to vertices and graph weights
abs.edges$from.x <- abs.fr.df$V1[match(abs.edges$from.Sample, abs.fr.df$Sample)]
abs.edges$from.y <- abs.fr.df$V2[match(abs.edges$from.Sample, abs.fr.df$Sample)]
abs.edges$to.x <- abs.fr.df$V1[match(abs.edges$to.Sample, abs.fr.df$Sample)]
abs.edges$to.y <- abs.fr.df$V2[match(abs.edges$to.Sample, abs.fr.df$Sample)]

message("Making ADT UMAP")
set.seed(42)
# can I plug in the graph co-ordinates as an initial layout?
abs.umap.map <- as.data.frame(umap(as.matrix(abs.pca$x[, c(1:30)]),
                                   n_neighbors=21,
                                   init='spectral',
                                   metric='cosine',
                                   n_components=2,
                                   min_dist=0.2,
                                   method='naive')$layout)

colnames(abs.umap.map) <- paste0("ADT.UMAP", c(1:2))
abs.umap.map$CellID <- colnames(abs.sce)

abs.meta.umap <- merge(abs.fr.merge, abs.umap.map, by='CellID')

write.table(abs.meta.umap,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/Covid_Abs-UMAP.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")


##############################################################################
## compute these without the doublets
##############################################################################
doublet.cells <- donor.id$CellID[grepl(donor.id$donor_id, pattern="doublet")]
keep.cells <- donor.id$CellID[!grepl(donor.id$donor_id, pattern="doublet")]

message(paste0("Computing after removing ", length(doublet.cells), " doublets"))
message("Find HVGs")
hvg.stats <- modelGeneVar(all.sce[, intersect(colnames(all.sce), keep.cells)])
table(hvg.stats$FDR <= 0.01)
# set NAs to 1
hvg.stats$FDR[is.na(hvg.stats$FDR)] <- 1

message("Do PCA")
all.pca <- prcomp_irlba(t(logcounts(all.sce[hvg.stats$FDR <= 0.01, intersect(colnames(all.sce), keep.cells)])), n=50)

message("Do graph building with k=21")
set.seed(42)
all.snn <- buildSNNGraph(all.pca$x[, c(1:30)], k=21,
                         d=NA, transposed=TRUE)
all.walktrap <- cluster_walktrap(all.snn, steps=4, membership=TRUE)
all.cluster <- data.frame("CellID"=colnames(all.sce[, intersect(colnames(all.sce), keep.cells)]), "Graph.Cluster"=as.character(all.walktrap$membership))
n.clusters <- length(unique(as.character(all.walktrap$membership)))
all.cluster$Graph.Cluster <- factor(all.cluster$Graph.Cluster,
                                    levels=c(1:n.clusters),
                                    labels=c(1:n.clusters))
all.meta.cluster <- merge(donor.id, all.cluster, by='CellID')
# all.meta.cluster <- all.cluster

message("Making force-directed layout")
all.fr.layout <- layout_with_fr(all.snn)

all.fr.df <- as.data.frame(all.fr.layout)
all.fr.df$CellID <- colnames(all.sce[, intersect(colnames(all.sce), keep.cells)])
all.fr.merge <- as.data.frame(merge(all.fr.df, all.meta.cluster, by='CellID'))

# get the edges of the layout
all.edges <- get.data.frame(all.snn)
all.edges$from.Sample <- all.fr.df$CellID[all.edges$from]
all.edges$to.Sample <- all.fr.df$CellID[all.edges$to]

# match edges to vertices and graph weights
all.edges$from.x <- all.fr.df$V1[match(all.edges$from.Sample, all.fr.df$Sample)]
all.edges$from.y <- all.fr.df$V2[match(all.edges$from.Sample, all.fr.df$Sample)]
all.edges$to.x <- all.fr.df$V1[match(all.edges$to.Sample, all.fr.df$Sample)]
all.edges$to.y <- all.fr.df$V2[match(all.edges$to.Sample, all.fr.df$Sample)]


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
all.umap.map$CellID <- colnames(all.sce[, intersect(colnames(all.sce), keep.cells)])

all.meta.umap <- merge(all.fr.merge, all.umap.map, by='CellID')

write.table(all.meta.umap,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/Covid-UMAP_sansDoublets.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

## just use the Ab information - is it concordant?
message("Making projections from protein-only")
message("Do ADT PCA")
abs.pca <- prcomp_irlba(t(logcounts(abs.sce[, intersect(colnames(abs.sce), keep.cells)])), n=50)

message("Do ADT graph building")
set.seed(42)
abs.snn <- buildSNNGraph(abs.pca$x[, c(1:30)], k=21,
                         d=NA, transposed=TRUE)
abs.walktrap <- cluster_walktrap(abs.snn, steps=4, membership=TRUE)
abs.cluster <- data.frame("CellID"=colnames(abs.sce[, intersect(colnames(abs.sce), keep.cells)]), "ADT.Cluster"=as.character(abs.walktrap$membership))
n.clusters <- length(unique(as.character(abs.walktrap$membership)))
abs.cluster$ADT.Cluster <- factor(abs.cluster$ADT.Cluster,
                                  levels=c(1:n.clusters),
                                  labels=c(1:n.clusters))
# abs.meta.cluster <- abs.cluster
abs.meta.cluster <- merge(donor.id, abs.cluster, by='CellID')

message("Making ADT force-directed layout")
abs.fr.layout <- layout_with_fr(abs.snn)

abs.fr.df <- as.data.frame(abs.fr.layout)
abs.fr.df$CellID <- colnames(abs.sce[, intersect(colnames(abs.sce), keep.cells)])
abs.fr.merge <- as.data.frame(merge(abs.fr.df, abs.meta.cluster, by='CellID'))

# get the edges of the layout
abs.edges <- get.data.frame(abs.snn)
abs.edges$from.Sample <- abs.fr.df$CellID[abs.edges$from]
abs.edges$to.Sample <- abs.fr.df$CellID[abs.edges$to]

# match edges to vertices and graph weights
abs.edges$from.x <- abs.fr.df$V1[match(abs.edges$from.Sample, abs.fr.df$Sample)]
abs.edges$from.y <- abs.fr.df$V2[match(abs.edges$from.Sample, abs.fr.df$Sample)]
abs.edges$to.x <- abs.fr.df$V1[match(abs.edges$to.Sample, abs.fr.df$Sample)]
abs.edges$to.y <- abs.fr.df$V2[match(abs.edges$to.Sample, abs.fr.df$Sample)]

message("Making ADT UMAP")
set.seed(42)
# can I plug in the graph co-ordinates as an initial layout?
abs.umap.map <- as.data.frame(umap(as.matrix(abs.pca$x[, c(1:30)]),
                                   n_neighbors=21,
                                   init='spectral',
                                   metric='cosine',
                                   n_components=2,
                                   min_dist=0.2,
                                   method='naive')$layout)

colnames(abs.umap.map) <- paste0("ADT.UMAP", c(1:2))
abs.umap.map$CellID <- colnames(abs.sce[, intersect(colnames(abs.sce), keep.cells)])

abs.meta.umap <- merge(abs.fr.merge, abs.umap.map, by='CellID')

write.table(abs.meta.umap,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/Covid_Abs-UMAP_sansDoublets.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")



