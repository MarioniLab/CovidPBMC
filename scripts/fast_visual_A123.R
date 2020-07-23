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

donor.id <- do.call(rbind.data.frame, list("SIGAA1"=sigaa1.donor.id, "SIGAA2"=sigaa2.donor.id, "SIGAA3"=sigaa3.donor.id))

donor.id <- donor.id[!grepl("doublet",donor.id$donor_id),]
donor.id <- donor.id[donor.id$CellID %in% colnames(all.sce),]
all.sce <- all.sce[,donor.id$CellID]

message("Find HVGs")
hvg.stats <- modelGeneVar(all.sce)
table(hvg.stats$FDR <= 0.01)
# set NAs to 1
hvg.stats$FDR[is.na(hvg.stats$FDR)] <- 1

message("Do PCA")
all.pca <- prcomp_irlba(t(logcounts(all.sce[hvg.stats$FDR <= 0.01 ,])), n=50)

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
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/SIGAA1_SIGAA2_SIGAA3-UMAPwoDoublets.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

#### Add graph construction from GEX ######
message("Constructing GEX SNN-graph")

broad.gex.snn <- buildSNNGraph(all.pca$x[, c(1:50)], d=NA, transposed=TRUE, k=21)
broad.gex.walktrap <- cluster_walktrap(broad.gex.snn, steps=4, membership=TRUE)
broad.clusters <- as.character(membership(broad.gex.walktrap))
n.gex.broad <- length(unique(broad.clusters))
broad.clust.df <- data.frame("CellID"=colnames(all.sce),
                             "GEX.Broad.Cluster"=broad.clusters)
broad.clust.df$GEX.Broad.Cluster <- ordered(broad.clust.df$GEX.Broad.Cluster,
                                            levels=c(1:n.gex.broad))

write.table(broad.clust.df,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/COVID_GEX_BroadCluster.tsv",
	    quote=FALSE, row.names=FALSE, sep="\t")

all.pcs <- as.data.frame(all.pca$x)
all.pcs$CellID <- colnames(all.sce)

write.table(all.pcs,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/COVID_GEX_PCA.tsv",
	    quote=FALSE, row.names=FALSE, sep="\t")

rm(list=c("all.pca", "broad.gex.snn"))
gc()

#### Add graph construction from ADT ######
message("Constructing ADT SNN-graph")

adt.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE_Ab.RDS")
adt.sce <- adt.sce[, intersect(colnames(adt.sce), colnames(all.sce))]
adt.cpm <- apply(counts(adt.sce), 2, FUN=function(X) log10((X+1)/((sum(X)+1)/(1e6+1))))
set.seed(42)
adt.pca <- prcomp_irlba(t(adt.cpm), n=50)

broad.adt.snn <- buildSNNGraph(adt.pca$x[, c(1:50)], d=NA, transposed=TRUE, k=21)
broad.adt.walktrap <- cluster_walktrap(broad.adt.snn, steps=4, membership=TRUE)
broad.clusters <- as.character(membership(broad.adt.walktrap))
n.adt.broad <- length(unique(broad.clusters))
broad.clust.df <- data.frame("CellID"=colnames(adt.sce),
                             "ADT.Broad.Cluster"=broad.clusters)
broad.clust.df$ADT.Broad.Cluster <- ordered(broad.clust.df$ADT.Broad.Cluster,
                                            levels=c(1:n.adt.broad))

write.table(broad.clust.df,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/COVID_ADT_BroadCluster.tsv",
	    quote=FALSE, row.names=FALSE, sep="\t")

adt.pcs <- as.data.frame(adt.pca$x)
adt.pcs$CellID <- colnames(adt.sce)

write.table(adt.pcs,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/COVID_ADT_PCA.tsv",
	    quote=FALSE, row.names=FALSE, sep="\t")
