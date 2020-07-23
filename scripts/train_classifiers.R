#! /usr/bin/env Rscript

## Train a series of random forests on different subsets of data from Catherine Blish's lab
library(igraph)
library(SingleCellExperiment)
library(scran)
library(randomForest)
library(caret)
library(biomaRt)

# I need to make sure that we are only using features present in both the training and testing sets
message("Fetching testing data features")
covid.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE.RDS")
covid.genes <- rownames(covid.sce)
rm(list=c("covid.sce"))
gc()

message("Retrieving gene annotation information")
biomaRt.connection <- useMart("ensembl", "hsapiens_gene_ensembl")

gene.df <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                 mart=biomaRt.connection)
rownames(gene.df) <- gene.df$ensembl_gene_id

message("Importing meta-data information")
covid.meta <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/external_data/Blish_meta.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)

message("Importing training data SCE")
# set genes to ensembl IDs
covid.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/external_data/Blish_SCE.RDS")
new.row.data <- rowData(covid.sce)

message("Converting gene names to Ensembl IDs")
new.row.data$external_gene_name <- rownames(new.row.data)
new.row.data <- merge(new.row.data, gene.df, by='external_gene_name', all.x=TRUE)
new.row.data <- new.row.data[!duplicated(new.row.data$external_gene_name), ]
rownames(new.row.data) <- new.row.data$external_gene_name
rowData(covid.sce) <- new.row.data[rownames(covid.sce), ]
covid.sce <- covid.sce[!is.na(rowData(covid.sce)$ensembl_gene_id), ]

rownames(covid.sce) <- rowData(covid.sce)$ensembl_gene_id

message(paste0("Subsetting to ", length(intersect(rownames(covid.sce), covid.genes)), " common features"))
covid.sce <- covid.sce[intersect(rownames(covid.sce), covid.genes), ]

message("Extracting healthy control cells")
healthy.meta <- covid.meta[covid.meta$Status %in% c("Healthy"), ]
rownames(healthy.meta) <-  healthy.meta$Sample

healthy.sce <- covid.sce[, colnames(covid.sce) %in% healthy.meta$Sample]
healthy.clusters <- as.character(healthy.meta[colnames(healthy.sce), ]$cell.type.coarse)

message(paste0("Finding marker genes for ", length(unique(healthy.clusters)), " cell type annotations"))
healthy.marker.genes <- findMarkers(healthy.sce, groups=healthy.clusters, pval.type="all")

message("Extracting top marker genes for all cell type annotations")
healthy.cell.types <- unique(healthy.clusters)
healthy.marker.gene.list <- list()

for(x in seq_along(healthy.cell.types)){
  x.cell <- healthy.cell.types[x]
  x.mark.genes <- healthy.marker.genes[[x.cell]]
  # select genes at 10% FDR
  x.top.genes <- x.mark.genes[x.mark.genes$FDR <= 1e-1, ]
  # save the top 100 for each cell type
  rank.genes <- x.top.genes[order(x.top.genes$FDR, decreasing=FALSE), ]
  healthy.marker.gene.list[[x.cell]] <- rownames(rank.genes)[1:100]
}

healthy.features <- unique(unlist(healthy.marker.gene.list))
healthy.features <- healthy.features[!is.na(healthy.features)]

message("Training healthy control random forest")
healthy.train <- randomForest(x=as.matrix(t(logcounts(healthy.sce[healthy.features, ]))),
                          y=as.factor(healthy.meta[colnames(healthy.sce), ]$cell.type.coarse),
                          ntree=500, importance=TRUE)

healthy.ofile <- "/mnt/scratchb/jmlab/morgan02/Covid/classifiers/Blish-Healthy_PBMC_RF_classifier.RDS"
message(paste0("Healthy control random forest saved to: ", healthy.ofile))
saveRDS(healthy.train, file=healthy.ofile)

message("Subsetting patient cells")
nonhealthy.meta <- covid.meta[!covid.meta$Status %in% c("Healthy"), ]
rownames(nonhealthy.meta) <-  nonhealthy.meta$Sample

nonhealthy.sce <- covid.sce[, colnames(covid.sce) %in% nonhealthy.meta$Sample]
nonhealthy.clusters <- as.character(nonhealthy.meta[colnames(nonhealthy.sce), ]$cell.type.coarse)

nonhealthy.marker.genes <- findMarkers(nonhealthy.sce, groups=nonhealthy.clusters, pval.type="all")

nonhealthy.cell.types <- unique(nonhealthy.clusters)
message(paste0("Finding marker genes for ", length(nonhealthy.cell.types), " cell type annotations"))

nonhealthy.marker.gene.list <- list()

for(x in seq_along(nonhealthy.cell.types)){
  x.cell <- nonhealthy.cell.types[x]
  x.mark.genes <- nonhealthy.marker.genes[[x.cell]]
  # select genes at 10% FDR
  x.top.genes <- x.mark.genes[x.mark.genes$FDR <= 1e-1, ]
  # save the top 100 for each cell type
  rank.genes <- x.top.genes[order(x.top.genes$FDR, decreasing=FALSE), ]
  nonhealthy.marker.gene.list[[x.cell]] <- rownames(rank.genes)[1:100]
}

nonhealthy.features <- unique(unlist(nonhealthy.marker.gene.list))
nonhealthy.features <- nonhealthy.features[!is.na(nonhealthy.features)]

message("Training paitent random forest")
nonhealthy.train <- randomForest(x=as.matrix(t(logcounts(nonhealthy.sce[nonhealthy.features, ]))),
                          y=as.factor(nonhealthy.meta[colnames(nonhealthy.sce), ]$cell.type.coarse),
                          ntree=500, importance=TRUE)

covid.ofile <-  "/mnt/scratchb/jmlab/morgan02/Covid/classifiers/Blish-COVID19_PBMC_RF_classifier.RDS"
message(paste0("Saving patient random forest to: ", covid.ofile))
saveRDS(nonhealthy.train, file=covid.ofile)

message("Building classifier across all cells")
rownames(covid.meta) <-  covid.meta$Sample
allpbmc.clusters <- as.character(covid.meta[colnames(covid.sce), ]$cell.type.coarse)

message(paste0("Finding marker genes across ", length(unique(allpmbc.clusters)), " cell type annotations"))
allpbmc.marker.genes <- findMarkers(covid.sce, groups=allpbmc.clusters, pval.type="all")

allpbmc.cell.types <- unique(allpbmc.clusters)
allpbmc.marker.gene.list <- list()

for(x in seq_along(allpbmc.cell.types)){
  x.cell <- allpbmc.cell.types[x]
  x.mark.genes <- allpbmc.marker.genes[[x.cell]]
  # select genes at 10% FDR
  x.top.genes <- x.mark.genes[x.mark.genes$FDR <= 1e-2, ]
  # save the top 100 for each cell type
  rank.genes <- x.top.genes[order(x.top.genes$FDR, decreasing=FALSE), ]
  allpbmc.marker.gene.list[[x.cell]] <- rownames(rank.genes)[1:100]
}

allpbmc.features <- unique(unlist(allpbmc.marker.gene.list))
allpbmc.features <- allpbmc.features[!is.na(allpbmc.features)]

message("Training random forest on all cells")
allpbmc.train <- randomForest(x=as.matrix(t(logcounts(covid.sce[allpbmc.features, ]))),
                          y=as.factor(covid.meta[colnames(covid.sce), ]$cell.type.coarse),
                          ntree=500, importance=TRUE)

all.ofile <- "/mnt/scratchb/jmlab/morgan02/Covid/classifiers/Blish-ALL_PBMC_RF_classifier.RDS"
message(paste0("Saving random forest classifier to: ", all.ofile))
saveRDS(allpbmc.train, file=all.ofile)