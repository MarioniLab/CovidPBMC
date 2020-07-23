#! /usr/bin/env Rscript

### quick and dirty visualisation - extract proteins of interest to overlay on UMAPs
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
abs.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE_ABS.RDS")
print(head(rowData(all.sce)))
print(head(colData(all.sce)))

goi.proteins <- c("CD3", "CD19", "CD14", "CD16", "CD19", "CD20", "CD11b", "CD79b", "CD80", "CD56",
                  "CD86", "CD8", "CD4", "CD11c", "IgM", "TCR", "HLA-DR", "CD1c", "CD64", "CD32",
                  "CX3CR1", "IgD", "CD45RA", "CD45RO", "CD197", "CD25", "IgG1_Ctrl", "CD21",
		  "CD127", "CD27",
                  "IgG2a_Ctrl", "IgG2b_Ctrl", "IgG2b_RatCtrl")
goi.prots.id <- rowData(abs.sce)[rowData(abs.sce)$Symbol %in% unique(goi.proteins) ,]$ID

goi.prot.exprs <- as.data.frame(as.matrix(t(logcounts(abs.sce[goi.prots.id, ]))))
colnames(goi.prot.exprs) <- rowData(abs.sce)[rowData(abs.sce)$Symbol %in% goi.proteins ,]$Symbol
goi.prot.exprs$CellID <- colnames(abs.sce)

write.table(goi.prot.exprs,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/SIGAA2_SIGAA3_ProteinGOI.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

# include the SARS-CoV-2 genes
sars.genes <- rowData(all.sce)[grepl(rowData(all.sce)$ID, pattern="SarsCov2"), ]$Symbol

goi.genes <- c("CD3E", "CD3G", "CD14", "FCGR3A", "FCGR3B", "FCGR2A", "FCGR2B", "FCGR2C",
               "MS4A1", "LYZ", "IL6", "IL6ST", "TNF", "IL1B", "CD4", "CD8A", "CD8B",
               "CD80", "CD86", "PTPRC", "CCR7", "RORC", "TBX21", "FOXP3", sars.genes,
               "IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNG", "IFNK", "IFNL1", "IFNL2", "IFNL3",
               "IFNB1")
goi.gene.id <- rowData(all.sce)[rowData(all.sce)$Symbol %in% goi.genes, ]$ID

goi.gene.exprs <-  as.data.frame(as.matrix(t(logcounts(all.sce[goi.gene.id, ]))))
colnames(goi.gene.exprs) <- rowData(all.sce)[rowData(all.sce)$Symbol %in% goi.genes ,]$Symbol
goi.gene.exprs$CellID <- colnames(all.sce)

write.table(goi.gene.exprs,
            file="/mnt/scratchb/jmlab/morgan02/Covid/reports/SIGAA2_SIGAA3_GeneGOI.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
