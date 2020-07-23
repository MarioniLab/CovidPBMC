#! /usr/bin/env Rscript

### Denoise ADT expression using 2 approaches:
### 1) Denoised PCA and regress out library sizes
### 2) John Tsangs background noise estimation and mixture model
# The order probably matters a lot here - so I'll remove the background signal first, 
# the de-sparse the resulting PCs

# ------- arg parsing ----------
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCElist"), type="character",
                     help="A set of comma-separated paths to the SCE objects, one per sample")

parser <- add_option(parser, c("-d", "--DonorList"), type="character",
                     help="Path to .txt containing Donor IDs of demultiplexed cells")

parser <- add_option(parser, c("-u", "--UMIthreshold"),  type="numeric", default=100,
                     help="A mimumum UMI threshold on which to select empty droplets for background estimation")

parser <- add_option(parser, c("-x", "--matrixlist"), type="character",
                     help="A set of comma-separated paths to the raw matrix.mtx.gz")

parser <- add_option(parser, c("-w", "--whitelist"), type="character",
                     help="A list of comma-separated paths to pre-QC barcode whitelists")

parser <- add_option(parser, c("-i", "--IgControls"), type="character",
                     help="A list of comma separated feature names representing IgG control features")

parser <- add_option(parser, c("-p", "--plots"), type="character",
                     help="Path to pdf for plots")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files denoised data and combined SCE object")

opt <- parse_args(parser)

library(igraph)
library(SingleCellExperiment)
library(scran)
library(igraph)
library(Matrix)
library(irlba)
library(umap)
library(dsb)
source("/mnt/scratchb/jmlab/morgan02/Covid/scripts/dsb_functions.R")

message(paste0("Readining in barcodes from: ", opt$whitelist))
barcode.list <- unlist(strsplit(opt$whitelist, split=",", fixed=TRUE))
samp.names <- lapply(barcode.list, 
                     FUN=function(P) gsub(unlist(lapply(strsplit(P, fixed=TRUE, split="/"), 
                                                        FUN=function(sP) paste0(sP[length(sP)]))), 
                                          pattern="_cells\\.txt", replacement=""))

samp.names <- as.factor(unlist(samp.names))

message("Reading in ADT SCE objects")
sce.file.list <- unlist(strsplit(opt$SCElist, split=",", fixed=TRUE))
sce.list <- lapply(sce.file.list, FUN=function(FB) readRDS(FB))
names(sce.list) <- levels(samp.names)

message(paste0("Readining in Donor information from: ", opt$DonorList))
donor.files <- unlist(strsplit(opt$DonorList, split=",", fixed=TRUE))
donor.list <- lapply(donor.files, FUN=function(DF) read.table(DF, stringsAsFactors=FALSE, header=TRUE, sep="\t"))
names(donor.list) <- levels(samp.names)

#######################################################################################
### -------- Estimate DSB (Denoised scaled by background) for each protein -------- ###
#######################################################################################
message("Estimating background staining from empty droplets")
# need both called barcodes (pre-QC) and all barcodes to select empty droplets with non-zero counts
# probs need empty droplets with ~100+ UMIs
# what effect does varying the background UMI threshold have?

message(paste0("Readining in non-empty droplets from: ", opt$whitelist))
white.list <- lapply(barcode.list, FUN=function(FB) read.table(FB, stringsAsFactors=FALSE, header=FALSE)[,1])
names(white.list) <- levels(samp.names)

message(paste0("Readining in counts matrices from: ", opt$matrixlist))
fldr.list <- unlist(lapply(strsplit(opt$matrixlist, split=",", fixed=TRUE), dirname))
matrix.list <- lapply(fldr.list, FUN=function(FM) readMM(paste0(FM, "/matrix.mtx.gz")))
names(matrix.list) <- samp.names

bcs.list <- lapply(fldr.list, FUN=function(FM) read.table(paste0(FM, "/barcodes.tsv.gz"), sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1])
names(bcs.list) <- levels(samp.names)

feature.list <- lapply(fldr.list, FUN=function(FM) read.table(paste0(FM, "/features.tsv.gz"), sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1])
names(feature.list) <- levels(samp.names)

message(paste0("Selecting empty droplets with sufficient UMIs: ", opt$UMIthreshold))
# select empty droplets for each sample
empty.list <- lapply(levels(samp.names), FUN=function(MX) colSums(matrix.list[[MX]]) > opt$UMIthreshold)
names(empty.list) <- levels(samp.names)

# remove non-empty droplets
message("Collating empty droplets")
empty.drops.list <- list()
for(x in seq_along(levels(samp.names))){
  x.samp <- levels(samp.names)[x]
  empty.mtx <- matrix.list[[x.samp]][, empty.list[[x.samp]]]
  empty.bcs <- bcs.list[[x.samp]][empty.list[[x.samp]]]
  print(dim(empty.mtx))
  # plot the distribution of all ADT library sizes
  adt.libsize <- log10(colSums(empty.mtx)+1)
  pdf(paste0(opt$plots, x.samp, "_ADT-LibSizes.pdf"), height=4.95, width=8.95, useDingbats=FALSE)
  hist(adt.libsize, xlab="Droplets", ylab="Counts")   
  dev.off()

  colnames(empty.mtx) <- empty.bcs
  rownames(empty.mtx) <- feature.list[[x.samp]]
  
  # remove valid droplet barcodes and subset to antibodies
  keep.rows <- feature.list[[x.samp]][!grepl(feature.list[[x.samp]], pattern="(ENSG)|(Sars)")]
  empty.mtx <- empty.mtx[keep.rows, setdiff(colnames(empty.mtx), white.list[[x.samp]])]
  
  message(paste0(ncol(empty.mtx), " empty droplets found with #UMI > ", opt$UMIthreshold))
  empty.drops.list[[x.samp]] <- empty.mtx
}

## normalise with IgG to estimate background staining
ig.controls <- unlist(strsplit(opt$IgControls, split=",", fixed=TRUE))

message(paste0("Performing DSB normalistion with controls: ", opt$IgControls))
out.sce.list <- list()
for(x in seq_along(levels(samp.names))){
  samp.x <- levels(samp.names)[x]
  
  donor.id <- donor.list[[samp.x]]
  print(head(paste(samp.x, donor.id$cell, sep="_")))
  donor.id$CellID <- paste(samp.x, donor.id$cell, sep="_")
  rownames(donor.id) <- donor.id$CellID
  remove.doublets <- donor.id$CellID[donor.id$donor_id %in% c("doublet", "unassigned")]
  keep.cells <- donor.id$CellID[!donor.id$donor_id %in% c("doublet", "unassigned")]
  
  message(paste0("Removing ", length(remove.doublets), " doublets and unassigned cells"))
  x.sce <- sce.list[[samp.x]][, intersect(colnames(sce.list[[samp.x]]), keep.cells)]
  print(dim(x.sce))
  
  # find the protein IDs that match the controls
  use.ig.controls <- rowData(x.sce)[rowData(x.sce)$Symbol %in% ig.controls, ]$ID
  message(paste0("DSB normalisation for sample: ", samp.x))
  print(head(empty.drops.list[[samp.x]][, 1:10]))
  # dsb.norm <- DSBNormalizeProtein(cell_protein_matrix=counts(x.sce),
  #                                 empty_drop_matrix=empty.drops.list[[samp.x]],
  #                                 pseudocount.use=1,
  #                                 use.isotype.control=TRUE,
  #                                 isotype.control.name.vec=use.ig.controls,
  #                                 define.pseudocount=TRUE)

  dsb.norm <- DSBNormalizeProtein_QC(cell_protein_matrix=counts(x.sce),
                                  empty_drop_matrix=empty.drops.list[[samp.x]],
                                  pseudocount.use=1,
				  mixture=FALSE,
                                  use.isotype.control=TRUE,
                                  isotype.control.name.vec=use.ig.controls,
                                  define.pseudocount=TRUE)

  
  # add this normalised matrix back to the SCE
  assay(x.sce, "DSB.counts") <- dsb.norm
  out.sce.list[[samp.x]] <- x.sce
}

message(paste0("Combining ADT SCE objects for ", sum(unlist(lapply(out.sce.list, ncol))), " cells"))
# combine all SCEs together
out.sce <- do.call(cbind, out.sce.list)

message("Cleaning up R environment")
# clean up the environment
sink(file="/dev/null")
rm(list=c("sce.list", "empty.drops.list", "matrix.list"))
gc()
sink(file=NULL)

message("Constructing ADT PCA after DSB normalisation")
abs.pca <- prcomp_irlba(t(assay(out.sce, "DSB.counts")), n=50)

message("Building SNN-graph with DSB normalised counts: d=30, k=21")
set.seed(42)
abs.snn <- buildSNNGraph(abs.pca$x[, c(1:30)], k=21,
                         d=NA, transposed=TRUE)
abs.walktrap <- cluster_walktrap(abs.snn, steps=4, membership=TRUE)
abs.cluster <- data.frame("CellID"=colnames(out.sce), "ADT.Cluster"=as.character(abs.walktrap$membership))
n.clusters <- length(unique(as.character(abs.walktrap$membership)))
abs.cluster$ADT.Cluster <- factor(abs.cluster$ADT.Cluster,
                                  levels=c(1:n.clusters),
                                  labels=c(1:n.clusters))

message("Making ADT force-directed layout from DSB SNN-graph")
abs.fr.layout <- layout_with_fr(abs.snn)

abs.fr.df <- as.data.frame(abs.fr.layout)
abs.fr.df$CellID <- colnames(out.sce)
abs.fr.merge <- as.data.frame(merge(abs.fr.df, abs.cluster, by='CellID'))

# get the edges of the layout
abs.edges <- get.data.frame(abs.snn)
abs.edges$from.Sample <- abs.fr.df$CellID[abs.edges$from]
abs.edges$to.Sample <- abs.fr.df$CellID[abs.edges$to]

# match edges to vertices and graph weights
abs.edges$from.x <- abs.fr.df$V1[match(abs.edges$from.Sample, abs.fr.df$Sample)]
abs.edges$from.y <- abs.fr.df$V2[match(abs.edges$from.Sample, abs.fr.df$Sample)]
abs.edges$to.x <- abs.fr.df$V1[match(abs.edges$to.Sample, abs.fr.df$Sample)]
abs.edges$to.y <- abs.fr.df$V2[match(abs.edges$to.Sample, abs.fr.df$Sample)]

message("Making DSB normalised UMAP")
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
abs.umap.map$CellID <- colnames(out.sce)

message("Adding projections and cluster to SCE object")
# add all projections, reduced dimensions and cluster info to the SCE
# also output as plain text meta-data
reducedDim(out.sce, "DSB.PCA") <- abs.pca$x
reducedDim(out.sce, "DSB.UMAP") <- abs.umap.map[, c(1,2)]
reducedDim(out.sce, "DSB.SNNgraph") <- abs.fr.df[, c(1,2)]
colData(out.sce)$DSB.Cluster <- abs.cluster$ADT.Cluster

sce.ofile <- paste0(opt$output, "_DSB_SCE.RDS")
message(paste0("SCE object written to: ", sce.ofile))
saveRDS(out.sce,
        file=sce.ofile)

# one big meta-data file
pca.df <- as.data.frame(abs.pca$x)
pca.df$CellID <- colnames(out.sce)

pca.merge <- merge(pca.df, abs.fr.merge, by='CellID')
pca.umap.merge <- merge(pca.merge, abs.umap.map, by='CellID')

out.meta.file <- paste0(opt$output, "_DSB_Projections.tsv")
message(paste0("Projections and meta-data written to: ", out.meta.file))
write.table(pca.umap.merge,
            file=out.meta.file, sep="\t", quote=FALSE, row.names=FALSE)


message("De-sparsing ADT data")
# estimate PCs on ADT data, remove technical noise PCs and re-compute reduced-rank ADT expression matrix
# do this with and without doublets
message("Estimating technical components of ADTs using trend var")
# Need more than 4 genes to estimate the technical variance
adt.trend <- modelGeneVar(out.sce)

trend.plot <- paste0(opt$plots, "ADT_TrendVar.pdf")
message(paste0("Mean-variance trend plot saved to : ", trend.plot))
pdf(trend.plot, height=2.95, width=3.95, useDingbats=FALSE)
hvadt <- adt.trend$FDR <= 0.1
plot(x=adt.trend[!hvadt, ]$mean, y=adt.trend[!hvadt, ]$total, xlab="Mean", ylab="Variance", pch=19)
points(x=adt.trend[hvadt, ]$mean, y=adt.trend[hvadt, ]$total, col='red', pch=19)
legend(4.5, 0.9, legend=c("Highly variable", "No Diff"),
       col=c("red", "black"), pch=19, box.lty=0)
dev.off()

message("Regressing out ADT library sizes from DSB PCs")
# desparse these PCs, then recompute the low-rank protein expression matrix.
# try regressing out the library size from each PC.
# Do I need to centre these at 0?
n.pcs <- ncol(reducedDim(out.sce, "DSB.PCA"))
resid.pcs <- matrix(ncol=n.pcs, nrow=ncol(out.sce))
colsums <- log10(colSums(counts(out.sce)))

for(x in seq_along(1:n.pcs)){
  x.pc <- reducedDim(out.sce, "DSB.PCA")[, x]
  x.df <- data.frame("PC"=x.pc, "SF"=colsums)
  x.fit <- lm(PC ~ SF, data=x.df)
  resid.pcs[, x] <- residuals(x.fit)
}

rownames(resid.pcs) <- colnames(out.sce)
colnames(resid.pcs) <- paste0("PC", 1:n.pcs)

message("Constructing de-sparsed low-rank approximation of protein expression")
abs.lowrank <- t(resid.pcs %*% t(abs.pca$rotation))
assay(out.sce, "lowrank") <- abs.lowrank

message("Creating desparsed UMAP")
desparse.umap.map <- as.data.frame(umap(resid.pcs,
                                        n_neighbors=21,
                                        init='spectral',
                                        metric='cosine',
                                        n_components=2,
                                        min_dist=0.1,
                                        method='naive')$layout)

colnames(desparse.umap.map) <- paste0("Desparsed.ADT.UMAP", c(1:2))
desparse.umap.map$CellID <- colnames(out.sce)

message("Creating desparsed SNN-graph")
set.seed(42)
# k will need to scale with the total sample size I expect
abs.snn <- buildSNNGraph(resid.pcs, k=21,
                         d=NA, transposed=TRUE)
abs.walktrap <- cluster_walktrap(abs.snn, steps=4, membership=TRUE) # walktrap is slow on big data sets
abs.cluster <- data.frame("CellID"=colnames(out.sce), "ADT.Cluster"=as.character(abs.walktrap$membership))
n.clusters <- length(unique(as.character(abs.walktrap$membership)))
abs.cluster$ADT.Cluster <- factor(abs.cluster$ADT.Cluster,
                                  levels=c(1:n.clusters),
                                  labels=c(1:n.clusters))

message("Making desparsed ADT force-directed layout")
abs.fr.layout <- layout_with_fr(abs.snn)

abs.fr.df <- as.data.frame(abs.fr.layout)
abs.fr.df$CellID <- colnames(out.sce)
abs.fr.merge <- as.data.frame(merge(abs.fr.df, abs.cluster, by='CellID'))

# get the edges of the layout
abs.edges <- get.data.frame(abs.snn)
abs.edges$from.Sample <- abs.fr.df$CellID[abs.edges$from]
abs.edges$to.Sample <- abs.fr.df$CellID[abs.edges$to]

# match edges to vertices and graph weights
abs.edges$from.x <- abs.fr.df$V1[match(abs.edges$from.Sample, abs.fr.df$Sample)]
abs.edges$from.y <- abs.fr.df$V2[match(abs.edges$from.Sample, abs.fr.df$Sample)]
abs.edges$to.x <- abs.fr.df$V1[match(abs.edges$to.Sample, abs.fr.df$Sample)]
abs.edges$to.y <- abs.fr.df$V2[match(abs.edges$to.Sample, abs.fr.df$Sample)]

desparse.pca <- as.data.frame(resid.pcs)
desparse.pca$CellID <- colnames(out.sce)
projection.df <- merge(desparse.umap.map, abs.fr.merge, by='CellID')
projection.merge <- merge(projection.df, desparse.pca, by='CellID')

desparse.ofile <- paste0(opt$output, "_Desparsed_Projections.tsv")
message(paste0("Writing de-sparsed projections to: ", desparse.ofile))
write.table(projection.merge,
            desparse.ofile,
            quote=FALSE, row.names=FALSE, sep="\t")

# write out De-sparsed SCE object
reducedDim(out.sce, "Desparse.PCA") <- abs.pca$x
reducedDim(out.sce, "Desparse.SNNgraph") <- projection.merge[, c("V1", "V2")]
reducedDim(out.sce, "Desparse.UMAP") <- projection.merge[, c("Desparsed.ADT.UMAP1", "Desparsed.ADT.UMAP2")]

desparse.ofile <- paste0(opt$output, "_Desparsed_SCE.RDS")
message(paste0("Saving de-sparsed SCE object to: ", desparse.ofile))
saveRDS(out.sce, file=desparse.ofile)
 
