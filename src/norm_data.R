####################
# Deconvolution size-factor normalisation
# for single-cell Covid19 PBMCs
####################

# ------- arg parsing ----------
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-x", "--matrixlist"), type="character",
       	             help="A set of comma-separated paths to the raw matrix.mtx.gz")

parser <- add_option(parser, c("-b", "--barcodeslist"), type="character",
       	             help="Path to .txt containing barcodes of called cells")

parser <- add_option(parser, c("-w", "--whitelists"), type="character",
                     help="Path to .txt file containing QC'd cells")

parser <- add_option(parser, c("-s", "--sparsity"), type="double", default=0.95,
                     help="threshold to remove genes with many zero counts prior to normalisation [default %default]", metavar="number")

parser <- add_option(parser, c("-l", "--logs"), type="character",
                     help="Path to csv for summary (optional)")

parser <- add_option(parser, c("-p", "--plots"), type="character",
       	             help="Path to pdf for plots")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Path to output SingleCellExperiment object")

parser <- add_option(parser, c("-q", "--outputqc"), type="logical", action="store_true",
                     help="Output normalised data for QC-failed cells as well")

parser <- add_option(parser, c("-f", "--sizefactors"), type="character",
       	  	     help="Output path for size factors")

parser <- add_option(parser, c("-n", "--samplenames"), type="character",
                     help="A list of sample names to prefix to cell barcodes - in the same order as the input data")

parser <- add_option(parser, c("-g", "--force"), type="logical", action="store_true", default=FALSE,
                     help="Force size factors to be positive")

opt <- parse_args(parser)

# ------ Load libraries ---------
library(SingleCellExperiment)
library(scran)
library(scater)
library(ggplot2)
library(cowplot)
library(DropletUtils)
theme_set(theme_cowplot())

# read in cell barcodes, features and counts matrix
# read in barcode whitelist to remove for estimating size factors
# estimate size factors with and without QC-fail cell barcodes

# this should be a comma-separated list of matrices
barcode.list <- unlist(strsplit(opt$barcodeslist, split=",", fixed=TRUE))
cells.list <- lapply(barcode.list, FUN=function(FB) read.table(FB, stringsAsFactors=FALSE, header=FALSE)[,1])
samp.names <- unlist(strsplit(opt$samplenames, split=",", fixed=TRUE))
names(cells.list) <- samp.names
message(paste0("Found ", length(unlist(cells.list)), " total cell barcodes across ", length(samp.names), " samples"))

# Read in raw counts
message("Making single cell object")
fldr.list <- unlist(lapply(strsplit(opt$matrixlist, split=",", fixed=TRUE), dirname))
sce.list <- lapply(fldr.list, FUN=function(FM) read10xCounts(FM, col.names=TRUE))
names(sce.list) <- samp.names

# samples may have the same barcodes - need to add a sample-specific prefix
pass_bcs <- list()
for(x in seq_along(samp.names)){
    samp.x <- samp.names[x]
    bcs.x <- paste(samp.x, cells.list[[samp.x]], sep="_")
    pass_bcs[[samp.x]] <- bcs.x
}

# Subset to called cells
for(x in seq_along(samp.names)){
    samp.x <- samp.names[x]
    sce.cols <- paste(samp.x, colnames(sce.list[[samp.x]]), sep="_")
    colnames(sce.list[[samp.x]]) <- sce.cols
}

message("Subsetting to common features")
# need to make sure the same genes are present in both SCE lists - what if the antibody data are different?
features_sce.list <- lapply(sce.list, FUN=function(Q) rowData(Q)$ID)
common.features <- Reduce(x=features_sce.list, f=function(x, y) intersect(x, y))

message(paste0(length(common.features), " common features found"))

message("Combining individual SCE objects together")
sce <- do.call(cbind, lapply(sce.list, FUN=function(SQ) SQ[common.features, ]))

# get the list of antibody features in each data set
ab.sce.list <- lapply(sce.list, FUN=function(Q) Q[rowData(Q)$Type == "Antibody Capture", ])
ab.sce.list <- lapply(samp.names, FUN=function(SP) ab.sce.list[[SP]][, pass_bcs[[SP]]])

message("Subsetting called cells")
sce <- sce[, unlist(pass_bcs)]
message(paste0(ncol(sce), " called cells in the SCE object"))

# read in barcode whitelist - a space-separated list of files
whitefile.list <- unlist(strsplit(opt$whitelists, split=",", fixed=TRUE))
white.list <- lapply(whitefile.list, FUN=function(FW) read.table(FW, sep="\t", stringsAsFactors=FALSE, header=FALSE)[, 1])
names(white.list) <- samp.names

white_bcs <- list()
for(x in seq_along(samp.names)){
    samp.x <- samp.names[x]
    white.x <- paste(samp.x, white.list[[samp.x]], sep="_")
    white_bcs[[samp.x]] <- white.x    
}

message(paste0("Found ", length(unlist(white_bcs)), " QC-passed cell barcodes"))
if(opt$outputqc){
    black.list <- setdiff(colnames(sce), unlist(white_bcs))
    message(paste0("Including ", length(black.list), " QC-failed cells for normalisation"))
    sce.fail <- sce
    sce <- sce[, unlist(white_bcs)]
    failab.sce.list <- ab.sce.list
    ab.sce.list <- lapply(samp.names, FUN=function(SP) ab.sce.list[[SP]][, white_bcs[[SP]]])
    print(lapply(ab.sce.list, FUN=dim))
} else{
    sce <- sce[, unlist(white_bcs)]
    ab.sce.list <- lapply(samp.names, FUN=function(SP) ab.sce.list[[SP]][, white_bcs[[SP]]])
}

# take an input gene X cell (barcode) dataframe of
# read counts/UMIs
# output the cell-specific size factor normalized log2 expression values

n.cells <- ncol(sce)
n.genes <- nrow(sce)

gene_sparsity <- apply(counts(sce) == 0, MARGIN = 1, sum)/n.cells
keep_genes <- gene_sparsity < opt$sparsity

message(paste0("Using ", sum(keep_genes), " genes for size factor estimation"))

# set cluster size to 5% of data
cluster.size <- ceiling(ncol(sce) * 0.05)
message(paste0("Cluster size set to ", cluster.size))

if(ncol(sce) > 300){
  clusters <- quickCluster(sce, min.size=cluster.size,
  	      		   subset.row=keep_genes,
                           method="igraph")
  max.size <- floor(cluster.size/2)
} else if(ncol(sce) < 100) {
  clusters <- quickCluster(sce, min.size=cluster.size, 
  	      		   subset.row=keep_genes,
                           method="igraph")
  max.size <- floor(cluster.size/2) + 1
} else{
  clusters <- quickCluster(sce, min.size=cluster.size,
  	      		   subset.row=keep_genes)
  max.size <- floor(cluster.size/2)
}

# change the window size in 50% increments
size.inc <- ceiling(max.size * 0.5)
message(paste0("Estimating size factors using ", size.inc, " cell increments"))
# how are there so many negative size factors?
sce <- computeSumFactors(sce,
                         max.cluster.size=max.size,
                         positive=opt$force,
                         assay.type='counts', clusters=clusters)
  
# the scater function is normalise, not normalize
# I appreciate the British spelling, but this might cause a clash
# with igraph for those not familiar
# count the number of negative size factors
neg.sf <- sum(sizeFactors(sce) < 0)
if(neg.sf > 0){
    message(paste0(neg.sf, " negative size factors estimated - consider a higher sparsity threshold"))
}

message("Normalising single-cell expression values")
sce <- normalize(sce)

message(paste0("Writing SCE object to: ", opt$output))
saveRDS(sce, file=opt$output)

message(paste0("Writing sum factors to: ", opt$sizefactors))
# output the sum factors in addition
s.factors <- data.frame("CellID"=colnames(sce), "SumFactor"=sizeFactors(sce)) 

# derive the output for size factors
write.table(s.factors, file=opt$sizefactors, quote=FALSE, sep="\t", row.names=FALSE)

message("Output antibody SCE for each sample")
for(x in seq_along(samp.names)){
    samp.x <- samp.names[x]
    ab.sce.x <- ab.sce.list[[samp.x]]
    gsub(opt$output, pattern=".RDS", replacement="_QCFail.RDS") 
    out.ab.sce <- gsub(opt$output, pattern=".RDS", replacement=paste0("-", samp.x, "_Ab.RDS"))
    saveRDS(ab.sce.x, file=out.ab.sce)
}

good.cell.libs <- log10(colSums(counts(sce)))
good.cell.df <- data.frame("CellID"=colnames(sce), "LibSize"=good.cell.libs, "SumFactor"=sizeFactors(sce),
	     		   "QC.Status"="Pass")
lapply(ab.sce.list, FUN=function(SP) print(SP))

ab.libs.list <- lapply(ab.sce.list, FUN=function(SP) data.frame("CellID"=colnames(SP), "Ab.LibSize"=colSums(counts(SP))))
good.ab.cell.df <- merge(good.cell.df, ab.libs, by='CellID')

# purge the old sce object
sink(file="/dev/null")
rm(list=c("sce"))
gc()
sink(file=NULL)

# estimate size factor for QC failed cells - if flag present
if(opt$outputqc){
	message("Estimating size factors for QC-failed cells")
	n.cells <- ncol(sce.fail)
	n.genes <- nrow(sce.fail)
  
	gene_sparsity <- (apply(counts(sce.fail) == 0, MARGIN = 1, sum)/n.cells)
	keep_genes <- gene_sparsity < opt$sparsity

	# set cluster size to 1% of data
	cluster.size <- ceiling(ncol(sce.fail) * 0.05)
  
	if(ncol(sce.fail)> 300){
  	    clusters <- quickCluster(sce.fail, min.size=cluster.size,
  	                             subset.row=keep_genes,
				     method="igraph")
	    max.size <- floor(cluster.size/2)
	} else if(ncol(sce.fail) < 100) {
	    clusters <- quickCluster(sce.fail, min.size=cluster.size, 
  	    	       		     subset.row=keep_genes,
                     		     method="igraph")
	    max.size <- floor(cluster.size/2) + 1
        } else{
	    clusters <- quickCluster(sce.fail, min.size=cluster.size,
  	           		     subset.row=keep_genes)
	    max.size <- floor(cluster.size/2)
	}

	# change the window size in 50% increments
	size.inc <- ceiling(max.size * 0.5)
	sce.fail <- computeSumFactors(sce.fail,
                         max.cluster.size=max.size,
                         positive=opt$force,
                         assay.type='counts', clusters=clusters)
  
	# the scater function is normalise, not normalize
	# I appreciate the British spelling, but this might cause a clash
	# with igraph for those not familiar
	sce.fail <- normalize(sce.fail)
	fail.out.sce <- gsub(opt$output, pattern=".RDS", replacement="_QCFail.RDS")
	message(paste0("Writing QC-failed SCE object to: ", fail.out.sce))  
	saveRDS(sce.fail, file=fail.out.sce)

	# output the sum factors in addition
	s.factors <- data.frame("CellID"=colnames(sce.fail), "SumFactor"=sizeFactors(sce.fail))
	out.fail.factors <- gsub(opt$sizefactors, pattern="\\.tsv", replacement="_QCFail.tsv")
	message(paste0("Writing QC-failed size factors to: ", out.fail.factors))  
	write.table(s.factors, file=out.fail.factors, quote=FALSE, sep="\t", row.names=FALSE)

	all.cell.libs <- log10(colSums(counts(sce)))
	all.cell.df <- data.frame("CellID"=colnames(sce.fail), "LibSize"=good.cell.libs, "SumFactor"=sizeFactors(sce.fail)) 

	all.ab.libs <- do.call(rbind.data.frame, lapply(failab.sce.list, FUN=function(SP) data.frame("CellID"=colnames(SP), "Ab.LibSize"=colSums(counts(SP)))))
	all.ab.cell.df <- merge(all.cell.df, all.ab.libs, by='CellID')
        all.ab.cell.df$QC <- "Fail"
	all.ab.cell.df$QC[all.ab.cell.df$CellID %in% good.cell.df$CellID] <- "Pass"
	fail.ab.cell.df <- all.ab.cell.df[all.ab.cell.df$QC %in% c("Fail"), ]
	out.cell.df <- do.call(rbind.data.frame, list("pass"=all.good.cell.df, "fail"=fail.ab.cell.df))
} else{
    out.cell.df <- all.good.cell.df
}

# what plots do we want to generate? Distributions of size factors with and without QC'd cells?
# distribution of gene sparsity for the same? Plot of sum factors vs. lib sizes for GEX and Abs
lib.plot <- ggplot(out.cell.df, aes(x=LibSize, y=SumFactor, colour=Pass)) +
	           geom_point() +
		   theme_bw()

ablib.plot <- ggplot(out.cell.df, aes(x=LibSize, y=Ab.LibSize, colour=Pass)) +
	             geom_point() +
		     theme_bw()

sfactor.plot <- ggplot(out.cell.df, aes(x=Ab.LibSize, y=SumFactor, colour=Pass)) +
	               geom_point() +
		       theme_bw()