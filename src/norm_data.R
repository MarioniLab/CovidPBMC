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

parser <- add_option(parser, c("-u", "--ADTthreshold"), type="numeric", default=100,
       	  	     help="Threshold for the total number of ADT UMIs on which to remove sparse ADT libraries")

parser <- add_option(parser, c("-l", "--logs"), type="character",
                     help="Path to csv for summary (optional)")

parser <- add_option(parser, c("-p", "--plots"), type="character",
       	             help="Path to pdf for plots")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Path to output SingleCellExperiment object")

parser <- add_option(parser, c("-q", "--outputqc"), type="logical", action="store_true", default=FALSE,
                     help="Output normalised data for QC-failed cells as well")

parser <- add_option(parser, c("-f", "--sizefactors"), type="character",
       	  	     help="Output path for size factors")

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
library(BiocParallel)
theme_set(theme_cowplot())

## Register BiocParallel params
ncores <- 4
mcparam <- MulticoreParam(workers=ncores)
register(mcparam)

# read in cell barcodes, features and counts matrix
# read in barcode whitelist to remove for estimating size factors
# estimate size factors with and without QC-fail cell barcodes

# this should be a comma-separated list of matrices
barcode.list <- unlist(strsplit(opt$barcodeslist, split=",", fixed=TRUE))
samp.names <- lapply(barcode.list, FUN=function(P) gsub(unlist(lapply(strsplit(P, fixed=TRUE, split="/"), 
       	  		       		       				      FUN=function(sP) paste0(sP[length(sP)]))), pattern="_cells\\.txt", replacement=""))
samp.names <- as.factor(unlist(samp.names))
cells.list <- lapply(barcode.list, FUN=function(FB) read.table(FB, stringsAsFactors=FALSE, header=FALSE)[,1])

# to keep a consistent ordering, this needs to become a factor - it also needs to be learnt from the filenames
# rather than strictly relying on the input order
names(cells.list) <- samp.names
message(paste0("Found ", length(unlist(cells.list)), " total cell barcodes across ", length(samp.names), " samples"))

# Read in raw counts
message("Making single cell object")
fldr.list <- unlist(lapply(strsplit(opt$matrixlist, split=",", fixed=TRUE), dirname))
sce.list <- lapply(fldr.list, FUN=function(FM) read10xCounts(FM, col.names=TRUE))
names(sce.list) <- samp.names

# samples may have the same barcodes - need to add a sample-specific prefix
pass_bcs <- list()
for(x in seq_along(levels(samp.names))){
    samp.x <- levels(samp.names)[x]
    bcs.x <- paste(samp.x, cells.list[[samp.x]], sep="_")
    pass_bcs[[samp.x]] <- bcs.x
    message(paste0(length(bcs.x), " passed barcodes found for sample: ", samp.x))
}

# Subset to called cells
for(x in seq_along(levels(samp.names))){
    samp.x <- levels(samp.names)[x]
    sce.cols <- paste(samp.x, colnames(sce.list[[samp.x]]), sep="_")
    colnames(sce.list[[samp.x]]) <- sce.cols
}

message("Subsetting to common features")
# need to make sure the same genes are present in both SCE lists - what if the antibody data are different?
# This will subset to just the common GEX features
features_sce.list <- lapply(sce.list, FUN=function(Q) rowData(Q)$ID[rowData(Q)$Type == "Gene Expression"])
common.features <- Reduce(x=features_sce.list, f=function(x, y) intersect(x, y))
message(paste0(length(common.features), " common features found"))

message("Combining individual SCE objects together")
sce <- do.call(cbind, lapply(samp.names, FUN=function(SQ) sce.list[[SQ]][common.features, ]))

# get the list of antibody features in each data set
ab.sce.list <- list()
for(x in seq_along(levels(samp.names))){
  samp.x <- levels(samp.names)[x]
  x.ab.sce <- sce.list[[samp.x]][rowData(sce.list[[samp.x]])$Type == "Antibody Capture",]
  ab.sce.list[[samp.x]] <- x.ab.sce
}

message("Subsetting protein data")
names(ab.sce.list) <- levels(samp.names)
ab.sce.list <- lapply(levels(samp.names), FUN=function(SP) ab.sce.list[[SP]][, intersect(colnames(ab.sce.list[[SP]]), pass_bcs[[SP]])])
names(ab.sce.list) <- levels(samp.names)

sink(file="/dev/null")
rm(list=c("sce.list"))
gc(reset=TRUE)
sink(file=NULL)

message("Subsetting called cells")
sce <- sce[, unlist(pass_bcs)]
message(paste0(ncol(sce), " called cells in the SCE object"))

# read in barcode whitelist - a space-separated list of files
whitefile.list <- unlist(strsplit(opt$whitelists, split=",", fixed=TRUE))
white.list <- lapply(whitefile.list, FUN=function(FW) read.table(FW, sep="\t", stringsAsFactors=FALSE, header=FALSE)[, 1])
names(white.list) <- samp.names

white_bcs <- list()
for(x in seq_along(levels(samp.names))){
    samp.x <- levels(samp.names)[x]
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
    # I think this step might switch the names around - grr!
    message("Subsetting ADT SCE objects")
    lapply(ab.sce.list, FUN=function(SO) print(dim(SO)))
    for(x in seq_along(levels(samp.names))){
       samp.x <- levels(samp.names)[x]
       x.ab.sce <- ab.sce.list[[samp.x]]
       x.ab.sce <- x.ab.sce[, colnames(x.ab.sce) %in% white_bcs[[samp.x]]]
       ab.sce.list[[samp.x]] <- x.ab.sce
       print(dim(x.ab.sce))
    }
    sink(file="/dev/null")
    gc(reset=TRUE)
    sink(file=NULL)
} else{
    sce <- sce[, unlist(white_bcs)]
    ab.sce.list <- lapply(samp.names, FUN=function(SP) ab.sce.list[[SP]][, white_bcs[[SP]]])
    names(ab.sce.list) <- samp.names
}

# take an input gene X cell (barcode) dataframe of
# read counts/UMIs
# output the cell-specific size factor normalized log2 expression values

n.cells <- ncol(sce)
n.genes <- nrow(sce)
message(paste0("Computing gene expression sparsity over ", n.cells, " droplets."))
gene_sparsity <- apply(counts(sce) == 0, MARGIN = 1, sum)/n.cells
keep_genes <- gene_sparsity < opt$sparsity

message(paste0("Using ", sum(keep_genes), " genes for size factor estimation"))
sink(file="/dev/null")
rm(list=c("gene_sparsity"))
gc(reset=TRUE)
sink(file=NULL)

# set cluster size to 5% of data
cluster.size <- ceiling(ncol(sce) * 0.05)
message(paste0("Cluster size set to ", cluster.size))

if(ncol(sce) > 300){
  clusters <- quickCluster(sce, min.size=cluster.size,
  	      		   subset.row=keep_genes,
			   #BPPARAM=mcparam,
                           method="igraph")
  max.size <- floor(cluster.size/2)
} else if(ncol(sce) < 100) {
  clusters <- quickCluster(sce, min.size=cluster.size, 
  	      		   subset.row=keep_genes,
                           method="igraph")
  max.size <- floor(cluster.size/2) + 1
} else{
  clusters <- quickCluster(sce, min.size=cluster.size,
                           #BPPARAM=mcparam,
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
			 subset.row=keep_genes,
                         #BPPARAM=mcparam,
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

sink(file="/dev/null")
rm(list=c("clusters"))
gc(reset=TRUE)
sink(file=NULL)

message(paste0("Writing SCE object to: ", opt$output))
saveRDS(sce, file=opt$output)

message(paste0("Writing sum factors to: ", opt$sizefactors))
# output the sum factors in addition
s.factors <- data.frame("CellID"=colnames(sce), "SumFactor"=sizeFactors(sce)) 

# derive the output for size factors
write.table(s.factors, file=opt$sizefactors, quote=FALSE, sep="\t", row.names=FALSE)

message("Output antibody SCE for each sample")
ab.out.list <- list()
for(x in seq_along(levels(samp.names))){
    samp.x <- levels(samp.names)[x]
    ab.sce.x <- ab.sce.list[[samp.x]]
    gsub(opt$output, pattern="\\.RDS", replacement="_QCFail.RDS")
    out.ab.sce <- gsub(opt$output, pattern="\\.RDS", replacement=paste0("-", samp.x, "_Ab.RDS"))
    # add size factors and normalise
    x.cells <- colnames(ab.sce.x)

    message(paste0("Removing cells with #UMI < ", opt$ADTthreshold))
    adt.libsize <- colSums(counts(ab.sce.x))
    keep.cells <- adt.libsize >= opt$ADTthreshold
    message(paste0("Keeping ", sum(keep.cells), " ADT libraries from sample: ", samp.x))
    sizeFactors(ab.sce.x) <- sizeFactors(sce[, x.cells])
    ab.sce.x <- ab.sce.x[, keep.cells]

    ab.sce.x <- scater::normalize(ab.sce.x)
    saveRDS(ab.sce.x, file=out.ab.sce)
    ab.out.list[[samp.x]] <- ab.sce.x
}

big.abs.sce <- do.call(cbind, ab.out.list)
abs.ofile <- gsub(opt$output, pattern="\\.RDS", replacement="_Ab.RDS")
saveRDS(big.abs.sce, file=abs.ofile)

sink(file="/dev/null")
rm(list=c("big.abs.sce"))
gc(reset=TRUE)
sink(file=NULL)

good.cell.libs <- log10(colSums(counts(sce)))
good.cell.df <- data.frame("CellID"=colnames(sce), "LibSize"=good.cell.libs, "SumFactor"=sizeFactors(sce),
	     		   "QC.Status"="Pass", "Sample"=gsub(colnames(sce), pattern="([a-zA-Z0-9]+)_([ATCG]+)-1", replacement="\\1"))

ab.libs.list <- lapply(ab.sce.list, FUN=function(SP) data.frame("CellID"=colnames(SP), "Ab.LibSize"=log10(colSums(counts(SP)))))
ab.libs <- do.call(rbind.data.frame, ab.libs.list)
good.ab.cell.df <- merge(good.cell.df, ab.libs, by='CellID')

# purge the old sce object
sink(file="/dev/null")
rm(list=c("sce"))
gc(reset=TRUE)
sink(file=NULL)

# estimate size factor for QC failed cells - if flag present
if(opt$outputqc){
	message("Estimating size factors for QC-failed cells")
	n.cells <- ncol(sce.fail)
	n.genes <- nrow(sce.fail)

	# set cluster size to 10% of data
	# keep the same genes as previous
	cluster.size <- ceiling(ncol(sce.fail) * 0.1)

	message("Quick clustering on all cells")
	# for some reason there is an error here as the internal quick-clustering 
	# does it normalised values, so a quick and dirty SF estimate gives -ve SFs
	# This puts a bit of a spanner in the works.
	# Can use.ranks get around this?
	# I'll remove the smallest 5% of QC-failed cells
	tryCatch(
  	clusters <- quickCluster(sce.fail, min.size=cluster.size,
  	                             subset.row=keep_genes,
				     use.ranks=TRUE,
                                     #BPPARAM=mcparam,
				     method="igraph"),
        error = function(c){
	   message("Negative size factors in quick cluster, removing smallest 5% of QC-failed cells")
	   col.sums <- colSums(counts(sce.fail))
	   tp.cells <- floor(n.cells * 0.05)
	   message(paste0("Removing ", tp.cells, " cells"))
	   tp.thresh <- max(col.sums[order(col.sums, decreasing=FALSE)][1:tp.cells])

	   sce.fail <- sce.fail[, col.sums > tp.thresh]
	   message(paste0(ncol(sce.fail), " cells remaining"))
	   clusters <- quickCluster(sce.fail, min.size=cluster.size,
                                     subset.row=keep_genes,
                                     use.ranks=FALSE,
                                     #BPPARAM=mcparam,
                                     method="igraph")
           return(clusters)
	   }, finally = {
	      if(exists("clusters")){
	          pass
	      } else{
	      message("Creating quick clusters")
	      col.sums <- colSums(counts(sce.fail))
	      tp.cells <- floor(n.cells * 0.05)
	      tp.thresh <- max(col.sums[order(col.sums, decreasing=FALSE)][1:tp.cells])
	      sce.fail <- sce.fail[, col.sums > tp.thresh]

	      clusters <- quickCluster(sce.fail, min.size=cluster.size,
                                     subset.row=keep_genes,
                                     use.ranks=TRUE,
                                     #BPPARAM=mcparam,
                                     method="igraph")
	      }
	   }
	)
	message("Quick clustering done")
	max.size <- floor(cluster.size/2)

	message("Computing sum factors")
	# change the window size in 50% increments
	# need to force positive size-factors for the cells that fail QC
	size.inc <- ceiling(max.size * 0.5)
	sce.fail <- computeSumFactors(sce.fail,
                         max.cluster.size=max.size,
                         positive=TRUE,
			 subset.row=keep_genes,
                         #BPPARAM=mcparam,
                         assay.type='counts', clusters=clusters)
  
	# the scater function is normalise, not normalize
	# I appreciate the British spelling, but this might cause a clash
	# with igraph for those not familiar
	message("Normalising QC-failed cells")
	sce.fail <- normalize(sce.fail)
	fail.out.sce <- gsub(opt$output, pattern=".RDS", replacement="_QCFail.RDS")
	message(paste0("Writing QC-failed SCE object to: ", fail.out.sce))  
	saveRDS(sce.fail, file=fail.out.sce)

	# output the sum factors in addition
	s.factors <- data.frame("CellID"=colnames(sce.fail), "SumFactor"=sizeFactors(sce.fail))
	out.fail.factors <- gsub(opt$sizefactors, pattern="\\.tsv", replacement="_QCFail.tsv")
	message(paste0("Writing QC-failed size factors to: ", out.fail.factors))
	write.table(s.factors, file=out.fail.factors, quote=FALSE, sep="\t", row.names=FALSE)

	all.cell.libs <- log10(colSums(counts(sce.fail)))
	all.cell.df <- data.frame("CellID"=colnames(sce.fail), "LibSize"=all.cell.libs, "SumFactor"=sizeFactors(sce.fail),
		       		  "Sample"=gsub(colnames(sce.fail), pattern="([a-zA-Z0-9]+)_([ATCG]+)-1", replacement="\\1")) 

	message("Collating Ab libraries on QC-failed cells")
	all.ab.libs.list <- lapply(failab.sce.list, FUN=function(SP) data.frame("CellID"=colnames(SP), "Ab.LibSize"=log10(colSums(counts(SP)))))
	all.ab.libs <- do.call(rbind.data.frame, all.ab.libs.list)
	all.ab.cell.df <- merge(all.cell.df, all.ab.libs, by='CellID')
        all.ab.cell.df$QC.Status <- "Fail"
	all.ab.cell.df$QC.Status[all.ab.cell.df$CellID %in% good.cell.df$CellID] <- "Pass"
	fail.ab.cell.df <- all.ab.cell.df[all.ab.cell.df$QC.Status %in% c("Fail"), ]
	out.cell.df <- do.call(rbind.data.frame, list("pass"=good.ab.cell.df, "fail"=fail.ab.cell.df))
} else{
    out.cell.df <- merge(good.cell.df, ab.libs, by='CellID')
    out.cell.df$QC.Status <- "Pass"
}

# what plots do we want to generate? Distributions of size factors with and without QC'd cells?
# distribution of gene sparsity for the same? Plot of sum factors vs. lib sizes for GEX and Abs
title <- ggdraw() +
            draw_label(paste0("Comparisons of library sizes across samples"),
	    		fontface="bold", x=0, hjust=0, size=20)

lib.plot <- ggplot(out.cell.df, aes(x=LibSize, y=SumFactor, colour=QC.Status)) +
	           geom_point() +
		   theme_bw() + facet_wrap(~Sample) +
		   labs(title="Library Size vs. Sum Factor")

ablib.plot <- ggplot(out.cell.df, aes(x=LibSize, y=Ab.LibSize, colour=QC.Status)) +
	             geom_point() +
		     theme_bw() + facet_wrap(~Sample) +
		     labs(title="Library Size vs. Ab Library Size")

sfactor.plot <- ggplot(out.cell.df, aes(x=Ab.LibSize, y=SumFactor, colour=QC.Status)) +
	               geom_point() +
		       theme_bw() + facet_wrap(~Sample) +
		       labs(title="Ab Library Size vs. Sum Factor")

out.plots <- plot_grid(lib.plot, ablib.plot, sfactor.plot,
	     	       ncol=2)
pout <- plot_grid(title, out.plots, ncol=1, rel_heights=c(0.1, 1))

ggsave(filename=opt$plots, plot=pout, width=12, height=12)