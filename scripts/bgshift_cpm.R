#! /usr/bin/env Rscript

## Extract background ADT signal from empty droplets 
# using empty droplets from GEX libraries
# subtract background estimated from a 2-component mixture model

# ------- arg parsing ----------
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-x", "--matrixlist"), type="character",
       	             help="A set of comma-separated paths to the raw matrix.mtx.gz")

parser <- add_option(parser, c("-f", "--featurelist"), type="character",
                     help="A set of comma-separated paths to the feature information")

parser <- add_option(parser, c("-b", "--barcodeslist"), type="character",
       	             help="Path to .txt containing barcodes of called cells")

parser <- add_option(parser, c("-w", "--whitelists"), type="character",
                     help="Path to .txt file containing QC'd cells")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files denoised data and combined SCE object")

parser <- add_option(parser, c("-p", "--plots"), type="character",
                     help="Path to directory for plotting")

opt <- parse_args(parser)

library(Matrix)
library(mclust)
library(ggplot2)
library(ggsci)
library(ggthemes)
library(reshape2)

# read in cell barcodes, features and counts matrix
# read in barcode whitelist to exclude QC-passed cells

# this should be a comma-separated list of matrices
barcode.list <- unlist(strsplit(opt$barcodeslist, split=",", fixed=TRUE))
samp.names <- lapply(barcode.list, FUN=function(P) gsub(unlist(lapply(strsplit(P, fixed=TRUE, split="/"), 
       	  		       		       				      FUN=function(sP) paste0(sP[length(sP)-3]))), pattern="_cells\\.txt", replacement=""))
samp.names <- gsub(samp.names, pattern="_CITE", replacement="")
samp.names <- as.factor(unlist(samp.names))
print(samp.names)

all.barcodes.list <- lapply(barcode.list, FUN=function(FB) read.table(FB, stringsAsFactors=FALSE, header=FALSE)[,1])
# to keep a consistent ordering, this needs to become a factor - it also needs to be learnt from the filenames
# rather than strictly relying on the input order
names(all.barcodes.list) <- samp.names

# Read in raw counts
message("Reading in raw counts matrices")
fldr.list <- unlist(strsplit(opt$matrixlist, split=",", fixed=TRUE))
matrix.list <- lapply(fldr.list, FUN=function(FM) readMM(FM))
names(matrix.list) <- samp.names

# Read in feature info
message("Reading in feature information")
feature.list <- unlist(strsplit(opt$featurelist, split=",", fixed=TRUE))
all.feature.list <- lapply(feature.list, FUN=function(FX) read.table(FX, stringsAsFactors=FALSE, header=FALSE, sep="\t"))
names(all.feature.list) <- samp.names

for(x in seq_along(levels(samp.names))){
    samp.x <- levels(samp.names)[x]
    x.mat <- matrix.list[[samp.x]]
    colnames(x.mat) <- all.barcodes.list[[samp.x]]
    matrix.list[[samp.x]] <- x.mat
}

rm(list=c("x.mat"))
gc()

# read in whitelist barcodes
white.list <- unlist(strsplit(opt$whitelists, fixed=TRUE, split=","))
keep.cell.list <- lapply(white.list, FUN=function(FW) read.table(FW, stringsAsFactors=FALSE, header=FALSE)[, 1])
names(keep.cell.list) <- samp.names
print(lapply(keep.cell.list, length))
print(names(keep.cell.list))

message(paste0("Found ", length(unlist(keep.cell.list)), " good cell barcodes to exclude"))

non.zero.list <- list()
for(x in seq_along(levels(samp.names))){
    samp.x <- levels(samp.names)[x]

    x.non <- colSums(matrix.list[[samp.x]][all.feature.list[[samp.x]]$V3 == "Gene Expression", ]) > 10 &
       colSums(matrix.list[[samp.x]][all.feature.list[[samp.x]]$V3 == "Gene Expression", ]) < 100

    x.bg.bcs <- setdiff(all.barcodes.list[[samp.x]][x.non], keep.cell.list[[samp.x]])
    non.zero.list[[samp.x]] <- x.bg.bcs
    print(length(x.bg.bcs))
}

message(paste0("Extracted ", length(unlist(non.zero.list)), " background barcodes"))

# extract the ADT counts for these background barcodes and log CPM normalize
adtcpm.bg.list <- list()
for(x in seq_along(levels(samp.names))){
    samp.x <- levels(samp.names)[x]
    message(paste("Extracting ADT counts and performing log CPM normalisation for sample: ", samp.x))
    x.adt <- matrix.list[[samp.x]][all.feature.list[[samp.x]]$V3 == "Antibody Capture", non.zero.list[[samp.x]]]
    x.adt <- apply(x.adt, 2, FUN=function(X) log10((X+1)/((sum(X)+1)/(1e6+1))))
    adtcpm.bg.list[[samp.x]] <- x.adt

    print(dim(x.adt))

    pdf(paste0(opt$plots, samp.x, "-ADT_logCPM_hist-badcells.pdf"), height=3.95, width=5.15, useDingbats=FALSE)
    hist(apply(x.adt, 2, mean), 100, main="ADT logCPM distribution for empty/bad cells")
    dev.off()

    x.adt <- as.data.frame(x.adt)
    x.adt$ADT <- all.feature.list[[samp.x]]$V2[all.feature.list[[samp.x]]$V3 == "Antibody Capture"]

    # write out these matrices for later use
    x.bgofile <- gzfile(paste0(opt$output, "/EmptyDroplets_", samp.x, "_bgCPM.txt.gz"), "w")
    write.table(x.adt, file=x.bgofile, sep="\t", quote=FALSE, row.names=FALSE)
    close(x.bgofile)

    x.counts <- as.data.frame(as.matrix(matrix.list[[samp.x]][all.feature.list[[samp.x]]$V3 == "Antibody Capture", non.zero.list[[samp.x]]]))
    x.counts$ADT <- all.feature.list[[samp.x]]$V2[all.feature.list[[samp.x]]$V3 == "Antibody Capture"]

    x.countofile <- gzfile(paste0(opt$output, "/EmptyDroplets_", samp.x, "_counts.txt.gz"), "w")
    write.table(x.counts, file=x.countofile, sep="\t", quote=FALSE, row.names=FALSE)
    close(x.countofile)

    sink(file="/dev/null")
    rm(list=c("x.counts", "x.adt"))
    gc()
    sink(file=NULL)
    
}

# fit a mixture model to these backgrounds
#######################################
#### Fitting a GMM to each protein ####
#######################################
message("Fitting a 2-component gaussian mixture model to each protein")
gmm.list <- list()
for(x in seq_along(levels(samp.names))){
    x.samp <- levels(samp.names)[x]
    x.mclust <- apply(adtcpm.bg.list[[x.samp]], 1, function(P) {
    	     g = mclust::Mclust(P, G=2, warn=FALSE , verbose=FALSE)
	     return(g$parameters$mean)})

    x.means <- do.call(rbind.data.frame, list(t(x.mclust)))
    colnames(x.means) <- paste0("Mean", 1:2)
    gmm.list[[x.samp]] <- x.means
    print(dim(x.means))

    # why is this model so bad at finding the 3 components?!
    # do I need to fit a separate model to each protein, rather than each cell??
    x.plot <- ggplot(melt(x.means), aes(x=value, fill=variable)) +
       geom_histogram(bins=100) +
       scale_fill_npg() +
       facet_wrap(~variable, ncol=1)
    ggsave(x.plot, filename=paste0(opt$plots, "/", samp.x, "-ADT_logCPM_hist-badcells.pdf"),
           height=4.95, width=4.95, useDingbats=FALSE)

}

message("Extracting ADT libraries to keep")
## The white list gives us the cell barcodes to keep
keep.adtcpm.list <- list()
for(x in seq_along(levels(samp.names))){
    samp.x <- levels(samp.names)[x]
    print(length(keep.cell.list[[samp.x]]))

    x.keep.mtx <- matrix.list[[samp.x]][all.feature.list[[samp.x]]$V3 == "Antibody Capture", keep.cell.list[[samp.x]]]
    x.keep.mtx <- apply(x.keep.mtx, 2, FUN=function(X) log10((X+1)/((sum(X)+1)/(1e6+1))))
    keep.adtcpm.list[[samp.x]] <- x.keep.mtx
}

message(paste0("Retained ", sum(unlist(lapply(keep.adtcpm.list, ncol))), " barcodes"))
message("Removing background signal per-protein")

for(x in seq_along(levels(samp.names))){
    x.samp <- levels(samp.names)[x]
    x.bgshift <- apply(keep.adtcpm.list[[x.samp]], 2,
                       FUN=function(X) X - gmm.list[[x.samp]][, 1])
    x.bgshift <- as.data.frame(x.bgshift)
    colnames(x.bgshift) <- paste0(x.samp, "_", colnames(x.bgshift))
    print(dim(x.bgshift))

    n.prots <- nrow(x.bgshift)
    pdf(paste0(opt$plots, "/", x.samp, "_logCPM_BgShift-histogram.pdf"),
	height=18.95, width=12.95, useDingbats=FALSE)
    par(mfrow=c(28, 7), mai=c(0, 0, 0, 0))
    for(i in seq_along(1:n.prots)){
        hist(unlist(x.bgshift[x, ]), xlab="", ylab="", main="", breaks=100)
     }
    dev.off()

    x.bgshift$ADT <- all.feature.list[[x.samp]]$V2[all.feature.list[[x.samp]]$V3 == "Antibody Capture"]

    x.ofile <- gzfile(paste0(opt$output, "/Covid_ADT_", x.samp, "_bgCPM.txt.gz"), "w")
    print(dim(x.bgshift))

    write.table(x.bgshift, file=x.ofile, quote=FALSE, row.names=FALSE, sep="\t")
    close(x.ofile)
}
