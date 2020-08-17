#! /usr/bin/env Rscript

## Prepare gzipped log CPM+1 matrices of ADT expression, 1 per experimental batch

# ------- arg parsing ----------
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCElist"), type="character",
                     help="A set of comma-separated paths to the SCE objects, one per sample")

parser <- add_option(parser, c("-d", "--DonorList"), type="character",
                     help="Path to .txt containing Donor IDs of demultiplexed cells")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files denoised data and combined SCE object")

opt <- parse_args(parser)

library(SingleCellExperiment)
library(scran)
library(Matrix)
library(reshape2)

message("Reading in ADT SCE objects")
barcode.list <- unlist(strsplit(opt$SCElist, split=",", fixed=TRUE))
samp.names <- lapply(barcode.list, 
                     FUN=function(P) gsub(unlist(lapply(strsplit(P, fixed=TRUE, split="/"), 
                                                        FUN=function(sP) paste0(sP[length(sP)]))), 
                                          pattern="_cells\\.txt", replacement=""))
samp.names <- gsub(gsub(samp.names, pattern="Covid_SCE-", replacement=""),
	           pattern="_Ab\\.RDS", replacement="")

samp.names <- as.factor(unlist(samp.names))

sce.file.list <- unlist(strsplit(opt$SCElist, split=",", fixed=TRUE))
sce.list <- lapply(sce.file.list, FUN=function(FB) readRDS(FB))
names(sce.list) <- levels(samp.names)

message(paste0("Reading in Donor information from: ", opt$DonorList))
donor.files <- unlist(strsplit(opt$DonorList, split=",", fixed=TRUE))
donor.list <- lapply(donor.files, FUN=function(DF) read.table(DF, stringsAsFactors=FALSE, header=TRUE, sep="\t"))
names(donor.list) <- levels(samp.names)

message(paste0(length(levels(samp.names)), " samples found: ", paste(levels(samp.names), collapse=", ")))

# subset to just the good cells
message("Selecting QC-passed cells only")

for(x in seq_along(levels(samp.names))){
      x.samp <- levels(samp.names)[x]
      x.donor <- donor.list[[x.samp]]
      x.donor$CellID <- paste0(x.samp, "_", x.donor$cell)
      rownames(x.donor) <- x.donor$CellID
      # remove doublets and unassigned cells
      x.donor <- x.donor[!x.donor$donor_id %in% c("doublet", "unassigned"), ]
      x.keep <- intersect(x.donor$CellID, colnames(sce.list[[x.samp]]))
      
      donor.list[[x.samp]] <- x.donor[x.keep, ]
      sce.list[[x.samp]] <- sce.list[[x.samp]][, x.keep]

      message(paste0("Keeping ", length(x.keep), " good cells in sample: ", x.samp))
}

#######################################################################################
### -------- log CPM+1 normalise the ADT counts -------- ###
#######################################################################################

message("Performing log CPM transform on ADT libraries")
## Use a log CPM + 1 to rescale counts
cpm.list <- list()
for(x in seq_along(levels(samp.names))){
      samp.x <- levels(samp.names)[x]
      
      x.cpm <- apply(counts(sce.list[[samp.x]]), 2,
                     FUN=function(X) (X+1)/((sum(X)+1)/(1e6+1)))
      x.cpm <- as.data.frame(as(x.cpm, "matrix"))

      colnames(x.cpm) <- colnames(sce.list[[samp.x]])
      rownames(x.cpm) <- rowData(sce.list[[samp.x]])$Symbol
      x.cpm$ADT <- rownames(x.cpm)

      # one file per group
      x.ofile <- gzfile(paste0(opt$output, "_", samp.x, "_cpm.txt.gz"), "w")

      message(paste0("Writing log CPM file to: ", paste0(opt$output, "_", samp.x, "_cpm.txt.gz")))
      write.table(x.cpm, file=x.ofile, quote=FALSE, row.names=FALSE, sep="\t")
      close(x.ofile)
      
      #print(head(x.cpm[, c((ncol(x.cpm)-10):ncol(x.cpm))]))
      x.melt <- melt(x.cpm, id.vars="ADT")
      x.melt$group <- samp.x
      x.melt$view <- "ADT"
      cpm.list[[samp.x]] <- x.melt
}

sink(file="/dev/null")
rm(list=c("x.cpm"))
gc()
sink(file=NULL)

# make a list of matrices
# and a vector group inclusion
#adt.mat <- do.call(rbind, cpm.list)
#message(paste0("Concatenated ADT matrix contains: ", length(unique(adt.mat$variable)), " single-cells"))

#long.ofile <- gzfile(paste0(opt$output, "_longCPM.txt.gz"), "w")
#write.table(adt.mat, file=long.ofile, sep="\t", quote=FALSE, row.names=FALSE)
#close(long.ofile)