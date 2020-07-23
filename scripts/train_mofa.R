#! /usr/bin/env Rscript

## Run MOFA on either a single sample or multiple samples, each of which is treated as a separate group

# ------- arg parsing ----------
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCElist"), type="character",
                     help="A set of comma-separated paths to the SCE objects, one per sample")

parser <- add_option(parser, c("-d", "--DonorList"), type="character",
                     help="Path to .txt containing Donor IDs of demultiplexed cells")

parser <- add_option(parser, c("-g", "--groups"), type="character",
       	  	     help="A comma-separated list of sample labels. If present then samples are run as groups")

parser <- add_option(parser, c("-f", "--factors"), type="numeric", default=20,
       	  	     help="The number of factors to learn")

parser <- add_option(parser, c("-r", "--runmode"), type="character", default="medium",
                     help="Mode to run MOFA in. Either fast, medium of slow")

parser <- add_option(parser, c("-p", "--plots"), type="character",
                     help="Path to pdf for plots")

parser <- add_option(parser, c("-m", "--MOFAoutput"), type="character",
       	  	     help="Path to HDF5 file to save trained MOFA model")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files denoised data and combined SCE object")

opt <- parse_args(parser)

library(SingleCellExperiment)
library(scran)
library(igraph)
library(Matrix)
library(irlba)
library(umap)
library(reticulate)
reticulate::use_python("~/Python-3.7.6/python", required=TRUE)
library(MOFA2)

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
      x.keep <- intersect(x.donor$CellID, colnames(sce.list[[x.samp]]))
      
      donor.list[[x.samp]] <- x.donor[x.keep, ]
      sce.list[[x.samp]] <- sce.list[[x.samp]][, x.keep]

      message(paste0("Keeping ", length(x.keep), " good cells in sample: ", x.samp))
}

#######################################################################################
### -------- Setup the MOFA model - individual samples are groups -------- ###
#######################################################################################

message("Performing log CPM transform on ADT libraries")
## Use a log CPM + 1 to rescale counts
cpm.list <- list()
for(x in seq_along(levels(samp.names))){
      samp.x <- levels(samp.names)[x]
      
      x.cpm <- apply(counts(sce.list[[samp.x]]), 2,
                     FUN=function(X) (X+1)/((sum(X)+1)/(1e6+1)))
      x.cpm <- as(x.cpm, "dgCMatrix")

      colnames(x.cpm) <- colnames(sce.list[[samp.x]])
      rownames(x.cpm) <- rownames(sce.list[[samp.x]])
      cpm.list[[samp.x]] <- x.cpm
}

sink(file="/dev/null")
rm(list=c("x.cpm"))
gc()
sink(file=NULL)

# make a list of matrices
# and a vector group inclusion
adt.mat.list <- list("ADT"=do.call(cbind, cpm.list))
message(paste0("Concatenated ADT matrix contains: ", ncol(adt.mat.list[[1]]), " single-cells"))

groups.vec <- unlist(lapply(levels(samp.names), FUN=function(X) rep(X, ncol(cpm.list[[X]]))))

# hopefully this works with one or many samples
mofa.obj <- create_mofa(adt.mat.list, groups=groups.vec)

# setup data options - scale and center views
data_opts <- get_default_data_options(mofa.obj)
data_opts$scale_views <- TRUE

# use default model options - stick with gaussian for the ADT and GEX views for now
# set the number of factors here
model_opts <- get_default_model_options(mofa.obj)
model_opts$num_factors <- opt$factors

# use default training options - this is just a test, so make it fast
train_opts <- get_default_training_options(mofa.obj)
train_opts$convergence_mode <- opt$runmode
message(paste0("Running MOFA in ", opt$runmode, " mode"))

message(paste0("Training MOFA with ", opt$factors, " factors"))
# now we can train the MOFA model
mofa.obj <- prepare_mofa(object=mofa.obj,
                         data_options=data_opts,
                         model_options=model_opts,
                         training_options=train_opts)

mofa.ofile <- paste0(opt$MOFAoutput, "_MOFA.hdf5")
mofa.trained <- run_mofa(mofa.obj, mofa.ofile)

message(paste0("Trained MOFA model saved to: ", mofa.ofile))

message("Extracting factors")
factor.list <- get_factors(mofa.trained, groups="all", factors="all")

message("Extracting feature weights")
weights.list <- get_weights(mofa.trained, views="all", factors="all")

message("Reconstructing low-rank approximation of ADT expression")
low.rank.list <- list()
for(x in seq_along(levels(samp.names))){
      x.samp <- levels(samp.names)[x]
      x.factor <- factor.list[[x.samp]]
      x.weight <- weights.list[[x.samp]]

      x.lowrank <- as.matrix(x.weight)  %*% as.matrix(t(x.factor))
      colnames(x.lowrank) <- colnames(adt.mat.list[[x.samp]])
      rownames(x.lowrank) <- rownames(adt.mat.list[[x.samp]])
      low.rank.list[[x.samp]] <- x.lowrank

      low.ofile <- paste0(opt$output, x.samp, "_lowrank.tsv")
      message(paste0("Saving low-rank matrix to: ", low.ofile))
      write.table(x.lowrank, file=low.ofile, quote=FALSE, sep="\t")
}

# make one big SCE with the lowrank as a new assay slot
big.sce <- do.call(cbind, sce.list)
big.low <- do.call(cbind, low.rank.list)

assay(big.sce, "lowrank") <- big.low[, colnames(big.sce)]

sce.ofile <- paste0(opt$output, "_Full_SCE.RDS")
message(paste0("Saving full ADT SCE object to: ", sce.ofile))
saveRDS(big.sce, file=sce.ofile)