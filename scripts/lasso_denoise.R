#! /usr/bin/env Rscript

## Denoise ADT expression.
## This takes a 3-step approach to denoising:
## 1) estimate background and signal components of expression for each single-cell
## using a mixture of gaussians
## 2) Using these latent factors we then estimate feature weights for each component 
## and protein using a LASSO regularized linear regression.
## 3) Finally, the feature weights and components are recombined to construct a lower rank
## approximation of the ADT expression with the background signal components removed.

## The LASSO uses 10-fold cross-validation to find the optimal lambda for the regularized model,
## so it can be a little slow and computationally heavy.

## No empty droplets were harmed in the making of this script.
library(mclust)
library(glmnet)

estimate_components <- function(x, debug=FALSE, plots=NULL, samp=NULL){
  #' Use a 3-component mixture of gaussians to find the per-cell mean for each component
  #' @param x a matrix of features X cells
  
  require(mclust)
  message(paste0(samp, ": Fitting a 3-component mixture model")  )
  cellwise_background_mclust <- apply(x, 2, function(P) {
    g = mclust::Mclust(P, G=3, warn=FALSE , verbose=FALSE)
    return(g$parameters$mean)})
  
  # extract mean components and plot the distributions.
  means <- t(cellwise_background_mclust)
  colnames(means) <- paste0("Mean", colnames(means))

  if(isTRUE(debug)){
      require(ggplot2)
      require(reshape2)
      ggplot(melt(means), aes(x=value, fill=Var1)) +
        geom_histogram(bins=50) +
	facet_wrap(~Var1, ncol=1)

      ggsave(paste0(plots, "_", samp, "_means-histogram.pdf"), 
      	     height=6.95, width=7.95, useDingbats=FALSE)
  }
  
  return(means)
}

denoise_adt <- function(A=NULL, background.means=NULL, adt.lib=NULL,
                        alpha=1, keep.model=FALSE){
  #' Use a LASSO model to remove background signal from ADT counts
  #' Makes use of k-fold cross validation to select the optimal 
  #' lambda parameter from the most regularised model within 1 standard
  #' error of the minimum lambda.
  #' This implicity assumes that 3 components were used in the model.
  #' It returns either a matrix of denoised ADT expression, or if keep.model=TRUE
  #' a list with slots: 'model' containing the CV LASSO model for inspection and plotting
  #' and the denoised matrix in 'x'.
  require(glmnet)
  
  # include ADT library size in Components matrix
  G <- as.matrix(do.call(cbind.data.frame,
                         list("mean"=background.means[, c(1:3)],
                              "lib"=adt.lib)))
  if(class(A) != "matrix"){
    A <- as.matrix(A)
  }
  
  # setup a LASSO model to derive the feature weights
  # we can now use k-fold cross validation to select the optimal value of lambda
  w.cvfit <- cv.glmnet(y=A, x=G, family="mgaussian", alpha=1)

  # optimal lambda gives the most regularized model that is within 1 SE
  optimal.lambda <- w.cvfit$lambda.1se

  # the most regularized model generally removes the 1st background component
  coeff.mat.1se <- t(do.call(cbind, coef(w.cvfit, s="lambda.1se")))
  rownames(coeff.mat.1se) <- names(coef(w.cvfit, s="lambda.1se"))
  resid.mat.1se <- t(A - t(coeff.mat.1se %*% t(cbind(1, G))))
  
  # I'll reconstruct the expression using the 2 different approaches
  # the coefficient/feature weights incldues the model intercept
  # implicitly assume the number of components is fixed at 3 + lib size
  C.1se <- (coeff.mat.1se[, c(2, 4)] %*% t(G[ , c(1, 3)])) + resid.mat.1se
  rownames(C.1se) <- rownames(A)
  colnames(C.1se) <- colnames(A)
  
  if(isTRUE(keep.model)){
    return(list("model"=w.cvfit, "x"=C.1se))
  } else{
    return(C.1se)
  }
}

# ------- arg parsing ----------
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCElist"), type="character",
                     help="A set of comma-separated paths to the SCE objects, one per sample")

parser <- add_option(parser, c("-d", "--DonorList"), type="character",
                     help="Path to .txt containing Donor IDs of demultiplexed cells")

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

# subset to just the good cells
message("Selecting QC-passed cells only")

message(paste0(length(levels(samp.names)), " samples found: ", paste(levels(samp.names), collapse=", ")))
#######################################################################################
### -------- Estimate mixture compoents for each cell -------- ###
#######################################################################################

## Use a log CPM + 1 to rescale counts
cpm.list <- list()
for(x in seq_along(levels(samp.names))){
      samp.x <- levels(samp.names)[x]
      
      x.cpm <- apply(counts(sce.list[[samp.x]]), 2,
                     FUN=function(X) (X+1)/((sum(X)+1)/(1e6+1)))

      rownames(x.cpm) <- rownames(sce.list[[samp.x]])
      cpm.list[[samp.x]] <- x.cpm
}

sink(file="/dev/null")
rm(list=c("x.cpm"))
gc()
sink(file=NULL)

# For now, don't include the IgG controls in the background estimation
ig.controls <- rowData(sce.list[[1]])[grep(rowData(sce.list[[1]])$Symbol, pattern="Ctrl"), ]$ID
protein.controls <- setdiff(rowData(sce.list[[1]])$ID, ig.controls)

message("Estimating signal components in QC-passed droplets")
means.list <- list()
for(x in seq_along(levels(samp.names))){
      x.samp <- levels(samp.names)[x]
      x.cpm <- cpm.list[[x.samp]][protein.controls, ]
      x.means <- estimate_components(x=x.cpm, debug=TRUE, plots=opt$plots, samp=x.samp)
      means.list[[x.samp]] <- x.means
}

#######################################################################################
### -------- Remove background and sequencing depth from ADT expression -------- ###
#######################################################################################
message("Calculating ADT library sizes")
adt.lib.list <- lapply(sce.list, FUN=function(X) log10(colSums(counts(X))))
names(adt.lib.list) <- names(sce.list)

message("Removing background and sequencing coverage from ADT expression")
denoise.list <- list()

# I need to decide about including IgG controls.
for(x in seq_along(levels(samp.names))){
    samp.x <- levels(samp.names)[x]
    x.A <- cpm.list[[samp.x]]
    x.adt <- adt.lib.list[[samp.x]]
    x.means <- means.list[[samp.x]]

    x.denoise <- denoise_adt(A=x.A, background.means=x.means,
    	      	 	     adt.lib=x.adt, keep.model=FALSE, alpha=1)
    message(paste0("Storing denoised matrix: ", dim(x.denoise)))
    denoise.list[[samp.x]] <- x.denoise
}

message("Combining denoised matrices")
print(lapply(denoise.list, dim))
denoise.adt <- do.call(cbind, denoise.list)

message("Combining SCE objects")
all.sce <- do.call(cbind, sce.list)

message("Adding denoise expression to SCE object")
assay(all.sce, "denoised") <- denoise.adt[, colnames(all.sce)]

sce.ofile <- paste0(opt$output, "-LASSOdenoise.SCE")
message(paste0("Saving full SCE object to: ", sce.ofile))
saveRDS(all.sce, file=sce.ofile)