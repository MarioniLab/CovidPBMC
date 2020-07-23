## Testing using GEX empty droplets to estimate background signal in CITE-seq data
library(Matrix)
library(mclust)
library(ggplot2)
library(ggsci)
library(ggthemes)
library(reshape2)
# read in the good barcodes and the empty droplet barcodes.
# Then select the non-zero, non-good cell barcodes to pull out the background signal from the ADT libraries.

keep.cells <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/QC/SIGAA1_barcode-whitelist.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")[, 1]
all.barcodes <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAA1_CITE/outs/raw_feature_bc_matrix/barcodes.tsv.gz",
                           sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1]
gex.features <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAA1_CITE/outs/raw_feature_bc_matrix/features.tsv.gz",
                           sep="\t", header=FALSE, stringsAsFactors=FALSE)
gex.matrix <- readMM("/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAA1_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz")
colnames(gex.matrix) <- all.barcodes

# anything with 10 <= UMIs <= 1000
non.zero <- colSums(gex.matrix[gex.features$V3 == "Gene Expression", ]) > 10 &
  colSums(gex.matrix[gex.features$V3 == "Gene Expression", ]) < 100
background.bcs <- setdiff(all.barcodes[non.zero], keep.cells)

# that gives a bucket load of bad droplets to hopefully estimate the background ADT
adt.matrix <- gex.matrix[gex.features$V3 == "Antibody Capture", background.bcs]

# log10 CPM + 1 this matrix
adt.bg.cpm <- apply(adt.matrix, 2, FUN=function(X) log10((X+1)/((sum(X)+1)/(1e6+1))))

# fit a mixture model to these backgrounds
# hist(log10(colSums(adt.matrix)+1), 100)
# fit a mixture model to these backgrounds
pdf("/mnt/scratchb/jmlab/morgan02/Covid/ADT_logCPM_hist-badcells.pdf", height=3.95, width=5.15, useDingbats=FALSE)
hist(apply(adt.bg.cpm, 2, mean), 100, main="ADT logCPM distribution for empty/bad cells")
dev.off()

####################################
#### Fitting a GMM to each cell ####
####################################
# cellwise_background_mclust <- apply(adt.bg.cpm, 2, function(P) {
#   g = mclust::Mclust(P, G=2, warn=FALSE , verbose=FALSE)
#   return(g$parameters$mean)})
# 
# # why does it return a list sometimes, and others a matrix?!
# means <- do.call(rbind.data.frame, cellwise_background_mclust)
# colnames(means) <- paste0("Mean", 1:2)
# 
# # why is this model so bad at finding the 3 components?!
# # do I need to fit a separate model to each protein, rather than each cell??
# ggplot(melt(means), aes(x=value, fill=variable)) +
#   geom_histogram(bins=100) +
#   scale_fill_npg() +
#   facet_wrap(~variable, ncol=1)

#######################################
#### Fitting a GMM to each protein ####
#######################################
# plot the distribution of each protein
n.prots <- nrow(adt.bg.cpm)
pdf("/mnt/scratchb/jmlab/morgan02/Covid/SIGAA1_logCPM_badcells-histogram.pdf",
    height=18.95, width=12.95, useDingbats=FALSE)
par(mfrow=c(28, 7), mai=c(0, 0, 0, 0))
for(x in seq_along(1:n.prots)){
  hist(adt.bg.cpm[x, ], xlab="", ylab="", main="", breaks=100)
}
dev.off()
dev.off()

proteinwise_background_mclust <- apply(adt.bg.cpm, 1, function(P) {
  g = mclust::Mclust(P, G=2, warn=FALSE , verbose=FALSE)
  return(g$parameters$mean)})

# why does it return a list sometimes, and others a matrix?!
means <- do.call(rbind.data.frame, list(t(proteinwise_background_mclust)))
colnames(means) <- paste0("Mean", 1:2)

# why is this model so bad at finding the 3 components?!
# do I need to fit a separate model to each protein, rather than each cell??
ggplot(melt(means), aes(x=value, fill=variable)) +
  geom_histogram(bins=100) +
  scale_fill_npg() +
  facet_wrap(~variable, ncol=1)

# for each protein we could subtract the value of the 1st component - but there still remains the same 
# contribution to the background signal...
keep.adt.matrix <- gex.matrix[gex.features$V3 == "Antibody Capture", setdiff(keep.cells, background.bcs)]

# log10 CPM + 1 this matrix
keep.adt.bg.cpm <- apply(keep.adt.matrix, 2, FUN=function(X) log10((X+1)/((sum(X)+1)/(1e6+1))))

# scale all of the logCPMs by this background distribution - should essentially be 0-centred now.
keep.adt.bg.shift <- apply(keep.adt.bg.cpm, 2,
                           FUN=function(X) X - means[, 1])

# what do the distributions look like now?
# plot the distribution of each protein
pdf("/mnt/scratchb/jmlab/morgan02/Covid/SIGAA1_logCPM_BgShift-histogram.pdf",
    height=18.95, width=12.95, useDingbats=FALSE)
par(mfrow=c(28, 7), mai=c(0, 0, 0, 0))
for(x in seq_along(1:n.prots)){
  hist(keep.adt.bg.shift[x, ], xlab="", ylab="", main="", breaks=100)
}
dev.off()
dev.off()

# can I use these shifted values as the input for MOFA??
keep.adt.bg.shift <- as.data.frame(keep.adt.bg.shift)
colnames(keep.adt.bg.shift) <- paste0("SIGAA1_", colnames(keep.adt.bg.shift))
keep.adt.bg.shift$ADT <- gex.features$V1[gex.features$V3 == "Antibody Capture"]

x.ofile <- gzfile("/mnt/scratchb/jmlab/morgan02/Covid/MOFA/Covid_ADT_SIGAA1_bgCPM.txt.gz", "w")
write.table(keep.adt.bg.shift, file=x.ofile, sep="\t", quote=FALSE, row.names=FALSE)
close(x.ofile)

############################
### Testing MOFA output ####
############################

# I've run MOFA using these inputs - let's see how it comes out - refer to the 
# MOFA_backgroundShift.Rmd for details - MOFA2 can't be installed on this Rstudio server











