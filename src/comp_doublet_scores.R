###########
#Identifying doublets
###########

# There are two types of possible doublets/multiplets
# Same cell type (A)
# Mixed cell types (B)

# We have the following per-cell information that we can use to flag and remove doublets
# 1. Genotype (as identified by vireo)
# 2. Doublet score (as computed by scran)
# 3. Antibody staining

# 1. Is the most reliable one (assuming we trust vireo) and will identify type A and B
# 2. Less reliable, i.e. wouldn't trust on a per cell basis, better at identifying B than A
# 3. Unexpected combinations of antibodies.Only useful for identifying B. Problem is we need to define the number of UMIs per antibody that we assume represents a sginal.

# For now I will leave out 3, as this still depends on denoising the ADT library and might need some more fiddeling. Will have to see if this is something that we only use post-pipeline when annotating cell types.

### Strategy
# Per Sample (Sample = lane on the 10x):
# 1. Estimate per-cell doublet score using scran
# 2. Highly resolved clustering 
# 3. Cells belonging to clusters with either high average cluster score or high number of vireo doublets are labeled as doublets
# 4. Finally combine all the samples and cluster together. Again clusters that have high number of previously flagged clusters are most likely doublets.
# 5. Remove Cells that were flagged as clusters by either, vireo, Step 3 or Step 4.

# This should remove the majority of doublet cells whilst being relatively conservative. We might have to remove additional doublets cluster during the cell type annotation but this step should clean up the dataset considerably.
# It's probably easiest to parallelize 1-2 within snakemake and then run the rest in a separate script.

# ---- Parsing ----
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character", 
		     help="Path to the normalized SCE object")
parser <- add_option(parser, c("-d", "--donorID"), type="character", 
		     help="Path to donor_ids.tsv containing the donor assignments")
parser <- add_option(parser, c("-t", "--Sample"), type="character", 
		     help="Sample name")
parser <- add_option(parser, c("-o", "--out"), type="character", 
		     help="Path to .csv to write doublet data")

opt <- parse_args(parser)

# ---- Load Data ----
library(DropletUtils)
library(scater)
library(scran)
library(irlba)
library(BiocSingular)
library(BiocNeighbors)
library(umap)

sce <- readRDS(opt$SCE)

donor.id <- read.table(opt$donorID, stringsAsFactors=FALSE, header=TRUE)

# Subset to Sample of interest
sce.full <- sce[,grepl(opt$Sample,colnames(sce))]

# Subset donor.id to QC pass
donor.id <- donor.id[donor.id$cell %in% sce.full$Barcode,]

# ---- DoubletScores ----

# This is only calculated on non-doublets 
nondoubs <- donor.id$cell[donor.id$donor_id!="doublet"]
nondoubs <- paste0(opt$Sample,"_",nondoubs)

sce <- sce[,nondoubs]

# Highly variable genes
dec.var <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec.var, prop=0.5)  # this is for the sake of speed

# Doublet Score
set.seed(42)
doublet.scores <- doubletCells(sce,BSPARAM=IrlbaParam(),
		     subset.row=hvgs)

# Quick SNN Graph 
set.seed(42)
igr <- buildSNNGraph(sce,
		     BSPARAM=IrlbaParam(),
		     assay.type="logcounts", 
		     subset.row=hvgs,
		     k=10,
		     BNPARAM=AnnoyParam()) 

# Clustering on graph
cl <- igraph::cluster_walktrap(igr, steps=3)
cluster <- cl$membership

# UMAP for visualisation in the next script
ump <- umap(as.matrix(t(logcounts(sce)[hvgs,])), random_state=42)
umap1 <- ump$layout[,1]
umap2 <- ump$layout[,2]

out <- data.frame("Barcode"=colnames(sce),
		  "Sample"=opt$Sample,
		  "Cluster"=cluster,
		  "UMAP1"=umap1,
		  "UMAP2"=umap2,
		  "DbltScore"=doublet.scores,
		  "UmiSums"=colSums(counts(sce)))

rownames(donor.id) <- paste0(opt$Sample,"_",donor.id$cell)
out$Donor <- donor.id[out$Barcode,"donor_id"]

# Save Data
write.csv(out,file=opt$out)
