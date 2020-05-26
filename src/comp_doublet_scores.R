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
# 3. Unexpected combinations of antibodies. Only useful for identifying B. 

# For now I will leave out 3. Will have to see if this is something that we only use post-pipeline when annotating cell types.

### Strategy
# Per Sample (Sample = lane on the 10x):
# 1. Estimate per-cell doublet score using scran (excluding the vireo doublets)
# 2. Highly resolved clustering (tried walktrap with small step size, two-step louvain seems better)
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
sce <- sce[,grepl(opt$Sample,colnames(sce))]

# Subset donor.id to QC pass
donor.id <- donor.id[donor.id$cell %in% sce$Barcode,]

# ---- DoubletScores ----

# This is only calculated on non-doublets 
nondoubs <- donor.id$cell[donor.id$donor_id!="doublet"]
nondoubs <- paste0(opt$Sample,"_",nondoubs)

sce.nodoub <- sce[,nondoubs]

# Highly variable genes
dec.var <- modelGeneVar(sce.nodoub)
hvgs <- getTopHVGs(dec.var)

# Doublet Score
set.seed(42)
doublet.scores <- doubletCells(sce.nodoub,BSPARAM=IrlbaParam(),
		     subset.row=hvgs)

# ---- Per-Sample-Clustering ----

# Quick SNN Graph note this is done including vireo doublets

dec.full <- modelGeneVar(sce)
hvgs.full <- getTopHVGs(dec.full)
set.seed(42)

# PCA
sce <- denoisePCA(sce, technical=dec.full,
		       subset.row=hvgs,
		       BSPARAM=BiocSingular::IrlbaParam())
# Graph with low k
igr <- buildSNNGraph(sce,
		     use.dimred="PCA",
		     k=10,
		     BNPARAM=AnnoyParam()) 

# Clustering on graph
cl <- igraph::cluster_louvain(igr)
sce$Cluster <- paste0("C",cl$membership)


# Sub Clustering
set.seed(42)
cluster.list <- quickSubCluster(sce, groups=sce$Cluster,
				  min.ncells=200,
    prepFUN=function(x) { 
        dec <- modelGeneVar(x)
        input <- denoisePCA(x, technical=dec,
            subset.row=getTopHVGs(dec),
            BSPARAM=BiocSingular::IrlbaParam())
    },
    clusterFUN=function(x) { 
        g <- buildSNNGraph(x, use.dimred="PCA", k=20,
			   BNPARAM=AnnoyParam())
        igraph::cluster_louvain(g)$membership
    }
)

# Extract barcode to cluster assignment
getSubClusters <- function(sce.obj){
    bcs <- colData(sce.obj)$Barcode
    clst <- colData(sce.obj)$subcluster
    out <- data.frame("Barcode"=bcs,
		      "Cluster"=clst)
    return(out)
}

clust.out <- lapply(cluster.list,getSubClusters)
out.cluster <- do.call(rbind,clust.out)



# UMAP for visualisation in the next script
pcs <- reducedDim(sce)
ump <- umap(pcs, random_state=42)

# Put all data in one DF
out.umap <- data.frame("UMAP1"=ump$layout[,1],
		       "UMAP2"=ump$layout[,2],
		       "Barcode"=sce$Barcode,
		       "UmiSums"=colSums(counts(sce)))

out.score <- data.frame("Barcode"=colData(sce.nodoub)$Barcode,
		  "DbltScore"=doublet.scores)

out <- dplyr::left_join(out.cluster,out.umap)
out <- dplyr::left_join(out,out.score)
out$Sample <- opt$Sample

# Add Donor information
rownames(donor.id) <- donor.id$cell
out$Donor <- donor.id[out$Barcode,"donor_id"]
out$Donor[is.na(out$Donor)] <- "doublet" # The doublets we removed earlier

out$Barcode <- rownames(out) <- paste0(out$Sample,"_",out$Barcode)

# Save Data
write.csv(out,file=opt$out)
