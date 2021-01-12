###########
#Identifying doublets Part 2
###########

# ---- Parsing ----
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-d", "--doubletdata"), type="character", 
		     help="Comma separated list of doublet data files")
parser <- add_option(parser, c("-o", "--out"), type="character", 
		     help="Path to .csv")
parser <- add_option(parser, c("-s", "--SCE"), type="character", 
		     help="Path to SCE object")
parser <- add_option(parser, c("-p", "--plot"), type="character", 
		     help="Path to QC plots")
opt <- parse_args(parser)

# ---- Load Data ----
library(ggplot2)
library(cowplot)
library(scran)
library(BiocSingular)
theme_set(theme_cowplot())

dd.files <- unlist(strsplit(opt$doubletdata,","))
ddata.list <- lapply(dd.files, function(DF) read.csv(DF,stringsAsFactors=FALSE,row.names=1))
message(paste0("Reading in ", length(dd.files), " doublet data files"))

# Sample names for list for plotting
list.names <- unlist(lapply(strsplit(dd.files,"/"),function(x) x[length(x)]))
list.names <- gsub("_doublet_data.csv","",list.names)
names(ddata.list) <- list.names

# ---- ClusterScores ----

ddata.df <- data.frame(do.call(rbind,ddata.list))

message("Computing per-cluster and per-sample scores")
# Log-Transform scores
# I decided against setting known GT doublets to high scores
# 1. This information should be captured independently when testing for percentage of GT doublets
# 2. The downstream effect this has is very much dependent on the exact value that you choose to set this on
#ddata.df$DbltScore[is.na(ddata.df$DbltScore)] <- quantile(ddata.df$DbltScore,.99, na.rm=TRUE)
ddata.df$DbltScore <- log10(ddata.df$DbltScore+1)

# Compute median scores per cluster and sample
cluster.scores <- aggregate(ddata.df$DbltScore, list(ddata.df$Sample, ddata.df$Cluster), function(x) median(x,na.rm=TRUE))
colnames(cluster.scores) = c("Sample", "Cluster", "Medscore")

# Number of cells per cluster
cluster.scores$ncells <- sapply(1:nrow(cluster.scores), function(RW){
			       sum(ddata.df$Cluster == cluster.scores$Cluster[RW] & ddata.df$Sample == cluster.scores$Sample[RW])
})

# Fraction of cluster vs total cells from that sample
cluster.scores$fracCells <- sapply(1:nrow(cluster.scores), function(RW){
				  cluster.scores$ncells[RW]/sum(cluster.scores$ncells[cluster.scores$Sample == cluster.scores$Sample[RW]])
})

#Fraction of genotype doublets within each cluster
cluster.scores$fracVireoDoublets = sapply(1:nrow(cluster.scores), function(RW){
				  ndoublet <- sum(ddata.df$Cluster == cluster.scores$Cluster[RW] & ddata.df$Sample == cluster.scores$Sample[RW] & ddata.df$Donor=="doublet")
				  ndoublet/cluster.scores$ncells[RW]})

# ---- DoubletID ----

message("Determining doublets")
# Define doublets based on fraction and median doublet score
# Will have to see whether it is sensible to identify outlier per sample rather than across all samples
cluster.scores$ScoreOutlier <- scater::isOutlier(cluster.scores$Medscore,nmad=2.5,type="higher",log=FALSE)
cluster.scores$ScoreIsDoublet <- cluster.scores$ScoreOutlier & cluster.scores$fracCells <= 0.07


# Define doublets based on fraction of cells within a cluster that are genotype doublets
cluster.scores$VireoOutlier <- scater::isOutlier(cluster.scores$fracVireo,nmad=2.5,type="higher",log=FALSE,
						 batch=cluster.scores$Sample)
cluster.scores$VireoIsDoublet <- cluster.scores$VireoOutlier & cluster.scores$fracCells <= 0.07


# Plot with median cluster scores
p0 <- ggplot(cluster.scores, aes(x=Sample, y=Medscore, colour=factor(Cluster), shape=VireoIsDoublet)) +
  geom_point(data=cluster.scores[cluster.scores$ScoreIsDoublet,],size=3) +
  geom_jitter(data=cluster.scores[!cluster.scores$ScoreIsDoublet,], col="darkgrey", size=2, width=.1) +
  guides(color="none") +
  labs(x="Sample", y="Median doublet score") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# Plot with fraction of vireo doublets per cluster
p1 <- ggplot(cluster.scores, aes(x=Sample, y=fracVireoDoublets, colour=factor(Cluster), shape=ScoreIsDoublet)) +
  geom_point(data=cluster.scores[cluster.scores$VireoIsDoublet,],size=2) +
  geom_jitter(data=cluster.scores[!cluster.scores$VireoIsDoublet,], colour="darkgrey", size=2, width=.1) +
  guides(color="none") +
  labs(x="Sample", y="Fraction of Vireo doublets") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
	legend.position="bottom")


message("Plotting")
# Combine plot
p_scores <- plot_grid(p0,p1,ncol=1)


# Visualize in UMAPs
cluster.scores.dbls <- cluster.scores[cluster.scores$VireoIsDoublet | cluster.scores$ScoreIsDoublet,]
smps <- names(ddata.list)[names(ddata.list) %in% cluster.scores$Sample]

doublet_plots <- lapply(smps, function(nm) {
		   x <- data.frame(ddata.list[[nm]])
		   p <- ggplot(x, aes(x=UMAP1, y=UMAP2)) +
			       geom_point(color="grey80",size=0.1) +
			       geom_point(data=x[x$Cluster %in% cluster.scores.dbls[cluster.scores.dbls$Sample==nm,"Cluster"],], aes(color=factor(Cluster)),size=.1) +
			       theme(legend.position="none", panel.background = element_rect(colour = "black", linetype = 1, size = 0.5),
				     axis.line = element_blank(),
				     axis.ticks = element_blank(),
				     axis.text = element_blank(),
				     axis.title = element_blank(),) +
			       ggtitle(paste0(nm," (n=",nrow(x),", k=",length(unique(x$Cluster)),", n_db=",sum(cluster.scores.dbls[cluster.scores.dbls$Sample==nm,"ncells"]),")")) 
			   return(p)
})
names(doublet_plots) <- smps
p_umaps <- plot_grid(plotlist=doublet_plots)


# Finalize per-sample doublet definition

# Cells identified by vireo
vir.dblts <- ddata.df[ddata.df$Donor=="doublet","Barcode"]

# Cluster with high doublet score 
scr.dblts.cl <- paste0(cluster.scores[cluster.scores$ScoreIsDoublet,"Sample"],
		       cluster.scores[cluster.scores$ScoreIsDoublet,"Cluster"])
scr.dblts <- ddata.df[paste0(ddata.df$Sample,ddata.df$Cluster) %in% scr.dblts.cl,"Barcode"]

# Cluster with high fraction of vireo cluster
virCluster.dblts.cl <- paste0(cluster.scores[cluster.scores$VireoIsDoublet,"Sample"],
		       cluster.scores[cluster.scores$VireoIsDoublet,"Cluster"])
virCluster.dblts <- ddata.df[paste0(ddata.df$Sample,ddata.df$Cluster) %in% virCluster.dblts.cl,"Barcode"]
virCluster.dblts <- virCluster.dblts[!virCluster.dblts %in% vir.dblts] # RM the cells that were called by vireo as doublets


# Cells identifed by one or more of the criteria
mult.dblts <- c(vir.dblts,scr.dblts,virCluster.dblts)[duplicated(c(vir.dblts,scr.dblts,virCluster.dblts))]
vir.dblts <- vir.dblts[!vir.dblts %in% mult.dblts]
scr.dblts <- scr.dblts[!scr.dblts %in% mult.dblts]
virCluster.dblts <- virCluster.dblts[!virCluster.dblts %in% mult.dblts]

per.smp.dblts <- c(vir.dblts,scr.dblts,virCluster.dblts,mult.dblts)

# ---- Cluster across whole dataset ----
message("Clustering whole dataset")
library(BiocNeighbors)
sce <- readRDS(opt$SCE)
dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec)
set.seed(42)

# PCA
sce <- denoisePCA(sce, technical=dec,
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
    bcs <- colnames(sce.obj)
    clst <- as.character(colData(sce.obj)$subcluster)
    out <- data.frame("Barcode"=bcs,
		      "Cluster"=clst)
    return(out)
}

# Compute a final UMAP of entire datasets
pcs <- reducedDim(sce)
library(umap)
ump <- umap(pcs, random_state=42)
ump <- data.frame("UMAP1"=ump$layout[,1],
		  "UMAP2"=ump$layout[,2],
		  "Barcode"=rownames(pcs))


clust.full <- lapply(cluster.list,getSubClusters)
cluster.ids <- do.call(rbind,clust.full)
colnames(cluster.ids)[2] <- "FullDataCluster"

# Add cluster assignment to ddatadf
ddata.df <- dplyr::left_join(ddata.df,cluster.ids)
ddata.df$PerSampleDoublet <- ddata.df$Barcode %in% per.smp.dblts

library(dplyr)
smry <- group_by(ddata.df, FullDataCluster, PerSampleDoublet) %>%
    summarize(NType=n()) %>%
    mutate(Prop=NType/sum(NType))

doubDef <- smry$Prop[smry$PerSampleDoublet] > 0.6
dblClusters <- smry$FullDataCluster[smry$PerSampleDoublet][doubDef]

#Define cells that are defined by guilty-by-association
ga.dblts <- ddata.df$Barcode[ddata.df$FullDataCluster %in% dblClusters]
ga.dblts <- ga.dblts[!ga.dblts %in% per.smp.dblts]

# ---- Save ----
message("Preparing output")

# Combine everything for output

all.dblts <- c(per.smp.dblts,ga.dblts)
ddata.df.dblts <- ddata.df[ddata.df$Barcode %in% all.dblts,]

out <- data.frame("Barcode"=ddata.df.dblts$Barcode)
out$DoubletType <- plyr::mapvalues(all.dblts,out$Barcode,
				   c(rep("Vireo",length(vir.dblts)),
				     rep("Score",length(scr.dblts)),
				     rep("VirCluster",length(virCluster.dblts)),
				     rep("Multiple",length(mult.dblts)),
				     rep("GBA",length(ga.dblts))))
message("Saving Output")
write.csv(out,opt$out)
name.scores <- paste0(opt$pl,"/DoubletScores.pdf")
name.umaps <- paste0(opt$pl,"/DoubletUMAPs.pdf")

ggsave(p_scores,file=name.scores)

# Save UMAPS on multiple pages of a single PDF
ln <- length(doublet_plots)
brks <- floor(ln/4)
brks <- 1:brks * 4
if (brks[length(brks)]!=ln) {
    brks <- c(brks,ln)
}
pdf(name.umaps, onefile=TRUE, width=10,height=10)
for (i in seq_along(brks)) {
    bgn <- 4 * (i-1) +1
    plot(plot_grid(plotlist=doublet_plots[bgn:brks[i]]))
}
# Add in Overall UMAP
ump <- left_join(ump,out)
ump$IsDoublet <- ump$Barcode %in% all.dblts
p <- ggplot(ump, aes(x=UMAP1, y=UMAP2)) +
	   geom_point(color="grey80",size=0.1) +
	   geom_point(data=ump[ump$IsDoublet,], aes(color=DoubletType),size=.1) +
	   theme(legend.position="none", panel.background = element_rect(colour = "black", linetype = 1, size = 0.5),
		 axis.line = element_blank(),
		 axis.ticks = element_blank(),
		 axis.text = element_blank(),
		 axis.title = element_blank(),) +
	   ggtitle(paste0("Full Dataset of ", nrow(ump), " with ", sum(ump$IsDoublet)," doublets")) +
	   facet_wrap(~DoubletType)
print(p)
dev.off()
