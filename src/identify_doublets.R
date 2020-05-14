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
parser <- add_option(parser, c("-p", "--plot"), type="character", 
		     help="Path to QC plots")

# ---- Load Data ----
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

dd.files <- unlist(strsplit(opt$doubletdata,","))
ddata.list <- lapply(dd.files, function(DF) read.csv(DF,stringsAsFactors=FALSE,row.names=1))

# Sample names for list for plotting
list.names <- unlist(lapply(strsplit(dd.files,"/"),function(x) x[length(x)]))
list.names <- gsub("_doublet_data.csv","",list.names)
names(ddata.list) <- list.names

# ---- ClusterScores ----

ddata.df <- data.frame(do.call(rbind,ddata.list))

# Log-Transform scores
ddata.df$DbltScore <- log10(ddata.df$DbltScore+1)

# Compute median scores per cluster and sample
cluster.scores <- aggregate(ddata.df$DbltScore, list(ddata.df$Sample, ddata.df$Cluster), median)
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

# Define doublets based on fraction and median doublet score
# Will have to see whether it is sensible to identify outlier per sample rather than across all samples
cluster.scores$ScoreOutlier <- scater::isOutlier(cluster.scores$Medscore,nmad=1.5,type="higher",log=FALSE)
cluster.scores$ScoreIsDoublet <- cluster.scores$ScoreOutlier & cluster.scores$fracCells <= 0.05

# Plot with median cluster scores
p0 <- ggplot(cluster.scores, aes(x=Sample, y=Medscore, colour=factor(Cluster))) +
  geom_point(data=cluster.scores[cluster.scores$ScoreIsDoublet,],size=2) +
  geom_jitter(data=cluster.scores[!cluster.scores$ScoreIsDoublet,], col="darkgrey", size=.5, width=.1) +
  guides(color="none") +
  labs(x="Sample", y="Median doublet score") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# Define doublets based on fraction of cells within a cluster that are genotype doublets
# This should be done per sample as this might vary depending on the amount of multiplexing
cluster.scores$VireoOutlier <- scater::isOutlier(cluster.scores$fracVireoDoublets,nmad=1.5,type="higher",log=FALSE, batch=cluster.scores$Sample)
cluster.scores$VireoIsDoublet <- cluster.scores$VireoOutlier & cluster.scores$fracCells <= 0.05

# Plot with fraction of vireo doublets per cluster
p1 <- ggplot(cluster.scores, aes(x=Sample, y=fracVireoDoublets, colour=factor(Cluster), shape=ScoreIsDoublet)) +
  geom_point(data=cluster.scores[cluster.scores$VireoIsDoublet,],size=2) +
  geom_jitter(data=cluster.scores[!cluster.scores$VireoIsDoublet,], colour="darkgrey", size=2, width=.1) +
  guides(color="none") +
  labs(x="Sample", y="Fraction of Vireo doublets") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
	legend.position="bottom")


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


#########
# Here will go the integration across datasets with identifying doublets that are guilty by association
# with per-sample doublets
#########


#### All of the below is just to label how each doublet was identified to make some checks later, might rm all of the below eventually

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

# Combine everything for output
all.dblts <- c(vir.dblts,scr.dblts,virCluster.dblts,mult.dblts)

ddata.df.dblts <- ddata.df[ddata.df$Barcode %in% all.dblts,]

out <- data.frame("Barcode"=ddata.df.dblts$Barcode)
out$DoubletType <- plyr::mapvalues(all.dblts,out$Barcode,
				   c(rep("Vireo",length(vir.dblts)),
				     rep("Score",length(scr.dblts)),
				     rep("VirCluster",length(virCluster.dblts)),
				     rep("Multiple",length(mult.dblts))))

# ---- Save ----
write.csv(out,opt$out)
name.scores <- paste0(opt$pl,"/DoubletScores.pdf")
name.umaps <- paste0(opt$pl,"/DoubletUMAPs.pdf")
ggsave(p_scores,file=name.scores)
ggsave(p_umaps,file=name.umaps)
