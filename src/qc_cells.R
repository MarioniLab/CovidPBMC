###################
# Basic QC for
# Covid19 PBMCs
###################

# ---- Parsing ----
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-x", "--matrix"), type="character", 
		     help="Path to the raw matrix.mtx.gz")
parser <- add_option(parser, c("-b", "--barcodes"), type="character", 
		     help="Path to .txt containing barcodes of called cells")
parser <- add_option(parser, c("-s", "--sparsitythreshold"), type="double", default=0.99, 
		     help="threshold to fraction of zeroes per cells [default %default]",
		     metavar="number")
parser <- add_option(parser, c("-u", "--umithreshold"), type="integer", default=1000, 
		     help="UMI threshold [default %default]",
		     metavar="number")
parser <- add_option(parser, c("-m", "--mtthreshold"), type="integer", default=2, 
		     help="MAD threshold for mitochondria [default %default]",
		     metavar="number")
parser <- add_option(parser, c("-l", "--logs"), type="character", 
		     help="Path to csv for summary (optional)")
parser <- add_option(parser, c("-p", "--plots"), type="character", 
		     help="Path to pdf for plots")
parser <- add_option(parser, c("-o", "--out"), type="character", 
		     help="Path to .txt to write QC pass cells")
opt <- parse_args(parser)

# ---- Load Data ----
library(DropletUtils)
library(scater)
library(ggplot2)
library(cowplot)
library(ggrepel)
theme_set(theme_cowplot())

# Read in called cell barcodes
cells <- read.csv(opt$barcodes,stringsAsFactors=FALSE,header=FALSE)[,1]
# Read in raw counts
fldr <- dirname(opt$matrix)
sce <- read10xCounts(fldr,col.names=TRUE)
# Subset to called cells
sce <- sce[,cells]

# ---- QC plots ----
# compute QC stats
qc <- data.frame(perCellQCMetrics(sce))

# Anything related to library size
gdthreshold <- (1-opt$sparsitythreshold)*nrow(sce)
qc$is.lib.fail <- qc$sum < opt$umithreshold | qc$detected < gdthreshold

lib.dist <- ggplot(qc, aes(x=sum)) +
    geom_histogram(bins=100) +
    scale_x_log10() +
    geom_vline(xintercept=opt$umithreshold) +
    xlab("Library Size") +
    ggtitle("Library Size")

gene.dist <- ggplot(qc, aes(x=detected)) +
    geom_histogram(bins=100) +
    scale_x_log10() +
    geom_vline(xintercept=gdthreshold) +
    xlab("# feautres") +
    ggtitle("Features Expressed")

det.trend <- ggplot(qc, aes(x=sum, y=detected, color=is.lib.fail)) +
    geom_point()  +
    scale_color_manual(values=c("black","grey80")) +
    xlab("Library Size") +
    ylab("# features") +
    ggtitle("Detection Tredn")

# Droplets with only GEX or Ab counts?
is.gex <- rowData(sce)$Type == "Gene Expression"
sum.expr <- colSums(counts(sce)[is.gex,])
sum.ab <- colSums(counts(sce)[!is.gex,])
fplot <- data.frame("GEXLib"=sum.expr,
		    "ABLib"=sum.ab,
		    "is.lib.fail"=qc$is.lib.fail)
gex.v.ab <- ggplot(fplot, aes(x=sum.expr, y=sum.ab, color=is.lib.fail)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values=c("black","grey80")) +
    xlab("GEX Sums") +
    ylab("Ab Sums") +
    coord_fixed(1) +
    ggtitle("Ab vs GEX counts")

# Mitochondiral counts
# Get genes on MT
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(sce), 
		   column="SEQNAME", keytype="GENEID")
is.mito <- location=="MT"
mitogenes <- names(location)[is.mito]
mitogenes <- mitogenes[!is.na(mitogenes)]


qc$mtsum <- colSums(counts(sce)[mitogenes,])
qc$mtrel <- qc$mtsum / qc$sum
qc$is.mt.fail <- isOutlier(qc$mtrel, type="higher", nmads=opt$mtthreshold)
relmtthreshold <- round(attr(qc$is.mt.fail,"thresholds")[2],3)

mt.dist <- ggplot(qc, aes(x=mtrel,y=sum, color=is.mt.fail)) +
    geom_point() +
    xlab("Fraction of MT counts") +
    ylab("Library Size") +
    scale_color_manual(values=c("black","grey80")) +
    ggtitle(paste0("Thresholded at ",relmtthreshold))

# ---- Define QC_Pass ----
keep <- !(qc$is.mt.fail | qc$is.lib.fail)

# ---- Difference between QC_Pass and QC_Fail ----

# MA-Plot between QC_Fail and QC_Pass
pass <- calculateAverage(counts(sce)[,!keep])
fail <- calculateAverage(counts(sce)[,keep])
library(edgeR)
logged <- cpm(cbind(pass, fail), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
names(logFC) <- uniquifyFeatureNames(rowData(sce)$ID,rowData(sce)$Symbol)
abundance <- rowMeans(logged)
fplot <- data.frame("abundance"=abundance,
		    "logFC"=logFC,
		    "Name"=names(logFC))
fplot <- dplyr::arrange(fplot,desc(abs(logFC)))
ma.plt <- ggplot(fplot, aes(x=abundance, y=logFC)) +
    geom_point() +
    geom_text_repel(data=fplot[abs(fplot$logFC)>1,], aes(label=Name)) +
    ggtitle("QC_Fail - QC_Pass") +
    ylab("logFC") +
    xlab("Average logExpr")

##################
# Classification of QC_Fail 
# tbc
##################


# ---- Summary ----
stat.names <- c("Total_Cells", "Lib_Fail", "Sparsity_Fail","Mt_Fail", "QC_Fail", "QC_Pass")
stat.values <- c(ncol(sce), 
		 sum(qc$sum < opt$umithreshold),
		 sum(qc$detected < gdthreshold),
		 sum(qc$is.mt.fail),
		 sum(keep),
		 sum(!keep))
qc.stats <- data.frame("Statistic"=stat.names,
		       "Value"=stat.values)


# ---- QC-Plot ----
ttle <- ggdraw() + 
  draw_label(paste0("Qualtiy Control of ", ncol(sce),
		    " cells, of which ", sum(keep),
		    " pass"), fontface='bold', x=0, hjust=0, size=20) 

plots <- plot_grid(lib.dist, gene.dist, 
	  det.trend, gex.v.ab,
	  mt.dist, ma.plt,ncol=2)
pout <- plot_grid(ttle, plots, ncol=1, rel_heights=c(0.1, 1))


# --- Save ----
# Barcodes
out <- data.frame("Barcode"=colnames(sce)[keep])
write.table(out, file=opt$out, row.names=FALSE,col.names=FALSE, sep=",")

#QC Plots
ggsave(filename=opt$plots,plot=pout,width=12,height=12)

#Statistics
if (!is.null(opt$logs)) {
    write.csv(x=qc.stats,file=opt$logs,row.names=FALSE)
}
