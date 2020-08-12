###################
# Demultiplex QC for
# Covid19 PBMCs
###################

# ---- Parsing ----
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-b", "--barcodes"), type="character", 
		     help="Path to .txt containing barcodes of called cells")
parser <- add_option(parser, c("-w", "--whitelist"), type="character", 
		     help="Path to .txt containing barcodes of QC-pass cells",
		     metavar="number")
parser <- add_option(parser, c("-d", "--donorIDs"), type="character", 
		     help="Path to donor_ids.tsv containing vireo output with donor assignments",
		     metavar="number")
parser <- add_option(parser, c("-x", "--matrix"), type="character", default=2, 
		     help="Path to raw feature_bc cellranger matrix",
		     metavar="number")
parser <- add_option(parser, c("-p", "--plots"), type="character", 
		     help="Path to pdf for output")
opt <- parse_args(parser)

# ---- Load Data ----
library(DropletUtils)
library(scater)
library(scran)
library(ggrepel)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Read in called cell barcodes
# Cells
message(paste0("Reading in barcode information: ", opt$barcodes))
cells <- read.csv(opt$barcodes,stringsAsFactors=FALSE,header=FALSE)[,1]

qcpass <- read.csv(opt$whitelist,stringsAsFactors=FALSE,header=FALSE)[,1]
message(paste0("Found ", length(qcpass), " QC-passed droplets"))
# Read in raw counts for all called cells
fldr <- dirname(opt$matrix)
sce <- read10xCounts(fldr, col.names=TRUE)
sce <- sce[,cells]

message(paste0("Extracting single-cell data for ", ncol(sce), " droplets"))
#Barcode-DonorID data
demux <- read.table(opt$donorIDs, stringsAsFactors=FALSE, header=TRUE)

message(paste0("Reading in demultiplexing info from: ", opt$donorIDs))
#Assignment Probabilites
demux.fldr <- dirname(opt$donorIDs)
singlet <- paste0(demux.fldr,"/prob_singlet.tsv.gz")
singlet <- gzfile(singlet,'rt')
doublet <- paste0(demux.fldr,"/prob_doublet.tsv.gz")
doublet <- gzfile(doublet,'rt')
probs <- read.table(singlet,stringsAsFactors=FALSE,header=TRUE)
probs.doublet <- read.table(doublet,stringsAsFactors=FALSE,skip=1)

message("Computing summary stats")
# ---- Summarize Data ----
# compute QC stats
qc <- data.frame(perCellQCMetrics(sce))
qc$cell <- colnames(sce)
add <- data.frame("cell"=cells,
		  "QCPass"=cells %in% qcpass,
		  stringsAsFactors=FALSE)
add <- dplyr::left_join(add,demux)
qc <- dplyr::left_join(add,qc)


tot <- sum(qc$QCPass)
dbl.rate <- round((sum(qc$QCPass & qc$donor_id=="doublet") / tot) * 100,2)
unassgnd.rate <- round((sum(qc$QCPass & qc$donor_id=="unassigned") / tot) * 100,2)

######## All of this is useless atm because the vireo output is buggy
######## 
# dbl.comb <- ggplot(qc[qc$donor_id=="doublet",], aes(x=best_doublet)) +
#     geom_bar(colour="black") +
#     ylab("# Cells") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# smry <- dplyr::group_by(qc[!qc$donor_id %in% c("unassigned", "doublet"),],donor_id)
# smry <- data.frame(dplyr::summarize(smry,n=n()))
# smry$nDoub <- sapply(smry$donor_id, function(donor) sum(grepl(donor,qc[qc$donor_id=="doublet","best_doublet"])))
# 
# ndobs <- ggplot(smry, aes(x=n, y=nDoub)) +
#     geom_point() +
#     geom_label_repel(aes(label=donor_id)) +
#     xlab("# Singlets of Donor") +
#     ylab("# Doublets of Donor")
######## 

# ---- QC plots ----
message("Generating QC plots")
# Donor distribution
fplot <- dplyr::group_by(qc, donor_id, QCPass) 
fplot <- dplyr::summarize(fplot,dplyr::n()) # what is 'n()'?
colnames(fplot)[3] <- "Freq"

don.dist <- ggplot(fplot, aes(x=donor_id, y=Freq, fill=QCPass)) +
    geom_bar(stat="identity", colour="black") +
    geom_text(aes(label=Freq),position=position_stack(vjust=.5)) +
    scale_fill_manual(values=c("grey80","grey30")) +
    ylab("# Cells") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Library size
lib.dist <- ggplot(qc, aes(y=sum, x=donor_id, fill=QCPass)) +
    geom_boxplot() +
    scale_fill_manual(values=c("grey80","grey30")) +
    scale_y_log10() +
    ylab("Library Size") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
	  legend.position="none")


# ---- Assignemnt -----

prob.hist <- ggplot(qc, aes(x=prob_max)) +
    geom_histogram() +
    ggtitle("Probabilty of assginment") 

vars.hist <- ggplot(qc, aes(x=n_vars)) +
    geom_histogram() +
    ggtitle("Number of variables") 


#Ordering should be identical anyways
rownames(probs) <- probs$cell
rownames(probs.doublet) <- probs.doublet$V1
probs.doublet <- probs.doublet[rownames(probs),]
probs$Doublet <- rowSums(probs.doublet[,-1])

library(viridis)
ord <- dplyr::arrange(qc, donor_id)[,"cell"]
probs.lt <- reshape2::melt(probs[ord,])
probs.lt$cell <- factor(probs.lt$cell,levels=ord)
colnames(probs.lt)[3] <- "Prob"
prob.heat <- ggplot(probs.lt, aes(x=variable, y=cell, fill=Prob)) +
    geom_tile() +
    scale_fill_viridis() +
    ggtitle("Assignemnt Probabilities") +
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
	axis.text.x = element_text(angle = 45, hjust = 1)) 

header <- paste0(tot, "QC-pass cells with ",dbl.rate,"% doublets and ", unassgnd.rate,"% unassigned")
ttle <- ggdraw() + 
  draw_label(header, fontface='bold', x=0, hjust=0, size=20) 

plots <- plot_grid(prob.heat, don.dist,
		   lib.dist, NULL,
		   prob.hist, vars.hist,
		   ncol=2)
pout <- plot_grid(ttle, plots, ncol=1, rel_heights=c(0.1, 1))

ggsave(opt$plots,width=10, height=15)
