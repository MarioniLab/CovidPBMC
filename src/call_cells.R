# Cell Calling with emptyDrops
# This scripts takes the following command line arguments

# ---- Parsing ----
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"), type="character", 
		     help="Path to the raw matrix.mtx.gz")
parser <- add_option(parser, c("-o", "--out"), type="character", 
		     help="Path to output file")
parser <- add_option(parser, c("-f", "--fdrthreshold"), type="double", default=0.001, 
		     help="FDR threshold to call cells [default %default]",
		     metavar="number")
parser <- add_option(parser, c("-u", "--umithreshold"), type="integer", default=100, 
		     help="UMI limit to call cells [default %default]",
		     metavar="number")
parser <- add_option(parser, c("-l", "--logs"), type="character", 
		     help="Path to pdf for plots")
opt <- parse_args(parser)

# ---- Load Data ----
library(DropletUtils)
library(ggplot2)
library(cowplot)
library(BiocParallel)
theme_set(theme_cowplot())

# Read in raw counts
# Using the mtx instead of h5 as emptyDrops doesn't support DelayedArray 
fldr <- dirname(opt$input)
sce <- read10xCounts(fldr,col.names=TRUE)

# Ensure to always set the seed before any call of empty drops
emptyDrops_seed <- function(...,rnd.seed=42) {
    set.seed(rnd.seed)
    out <- emptyDrops(...)
    return(out)
}


# ---- EmptyDrops ----

# Test for empty droplets
niter <- 10000
e.out <- emptyDrops_seed(counts(sce), lower=opt$umithreshold, niters=niter)

# Are any cells above the thresold due to limited precision of estimated p-values
tbl <- table(Sig=e.out$FDR <= opt$fdrthreshold, Limited=e.out$Limited)
is.lim <- tbl[1,2] > 0

# If true increase niter by one order of magnitude and recompute
# This might need to be adopted if there are still limited pvalues
if (is.lim) {
    n.iter <- niter * 10
    e.out <- emptyDrops_seed(counts(sce), lower=opt$umithreshold, niters=niter)
    # Throw warning if unsignificant cells still have limited p-values
    tbl <- table(Sig=e.out$FDR <= opt$fdrthreshold, Limited=e.out$Limited)
    still.lim <- tbl[1,2] > 0
    if (still.lim) {
	warning("Limited P-Values above FDR threshold, some cells might be wrongly labelled as ambient")
    }
}

# ---- QCPlots ----

# P-Value Histogram

# Run emptyDrops with ambient true for P-Val Histogram
all.out <- emptyDrops_seed(counts(sce), lower=opt$umithreshold, test.ambient=TRUE, niter=niter)

# Plot
fplot <- data.frame("X"=all.out$PValue[all.out$Total <= opt$umithreshold & all.out$Total > 0])
pval.hist <- ggplot(fplot, aes(x=X)) +
    geom_histogram(bins=20) +
    xlab("P-Value") +
    ggtitle("P-Value Histogram")

# Knee Plot

# Compute barcode ranks
bcrank <- barcodeRanks(counts(sce))

# Plot
uniq <- !duplicated(bcrank$rank)
fplot <- data.frame("rank"=bcrank$rank[uniq],
		    "total"=bcrank$total[uniq]) 
bcplt <- ggplot(fplot, aes(x=rank, y=total)) +
	    geom_point() +
	    scale_x_log10() +
	    scale_y_log10() +
	    geom_hline(yintercept=metadata(bcrank)$inflection,
		       lty="dashed",color="blue") +
	    geom_hline(yintercept=metadata(bcrank)$knee,
		       lty="dashed") +
	    xlab("Barcode Rank") +
	    ylab("Total UMIs") +
	    ggtitle("Knee plot")


# Combine the two plots
plots <- plot_grid(bcplt, pval.hist, nrow=1)


# ---- Save ----

# Barcodes associated with cells
out <- data.frame("Barcode"=sce$Barcode[which(e.out$FDR<=0.001)])
write.csv(out, file=opt$out, row.names=FALSE)

# QCplots as pdf
ggsave(filename=opt$logs, plot=plots, width=10, height=5)
