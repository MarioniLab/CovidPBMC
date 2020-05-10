#! /usr/bin/env Rscript

### ADT QC testing
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(DropletUtils)
library(Matrix)
library(reldist)
library(optparse)
library(reshape2)

# rather than pull in all of edgeR we only need this one function
# it operates over a matrix and returns a vector
parser <- OptionParser() 
parser <- add_option(parser, c("-p", "--plots"), type="character",
                     help="Path to pdf for plots")

parser <- add_option(parser, c("-t", "--CITEdirectory"), type="character",
                     help="Paths to CellRanger CITE-seq output directories, comma-separated list")

parser <- add_option(parser, c("-t", "--ADTSCE"), type="character",
                     help="Paths to ADT-specific SCE object, comma-separated list")


message("Extracting CITE-seq CellRanger metrics")
sample.list <- unlist(strsplit(opt$CITEdirectory, split=",", fixed=TRUE))

# sample names are in the directory name
samp.names <- lapply(sample.list, FUN=function(P) unlist(lapply(strsplit(P, fixed=TRUE, split="/"),
                                                          FUN=function(sP) paste0(sP[length(sP)-2]))))
samp.names <- as.factor(unlist(samp.names))
names(sample.list) <- samp.names

# assume ordering is the same
sce.sample.list <- unlist(strsplit(opt$ADTSCE, split=",", fixed=TRUE))
names(sce.sample.list) <- samp.names
cite.metric.list <- list()
cite.sce.list <- list()

for(x in levels(samp.names)){
      x.samp <- levels(samp.names)[x]
      x.dir <- sample.list[[x.samp]]
      x.metric.file <- paste0(x.dir, "outs/metrics_summary.csv")
      cite.metrics <- read.table(x.metric.file,
                               sep=",", header=TRUE, stringsAsFactors=FALSE)
      cite.metrics <- as.data.frame(t(apply(cite.metrics, 2, FUN=function(X) as.numeric(gsub(X, pattern="[\\,\\%]", replacement="")))))
      cite.metrics$SampID <- x.samp

      percent.colnames <- c("Sequencing.Saturation", "Q30.Bases.in.Barcode", "Q30.Bases.in.RNA.Read", "Q30.Bases.in.Sample.Index",
                            "Q30.Bases.in.UMI", "Fraction.Reads.in.Cells", "Antibody:.Mean.Reads.per.Cell", "Antibody:.Sequencing.Saturation")
      cell.colnames <- c("Estimated.Number.of.Cells")
      adt.colnames <- c("Antibody:.Q30.Bases.in.Barcode", "Antibody:.Q30.Bases.in.Antibody.Read", "Antibody:.Q30.Bases.in.Sample.Index",
                        "Antibody:.Q30.Bases.in.UMI", "Antibody:.Fraction.Antibody.Reads", "Antibody:.Fraction.Antibody.Reads.Usable",
                        "Antibody:.Fraction.Reads.in.Barcodes.with.High.UMI.Counts", "Antibody:.Fraction.Unrecognized.Antibody")

     x.cite.reads.melt <- melt(cite.metrics[, c(setdiff(colnames(cite.metrics), c(cell.colnames, percent.colnames, adt.colnames)))], 
                       id.vars="SampID")
     x.cite.percent.melt <- melt(cite.metrics[, c("SampID", percent.colnames)], id.vars="SampID")
     x.cite.cells.melt <- melt(cite.metrics[, c("SampID", cell.colnames)], id.vars="SampID")
     x.cite.adt.melt <- melt(cite.metrics[, c("SampID", adt.colnames)], id.vars="SampID")

     cite.metric.list[[x.samp]] <- list("reads"=x.cite.reads.melt, "percent"=x.cite.percent.melt, "cells"=x.cite.cells.melt, "adt"=x.cite.adt.melt)
     
     ## read in ADT SCE objects
     sce.file <- sce.sample.list[[x.samp]]
     x.sce <- readRDS(sce.file)
     cite.sce.list[[x.samp]] <- x.sce
}

cite.reads.melt <- do.call(rbind.data.frame, list=lapply(cite.metric.list, FUN=function(MB) MB[["reads"]]))
cite.cells.melt <- do.call(rbind.data.frame, list=lapply(cite.metric.list, FUN=function(MB) MB[["cells"]]))
cite.percent.melt <- do.call(rbind.data.frame, list=lapply(cite.metric.list, FUN=function(MB) MB[["percent"]]))
cite.adt.melt <- do.call(rbind.data.frame, list=lapply(cite.metric.list, FUN=function(MB) MB[["adt"]]))

reads.plot <- ggplot(cite.reads.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_y_log10() +
  labs(x="Metric", y="Value")

cells.plot <- ggplot(cite.cells.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_y_log10() +
  labs(x="Metric", y="Value")

percent.plot <- ggplot(cite.percent.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Metric", y="Value")

adt.plot <- ggplot(cite.adt.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Metric", y="Value")

message("Iterating over Antibody targets")
## Is there enough coverage of the antibody targets given the sequencing?
# make a cell x ADT matrix

# pull in the ADT matrix for the sample

cite.expmat.list <- list()
cite.adt.list <- list()
for(x in levels(samp.names)){
  samp.x <- levels(samp.names)[x]
  
  cite.sce <- cite.sce.list[[samp.x]]
  adt.mat <- counts(cite.sce)
  
  # I can randomly down sample this matrix at different proportions to see how bad the drop in ADT coverage is
  # the y-axis will be the Gini index as the measure of evenness of coverage
  down.props <- c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01)
  adt.gini.vec <- c()
  adt.n.targets <- c()
  
  for(x in seq_along(down.props)){
    x.prop <- down.props[x]
    adt.down.mat <- downsampleMatrix(adt.mat, prop=x.prop, bycol = FALSE)
       
    # count non-zero antibody targets
    # gini index is a measure of how even the coverage is across targets
    adt.freq <- apply(adt.down.mat, 1, function(S) sum(S > 0))
    x.gini.adt <- gini(adt.freq[adt.freq > 0])
    
    adt.gini.vec <- c(adt.gini.vec, x.gini.adt)
    adt.n.targets <- c(adt.n.targets, sum(adt.freq > 0))
  }
  
  x.adt.gini.df <- data.frame("Props"=down.props, "Gini"=adt.gini.vec, "NTargets"=adt.n.targets, "SampID"=samp.x)
  cite.adt.list[[samp.x]] <- x.adt.gini.df
  
  # calculate for each target the percentage of ADT counts for that protein
  col.sums <- colSums(adt.mat)
  adt.pc.exprs <- apply(adt.mat, 1, FUN=function(A) A/col.sums)
  # in each cell rank the proteins with their % expression
  adt.rank.exprs <- apply(adt.pc.exprs, 2, FUN=function(AT) rank(AT))
  ave.ranks <- rowMeans(adt.rank.exprs)
  names(ave.ranks) <- rowData(cite.sce)$Symbol
  
  # output the top 50 ADTs that contribute most to the protein expression in each single cell
  top.50 <- names(ave.ranks[order(ave.ranks, decreasing=TRUE)][1:50])
  top.50.expr.df <- as.data.frame(as.matrix(adt.pc.exprs[top.50, ]))
  top.50.expr.df$CellID <- colnames(cite.sce)
  
  top.50.expr.melt <- melt(top.50.expr.df, id.vars="CellID")
  top.50.expr.melt$SampID <- samp.x
  
  cite.expmat.list[[samp.x]] <- top.50.expr.melt
}

cite.gini.df <- do.call(rbind.data.frame, cite.adt.list)
cite.pcexprs.df <- do.call(rbind.data.frame, cite.expmat.list)

cite.gini.plot <- ggplot(cite.gini.df, aes(x=Props, y=Gini, fill=SampID)) +
  geom_point(shape=21, size=3) +
  theme_cowplot() +
  #scale_fill_colorblind() +
  expand_limits(y=c(0, 1)) +
  labs(x="Proporiton of data", y="Gini index")

cite.adt.plot <- ggplot(cite.gini.df, aes(x=Props, y=NTargets, fill=SampID)) +
  geom_point(shape=21, size=3) +
  theme_cowplot() +
  #scale_fill_colorblind() +
  expand_limits(y=c(0, 1)) +
  labs(x="Proporiton of data", y="#Target")

cite.pc.plot <- ggplot(cite.pcexprs.df, aes(x=variable, y=value)) +
  geom_point(size=1) +
  theme_cowplot() +
  facet_wrap(~SampID) +
  coord_flip() +
  labs(x="Percent ADT library", y="Ab Target")


plt.title <- ggdraw() + draw_label("ADT QC metrics",
                                   fontface="bold", x=0, y=0.5, hjust=0, size=20)

plots <- plot_grid(reads.plot, cells.plot, percent.plot, adt.plot,
                   cite.gini.plot, cite.adt.plot, cite.pc.plot,
                   ncol=2) 
plot.out <- plot_grid(plt.title, plots, ncol=1, rel_heights=c(0.1, 0.95))

ggsave(plot=plot.out, file=opt$plots, height=16.95, width=16.95)

sink(file="/dev/null")
rm(list=ls())
gc()
sink(file=NULL)




