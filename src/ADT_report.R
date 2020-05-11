#! /usr/bin/env Rscript

### ADT QC testing
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(DropletUtils)
library(Matrix)
library(reldist)
library(ggthemes)
library(optparse)
library(reshape2)

# rather than pull in all of edgeR we only need this one function
# it operates over a matrix and returns a vector
parser <- OptionParser() 
parser <- add_option(parser, c("-p", "--plots"), type="character",
                     help="Path to pdf for plots")

parser <- add_option(parser, c("-t", "--CITEdirectory"), type="character",
                     help="Paths to CellRanger CITE-seq output directories, comma-separated list")

parser <- add_option(parser, c("-a", "--ADTSCE"), type="character",
                     help="Paths to complete SCE object - individual ADT SCE objects will be derived from this")

opt <- parse_args(parser)

message("Extracting CITE-seq CellRanger metrics")
sample.list <- strsplit(opt$CITEdirectory, split=",", fixed=TRUE)
# sample names are in the directory name
samp.names <- lapply(sample.list, FUN=function(P) unlist(lapply(strsplit(P, fixed=TRUE, split="/"),
                                                          FUN=function(sP) paste0(sP[length(sP)]))))
samp.names <- as.factor(unlist(samp.names))
names(sample.list) <- samp.names

# given an input SCE assum ADT SCE have the structure <name>-<sample>_Abs.RDS
sce.sample.list <- list()
for(x in seq_along(levels(samp.names))){
  x.samp <- levels(samp.names)[x]
  x.dir <- dirname(opt$ADTSCE)
  x.sce.file <- list.files(x.dir, pattern="Ab\\.RDS")
  x.sce.file <- x.sce.file[grepl(x.sce.file, pattern=gsub(x.samp, pattern="_CITE", replacement=""))]
  print(x.sce.file)
  sce.sample.list[[x.samp]] <- x.sce.file
}

cite.metric.list <- list()
cite.sce.list <- list()
print(sample.list)
print(names(sample.list))
print(samp.names)

for(x in seq_along(levels(samp.names))){
  x.samp <- levels(samp.names)[x]
  print(x.samp)
  x.dir <- sample.list[[x.samp]]
  x.metric.file <- paste0(x.dir, "/outs/metrics_summary.csv")
  cite.metrics <- read.table(x.metric.file,
                             sep=",", header=TRUE, stringsAsFactors=FALSE)
  cite.metrics <- as.data.frame(t(apply(cite.metrics, 2, FUN=function(X) as.numeric(gsub(X, pattern="[\\,\\%]", replacement="")))))
  cite.metrics$SampID <- x.samp
  
  percent.colnames <- c("Sequencing.Saturation", "Q30.Bases.in.Barcode", "Q30.Bases.in.RNA.Read", "Q30.Bases.in.Sample.Index",
                        "Q30.Bases.in.UMI", "Fraction.Reads.in.Cells", "Antibody..Sequencing.Saturation")
  cell.colnames <- c("Estimated.Number.of.Cells")
  adt.colnames <- c("Antibody..Q30.Bases.in.Barcode", "Antibody..Q30.Bases.in.Antibody.Read", "Antibody..Q30.Bases.in.Sample.Index",
                    "Antibody..Q30.Bases.in.UMI", "Antibody..Fraction.Antibody.Reads", "Antibody..Fraction.Antibody.Reads.Usable",
                    "Antibody..Fraction.Reads.in.Barcodes.with.High.UMI.Counts", "Antibody..Fraction.Unrecognized.Antibody")
  exclude.column <- c("Antibody..Median.UMIs.per.Cell..summed.over.all.recognized.antibody.barcodes.")
  
  # change Antibody for Ab to save on report space later
  x.cite.reads.melt <- melt(cite.metrics[, c(setdiff(colnames(cite.metrics), c(exclude.column, cell.colnames, percent.colnames, adt.colnames)))], 
                            id.vars="SampID")
  x.cite.reads.melt$variable <- as.character(x.cite.reads.melt$variable)
  x.cite.reads.melt$variable <- gsub(x.cite.reads.melt$variable, pattern="Antibody", replacement="Ab")
  
  x.cite.percent.melt <- melt(cite.metrics[, c("SampID", percent.colnames)], id.vars="SampID")
  x.cite.percent.melt$variable <- as.character(x.cite.percent.melt$variable)
  x.cite.percent.melt$variable <- gsub(x.cite.percent.melt$variable, pattern="Antibody", replacement="Ab")
  
  x.cite.cells.melt <- melt(cite.metrics[, c("SampID", cell.colnames)], id.vars="SampID")
  x.cite.cells.melt$variable <- as.character(x.cite.cells.melt$variable)
  x.cite.cells.melt$variable <- gsub(x.cite.cells.melt$variable, pattern="Antibody", replacement="Ab")
  
  x.cite.adt.melt <- melt(cite.metrics[, c("SampID", adt.colnames)], id.vars="SampID")
  x.cite.adt.melt$variable <- as.character(x.cite.adt.melt$variable)
  x.cite.adt.melt$variable <- gsub(x.cite.adt.melt$variable, pattern="Antibody", replacement="Ab")
  
  cite.metric.list[[x.samp]] <- list("reads"=x.cite.reads.melt, "percent"=x.cite.percent.melt, "cells"=x.cite.cells.melt, "adt"=x.cite.adt.melt)
     
  ## read in ADT SCE objects
  sce.file <- sce.sample.list[[x.samp]]
  sce.dir <- lapply(strsplit(opt$ADTSCE, fixed=TRUE, split="/"),
                    FUN=function(sP) paste(sP[1:(length(sP)-1)], collapse="/"))
  print(sce.dir)
  print(paste(sce.dir, sce.file, sep="/"))
  x.sce <- readRDS(paste(sce.dir, sce.file, sep="/"))
  cite.sce.list[[x.samp]] <- x.sce
}


cite.reads.melt <- do.call(rbind.data.frame, lapply(cite.metric.list, FUN=function(MB) MB[["reads"]]))
cite.cells.melt <- do.call(rbind.data.frame, lapply(cite.metric.list, FUN=function(MB) MB[["cells"]]))
cite.percent.melt <- do.call(rbind.data.frame, lapply(cite.metric.list, FUN=function(MB) MB[["percent"]]))
cite.adt.melt <- do.call(rbind.data.frame, lapply(cite.metric.list, FUN=function(MB) MB[["adt"]]))

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
  coord_flip() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Metric", y="Percent")

adt.plot <- ggplot(cite.adt.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity') +
  theme_cowplot() +
  scale_fill_colorblind() +
  coord_flip() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Metric", y="Value")

message("Iterating over Antibody targets")
## Is there enough coverage of the antibody targets given the sequencing?
# make a cell x ADT matrix

# pull in the ADT matrix for the sample
cite.expmat.list <- list()
cite.adt.list <- list()
for(x in seq_along(levels(samp.names))){
  samp.x <- levels(samp.names)[x]
  
  cite.sce <- cite.sce.list[[samp.x]]
  adt.mat <- counts(cite.sce)
  
  # I can randomly down sample this matrix at different proportions to see how bad the drop in ADT coverage is
  # the y-axis will be the Gini index as the measure of evenness of coverage
  down.props <- c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01)
  adt.gini.vec <- c()
  adt.n.targets <- c()
  message(paste0("Downsampling ADT counts for sample: ", samp.x))
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
  
  message(paste0("Calculating ADT library coverage for sample: ", samp.x))
  # calculate for each target the percentage of ADT counts for that protein
  col.sums <- colSums(adt.mat)
  adt.pc.exprs <- t(apply(adt.mat, 1, FUN=function(A) A/col.sums))
  rownames(adt.pc.exprs) <- rowData(cite.sce)$Symbol

  # in each cell rank the proteins with their % expression
  adt.rank.exprs <- apply(adt.pc.exprs, 2, FUN=function(AT) rank(AT))
  ave.ranks <- rowMeans(adt.rank.exprs)
  names(ave.ranks) <- rowData(cite.sce)$Symbol
  
  # output the top 25 ADTs that contribute most to the protein expression in each single cell
  top.50 <- names(ave.ranks[order(ave.ranks, decreasing=TRUE)][1:25])

  # Do I really want to transpose this matrix?!
  top.50.expr.df <- as.data.frame(t(as.matrix(adt.pc.exprs[top.50, ])))
  top.50.expr.df$CellID <- colnames(cite.sce)
  
  top.50.expr.melt <- melt(top.50.expr.df, id.vars="CellID")
  top.50.expr.melt$SampID <- samp.x
  
  cite.expmat.list[[samp.x]] <- top.50.expr.melt
}

cite.gini.df <- do.call(rbind.data.frame, cite.adt.list)
cite.pcexprs.df <- do.call(rbind.data.frame, cite.expmat.list)

message("Making plots")
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
  geom_boxplot() +
  theme_cowplot() +
  facet_wrap(~SampID) +
  coord_flip() +
  labs(y="Percent ADT library", x="Ab Target")

plt.title <- ggdraw() + draw_label("ADT QC metrics",
                                   fontface="bold", x=0, y=0.5, hjust=0, size=20)

panel_a <- plot_grid(reads.plot, cells.plot, nrow=1, rel_widths=c(1, 0.5))
panel_b <- plot_grid(percent.plot, adt.plot, nrow=1, rel_widths=c(1, 1))
panel_c.1 <- plot_grid(cite.gini.plot, cite.adt.plot, ncol=1, rel_widths=c(1, 1))
panel_c.2 <- plot_grid(panel_c.1, cite.pc.plot, nrow=1, rel_widths=c(1, 1))

plots <- plot_grid(panel_a, panel_b, panel_c.2,
                   nrow=3, rel_heights=c(1, 1, 1))
plot.out <- plot_grid(plt.title, plots, ncol=1, rel_heights=c(0.1, 0.95))

message(paste0("Saving plots to: ", opt$plots))
ggsave(plot=plot.out, file=opt$plots, height=16.95, width=16.95)

sink(file="/dev/null")
rm(list=ls())
gc()
sink(file=NULL)




