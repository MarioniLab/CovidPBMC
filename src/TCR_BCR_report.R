#! /usr/bin/env Rscript

### BCR, TCR QC testing
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(DropletUtils)
library(Matrix)
library(reldist)
library(optparse)
library(ggthemes)

# rather than pull in all of edgeR we only need this one function
# it operates over a matrix and returns a vector

parser <- OptionParser() 
parser <- add_option(parser, c("-l", "--logs"), type="character",
   	  	     help="Path to csv for summary (optional)")

parser <- add_option(parser, c("-p", "--plots"), type="character",
                     help="Path to pdf for plots")

parser <- add_option(parser, c("-t", "--TCRdirectory"), type="character",
                     help="Paths to CellRanger TCR output directories, comma-separated list")

parser <- add_option(parser, c("-b", "--BCRdirectory"), type="character",
                     help="Paths to CellRanger BCR output directories, comma-separated list")

opt <- parse_args(parser)

message("Extracting BCR CellRanger metrics")
sample.list <- unlist(strsplit(opt$BCRdirectory, split=",", fixed=TRUE))
# sample names are in the directory name
samp.names <- lapply(sample.list, FUN=function(P) unlist(lapply(strsplit(P, fixed=TRUE, split="/"),
                                                          FUN=function(sP) paste0(sP[length(sP)]))))
samp.names <- as.factor(unlist(samp.names))
print(samp.names)
names(sample.list) <- levels(samp.names)
print(names(sample.list))

bcr.metric.list <- list()
bcr.annot.list <- list()
for(x in seq_along(levels(samp.names))){
      x.samp <- levels(samp.names)[x]
      x.dir <- sample.list[[x.samp]]
      x.metric.file <- paste0(x.dir, "/outs/metrics_summary.csv")
      bcr.metrics <- read.table(x.metric.file,
                               sep=",", header=TRUE, stringsAsFactors=FALSE)
      bcr.metrics <- as.data.frame(t(apply(bcr.metrics, 2, FUN=function(X) as.numeric(gsub(X, pattern="[\\,\\%]", replacement="")))))
      bcr.metrics$SampID <- x.samp

      x.annot.file <- paste0(x.dir, "/outs/filtered_contig_annotations.csv")
      bcr.annots <- read.table(x.annot.file,
                               sep=",", header=TRUE, stringsAsFactors=FALSE)
      bcr.annots$SampID <- x.samp
      print(head(bcr.annots))
      bcr.annot.list[[x.samp]] <- bcr.annots

      percent.colnames <- c("Valid.Barcodes", "Q30.Bases.in.Barcode", "Q30.Bases.in.RNA.Read.1", "Q30.Bases.in.Sample.Index",
                      "Q30.Bases.in.UMI", "Reads.Mapped.to.Any.V.D.J.Gene", "Reads.Mapped.to.IGH", "Reads.Mapped.to.IGK",
                      "Reads.Mapped.to.IGL", "Fraction.Reads.in.Cells", "Cells.With.Productive.V.J.Spanning.Pair",
                      "Cells.With.Productive.V.J.Spanning..IGK..IGH..Pair", "Cells.With.Productive.V.J.Spanning..IGL..IGH..Pair",
                      "Cells.With.IGH.Contig", "Cells.With.IGK.Contig", "Cells.With.IGL.Contig",
                      "Cells.With.CDR3.annotated.IGH.Contig", "Cells.With.CDR3.annotated.IGK.Contig", "Cells.With.CDR3.annotated.IGL.Contig",
                      "Cells.With.V.J.Spanning.IGH.Contig", "Cells.With.V.J.Spanning.IGK.Contig", "Cells.With.V.J.Spanning.IGL.Contig",
                      "Cells.With.Productive.IGH.Contig", "Cells.With.Productive.IGK.Contig", "Cells.With.Productive.IGL.Contig")
     cell.colnames <- c("Estimated.Number.of.Cells", "Number.of.Cells.With.Productive.V.J.Spanning.Pair")
     div.colnames <- c("Paired.Clonotype.Diversity")

     x.bcr.reads.melt <- melt(bcr.metrics[, c(setdiff(colnames(bcr.metrics), c(cell.colnames, percent.colnames, div.colnames)))], 
                       id.vars="SampID")
     x.bcr.percent.melt <- melt(bcr.metrics[, c("SampID", percent.colnames)], id.vars="SampID")
     x.bcr.cells.melt <- melt(bcr.metrics[, c("SampID", cell.colnames)], id.vars="SampID")
     x.bcr.diversity.melt <- melt(bcr.metrics[, c("SampID", div.colnames)], id.vars="SampID")

     bcr.metric.list[[x.samp]] <- list("reads"=x.bcr.reads.melt, "percent"=x.bcr.percent.melt, "cells"=x.bcr.cells.melt, "diversity"=x.bcr.diversity.melt)
}

bcr.reads.melt <- do.call(rbind.data.frame, lapply(bcr.metric.list, FUN=function(MB) MB[["reads"]]))
bcr.cells.melt <- do.call(rbind.data.frame, lapply(bcr.metric.list, FUN=function(MB) MB[["cells"]]))
bcr.percent.melt <- do.call(rbind.data.frame, lapply(bcr.metric.list, FUN=function(MB) MB[["percent"]]))
bcr.diversity.melt <- do.call(rbind.data.frame, lapply(bcr.metric.list, FUN=function(MB) MB[["diversity"]]))

reads.plot <- ggplot(bcr.reads.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity', position='dodge') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_y_log10() +
  labs(x="Metric", y="Value")

cells.plot <- ggplot(bcr.cells.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity', position='dodge') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_y_log10() +
  labs(x="Metric", y="Value")

percent.plot <- ggplot(bcr.percent.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity', position='dodge') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Metric", y="Value")

diversity.plot <- ggplot(bcr.diversity.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity', position='dodge') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Metric", y="Value")

message("Iterating over BCR annotations")
## Is there enough coverage of the clonotypes given the sequencing?
# make a cell x clonotype matrix
# this should be split by light and heavy chain
## Is there enough coverage of the clonotypes given the sequencing?
# make a cell x clonotype matrix
# this should be split by light and heavy chain
bcr.clono.list <- list()
for(x in seq_along(levels(samp.names))){
      samp.x <- levels(samp.names)[x]

      bcr.annots <- bcr.annot.list[[samp.x]]
      print(head(bcr.annots))
      bcrL.clono.mat <- acast(formula=cdr3 ~ barcode,
                       value.var='umis', 
                       fill=0, # need to fill in the missing values with 0's
                       fun.aggregate = median, na.rm=TRUE,
                       data=bcr.annots[bcr.annots$chain %in% c("IGH"), ])
     bcrL.clono.mat <- as(bcrL.clono.mat, "dgCMatrix")

     # heavy chain
     bcrH.clono.mat <- acast(formula=cdr3 ~ barcode, 
                        value.var='umis', 
                        fill=0, # need to fill in the missing values with 0's
                        fun.aggregate = median, na.rm=TRUE,
                        data=bcr.annots[bcr.annots$chain %in% c("IGL"), ])
     bcrH.clono.mat <- as(bcrH.clono.mat, "dgCMatrix")

     # kappa chain
     # heavy chain
     bcrK.clono.mat <- acast(formula=cdr3 ~ barcode, 
                        value.var='umis', 
                        fill=0, # need to fill in the missing values with 0's
                        fun.aggregate = median, na.rm=TRUE,
                        data=bcr.annots[bcr.annots$chain %in% c("IGK"), ])
     bcrK.clono.mat <- as(bcrK.clono.mat, "dgCMatrix")

     # I can randomly down sample this matrix at different proportions to see how bad the drop in clonotype diversity is
     # the y-axis will be the Gini index as the measure of evenness of coverage
     down.props <- c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01)
     heavy.gini.vec <- c()
     light.gini.vec <- c()
     kappa.gini.vec <- c()

     kappa.n.clones <- c()
     heavy.n.clones <- c()
     light.n.clones <- c()
     for(x in seq_along(down.props)){
           x.prop <- down.props[x]
	   light.down.mat <- downsampleMatrix(bcrL.clono.mat, prop=x.prop, bycol = FALSE)
	   heavy.down.mat <- downsampleMatrix(bcrH.clono.mat, prop=x.prop, bycol = FALSE)
  	   kappa.down.mat <- downsampleMatrix(bcrK.clono.mat, prop=x.prop, bycol = FALSE)
  
           # count non-zero clonotypes
	   light.clon.freq <- apply(light.down.mat, 1, function(S) sum(S > 0))
  	   heavy.clon.freq <- apply(heavy.down.mat, 1, function(S) sum(S > 0))
  	   kappa.clon.freq <- apply(kappa.down.mat, 1, function(S) sum(S > 0))
  
	   x.gini.light <- gini(light.clon.freq[light.clon.freq > 0])
  	   x.gini.heavy <- gini(heavy.clon.freq[heavy.clon.freq > 0])
  	   x.gini.kappa <- gini(kappa.clon.freq[kappa.clon.freq > 0])
  
	   light.gini.vec <- c(light.gini.vec, x.gini.light)
  	   heavy.gini.vec <- c(heavy.gini.vec, x.gini.heavy)
  	   kappa.gini.vec <- c(kappa.gini.vec, x.gini.kappa)
  
	   light.n.clones <- c(light.n.clones, sum(light.clon.freq > 0))
           heavy.n.clones <- c(heavy.n.clones, sum(heavy.clon.freq > 0))
           kappa.n.clones <- c(kappa.n.clones, sum(kappa.clon.freq > 0))
    }

    x.bcrL.gini.df <- data.frame("Props"=down.props, "Gini"=light.gini.vec, "NClones"=light.n.clones, "Chain"="Lambda", "SampID"=samp.x)
    x.bcrK.gini.df <- data.frame("Props"=down.props, "Gini"=kappa.gini.vec, "NClones"=kappa.n.clones, "Chain"="Kappa", "SampID"=samp.x)
    x.bcrH.gini.df <- data.frame("Props"=down.props, "Gini"=heavy.gini.vec, "NClones"=heavy.n.clones, "Chain"="Heavy", "SampID"=samp.x)
    x.bcr.gini.df <- do.call(rbind.data.frame, list("light"=x.bcrL.gini.df, "heavy"=x.bcrH.gini.df, "kappa"=x.bcrK.gini.df))

    bcr.clono.list[[samp.x]] <- x.bcr.gini.df    
}

bcr.gini.df <- do.call(rbind.data.frame, bcr.clono.list)

bcr.gini.plot <- ggplot(bcr.gini.df, aes(x=Props, y=Gini, fill=SampID)) +
  geom_point(shape=21, size=3) +
  theme_cowplot() +
  scale_fill_colorblind() +
  expand_limits(y=c(0, 1)) +
  facet_wrap(~Chain, nrow=1) + 
  labs(x="Proporiton of data", y="Gini index")

bcr.clono.plot <- ggplot(bcr.gini.df, aes(x=Props, y=NClones, fill=SampID)) +
  geom_point(shape=21, size=3) +
  theme_cowplot() +
  scale_fill_colorblind() +
  expand_limits(y=c(0, 1)) +
  facet_wrap(~Chain, nrow=1) +
  labs(x="Proporiton of data", y="#Clonotypes")

plt.title <- ggdraw() + draw_label("BCR QC metrics",
                                   fontface="bold", x=0, y=0.5, hjust=0, size=20)

plots <- plot_grid(reads.plot, cells.plot, percent.plot, diversity.plot,
                   bcr.gini.plot, bcr.clono.plot,
                   ncol=2) 
plot.out <- plot_grid(plt.title, plots, ncol=1, rel_heights=c(0.1, 0.95))

ggsave(plot=plot.out, file=gsub(opt$plots, pattern="\\.pdf", replacement="_BCR.pdf"), height=16.95, width=16.95)

sink(file="/dev/null")
rm(list=c("bcr.clono.list", "x.bcr.gini.df", "bcr.annot.list"))
gc()
sink(file=NULL)


### ---------- TCR Seq results ------------- ###
message("Extracting TCR CellRanger metrics")
print(opt$TCRdirectory)
sample.list <- unlist(strsplit(opt$TCRdirectory, split=",", fixed=TRUE))
# sample names are in the directory name
samp.names <- lapply(sample.list, FUN=function(P) unlist(lapply(strsplit(P, fixed=TRUE, split="/"),
                                                          FUN=function(sP) paste0(sP[length(sP)]))))
samp.names <- as.factor(unlist(samp.names))
names(sample.list) <- samp.names
tcr.metric.list <- list() 
tcr.annot.list <- list()

for(x in seq_along(levels(samp.names))){
      x.samp <- levels(samp.names)[x]
      x.dir <- sample.list[[x.samp]]
      x.metric.file <- paste0(x.dir, "/outs/metrics_summary.csv")

      tcr.metrics <- read.table(x.metric.file,
                          sep=",", header=TRUE, stringsAsFactors=FALSE)
      tcr.metrics <- as.data.frame(t(apply(tcr.metrics, 2, FUN=function(X) as.numeric(gsub(X, pattern="[\\,\\%]", replacement="")))))
      tcr.metrics$SampID <- x.samp

      x.annots.file <- paste0(x.dir, "/outs/filtered_contig_annotations.csv")
      tcr.annots <- read.table(x.annots.file,
                         sep=",", header=TRUE, stringsAsFactors=FALSE)
      tcr.annots$SampID <- x.samp
      tcr.annot.list[[x.samp]] <- tcr.annots

      percent.colnames <- c("Valid.Barcodes", "Q30.Bases.in.Barcode", "Q30.Bases.in.RNA.Read.1", "Q30.Bases.in.Sample.Index",
                      "Q30.Bases.in.UMI", "Reads.Mapped.to.Any.V.D.J.Gene", "Reads.Mapped.to.TRA", "Reads.Mapped.to.TRB",
                      "Fraction.Reads.in.Cells",
                      "Cells.With.Productive.V.J.Spanning.Pair", "Cells.With.Productive.V.J.Spanning..TRA..TRB..Pair",
                      "Cells.With.TRA.Contig", "Cells.With.TRB.Contig", "Cells.With.CDR3.annotated.TRA.Contig",
                      "Cells.With.CDR3.annotated.TRB.Contig", "Cells.With.V.J.Spanning.TRA.Contig", "Cells.With.V.J.Spanning.TRB.Contig",
                      "Cells.With.Productive.TRA.Contig", "Cells.With.Productive.TRB.Contig")
      cell.colnames <- c("Estimated.Number.of.Cells", "Number.of.Cells.With.Productive.V.J.Spanning.Pair")
      div.colnames <- c("Paired.Clonotype.Diversity")

      x.tcr.reads.melt <- melt(tcr.metrics[, c(setdiff(colnames(tcr.metrics), c(cell.colnames, percent.colnames, div.colnames)))], 
                       id.vars="SampID")
      x.tcr.percent.melt <- melt(tcr.metrics[, c("SampID", percent.colnames)], id.vars="SampID")
      x.tcr.cells.melt <- melt(tcr.metrics[, c("SampID", cell.colnames)], id.vars="SampID")
      x.tcr.diversity.melt <- melt(tcr.metrics[, c("SampID", div.colnames)], id.vars="SampID")

      tcr.metric.list[[x.samp]] <- list("reads"=x.tcr.reads.melt, "percent"=x.tcr.percent.melt, "cells"=x.tcr.cells.melt, "diversity"=x.tcr.diversity.melt)
}

tcr.reads.melt <- do.call(rbind.data.frame, lapply(tcr.metric.list, FUN=function(MB) MB[["reads"]]))
tcr.cells.melt <- do.call(rbind.data.frame, lapply(tcr.metric.list, FUN=function(MB) MB[["cells"]]))
tcr.percent.melt <- do.call(rbind.data.frame, lapply(tcr.metric.list, FUN=function(MB) MB[["percent"]]))
tcr.diversity.melt <- do.call(rbind.data.frame, lapply(tcr.metric.list, FUN=function(MB) MB[["diversity"]]))

reads.plot <- ggplot(tcr.reads.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity', position='dodge') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_y_log10() +
  labs(x="Metric", y="Value")

cells.plot <- ggplot(tcr.cells.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity', position='dodge') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_y_log10() +
  labs(x="Metric", y="Value")

percent.plot <- ggplot(tcr.percent.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity', position='dodge') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Metric", y="Value")

diversity.plot <- ggplot(tcr.diversity.melt, aes(x=variable, y=value, fill=SampID)) +
  geom_bar(stat='identity', position='dodge') +
  theme_cowplot() +
  scale_fill_colorblind() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Metric", y="Value")


## Is there enough coverage of the clonotypes given the sequencing?
# make a cell x clonotype matrix
# I think I need a separate matrix for the TRA and TRB chains
# why is htis so memory intensive????? would acast be better?
message("Iterating over TCR clonotype results")
tcr.clono.list <- list()
for(x in seq_along(levels(samp.names))){
      samp.x <- levels(samp.names)[x]
      tcr.annots <- tcr.annot.list[[samp.x]]

      tcra.clono.mat <- acast(formula=cdr3 ~ barcode, 
                       value.var='umis', 
                       fill=0, # need to fill in the missing values with 0's
                       fun.aggregate = median, na.rm=TRUE,
                       data=tcr.annots[tcr.annots$chain %in% c("TRA"), ])
       tcra.clono.mat <- as(tcra.clono.mat, "dgCMatrix")

       tcrb.clono.mat <- acast(formula=cdr3 ~ barcode,
                        value.var='umis', 
                        fill=0, # need to fill in the missing values with 0's
                        fun.aggregate = median, na.rm=TRUE,
                        data=tcr.annots[tcr.annots$chain %in% c("TRB"), ])

       tcrb.clono.mat <- as(tcrb.clono.mat, "dgCMatrix")

       multi.clono.mat <- acast(formula=cdr3 ~ barcode, 
                        value.var='umis', 
                        fill=0, # need to fill in the missing values with 0's
                        fun.aggregate = median, na.rm=TRUE,
                        data=tcr.annots[tcr.annots$chain %in% c("Multi"), ])
       multi.clono.mat <- as(multi.clono.mat, "dgCMatrix")

       # I can randomly down sample this matrix at different proportions to see how bad the drop in clonotype diversity is
       # the y-axis will be the Gini index as the measure of evenness of coverage
       down.props <- c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01)
       tra.gini.vec <- c()
       trb.gini.vec <- c()
       multi.gini.vec <- c()

       tra.n.clones <- c()
       trb.n.clones <- c()
       multi.n.clones <- c()

       for(x in seq_along(down.props)){
           x.prop <- down.props[x]
	   # split each matrix into TCRA,TCRB, and multi
  	   tra.down.mat <- downsampleMatrix(tcra.clono.mat,
                                   prop=x.prop, bycol = FALSE)
  	   trb.down.mat <- downsampleMatrix(tcrb.clono.mat,
                                   prop=x.prop, bycol = FALSE)
           multi.down.mat <- downsampleMatrix(multi.clono.mat,
                                     prop=x.prop, bycol = FALSE)
           # count non-zero clonotypes
  	   tra.clon.freq <- apply(tra.down.mat, 1, function(S) sum(S > 0))
  	   trb.clon.freq <- apply(trb.down.mat, 1, function(S) sum(S > 0))
  	   multi.clon.freq <- apply(multi.down.mat, 1, function(S) sum(S > 0))
  
	   x.gini.tra <- gini(tra.clon.freq[tra.clon.freq > 0])
  	   tra.gini.vec <- c(tra.gini.vec, x.gini.tra)
  	   x.gini.trb <- gini(trb.clon.freq[trb.clon.freq > 0])
  	   trb.gini.vec <- c(trb.gini.vec, x.gini.trb)
  	   x.gini.multi <- gini(multi.clon.freq[multi.clon.freq > 0])
  	   multi.gini.vec <- c(multi.gini.vec, x.gini.multi)
  
	   tra.n.clones <- c(tra.n.clones, sum(tra.clon.freq > 0))
  	   trb.n.clones <- c(trb.n.clones, sum(trb.clon.freq > 0))
  	   multi.n.clones <- c(multi.n.clones, sum(multi.clon.freq > 0))
       }

       x.tcra.gini.df <- data.frame("Props"=down.props, "Gini"=tra.gini.vec, "NClones"=tra.n.clones, "Chain"="TRA", "SampID"=samp.x)
       x.tcrb.gini.df <- data.frame("Props"=down.props, "Gini"=trb.gini.vec, "NClones"=trb.n.clones, "Chain"="TRB", "SampID"=samp.x)
       x.multi.gini.df <- data.frame("Props"=down.props, "Gini"=multi.gini.vec, "NClones"=multi.n.clones, "Chain"="Multi", "SampID"=samp.x)
       x.tcr.gini.df <- do.call(rbind.data.frame, list("tra"=x.tcra.gini.df, "trb"=x.tcrb.gini.df, "multi"=x.multi.gini.df))
       tcr.clono.list[[samp.x]] <- x.tcr.gini.df
}

tcr.gini.df <- do.call(rbind.data.frame, tcr.clono.list)

tcr.gini.plot <- ggplot(tcr.gini.df, aes(x=Props, y=Gini, fill=SampID)) +
  geom_point(shape=21, size=3) +
  theme_cowplot() +
  scale_fill_colorblind() +
  facet_wrap(~Chain, nrow=1) +
  expand_limits(y=c(0, 1)) +
  labs(x="Proporiton of data", y="Gini index")

tcr.clono.plot <- ggplot(tcr.gini.df, aes(x=Props, y=NClones, fill=SampID)) +
  geom_point(shape=21, size=3) +
  theme_cowplot() +
  scale_fill_colorblind() +
  facet_wrap(~Chain, nrow=1) +
  expand_limits(y=c(0, 1)) +
  labs(x="Proporiton of data", y="#Clonotypes")

plt.title <- ggdraw() + draw_label("TCR CellRanger metrics",
                                   fontface="bold", x=0, y=0.5, hjust=0, size=20)

plots <- plot_grid(reads.plot, cells.plot, percent.plot, diversity.plot,
      	 	   tcr.gini.plot, tcr.clono.plot,
                   ncol=2)
plot.out <- plot_grid(plt.title, plots, ncol=1, rel_heights=c(0.1, 0.95))

ggsave(plot=plot.out, file=gsub(opt$plots, pattern="\\.pdf", replacement="_TCR.pdf"), height=16.95, width=15.5)



