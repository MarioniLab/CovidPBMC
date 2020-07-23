#! /usr/bin/env Rscript

library(SingleCellExperiment)
library(scater)
library(scran)

covid.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE.RDS")
covid.fail.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE_QCFail.RDS")


sigac2.tcr.contigs <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAC2_TCR/outs/filtered_contig_annotations.csv",
		                sep=",", header=TRUE, stringsAsFactors=FALSE)
sigac2.tcr.contigs$Sample <- paste0("SIGAA2_", sigac2.tcr.contigs$barcode)
sigac2.tcr.contigs$ExpSample <- "SIGAC2"

sigac3.tcr.contigs <- read.table("/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAC3_TCR/outs/filtered_contig_annotations.csv",
		                sep=",", header=TRUE, stringsAsFactors=FALSE)
sigac3.tcr.contigs$Sample <- paste0("SIGAA3_", sigac3.tcr.contigs$barcode)
sigac3.tcr.contigs$ExpSample <- "SIGAC3"

qc.pass.cells <- colnames(covid.sce)
qc.fail.cells <- setdiff(colnames(covid.fail.sce), colnames(covid.sce))

tcr.contigs <- do.call(rbind.data.frame, list(sigac2.tcr.contigs, sigac3.tcr.contigs)) 
tcr.contigs$QCFail <- "Pass"
tcr.contigs$QCFail[tcr.contigs$Sample %in% qc.fail.cells] <- "Fail"

write.table(tcr.contigs, file="/mnt/scratchb/jmlab/morgan02/Covid/SCE/TCR_contigs_ALL.tsv",
			 sep="\t", quote=FALSE, row.names=FALSE)



