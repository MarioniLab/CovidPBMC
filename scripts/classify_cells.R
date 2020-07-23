#! /usr/bin/env Rscript

# Classify input cells using the Random Forest classifer built using the Blish lab annotations
### COVID-19 predicting cell types
library(SingleCellExperiment)
library(scran)
library(scater)
library(ggplot2)
library(ggthemes)
library(scales)
library(cowplot)
library(biomaRt)
library(randomForest)

# read in the Random Forest classifier trained on the 10X Genomics GEX data set
# need a gene.df
biomaRt.connection <- useMart("ensembl", "hsapiens_gene_ensembl")

gene.df <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                 mart=biomaRt.connection)
rownames(gene.df) <- gene.df$ensembl_gene_id

## predict cell types
covid.sce <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE.RDS")

########## ------------ PBMC Random Forest classifier ------------ ##########
pbmc.rf <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/classifiers/PBMC_RF_classifier.RDS")

tier1.importance <- as.data.frame(pbmc.rf$importance)
tier1.importance$gene_id <- rownames(pbmc.rf$importance)
tier1.importance <- merge(tier1.importance, gene.df, by.x='gene_id', by.y='ensembl_gene_id')
accuracy <- tier1.importance[order(tier1.importance$MeanDecreaseAccuracy, decreasing=TRUE), ]

mean.decrease.p <- ggplot(accuracy[1:30, ], 
                          aes(x=reorder(external_gene_name,
                                        MeanDecreaseAccuracy),
                              y=MeanDecreaseAccuracy)) +
  geom_point() +
  coord_flip() +
  theme_cowplot() +
  labs(x="Gene", y="Mean Decrease Accuracy")

accuracy <- accuracy[order(accuracy$MeanDecreaseGini, decreasing=TRUE), ]

gini.decrease.p <- ggplot(accuracy[1:30, ], 
                          aes(x=reorder(external_gene_name,
                                        MeanDecreaseGini),
                              y=MeanDecreaseGini)) +
  geom_point() +
  coord_flip() +
  theme_cowplot() +
  labs(x="Gene", y="Mean Decrease Gini")

plot_grid(mean.decrease.p, gini.decrease.p, ncol=2)

ggsave(filename="/mnt/scratchb/jmlab/morgan02/Covid/reports/PBMC_RF_topFeatures.pdf",
       height=5.95, width=6.95, useDingbats=FALSE)

pdf("/mnt/scratchb/jmlab/morgan02/Covid/reports/PBMC_RF_ClassErrors.pdf", height=4.95, width=9.75, useDingbats=FALSE)
par(mfrow=c(3, 3))
plot(pbmc.rf$err.rate[ ,1], main="OOB error")
plot(pbmc.rf$err.rate[, 2], main="B cells")
plot(pbmc.rf$err.rate[, 3], main="Dendritic cells")
plot(pbmc.rf$err.rate[, 4], main="Monocytes")
plot(pbmc.rf$err.rate[, 5], main="NK cells")
plot(pbmc.rf$err.rate[, 6], main="NKT cells")
plot(pbmc.rf$err.rate[, 7], main="Platelets")
plot(pbmc.rf$err.rate[, 8], main="T cells")
dev.off()

pbmc.features <- rownames(pbmc.rf$importance)
pbmc.train <- predict(pbmc.rf, newdata=as.matrix(t(logcounts(covid.sce[pbmc.features, ]))))

# compare classifiers
pbmc.class.df <- data.frame("CellID"=names(pbmc.train), "PBMC.Class"=pbmc.train)
write.table(pbmc.class.df,
            file="/mnt/scratchb/jmlab/morgan02/Covid/classifiers/PBMC_classifed.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

########## ------------ Blish Healthy control Random Forest classifier ------------ ##########
healthy.rf <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/classifiers/Blish-Healthy_PBMC_RF_classifier.RDS")

tier1.importance <- as.data.frame(healthy.rf$importance)
tier1.importance$gene_id <- rownames(healthy.rf$importance)
tier1.importance <- merge(tier1.importance, gene.df, by.x='gene_id', by.y='ensembl_gene_id')
accuracy <- tier1.importance[order(tier1.importance$MeanDecreaseAccuracy, decreasing=TRUE), ]

mean.decrease.p <- ggplot(accuracy[1:30, ], 
                          aes(x=reorder(external_gene_name,
                                        MeanDecreaseAccuracy),
                              y=MeanDecreaseAccuracy)) +
  geom_point() +
  coord_flip() +
  theme_cowplot() +
  labs(x="Gene", y="Mean Decrease Accuracy")

accuracy <- accuracy[order(accuracy$MeanDecreaseGini, decreasing=TRUE), ]

gini.decrease.p <- ggplot(accuracy[1:30, ], 
                          aes(x=reorder(external_gene_name,
                                        MeanDecreaseGini),
                              y=MeanDecreaseGini)) +
  geom_point() +
  coord_flip() +
  theme_cowplot() +
  labs(x="Gene", y="Mean Decrease Gini")

plot_grid(mean.decrease.p, gini.decrease.p, ncol=2)

ggsave(filename="/mnt/scratchb/jmlab/morgan02/Covid/reports/Blish_healthy_RF_topFeatures.pdf",
       height=5.95, width=6.95, useDingbats=FALSE)

pdf("/mnt/scratchb/jmlab/morgan02/Covid/reports/Blish_healthy_RF_ClassErrors.pdf", height=9.95, width=9.75, useDingbats=FALSE)
par(mfrow=c(5, 5))
healthy.class <- colnames(healthy.rf$err.rate)
for(x in seq_along(healthy.class)){
  plot(healthy.rf$err.rate[ ,x], main=healthy.class[x])
}
dev.off()

healthy.features <- rownames(healthy.rf$importance)
covid.symbol <- as.matrix(t(logcounts(covid.sce[intersect(rownames(covid.sce), healthy.features), ])))

healthy.train <- predict(healthy.rf,
                         newdata=covid.symbol)

healthy.class.df <- data.frame("CellID"=names(healthy.train), "Healthy.Class"=healthy.train)
write.table(healthy.class.df,
            file="/mnt/scratchb/jmlab/morgan02/Covid/classifiers/Blish_healthy_classifed.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

########## ------------ Blish COVID19 control Random Forest classifier ------------ ##########
covid.rf <- readRDS("/mnt/scratchb/jmlab/morgan02/Covid/classifiers/Blish-COVID19_PBMC_RF_classifier.RDS")

tier1.importance <- as.data.frame(covid.rf$importance)
tier1.importance$gene_id <- rownames(covid.rf$importance)
tier1.importance <- merge(tier1.importance, gene.df, by.x='gene_id', by.y='ensembl_gene_id')
accuracy <- tier1.importance[order(tier1.importance$MeanDecreaseAccuracy, decreasing=TRUE), ]

mean.decrease.p <- ggplot(accuracy[1:30, ], 
                          aes(x=reorder(external_gene_name,
                                        MeanDecreaseAccuracy),
                              y=MeanDecreaseAccuracy)) +
  geom_point() +
  coord_flip() +
  theme_cowplot() +
  labs(x="Gene", y="Mean Decrease Accuracy")

accuracy <- accuracy[order(accuracy$MeanDecreaseGini, decreasing=TRUE), ]

gini.decrease.p <- ggplot(accuracy[1:30, ], 
                          aes(x=reorder(external_gene_name,
                                        MeanDecreaseGini),
                              y=MeanDecreaseGini)) +
  geom_point() +
  coord_flip() +
  theme_cowplot() +
  labs(x="Gene", y="Mean Decrease Gini")

plot_grid(mean.decrease.p, gini.decrease.p, ncol=2)

ggsave(filename="/mnt/scratchb/jmlab/morgan02/Covid/reports/Blish_covid_RF_topFeatures.pdf",
       height=5.95, width=6.95, useDingbats=FALSE)

pdf("/mnt/scratchb/jmlab/morgan02/Covid/reports/Blish_covid_RF_ClassErrors.pdf", height=9.95, width=9.75, useDingbats=FALSE)
par(mfrow=c(5, 5))
covid.class <- colnames(covid.rf$err.rate)
for(x in seq_along(covid.class)){
  plot(covid.rf$err.rate[ ,x], main=covid.class[x])
}
dev.off()

covid.features <-rownames(covid.rf$importance)
covid.symbol <- as.matrix(t(logcounts(covid.sce[intersect(rownames(covid.sce), covid.features), ])))

covid.train <- predict(covid.rf,
                         newdata=covid.symbol)

covid.class.df <- data.frame("CellID"=names(covid.train), "Covid.Class"=covid.train)
write.table(covid.class.df,
            file="/mnt/scratchb/jmlab/morgan02/Covid/classifiers/Blish_covid_classifed.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

