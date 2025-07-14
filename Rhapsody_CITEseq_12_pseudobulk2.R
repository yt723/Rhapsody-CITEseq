#### Comparison of gene expression by pseudobulk 

library(Seurat)
library(ggplot2)
library(dplyr)
library(edgeR)
library(tidyr)

rm(list=ls())

## IFNG expression by CD8T, CD4T, NK clusters 
CD8T <- readRDS(file = "output/4_CD8T_clusters/CD8T_clusters.rds")
CD4T <- readRDS(file = "output/5_CD4T_clusters/CD4T_clusters.rds")
NK <- readRDS(file = "output/7_NK_clusters/NK_clusters.rds")

CD8T$cluster <- paste0("CD8T_C",CD8T$cluster_0.3)
cts.cd8 <- AggregateExpression(CD8T, assays = "RNA", group.by = c("cluster","SampleName"), return.seurat = F)
cts.cd8 <- cts.cd8$RNA

CD4T$cluster <- paste0("CD4T_C",CD8T$cluster_0.3)
cts.cd4 <- AggregateExpression(CD4T, assays = "RNA", group.by = c("cluster","SampleName"), return.seurat = F)
cts.cd4 <- cts.cd4$RNA

NK$cluster <- paste0("NK_C",CD8T$cluster_0.3)
cts.nk <- AggregateExpression(NK, assays = "RNA", group.by = c("cluster","SampleName"), return.seurat = F)
cts.nk <- cts.nk$RNA

cts <- cbind(cts.cd8, cts.cd4, cts.nk)
colData <- data.frame(sample = colnames(cts))
col.order <- colnames(cts)
rownames(colData) <- colData$sample
colData <- separate(colData, col = "sample", into = c("cluster","SampleName"), sep = "_",
                    remove = F)

meta <- distinct(CD8T@meta.data[,c("SampleName","SampleType","Library") ])
meta$SampleName <- gsub("_","-", meta$SampleName)
colData <- colData %>% merge(meta, by.x = "SampleName", by.y = "SampleName")
rownames(colData) <- colData$sample
colData <- colData[order(factor(rownames(colData), levels = col.order)),]

obj <- DGEList(counts = cts, samples = colData)
#head(obj$samples, 25)
obj <- calcNormFactors(obj)
log2CPM <- cpm(obj, log=TRUE)
out <- removeBatchEffect(log2CPM, batch = obj$samples$Library)
out <- as.data.frame(t(out))

out.z <- sapply(out, function(dat) (dat - mean(dat))/sd(dat))
rownames(out.z) <- rownames(out)
out.z <- out.z %>% merge(colData, by = 0)

dat.summary <- out.z %>% group_by(cluster) %>% 
  summarise(IFNG = mean(IFNG), TNF = mean(TNF), GZMK = mean(GZMK), GZMB = mean(GZMB), PRF1 = mean(PRF1))
dat.summary <- dat.summary %>% pivot_longer(cols = 2:6, names_to = "marker", values_to = "average")
dat.summary$cluster <- gsub("-","_", dat.summary$cluster)
dat.summary$marker <- factor(dat.summary$marker,
                             levels = c("IFNG","TNF","GZMK","GZMB","PRF1"))

# *Figure S8B
ggplot(dat.summary, aes(y=marker, x=cluster)) +
  geom_point(aes(colour = average), size=3) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size = 9, face='bold',,angle = 270, hjust = 0, vjust=0.5),
        axis.text.y = element_text(size = 9, face='bold'),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 9), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray50", linetype = 3))+
  scale_color_gradient(low = 'blue', high = 'yellow', name = 'Pseudobulk \nexpression \nzscore') +
  xlab("") + ylab("")
ggsave(file = "output/12_pseudobulk2/dotplot_cytotoxic_CD8_CD4_NK_clusters.png",
       width = 5, height = 3, units = "in", bg = 'white')

# *** end of pseudobulk analysis2 ***