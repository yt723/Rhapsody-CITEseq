#### Comparison of IFNG expression by major cell types 

library(Seurat)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(cowplot)
library(edgeR)
library(rstatix)
library(ggpubr)
library(tidyr)

rm(list=ls())

path = "~/Documents/UCSD_CAV/CITEseq"

#load the data
pbmc <- readRDS(file = paste0(path, "/output/2_major_clusters/pbmc_major_celltype_annotated.rds"))

pbmc$major_celltype2 <- as.character(pbmc$major_celltype)
pbmc$major_celltype2 <- ifelse(pbmc$major_celltype2 %in% c("Naive_CD8T", "Memory_CD8T"), "CD8T", pbmc$major_celltype2)
pbmc$major_celltype2 <- ifelse(pbmc$major_celltype2 %in% c("Naive_CD4T", "Memory_CD4T"), "CD4T", pbmc$major_celltype2)
table(pbmc$major_celltype2)

# obtain aggregation matrix (sum)
cts <- AggregateExpression(pbmc, assays = "RNA", group.by = c("major_celltype2","SampleName"), return.seurat = F)
cts <- cts$RNA
cts.t <- as.data.frame(t(cts))
#cts.t[1:10, 1:10]

# modify labels 
cell.label <- gsub('_.*',"", rownames(cts.t)) 
cts.split <- split.data.frame(cts.t, f = factor(cell.label))

cts.split <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)','\\1', rownames(x))
  rownames(x) <- gsub('-','_', rownames(x))
  t(x)
})

# prepare metadata
colData <- unique(pbmc@meta.data[,c("SampleName","SampleType","Library") ])
colData <- colData[order(colData$SampleName),] 
rownames(colData) <- colData$SampleName

## CD8T 
cts.celltype <- cts.split$`CD8T`
obj <- DGEList(counts = cts.celltype, samples = colData)
#head(obj$samples,25)
obj <- calcNormFactors(obj)
log2CPM <- cpm(obj, log=TRUE)
exp_CD8T <- removeBatchEffect(log2CPM, batch = obj$samples$Library)
exp_CD8T <- as.data.frame(t(exp_CD8T))
exp_CD8T$SampleType <- colData$SampleType

stat_test <- exp_CD8T %>% wilcox_test(IFNG~SampleType) %>% add_xy_position() # Wilcoxon rank sum test
stat_test$p <- round(stat_test$p, 3)
y.pos <- stat_test$y.position 

gg1 <- ggplot(exp_CD8T, aes(x=SampleType, y=IFNG)) +
  geom_boxplot(fill= c("maroon","grey"), width=0.5) +
  geom_point(aes(group = SampleType), size=0.5) +
  theme_cowplot()+
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, face = 'bold'),
        axis.title.y = element_text(size = 9, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  stat_pvalue_manual(stat_test, label = 'p={p}', hide.ns = FALSE, label.size = 3) +
  ylim(0,10) + xlab("") +ylab("Pseudobulk IFNG expression")+ ggtitle("CD8T"); gg1

## CD4T 
cts.celltype <- cts.split$`CD4T`
obj <- DGEList(counts = cts.celltype, samples = colData)
head(obj$samples,25)
obj <- calcNormFactors(obj)
log2CPM <- cpm(obj, log=TRUE)
exp_CD4T <- removeBatchEffect(log2CPM, batch = obj$samples$Library)
exp_CD4T <- as.data.frame(t(exp_CD4T))
exp_CD4T$SampleType <- colData$SampleType

stat_test <- exp_CD4T %>% wilcox_test(IFNG~SampleType) %>% add_xy_position() # Wilcoxon rank sum test 
stat_test$p <- round(stat_test$p, 3)
stat_test$y.position <- y.pos

gg2 <- ggplot(exp_CD4T, aes(x=SampleType, y=IFNG)) +
  geom_boxplot(fill= c("maroon","grey"), width=0.5) +
  geom_point(aes(group = SampleType), size=0.5) +
  theme_cowplot()+
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, face = 'bold'),
        axis.title.y = element_text(size = 9, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  stat_pvalue_manual(stat_test, label = 'p={p}', hide.ns = FALSE, label.size = 3) +
  ylim(0,10) + xlab("") +ylab("Pseudobulk IFNG expression")+ ggtitle("CD4T"); gg2

## NK
cts.celltype <- cts.split$`Natural-killer`
obj <- DGEList(counts = cts.celltype, samples = colData)
#head(obj$samples,25)
obj <- calcNormFactors(obj)
log2CPM <- cpm(obj, log=TRUE)
exp_NK <- removeBatchEffect(log2CPM, batch = obj$samples$Library)
exp_NK <- as.data.frame(t(exp_NK))
exp_NK$SampleType <- colData$SampleType

stat_test <- exp_NK %>% wilcox_test(IFNG~SampleType) %>% add_xy_position() # Wilcoxon rank sum test 
stat_test$p <- round(stat_test$p, 3)
stat_test$y.position <- y.pos

gg3 <- ggplot(exp_NK, aes(x=SampleType, y=IFNG)) +
  geom_boxplot(fill= c("maroon","grey"), width=0.5) +
  geom_point(aes(group = SampleType), size=0.5) +
  theme_cowplot()+
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, face = 'bold'),
        axis.title.y = element_text(size = 9, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  stat_pvalue_manual(stat_test, label = 'p={p}', hide.ns = FALSE, label.size = 3) +
  ylim(0,10) + xlab("") +ylab("Pseudobulk IFNG expression")+ ggtitle("NK"); gg3

## *Figure 1D*
ggarrange(gg1, gg2, gg3, ncol = 3)
ggsave(file = paste0(path, "/output/3_pseudobulk1/boxplot_IFNG_ggarrange_pseudobulk_cpm.png"),
       width = 4.8, height = 2.3, units = "in", bg = 'white')


## Correlation between CD4T and CD8T IFNG expression  
IFNG_exp <- data.frame(CD8T = exp_CD8T$IFNG, CD4T = exp_CD4T$IFNG, 
                       SampleType=c(rep("CAV",6), rep("Normal_HTx",12)))

cor <- cor.test(IFNG_exp$CD8T, IFNG_exp$CD4T)
cor[["p.value"]] # p=0.0007098253
cor[["estimate"]] # r=0.722361 

## *Figure 1E*
ggplot(IFNG_exp, aes(x=CD8T, y=CD4T)) +
  geom_point(aes(shape = SampleType, color = SampleType), size=2)+
  scale_color_manual(values = c('maroon','grey30')) +
  geom_smooth(method = "lm", color = 'black', linetype = "dashed") +
  theme_classic() +
  theme(legend.text = element_text(size = 9, face = 'bold'), 
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9, face = 'bold')) + 
  labs(color = "", shape ="") +
  xlab("IFNG expression by CD8+ T cell") + ylab("IFNG expression by CD4+ T cell")
ggsave(file = paste0(path, "/output/3_pseudobulk1/cor_scatter_IFNG_CD4T_vs_CD8T.png"),
       width = 3.7, height = 2.3, units = "in", bg = 'white')

## major_celltype cytotoxic marker expression 
cts <- AggregateExpression(pbmc, assays = "RNA", group.by = c("major_celltype","SampleName"), return.seurat = F)
cts <- cts$RNA
col.order <- colnames(cts)
colData <- data.frame(sample = colnames(cts))
rownames(colData) <- colData$sample
colData <- separate(colData, col = "sample", into = c("major_celltype", "SampleName"), sep = "_", remove = F)
meta <- unique(pbmc@meta.data[,c("SampleName", "Library", "SampleType")])
meta$SampleName <- gsub("_","-", meta$SampleName)

colData <- colData %>% merge(meta, by.x = c("SampleName"), by.y = "SampleName")
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
dat.summary <- out.z %>% group_by(major_celltype) %>% 
  summarise(IFNG = mean(IFNG), TNF = mean(TNF), GZMK = mean(GZMK), GZMB = mean(GZMB), PRF1 = mean(PRF1))
dat.summary <- dat.summary %>% pivot_longer(cols = 2:6, names_to = "marker", values_to = "average")
dat.summary$major_celltype <- gsub("-","_", dat.summary$major_celltype)
dat.summary$major_celltype <- factor(dat.summary$major_celltype, 
                                     levels = c("Naive_CD4T","Memory_CD4T","Naive_CD8T","Memory_CD8T",
                                                "Natural_killer", "B", "Classical_mono","Nonclassical_mono","Dendritic" ))
dat.summary$marker <- factor(dat.summary$marker,
                             levels = c("IFNG","TNF","GZMK","GZMB","PRF1"))

## *Figure 1C*
ggplot(dat.summary, aes(x=marker, y=major_celltype)) +
  geom_point(aes(colour = average), size=5) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size = 12, face='bold',,angle = 270, hjust = 0, vjust=0.5),
        axis.text.y = element_text(size = 12, face='bold'),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 11), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray50", linetype = 2))+
  scale_color_gradient(low = 'blue', high = 'yellow', name = 'Pseudobulk \nexpression \nzscore') +
  xlab("") + ylab("")
ggsave(file = paste0(path, "/output/3_pseudobulk1/dotplot_cytotoxic_major_celltype.png"),
       width = 5.3, height = 3, units = "in", bg = 'white')


# *** end of pseudobulk 1 ***