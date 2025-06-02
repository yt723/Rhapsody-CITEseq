#### CD8 T cell sub-clusters 

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

rm(list=ls())

#load the data
pbmc <- readRDS(file ="output/2_major_clusters/pbmc_major_celltype_annotated.rds")

# subset CD8T cells
CD8T <- subset(pbmc, major_celltype %in% c("Naive_CD8T", "Memory_CD8T"))

## Cluster double positive and double negative T cells  
DefaultAssay(CD8T) <- "ADT18"
CD8T <- CD8T %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
#ElbowPlot(CD8T, ndims = 25)
CD8T<- CD8T %>% RunUMAP(dims=1:10) %>% FindNeighbors() %>% FindClusters(resolution=0.2) %>% identity()
CD8T$seurat_clusters <- as.factor(as.numeric(as.character(CD8T$seurat_clusters)) + 1)

## *Figure S3A*
DimPlot(CD8T, group.by ="seurat_clusters",pt.size = 0.3, label = T, repel=F, label.size = 6) + 
  ggsci::scale_color_d3(palette="category20") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 11, face = 'bold'),
        legend.key.height = unit(0.15, "in"),
        legend.key.spacing = unit(0.05, "in"))+
  guides(color = guide_legend(override.aes = list(size=3))) +ggtitle("")
ggsave(file = "output/4_CD8T_clusters/Dimplot_CD8T_clusters_initial.png", 
       width=3.5, height=3, units = "in", bg = 'white')

## *Figure S3A* 
Idents(CD8T) <- "seurat_clusters"
markers <- c("CD3", "CD4", "CD8")
DefaultAssay(CD8T) <- "ADT18"
for(i in markers){
  FeaturePlot(CD8T, features =i, cols =  c(alpha("lightgrey",0.4), "red"), 
                 pt.size = 0.3, min.cutoff =0 ,max.cutoff =5, label = T, label.size = 5)+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 16), 
          legend.text = element_text(size = 9),
          legend.position = "inside",
          legend.position.inside = c(0.9, 0.2),
          legend.key.height = unit(0.15, "in"),
          legend.key.width = unit(0.15, "in"))
    
  ggsave(file = paste0("output/4_CD8T_clusters/FeaturePlot_1/FeaturePlot_", i, ".png"), 
         width=3, height=3, units = "in", bg = 'white')
}

# remove double negative & double positive T cell clusters 
# seurat_clusters == 4: double negative T
# seurat_clusters == 5: double positive T
CD8T <- subset(CD8T, seurat_clusters %in% c("1","2","3"))

# save RDS file 
saveRDS(CD8T, file = "output/4_CD8T_clusters/CD8T_cleaned.rds")

## CD8T cell unsupervised clustering (WNN)
DefaultAssay(CD8T) <- 'integrated'
CD8T <- CD8T %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'pca')
#ElbowPlot(CD8T)
DefaultAssay(CD8T) <- 'integratedADT'
CD8T <- CD8T %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'apca')
CD8T <- FindMultiModalNeighbors(CD8T, reduction.list = list("pca", "apca"), dims.list = list(1:10, 1:10),
                                modality.weight.name = "RNA.weight")
CD8T <- RunUMAP(CD8T, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
  FindClusters(graph.name="wsnn", 
               resolution = c(0.25, 0.3, 0.35), 
               cluster.name = c("cluster_0.25", "cluster_0.3","cluster_0.35"), 
               random.seed=42) %>% identity()
CD8T$cluster_0.25<- as.factor(as.numeric(as.character(CD8T$cluster_0.25)) + 1)
CD8T$cluster_0.3<- as.factor(as.numeric(as.character(CD8T$cluster_0.3)) + 1)
CD8T$cluster_0.35<- as.factor(as.numeric(as.character(CD8T$cluster_0.35)) + 1)

cluster.name = c("cluster_0.25", "cluster_0.3","cluster_0.35")

## *Figure S3B*
for (i in cluster.name){
  DimPlot(CD8T, group.by = i, label = T, repel = F, pt.size = 0.3, label.size = 5) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 16), 
          legend.text = element_text(size = 9),
          legend.key.height = unit(0.15, "in"),
          legend.key.width = unit(0.15, "in")) +
    guides(color = guide_legend(override.aes = list(size=4)))
  
  ggsave(file = paste0("output/4_CD8T_clusters/umap_comparison/DimPlot_CD8T_", i, ".png"),
         height = 3, width = 3.5, units = "in", bg = 'white')
}

# Use "cluster_0.3"
table(CD8T@meta.data$cluster_0.3)
#C8 is too small for characterization (n = 48) 
CD8T <- subset(CD8T, cluster_0.3 == "8", invert = T)

## comparison with SingleR annotations 
library(SingleR)
library(celldex)

monaco.ref <- MonacoImmuneData()
celltype <- c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells",
              "Terminal effector CD8 T cells", "MAIT cells", "Vd2 gd T cells", "Non-Vd2 gd T cells")
ref.subset <- monaco.ref[, colData(monaco.ref)$label.fine %in% celltype]

# convert Seurat to SCE object
sce <- as.SingleCellExperiment(DietSeurat(CD8T), assay = 'RNA')

# run single R
pred.fine <- SingleR(sce, ref = ref.subset, labels = ref.subset$label.fine, de.method = 'classic')

#table(pred.fine$pruned.labels)
df_monaco.fine <- data.frame(monaco.fine =pred.fine$pruned.labels)
rownames(df_monaco.fine) <- pred.fine@rownames
CD8T <- AddMetaData(CD8T, metadata = df_monaco.fine)

## *Figure S3C* 
DimPlot(CD8T, group.by = "monaco.fine", pt.size = 0.3) + 
  scale_color_brewer(palette = "Set2") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 14),
        legend.key.height = unit(0.15, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.key.spacing = unit(0.05, "in")) +
  guides(color = guide_legend(override.aes = list(size=4)))+ ggtitle("")
ggsave(file = "output/4_CD8T_clusters/DimPlot_CD8T_monaco.png",
       height = 3, width = 5.5, units = "in", bg = 'white')

# add cluster names as meta data  
cluster_name = c("1:EMRA","2:EM","3:naive","4:CD38+HLADR+EM","5:CD56+EMRA","6:MAIT","7:CD25+EM")
cluster_df <- data.frame(cluster_0.3 = as.character(1:7), cluster_name0.3 = cluster_name)
df <- CD8T@meta.data 
df$cell_id <- rownames(df)
df <- df %>% merge(cluster_df, by.x = "cluster_0.3", by.y = "cluster_0.3")
rownames(df) <- df$cell_id
CD8T <- AddMetaData(CD8T, df)

# Save RDS file 
saveRDS(CD8T, file = "output/4_CD8T_clusters/CD8T_clusters.rds")


# *Figure 2A* 
col <- ggsci::pal_d3("category10")(7)
colset <- c(col[1:3],col[6],col[4],col[7],col[5])

DimPlot(CD8T, group.by ="cluster_name0.3", label = F, pt.size =0.3) + 
  scale_color_manual(values = colset)+
  NoAxes() + ggtitle("") +
  theme(legend.text = element_text(size = 12, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(file = "output/4_CD8T_clusters/DimPlot_CD8T_cluster_name0.3.png",
       height = 4, width = 6, units = "in", bg = 'white')

# *Figure 2A* 
DimPlot(CD8T, group.by ="SampleType", label = F, pt.size =0.05) + 
  scale_color_manual(values = c("maroon", "gray80")) +
  NoLegend() +NoAxes() + ggtitle("")
ggsave(file = "output/4_CD8T_clusters/DimPlot_SampleType.png",
       height = 2, width = 2, units = "in", bg = 'white')

# *Figure 3B*
features <- rownames(CD8T@assays$ADT)
Idents(CD8T) <- "cluster_0.3"
for(i in features){
  FeaturePlot(CD8T, features =i, cols =  c(alpha("lightgrey",0.4), "blue"),label = T,  
                 label.size = 8, pt.size =0.3, min.cutoff =1 ,max.cutoff =5) +
    theme(plot.title = element_text(size = 28),
          legend.position = "inside",
          legend.position.inside = c(0.15, 0.3)) +
    NoAxes()
  ggsave(file = paste0("output/4_CD8T_clusters/FeaturePlot_2/FeaturePlot_", i, ".png"), 
         width=4, height=3.7, units = "in", bg = 'white')
}

## Cluster cell counts per sample (*Excel file 2*)
df_count <- CD8T@meta.data %>% group_by(cluster_name0.3, SampleName) %>% 
  summarise(num=n()) %>% pivot_wider(names_from = "cluster_name0.3", values_from = "num", values_fill = 0)
write.csv(df_count, file = "output/4_CD8T_clusters/CD8T_cluster_counts_patient.csv", row.names = F) 

## Cluster proportion per patient: CAV vs Normal HTx 
df_prop <- CD8T@meta.data %>% group_by(cluster_name0.3, SampleName, SampleType) %>% summarise(num=n())
df_prop <- df_prop %>% pivot_wider(names_from = 'cluster_name0.3', values_from = 'num', values_fill = 0)

mat_prop <- df_prop[, 3:ncol(df_prop)] %>% as.matrix()
mat_prop <- round(prop.table(mat_prop, margin = 1) *100, 2)
df_prop[,3:ncol(df_prop)] <- mat_prop 
write.csv(df_prop, file = "output/4_CD8T_clusters/CD8T_cluster_proportion_patient.csv", row.names = F)

cav <- df_prop[df_prop$SampleType == "CAV",]
summary(cav$'4:CD38+HLADR+EM') 
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.530   8.665  11.260  20.162  15.152  66.670 
control <- df_prop[df_prop$SampleType == "Normal_HTx",]
summary(control$'4:CD38+HLADR+EM')
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.760   3.180   5.340   5.003   5.945  12.360 

## DGE (MAST)
Idents(CD8T) <- "cluster_0.3"
one_vs_all_rna <- FindAllMarkers(CD8T, assay = "RNA", test.use = "MAST", only.pos = F, latent.vars = "SampleName", min.pct = 0.1)
write.csv(one_vs_all_rna, file = "output/4_CD8T_clusters/dge/CD8T_one_vs_all_rna_MAST.csv")

## DGE (LR)
Idents(CD8T) <- "cluster_0.3"
one_vs_all_rna <- FindAllMarkers(CD8T, assay = "RNA", test.use = "LR", only.pos = F, latent.vars = "SampleName", min.pct = 0.1)

# *Excel file 3*
Idents(CD8T) <- "cluster_name0.3"
one_vs_all_rna <- FindAllMarkers(CD8T, assay = "RNA", test.use = "MAST", only.pos = T, latent.vars = "SampleName", min.pct = 0.1)
one_vs_all_rna <- one_vs_all_rna %>%
  filter(p_val_adj < 0.05)
write.csv(one_vs_all_rna, file = "output/4_CD8T_clusters/dge/CD8T_one_vs_all_rna_cluster_name0.3_MAST.csv")

## Dotplot for surface marker expression
markers <- c("CCR7", "CD45RA", "CD45RO", "CD25", "CD27", "CD28", "CD38", "HLA.DR", "CD56", "CD161")
dotp <- DotPlot(CD8T, features = markers, assay = "ADT", scale = T, group.by = "cluster_name0.3") ;dotp
df.dot <- dotp$data

# *Figure 3A*
ggplot(df.dot, aes(y = id, x = features.plot)) +
  geom_point(aes(color = avg.exp.scaled, size= pct.exp)) +
  scale_color_gradient2(low = "navy", mid = 'white', high = 'red', name ="Average\nExpression" )  +
  scale_radius(limits = c(0,100), range = c(0,5), name = "Percent\nExpression") +
  theme(axis.text.x = element_text(size = 14, face='bold', angle = 270, hjust = 0, vjust=0.5),
        axis.text.y = element_text(size = 14, face='bold'),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray", linetype = 2),
        legend.key.height = unit(0.1, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.key.spacing = unit(0.05, "in"))+
  xlab("") +ylab("")
ggsave(file = "output/4_CD8T_clusters/DotPlot_CD8T_cluster_ADT_markers.png",
       height = 3.5, width = 6, units = 'in', bg = 'white')

## Heatmap for DEGs
dge <- read.csv(file = "output/4_CD8T_clusters/dge/CD8T_one_vs_all_rna_MAST.csv", row.names = 1)

dge$p_val_adj[which(dge$p_val_adj == 0)] <- 1e-300
dge$rank <- dge$avg_log2FC * (-1) * log10(dge$p_val_adj)

markers <- dge %>%  group_by(cluster) %>% top_n(8, rank)  %>% dplyr::arrange(desc(rank), .by_group = T)
markers2 <- markers[!duplicated(markers$gene),]
top.genes <- markers2$gene

dotp <- DotPlot(CD8T, features = top.genes, assay = "RNA", scale = T, group.by = "cluster_name0.3");dotp
df_mat <- dotp$data
df_mat <- df_mat %>% select(features.plot, id, avg.exp.scaled)
df_mat <- df_mat %>% pivot_wider(names_from = "features.plot", values_from = "avg.exp.scaled")
annotation_col <- data.frame(cluster = df_mat$id)
annotation_col <- data.frame(cluster = factor(sort(annotation_col$cluster), labels = sort(annotation_col$cluster)))
rownames(annotation_col) <- annotation_col$cluster

col <- ggsci::pal_d3("category10")(7)
cluster <- c(col[1:3],col[6],col[4],col[7],col[5])
names(cluster) <- sort(annotation_col$cluster)
annotation_colors <- list(cluster = cluster)

hm.mat <- df_mat[,-1] %>% as.matrix() %>% t()
colnames(hm.mat) <- df_mat$id
pal_cols <- colorRampPalette(c("navy","white","red")) 
newnames_row <- lapply(
  rownames(hm.mat) ,
  function(x) bquote(bold(.(x))))

newnames_cols <- lapply(
  colnames(hm.mat),
  function(x) bquote(bold(.(x))))

# *Figure 3C* 
png(filename = "output/4_CD8T_clusters/Heatmap_CD8T_cluster_name0.3_MAST.png", height =15, width =7.5, res=500, units = "in") 
pheatmap::pheatmap(hm.mat, 
                   color = pal_cols(100), 
                   fontsize = 15, fontsize_row = 15, fontsize_number=15, 
                   labels_row = as.expression(newnames_row),
                   annotation_names_col = F,
                   show_colnames = F,
                   border_color = "white",
                   annotation_col = annotation_col,
                   annotation_colors = annotation_colors) #cluster_rows = FALSE
dev.off()


## Comparison of expression of cytotoxic genes between clusters (pseudobulk) 
library(edgeR)

cts <- AggregateExpression(CD8T, assays = "RNA", group.by = c("cluster_0.3","SampleName"), return.seurat = F)
cts <- cts$RNA

col.order <- colnames(cts)
colData <- data.frame(sample = colnames(cts))
rownames(colData) <- colData$sample
colData <- separate(colData, col = "sample", into = c("cluster", "SampleName"), sep = "_", remove = F)
meta <- unique(CD8T@meta.data[,c("SampleName", "Library", "SampleType")])
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

clusters <- distinct(data.frame(cluster = CD8T$cluster_0.3, cluster_name = CD8T$cluster_name0.3))
clusters$cluster <- paste0("g",clusters$cluster)
out.z <- out.z %>% merge(clusters, by.x = "cluster", by.y = "cluster")

dat.summary <- out.z %>% group_by(cluster_name) %>% 
  summarise(IFNG = mean(IFNG), TNF = mean(TNF), GZMK = mean(GZMK), GZMB = mean(GZMB), PRF1 = mean(PRF1))
dat.summary <- dat.summary %>% pivot_longer(cols = 2:6, names_to = "marker", values_to = "average")
dat.summary$marker <- factor(dat.summary$marker,
                             levels = c("IFNG","TNF","GZMK","GZMB","PRF1"))

# *Figure 3F* 
ggplot(dat.summary, aes(x=marker, y=cluster_name)) +
  geom_point(aes(colour = average), size=5) +
  cowplot::theme_cowplot() +
  scale_color_gradient(low = 'blue', high = 'yellow', name = 'Pseudobulk \nexpression \nzscore') +
  theme(axis.text.x = element_text(size = 12, face='bold',,angle = 270, hjust = 0, vjust=0.5),
        axis.text.y = element_text(size = 12, face='bold'),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 11), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray50", linetype = 2))+
  xlab("") +ylab("")
ggsave(file = "output/4_CD8T_clusters/dotplot_cytotoxic_CD8T_cluster_pseudobulk.png",
       width = 5.3, height = 3, units = "in", bg = 'white')

## GO
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

dge <- read.csv(file = "output/4_CD8T_clusters/dge/CD8T_one_vs_all_rna_MAST.csv", row.names = 1)
dge <- dge %>% filter(avg_log2FC>0.3, p_val_adj <0.05)

for (id in as.character(1:7)){
  genes.test <- dge$gene[dge$cluster== id] 
  GO.test <- enrichGO(gene = genes.test, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL")
  GO <- as.data.frame(GO.test)
  write.csv(GO, file = paste0("output/4_CD8T_clusters/GO/CD8T_enrichGO_C",id, ".csv"))
}

go <- paste0("GO_C", as.character(1:7))
for (x in go){
  assign(x, read.csv(file = paste0("output/4_CD8T_clusters/GO/CD8T_enrich",x, ".csv")))
}
topGO <- GO_C4[1:12,]$Description
duplicate <- c("viral process", "cell activation involved in immune response",
               "leukocyte cell-cell adhesion" , "regulation of leukocyte cell-cell adhesion" )
topGO <- topGO[!topGO %in% duplicate]

for (name in go) {
  x <- get(name)
  x <- x %>% filter(Description %in% topGO)
  x$cluster <- name
  assign(name, x)
}

GO <- rbind(GO_C1, GO_C2, GO_C3, GO_C4, GO_C5, GO_C6, GO_C7)
GO.sub <- GO[,c("Description","p.adjust","cluster")]
GO.sub$cluster <- gsub("GO_","",GO.sub$cluster)
GO.sub <- GO.sub %>% pivot_wider(names_from = "cluster", values_from = "p.adjust", values_fill = 1)
GO.sub <- GO.sub %>% pivot_longer(cols = c(2:8), names_to = "cluster", values_to = "p.adjust")
GO.sub$log10.p_adj <- (-1)*log10(GO.sub$p.adjust) 

GO.sub2 <- GO[,c("Description","zScore","cluster")]
GO.sub2$cluster <- gsub("GO_","",GO.sub2$cluster)
GO.sub2 <- GO.sub2 %>% pivot_wider(names_from = "cluster", values_from = "zScore", values_fill = 0)
GO.sub2 <- GO.sub2 %>% pivot_longer(cols = c(2:8), names_to = "cluster", values_to = "zScore")
GO.sub <- GO.sub %>% merge(GO.sub2, by.x = c("cluster", "Description"), by.y = c("cluster", "Description"))

# *Figure 3E* 
ggplot(GO.sub, aes(x = cluster,  y =str_wrap(Description, width = 30), colour  = zScore, size = log10.p_adj)) +
  geom_point() +
  scale_color_viridis_b() +
  scale_radius(range = c(0,5), name = "-log10(p_adj)") +
  theme(axis.text.x = element_text(size = 9, face='bold'),
        axis.text.y = element_text(size = 10, face='bold'),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray", linetype = 3),
        legend.key.height = unit(0.1, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.key.spacing = unit(0.05, "in"))+
  xlab("") +ylab("")
ggsave(file = "output/4_CD8T_clusters/DotPlot_CD8T_cluster_GO.png",
       height = 3, width = 5.5, units = 'in', bg = 'white')
  

## Volcano plot
library(ggrepel)

dge <- read.csv(file = "output/4_CD8T_clusters/dge/CD8T_one_vs_all_rna_MAST.csv", row.names = 1)

dge$p_val_adj[which(dge$p_val_adj == 0)] <- 1e-300
dge <- dge %>% filter(avg_log2FC >0)

for (id in as.character(1:7)){
  res <- dge %>% filter(cluster == id)
  res$diff <- ""
  res$diff[res$avg_log2FC>0.3 & res$p_val_adj<0.05] <- "up"
  
  for (x in 1:nrow(res)){
    if (abs(res$avg_log2FC[x]) >0.3 & res$p_val_adj[x]<0.05){
      res$gene[x] <- res$gene[x]
    }else {
      res$gene[x] <- NA
    }
  } 
  
  ggplot(res, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diff, label=gene))+
    geom_point()+
    theme_minimal()+
    scale_color_manual(values=c("gray","red"))+
    geom_text_repel()+
    xlab("average log2 fold change") + ylab("-log10(adj. p-value)")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10, face = 'bold'),
          legend.title=element_blank(), legend.position = 'none',
          axis.line = element_line(color = 'black'),
          axis.ticks = element_line(color = 'black', linewidth = 0.2),
    )+
    geom_hline(aes(yintercept = 1.3), linetype = 'dashed', color='black')+
    geom_vline(xintercept = 0.3, linetype = 'dashed', color ='black') +
    ggtitle("")
  ggsave(file = paste0("output/4_CD8T_clusters/volcano/CD8T_volcano_C",id,"_pos_MAST.png"),
         width=3.5, height=3, units="in",bg = 'white')
  }

## Differential composition analysis: sccomp 
library(sccomp)

CD8T$SampleType <- factor(CD8T$SampleType, levels = c("Normal_HTx","CAV"))
set.seed(47)
sccomp_result <- CD8T %>% 
  sccomp_estimate(formula_composition = ~SampleType,
                  formula_variability = ~SampleType,
                  .sample = SampleName, 
                  .cell_group = cluster_name0.3,
                  variational_inference = F,
                  bimodal_mean_variability_association = TRUE)

plots <- sccomp_result %>% sccomp_test() %>% plot()

col <- c("gray50","red")
names(col) <- c(FALSE, TRUE)

# *Figure 2A*
plots$credible_intervals_1D[[1]] +
  theme_classic() +
  theme(axis.text.y = element_text(size=9, face = 'bold'),
        axis.text.x = element_text(size = 7, face = 'bold'),
        axis.title.x = element_text(size = 9, face = 'bold'),
        legend.text = element_text(size = 7),
        legend.position = 'none') +
  scale_color_manual(values = col) +
  ylab("")+ xlab("Log odds") +
  guides(color = guide_legend(title  = "FDR < significance threshold (0.025)")) +
  ggtitle("") 
ggsave(file = "output/4_CD8T_clusters/CD8T_cell_composition_sccomp.png", 
       width=3.5, height=2, units="in",bg = 'white')

## cell composition analysis for comparison between the timepoints post-HTx 
CD8T$Post_HTx_5y <- ifelse(CD8T$Post_HTx_Interval_y <=5, "<=5y", ">5y")
sccomp_result <- CD8T %>% 
  sccomp_estimate(formula_composition = ~Post_HTx_5y,
                  .sample = SampleName, 
                  .cell_group = cluster_name0.3,
                  bimodal_mean_variability_association = TRUE,
                  variational_inference = F)

plots <- sccomp_result %>% sccomp_test() %>% plot()

col <- c("gray50","red")
names(col) <- c(FALSE, TRUE)

# *Figure S11C*
plots$credible_intervals_1D[[1]] +
  theme_classic() +
  theme(axis.text.y = element_text(size=9, face = 'bold'),
        axis.text.x = element_text(size = 7, face = 'bold'),
        axis.title.x = element_text(size = 9, face = 'bold'),
        legend.text = element_text(size = 7),
        legend.position = 'none') +
  scale_color_manual(values = col) +
  ylab("")+ xlab("Log odds") +
  guides(color = guide_legend(title  = "FDR < significance threshold (0.025)")) +
  ggtitle("") 
ggsave(file = "output/4_CD8T_clusters/CD8T_cell_composition_sccomp_Post_HTx_5y.png", 
       width=3.5, height=2, units="in",bg = 'white')

# *** end of CD8T clusters ***