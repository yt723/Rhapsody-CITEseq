#### B cell sub-clusters 

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

rm(list=ls())

#load the data
pbmc <- readRDS(file = "output/2_major_clusters/pbmc_major_celltype_annotated.rds")

# subset B cells
Bcell <- subset(pbmc, major_celltype == "B")

## Bcell cell unsupervised clustering (WNN)
DefaultAssay(Bcell) <- 'integrated'
Bcell <- Bcell %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'pca')

DefaultAssay(Bcell) <- 'integratedADT'
Bcell <- Bcell %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'apca')
Bcell <- FindMultiModalNeighbors(Bcell, reduction.list = list("pca", "apca"), dims.list = list(1:10, 1:10), 
                                modality.weight.name = "RNA.weight")

Bcell <- RunUMAP(Bcell, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
  FindClusters(graph.name="wsnn", 
               resolution = c(0.25, 0.3, 0.35), 
               cluster.name = c("cluster_0.25", "cluster_0.3","cluster_0.35"), 
               random.seed=42) %>% identity()
Bcell$cluster_0.25<- as.factor(as.numeric(as.character(Bcell$cluster_0.25)) + 1)
Bcell$cluster_0.3<- as.factor(as.numeric(as.character(Bcell$cluster_0.3)) + 1)
Bcell$cluster_0.35<- as.factor(as.numeric(as.character(Bcell$cluster_0.35)) + 1)

cluster.name = c("cluster_0.25", "cluster_0.3","cluster_0.35")

for (i in cluster.name){
  dimp <- DimPlot(Bcell, group.by = i, label = T, repel = F, pt.size = 0.3, label.size = 4) +
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 8, face = 'bold'),
          title = element_text(size = 8))+
    guides(color = guide_legend(override.aes = list(size=3)))
  ggsave(dimp, file = paste0("output/6_Bcell_clusters/umap_comparison/DimPlot_Bcell_", i, ".png"),
         height = 3, width = 3.5, units = "in", bg = 'white')
}

features <- rownames(Bcell@assays$ADT)
for(i in features){
  FeaturePlot(Bcell, features =i, cols =  c(alpha("lightgrey",0.4), "blue"), pt.size =1, min.cutoff =1 ,max.cutoff =5)
  ggsave(file = paste0("output/6_Bcell_clusters/FeaturePlot_1/FeaturePlot_", i, ".png"), 
         width=7, height=6, units = "in", bg = 'white')
}

# remove doublet cluster 
# cluster_0.3 = 4 -->doublet (CD3+CD19+)
Bcell <- subset(Bcell, cluster_0.3 == "4", invert = T)

## comparison with SingleR annotations  
library(SingleR)
library(celldex)

monaco.ref <- MonacoImmuneData()
celltype <- c("Naive B cells", "Non-switched memory B cells", "Exhausted B cells", "Switched memory B cells")
ref.subset <- monaco.ref[, colData(monaco.ref)$label.fine %in% celltype]

# convert Seurat to SCE object
sce <- as.SingleCellExperiment(DietSeurat(Bcell), assay = 'RNA')

# run single R
pred.fine <- SingleR(sce, ref = ref.subset, labels = ref.subset$label.fine, de.method = 'classic')

table(pred.fine$pruned.labels)
df_monaco.fine <- data.frame(monaco.fine =pred.fine$pruned.labels)
rownames(df_monaco.fine) <- pred.fine@rownames
Bcell <- AddMetaData(Bcell, metadata = df_monaco.fine)

DimPlot(Bcell, group.by = "monaco.fine", pt.size = 0.3) + 
  scale_color_brewer(palette = "Set2") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 14),
        legend.key.height = unit(0.15, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.key.spacing = unit(0.05, "in")) +
  guides(color = guide_legend(override.aes = list(size=4)))+ ggtitle("")
ggsave(file = "output/6_Bcell_clusters/DimPlot_Bcell_monaco.png",
       height = 3, width = 5.5, units = "in", bg = 'white')


# Use "cluster_0.3"
table(Bcell@meta.data$cluster_0.3)

# add cluster names as meta data  
cluster_name = c("1:naive","2:memory CD11c+","3:memory")
cluster_df <- data.frame(cluster_0.3 = as.character(1:3), cluster_name0.3 = cluster_name)
df <- Bcell@meta.data 
df$cell_id <- rownames(df)
df <- df %>% merge(cluster_df, by.x = "cluster_0.3", by.y = "cluster_0.3")
rownames(df) <- df$cell_id
Bcell <- AddMetaData(Bcell, df)

# *Save RDS file* 
saveRDS(Bcell, file = "output/6_Bcell_clusters/Bcell_clusters.rds")

# *Figure 2C*
DimPlot(Bcell, group.by ="cluster_name0.3", label = F, pt.size =0.3) + 
  scale_color_brewer(palette ="Set1")+
  theme(legend.text = element_text(size = 12, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 5))) + ggtitle("") +NoAxes()
ggsave(file = "output/6_Bcell_clusters/DimPlot_Bcell_cluster_name0.3.png",
       height = 4, width = 6, units = "in", bg = 'white')

DimPlot(Bcell, group.by ="SampleType", label = F, pt.size =0.05) + 
  scale_color_manual(values = c("maroon", "gray80")) +
  theme(legend.position = 'top',
        legend.spacing.x = unit(0.1, "in"),
        legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=5), nrow = 1, title.hjust=0.5)) +
  NoAxes() + ggtitle("")
ggsave(file = "output/6_Bcell_clusters/DimPlot_Bcell_SampleType.png",
       height = 3.4, width = 3.1, units = "in", bg = 'white')


# FeaturePlot (*Fig S6B*)
features <- rownames(Bcell@assays$ADT)
Idents(Bcell) <- "cluster_0.3"
for(i in features){
  FeaturePlot(Bcell, features =i, cols =  c(alpha("lightgrey",0.4), "blue"),label = T,  
                 label.size = 8, pt.size =0.3, min.cutoff =1 ,max.cutoff =5) +
    theme(plot.title = element_text(size = 28),
          legend.position = "inside",
          legend.position.inside = c(0.05, 0.8)) +
    NoAxes()
  ggsave(file = paste0("output/6_Bcell_clusters/FeaturePlot_2/FeaturePlot_", i, ".png"), 
         width=4.3, height=3.7, units = "in", bg = 'white')
}


## Cluster cell counts per sample (*Excel file 2*) 
df_count <- Bcell@meta.data %>% group_by(cluster_name0.3, SampleName) %>% 
  summarise(num=n()) %>% pivot_wider(names_from = "cluster_name0.3", values_from = "num", values_fill = 0)
write.csv(df_count, file = "output/6_Bcell_clusters/Bcell_cluster_counts_patient.csv", row.names = F) 

## Cluster proportion per patient: CAV vs Normal HTx 
df_prop <- Bcell@meta.data %>% group_by(cluster_name0.3, SampleName, SampleType) %>% dplyr::summarise(num=n())
df_prop <- df_prop %>% pivot_wider(names_from = 'cluster_name0.3', values_from = 'num', values_fill = 0)

mat_prop <- df_prop[, 3:ncol(df_prop)] %>% as.matrix()
mat_prop <- round(prop.table(mat_prop, margin = 1) *100, 2)
df_prop[,3:ncol(df_prop)] <- mat_prop 

write.csv(df_prop, file = "output/6_Bcell_clusters/Bcell_cluster_proportion_patient.csv", row.names = F) 

## Differential composition analysis (*Figure 2C*) 
library(sccomp)

Bcell$SampleType <- factor(Bcell$SampleType, levels = c("Normal HTx","CAV"))
set.seed(47)
sccomp_result <- Bcell %>% 
  sccomp_estimate(formula_composition = ~SampleType,
                  formula_variability = ~SampleType,
                  .sample = SampleName, 
                  .cell_group = cluster_name0.3,
                  variational_inference = F,
                  bimodal_mean_variability_association = TRUE)

plots <- sccomp_result %>% sccomp_test() %>% plot()

col <- c("gray50","red")
names(col) <- c(FALSE, TRUE)
plots$credible_intervals_1D[[1]] +
  theme_classic() +
  theme(axis.text.y = element_text(size=9, face = 'bold'),
        axis.text.x = element_text(size = 7, face = 'bold'),
        axis.title.x = element_text(size = 9, face = 'bold'),
        legend.text = element_text(size = 7),
        legend.position = '') +
  scale_color_manual(values = col) +
  ylab("")+ xlab("Log odds") +
  guides(color = guide_legend(title  = "FDR < significance threshold (0.025)")) +
  ggtitle("") 
ggsave(file = "output/6_Bcell_clusters/Bcell_cell_composition_sccomp.png", 
       width=3.5, height=1.5, units="in",bg = 'white')

## DGE: MAST
Idents(Bcell) <- "cluster_0.3"
one_vs_all_rna <- FindAllMarkers(Bcell, assay = "RNA", test.use = "MAST", only.pos = F, latent.vars = "SampleName", min.pct = 0.1)
write.csv(one_vs_all_rna, file = "output/6_Bcell_clusters/dge/Bcell_one_vs_all_rna_MAST.csv")

# *Excel file 3*
Idents(Bcell) <- "cluster_name0.3"
one_vs_all_rna <- FindAllMarkers(Bcell, assay = "RNA", test.use = "MAST", only.pos = T, latent.vars = "SampleName", min.pct = 0.1)
one_vs_all_rna <- one_vs_all_rna %>%
  filter(p_val_adj < 0.05)
write.csv(one_vs_all_rna, file = "output/6_Bcell_clusters/dge/Bcell_one_vs_all_rna_cluster_name0.3_MAST.csv")

# Dotplot for surface marker expression
markers <- c("CD19","CD20","CD27","CD24","CD11c","CD38","CXCR5","IgD","IgM")

dotp <- DotPlot(Bcell, features = markers, assay = "ADT", scale = T, group.by = "cluster_name0.3"); dotp
df.dot <- dotp$data
# *Figure S6A* 
ggplot(df.dot, aes(y = id, x = features.plot)) +
  geom_point(aes(color = avg.exp.scaled, size= pct.exp)) +
  scale_color_gradient2(low = "navy", mid = 'white', high = 'red', name ="Average\nExpression" )  +
  scale_radius(limits = c(0,100), range = c(0,5), name = "Percent\nExpression") +
  theme(axis.text.x = element_text(size = 10, face='bold', angle = 270, hjust = 0, vjust=0.5),
        axis.text.y = element_text(size = 10, face='bold'),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray", linetype = 2),
        legend.key.height = unit(0.1, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.box.spacing = unit(0.05, "in"),
        legend.key.spacing = unit(0.01, "in"))+
  xlab("") +ylab("")
ggsave(file = "output/6_Bcell_clusters/DotPlot_Bcell_cluster_ADT_markers.png",
       height = 3, width = 4.5, units = 'in', bg = 'white')

## Heatmap for DEGs
dge <- read.csv(file = "output/6_Bcell_clusters/dge/Bcell_one_vs_all_rna_MAST.csv", row.names = 1)

dge$p_val_adj[which(dge$p_val_adj == 0)] <- 1e-300
dge$rank <- dge$avg_log2FC * (-1) * log10(dge$p_val_adj)

markers <- dge %>%  group_by(cluster) %>% top_n(10, rank)  %>% dplyr::arrange(desc(rank), .by_group = T)
markers2 <- markers[!duplicated(markers$gene),]
top.genes <- markers2$gene

dotp <- DotPlot(Bcell, features = top.genes, assay = "RNA", scale = T, group.by = "cluster_name0.3");dotp
df_mat <- dotp$data
df_mat <- df_mat %>% select(features.plot, id, avg.exp.scaled)
df_mat <- df_mat %>% pivot_wider(names_from = "features.plot", values_from = "avg.exp.scaled")

hm.mat <- df_mat[,-1] %>% as.matrix() 
rownames(hm.mat) <- df_mat$id
pal_cols <- colorRampPalette(c("navy","white","red")) 
newnames_row <- lapply(
  rownames(hm.mat) ,
  function(x) bquote(bold(.(x))))

newnames_cols <- lapply(
  colnames(hm.mat),
  function(x) bquote(bold(.(x))))

# *Figure S6C* 
png(filename = "output/6_Bcell_clusters/Heatmap_Bcell_cluster_name0.3.png", height =3.5, width =7.5, res=500, units = "in") 
pheatmap::pheatmap(hm.mat, 
                   color = pal_cols(100), 
                   fontsize = 12, fontsize_row = 12, fontsize_number=12, 
                   labels_row = as.expression(newnames_row),
                   labels_col = as.expression(newnames_cols),
                   border_color = "white") 
dev.off()





# *** end of Bcell clusters ***