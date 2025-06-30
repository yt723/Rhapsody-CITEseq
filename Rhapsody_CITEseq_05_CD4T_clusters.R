#### CD4 T cell sub-clusters 

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

rm(list=ls())

#load the data
pbmc <- readRDS(file = "output/2_major_clusters/pbmc_major_celltype_annotated.rds")

# subset CD4T cells
CD4T <- subset(pbmc, major_celltype %in% c("Naive CD4T", "Memory CD4T"))

## CD4T cell unsupervised clustering (WNN)
DefaultAssay(CD4T) <- 'integrated'
CD4T <- CD4T %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'pca')

DefaultAssay(CD4T) <- 'integratedADT'
CD4T <- CD4T %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'apca')
CD4T <- FindMultiModalNeighbors(CD4T, reduction.list = list("pca", "apca"), dims.list = list(1:10, 1:10), 
                                modality.weight.name = "RNA.weight")

CD4T <- RunUMAP(CD4T, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
  FindClusters(graph.name="wsnn", 
               resolution = c(0.25, 0.3, 0.35), 
               cluster.name = c("cluster_0.25", "cluster_0.3","cluster_0.35"), 
               random.seed=42) %>% identity()
CD4T$cluster_0.25<- as.factor(as.numeric(as.character(CD4T$cluster_0.25)) + 1)
CD4T$cluster_0.3<- as.factor(as.numeric(as.character(CD4T$cluster_0.3)) + 1)
CD4T$cluster_0.35<- as.factor(as.numeric(as.character(CD4T$cluster_0.35)) + 1)

cluster.name = c("cluster_0.25", "cluster_0.3","cluster_0.35")

for (i in cluster.name){
  dimp <- DimPlot(CD4T, group.by = i, label = T, repel = F, pt.size = 0.3, label.size = 4) +
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 8, face = 'bold'),
          title = element_text(size = 8))+
    guides(color = guide_legend(override.aes = list(size=3)))
  ggsave(dimp, file = paste0("output/5_CD4T_clusters/umap_comparison/DimPlot_CD4T_", i, ".png"),
         height = 3, width = 3.5, units = "in", bg = 'white')
}

# Use "cluster_0.3"
table(CD4T@meta.data$cluster_0.3)
# C10 is too small (n=29) for analysis
CD4T <- subset(CD4T, cluster_0.3 == "10", invert = T)

## comparison with SingleR annotations 
library(SingleR)
library(celldex)

monaco.ref <- MonacoImmuneData()
celltype <- c("Follicular helper T cells", "T regulatory cells", "Th1 cells", "Th1/Th17 cells", "Th17 cells",
              "Th2 cells", "Naive CD4 T cells", "Terminal effector CD4 T cells")
ref.subset <- monaco.ref[, colData(monaco.ref)$label.fine %in% celltype]

# convert Seurat to SCE object
sce <- as.SingleCellExperiment(DietSeurat(CD4T), assay = 'RNA')

# run single R
pred.fine <- SingleR(sce, ref = ref.subset, labels = ref.subset$label.fine, de.method = 'classic')

table(pred.fine$pruned.labels)
df_monaco.fine <- data.frame(monaco.fine =pred.fine$pruned.labels)
rownames(df_monaco.fine) <- pred.fine@rownames
CD4T <- AddMetaData(CD4T, metadata = df_monaco.fine)

# *Figure S5D*
DimPlot(CD4T, group.by = "monaco.fine", pt.size = 0.3) + 
  scale_color_brewer(palette = "Set2") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 14),
        legend.key.height = unit(0.15, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.key.spacing = unit(0.05, "in")) +
  guides(color = guide_legend(override.aes = list(size=4)))+ ggtitle("")
ggsave(file = "output/5_CD4T_clusters/DimPlot_CD4T_monaco.png",
       height = 3, width = 5.5, units = "in", bg = 'white')

# add cluster names as meta data  
cluster_name = c("1:naive","2:EM1","3:Tfh","4:cytotoxic","5:CM",
                 "6:EM2","7:Treg", "8:NK like", "9:CD38+HLADR+EM")
cluster_df <- data.frame(cluster_0.3 = as.character(1:9), cluster_name0.3 = cluster_name)
df <- CD4T@meta.data 
df$cell_id <- rownames(df)
df <- df %>% merge(cluster_df, by.x = "cluster_0.3", by.y = "cluster_0.3")
rownames(df) <- df$cell_id
CD4T <- AddMetaData(CD4T, df)

# *Save RDS file* 
saveRDS(CD4T, file = "output/5_CD4T_clusters/CD4T_clusters.rds")

# *Figure 2B*
colset <- ggsci::pal_d3("category20")(20)[8:16]
DimPlot(CD4T, group.by ="cluster_name0.3", label = F, pt.size =0.3) + 
  scale_color_manual(values = colset)+
  NoAxes() + ggtitle("") +
  theme(legend.text = element_text(size = 12, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(file = "output/5_CD4T_clusters/DimPlot_CD4T_cluster_name0.3.png",
       height = 4, width = 6, units = "in", bg = 'white')

DimPlot(CD4T, group.by ="SampleType", label = F, pt.size =0.05) + 
  scale_color_manual(values = c("maroon", "gray80")) +
  theme(legend.position = 'top',
        legend.spacing.x = unit(0.1, "in"),
        legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=5), nrow = 1, title.hjust=0.5)) +
  NoAxes() + ggtitle("")
ggsave(file = "output/5_CD4T_clusters/DimPlot_CD4T_SampleType.png",
       height = 3.5, width = 3.2, units = "in", bg = 'white')

# FeaturePlot (*Figure S5B*)
features <- rownames(CD4T@assays$ADT)
Idents(CD4T) <- "cluster_0.3"
for(i in features){
  FeaturePlot(CD4T, features =i, cols =  c(alpha("lightgrey",0.4), "blue"),label = T,  
                 label.size = 8, pt.size =0.3, min.cutoff =1 ,max.cutoff =5) +
    theme(plot.title = element_text(size = 28),
          legend.position = "inside",
          legend.position.inside = c(0.05, 0.8)) +
    NoAxes()
  ggsave(file = paste0("output/5_CD4T_clusters/FeaturePlot_1/FeaturePlot_", i, ".png"), 
         width=4.3, height=3.7, units = "in", bg = 'white')
}

## Cluster cell counts per sample (*Excel file 2*)
df_count <- CD4T@meta.data %>% group_by(cluster_name0.3, SampleName) %>% 
  summarise(num=n()) %>% pivot_wider(names_from = "cluster_name0.3", values_from = "num", values_fill = 0)
write.csv(df_count, file = "output/5_CD4T_clusters/CD4T_cluster_counts_patient.csv", row.names = F) 

## Cluster proportion per patient: CAV vs Normal HTx 
df_prop <- CD4T@meta.data %>% group_by(cluster_name0.3, SampleName, SampleType) %>% dplyr::summarise(num=n())
df_prop <- df_prop %>% pivot_wider(names_from = 'cluster_name0.3', values_from = 'num', values_fill = 0)

mat_prop <- df_prop[, 3:ncol(df_prop)] %>% as.matrix()
mat_prop <- round(prop.table(mat_prop, margin = 1) *100, 2)
df_prop[,3:ncol(df_prop)] <- mat_prop 

write.csv(df_prop, file = "output/5_CD4T_clusters/CD4T_cluster_proportion_patient.csv", row.names = F) 

## Differential composition analysis
library(sccomp)

CD4T$SampleType <- factor(CD4T$SampleType, levels = c("Normal HTx","CAV"))
set.seed(47)
sccomp_result <- CD4T %>% 
  sccomp_estimate(formula_composition = ~SampleType,
                  formula_variability = ~SampleType,
                  .sample = SampleName, 
                  .cell_group = cluster_name0.3,
                  variational_inference = F,
                  bimodal_mean_variability_association = TRUE)

plots <- sccomp_result %>% sccomp_test() %>% plot()

col <- c("gray50","red")
names(col) <- c(FALSE, TRUE)
# *Figure 2B*
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
ggsave(file = "output/5_CD4T_clusters/CD4T_cell_composition_sccomp.png", 
       width=3.5, height=2, units="in",bg = 'white')


## DGE: MAST
Idents(CD4T) <- "cluster_0.3"
one_vs_all_rna <- FindAllMarkers(CD4T, assay = "RNA", test.use = "MAST", only.pos = F, latent.vars = "SampleName", min.pct = 0.1)
write.csv(one_vs_all_rna, file = "output/5_CD4T_clusters/dge/CD4T_one_vs_all_rna_MAST.csv")

# *Excel file 3*
Idents(CD4T) <- "cluster_name0.3"
one_vs_all_rna <- FindAllMarkers(CD4T, assay = "RNA", test.use = "MAST", only.pos = T, latent.vars = "SampleName", min.pct = 0.1)
one_vs_all_rna <- one_vs_all_rna %>%
  filter(p_val_adj < 0.05)
write.csv(one_vs_all_rna, file = "output/5_CD4T_clusters/dge/CD4T_one_vs_all_rna_cluster_name0.3_MAST.csv")

## Dotplot for surface marker expression
markers <- c("CCR7", "CD62L","CD45RA", "CD45RO","CD192","CXCR5", "CD69","CD25", "CD27","CD183","CD16", "CD56", "CD161", "CD38", "HLA.DR")

dotp <- DotPlot(CD4T, features = markers, assay = "ADT", scale = T, group.by = "cluster_name0.3"); dotp 
df.dot <- dotp$data

# *Figure S5A*
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
ggsave(file = "output/5_CD4T_clusters/DotPlot_CD4T_cluster_ADT_markers.png",
       height = 3.5, width = 6, units = 'in', bg = 'white')

## Heatmap for DEGs
dge <- read.csv(file = "output/5_CD4T_clusters/dge/CD4T_one_vs_all_rna_MAST.csv", row.names = 1)

dge$p_val_adj[which(dge$p_val_adj == 0)] <- 1e-300
dge$rank <- dge$avg_log2FC * (-1) * log10(dge$p_val_adj)

markers <- dge %>%  group_by(cluster) %>% top_n(8, rank)  %>% dplyr::arrange(desc(rank), .by_group = T)
markers2 <- markers[!duplicated(markers$gene),]
top.genes <- markers2$gene

dotp <- DotPlot(CD4T, features = top.genes, assay = "RNA", scale = T, group.by = "cluster_name0.3");dotp
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

# *Figure S5C* 
png(filename = "output/5_CD4T_clusters/Heatmap_CD4T_cluster_name0.3.png", height =4, width =14, res=500, units = "in") 
pheatmap::pheatmap(hm.mat, 
                   color = pal_cols(100), 
                   fontsize = 15, fontsize_row = 15, fontsize_number=15, 
                   labels_row = as.expression(newnames_row),
                   labels_col = as.expression(newnames_cols),
                   border_color = "white") 
dev.off()



# *** end of CD4T clusters ***

