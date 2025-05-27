#### Natural killer cell sub-clusters 

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

rm(list=ls())

path = "~/Documents/UCSD_CAV/CITEseq"

#load the data
pbmc <- readRDS(file = paste0(path, "/output/2_major_clusters/pbmc_major_celltype_annotated.rds"))

# subset NK cells
NK <- subset(pbmc, major_celltype %in% c("Natural_killer"))

## NK cell unsupervised clustering (WNN)
DefaultAssay(NK) <- 'integrated'
NK <- NK %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'pca')

DefaultAssay(NK) <- 'integratedADT'
NK <- NK %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'apca')
NK <- FindMultiModalNeighbors(NK, reduction.list = list("pca", "apca"), dims.list = list(1:10, 1:10), 
                                modality.weight.name = "RNA.weight")
NK <- RunUMAP(NK, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
  FindClusters(graph.name="wsnn", 
               resolution = c(0.25, 0.3, 0.35), 
               cluster.name = c("cluster_0.25", "cluster_0.3","cluster_0.35"),  
               random.seed=42) %>% identity()
NK$cluster_0.25<- as.factor(as.numeric(as.character(NK$cluster_0.25)) + 1)
NK$cluster_0.3<- as.factor(as.numeric(as.character(NK$cluster_0.3)) + 1)
NK$cluster_0.35<- as.factor(as.numeric(as.character(NK$cluster_0.35)) + 1)

cluster.name = c("cluster_0.25", "cluster_0.3","cluster_0.35")

for (i in cluster.name){
  DimPlot(NK, group.by = i, label = T, repel = F, pt.size = 0.3, label.size = 4) +
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 8, face = 'bold'),
          title = element_text(size = 8))+
    guides(color = guide_legend(override.aes = list(size=3)))
  ggsave(file = paste0(path,"/output/7_NK_clusters/umap_comparison/DimPlot_NK_", i, ".png"),
         height = 3, width = 3.5, units = "in", bg = 'white')
}

# FeaturePlot 
features <- rownames(NK@assays$ADT)
for(i in features){
  FeaturePlot(NK, features =i, cols =  c(alpha("lightgrey",0.4), "blue"), pt.size =1, min.cutoff =1 ,max.cutoff =5)
  ggsave(file = paste0(path, "/output/7_NK_clusters/FeaturePlot_1/FeaturePlot_", i, ".png"), 
         width=7, height=6, units = "in", bg = 'white')
}

# Use "cluster_0.3"
# C5 is a suspected doublet (CD3+) 
# C7 is unidentifiable  
table(NK@meta.data$cluster_0.3)
NK <- subset(NK, cluster_0.3 == "5", invert = T)
NK <- subset(NK, cluster_0.3 == "7", invert = T)
NK$cluster_0.3 <- gsub("6","5", NK$cluster_0.3)

# add cluster names as meta data  
cluster_name = c("1:CD161+mature","2:mature","3:immature","4:CD14+mature","5:proliferating")
cluster_df <- data.frame(cluster_0.3 = as.character(1:5), cluster_name0.3 = cluster_name)
df <- NK@meta.data 
df$cell_id <- rownames(df)
df <- df %>% merge(cluster_df, by.x = "cluster_0.3", by.y = "cluster_0.3")
rownames(df) <- df$cell_id
NK <- AddMetaData(NK, df)

# *Save RDS file* 
saveRDS(NK, file = paste0(path, "/output/7_NK_clusters/NK_clusters.rds"))

# *Figure 2D*
DimPlot(NK, group.by ="cluster_name0.3", label = F, repel = F, pt.size =0.03, 
        label.size = 3) + 
  scale_color_brewer(palette ="Set3")+
  NoLegend() +NoAxes() + ggtitle("")
ggsave(file = paste0(path,"/output/7_NK_clusters/DimPlot_NK_cluster_name0.3_nolegend.png"),
       height = 1.7, width = 1.7, units = "in", bg = 'white')

DimPlot(NK, group.by ="cluster_name0.3", label = F, pt.size =0.3) + 
  scale_color_brewer(palette ="Set3")+
  NoAxes() + ggtitle("") +
  theme(legend.text = element_text(size = 12, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(file = paste0(path,"/output/7_NK_clusters/DimPlot_NK_cluster_name0.3.png"),
       height = 4, width = 6, units = "in", bg = 'white')

# *Figure 2D*
DimPlot(NK, group.by ="SampleType", label = F, pt.size =0.05) + 
  scale_color_manual(values = c("maroon", "gray80")) +
  NoLegend() +NoAxes() + ggtitle("")
ggsave(file = paste0(path,"/output/7_NK_clusters/DimPlot_NK_SampleType.png"),
       height = 1.7, width = 1.7, units = "in", bg = 'white')

# *Fig S6B*
features <- rownames(NK@assays$ADT)
Idents(NK) <- "cluster_0.3"
for(i in features){
  FeaturePlot(NK, features =i, cols =  c(alpha("lightgrey",0.4), "blue"),label = T,  
                 label.size = 8, pt.size =0.3, min.cutoff =1 ,max.cutoff =5) +
    theme(plot.title = element_text(size = 28),
          legend.position = "inside",
          legend.position.inside = c(0.05, 0.8)) +
    NoAxes()
  ggsave(file = paste0(path, "/output/7_NK_clusters/FeaturePlot_2/FeaturePlot_", i, ".png"), 
         width=4.3, height=3.7, units = "in", bg = 'white')
}

## Cluster cell counts per sample (*Excel file 2*)
df_count <- NK@meta.data %>% group_by(cluster_name0.3, SampleName) %>% 
  summarise(num=n()) %>% pivot_wider(names_from = "cluster_name0.3", values_from = "num", values_fill = 0)
write.csv(df_count, file = paste0(path,"/output/7_NK_clusters/NK_cluster_counts_patient.csv"), row.names = F) 

## Cluster proportion per patient: CAV vs Normal HTx 
df_prop <- NK@meta.data %>% group_by(cluster_name0.3, SampleName, SampleType) %>% summarise(num=n())
df_prop <- df_prop %>% pivot_wider(names_from = 'cluster_name0.3', values_from = 'num', values_fill = 0)

mat_prop <- df_prop[, 3:ncol(df_prop)] %>% as.matrix()
mat_prop <- round(prop.table(mat_prop, margin = 1) *100, 2)
df_prop[,3:ncol(df_prop)] <- mat_prop 

write.csv(df_prop, file = paste0(path,"/output/7_NK_clusters/NK_cluster_proportion_patient.csv"), row.names = F) 

## DGE: MAST
Idents(NK) <- "cluster_0.3"
one_vs_all_rna <- FindAllMarkers(NK, assay = "RNA", test.use = "MAST", only.pos = F, latent.vars = "SampleName", min.pct = 0.1)
write.csv(one_vs_all_rna, file = paste0(path,"/output/7_NK_clusters/dge/NK_one_vs_all_rna_MAST.csv"))

## DGE: MAST *Excel file 3*
Idents(NK) <- "cluster_name0.3"
one_vs_all_rna <- FindAllMarkers(NK, assay = "RNA", test.use = "MAST", only.pos = T, latent.vars = "SampleName", min.pct = 0.1)
one_vs_all_rna <- one_vs_all_rna %>%
  filter(p_val_adj < 0.05)
write.csv(one_vs_all_rna, file = paste0(path,"/output/7_NK_clusters/dge/NK_one_vs_all_rna_cluster_name0.3_MAST.csv"))

# Dotplot for surface marker expression
markers <- c("CD16", "CD56","CCR7", "CD62L","CD25","CD27", "CD45RA", "CD45RO",
             "CD14","CD36", "CD192","CD161", "CD38", "HLA.DR", "CD69")

dotp <- DotPlot(NK, features = markers, assay = "ADT", scale = T, 
                group.by = "cluster_name0.3"); dotp

df.dot <- dotp$data

# *Figure S6A*
ggplot(df.dot, aes(x = id, y = features.plot)) +
  geom_point(aes(color = avg.exp.scaled, size= pct.exp)) +
  scale_color_gradient2(low = "navy", mid = 'white', high = 'red', name ="Average\nExpression" )  +
  scale_radius(limits = c(0,100), range = c(0,5), name = "Percent\nExpression") +
  theme(axis.text.x = element_text(size = 12, face='bold', angle = 270, hjust = 0, vjust=0.5),
        axis.text.y = element_text(size = 12, face='bold'),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray", linetype = 2),
        legend.position = "top", legend.box = "vertical",
        legend.key.height = unit(0.1, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.box.spacing = unit(0.05, "in"),
        legend.key.spacing = unit(0.01, "in"))+
  xlab("") +ylab("")
ggsave(file = paste0(path,"/output/7_NK_clusters/DotPlot_NK_cluster_ADT_markers.png"),
       height = 5.5, width = 3.3, units = 'in', bg = 'white')

## Heatmap for DGE
dge <- read.csv(file = paste0(path,"/output/7_NK_clusters/dge/NK_one_vs_all_rna_MAST.csv"), row.names = 1)

dge$p_val_adj[which(dge$p_val_adj == 0)] <- 1e-300
dge$rank <- dge$avg_log2FC * (-1) * log10(dge$p_val_adj)

markers <- dge %>%  group_by(cluster) %>% top_n(8, rank)  %>% dplyr::arrange(desc(rank), .by_group = T)
markers2 <- markers[!duplicated(markers$gene),]
top.genes <- markers2$gene

dotp <- DotPlot(NK, features = top.genes, assay = "RNA", scale = T, group.by = "cluster_name0.3");dotp

df_mat <- dotp$data %>% select(features.plot, id, avg.exp.scaled)
df_mat <- df_mat %>% pivot_wider(names_from = "features.plot", values_from = "avg.exp.scaled")
annotation_col <- data.frame(cluster = df_mat$id)
annotation_col <- data.frame(cluster = factor(sort(annotation_col$cluster), labels = sort(annotation_col$cluster)))
rownames(annotation_col) <- annotation_col$cluster
col <-scales::pal_brewer(palette ="Set3")(5)
cluster <-col
names(cluster) <- sort(annotation_col$cluster)
annotation_colors <- list(cluster = cluster)

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
png(filename = paste0(path, "/output/7_NK_clusters/Heatmap_NK_cluster_name0.3_MAST.png"), height =3.5, width =11, res=500, units = "in") 
pheatmap::pheatmap(hm.mat, 
                   color = pal_cols(100), cellwidth = 13, 
                   fontsize = 14, fontsize_row = 14, fontsize_number=14, 
                   labels_col = as.expression(newnames_cols),
                   labels_row = as.expression(newnames_row),
                   annotation_names_row = F,
                   show_rownames = FALSE,
                   annotation_legend = T,
                   border_color = "white",
                   annotation_row = annotation_col,
                   annotation_colors = annotation_colors) 
dev.off()

## Differential composition analysis: sccomp 
library(sccomp)

NK$SampleType <- factor(NK$SampleType, levels = c("Normal_HTx","CAV"))
set.seed(47)
sccomp_result <- NK %>% 
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
        legend.position = 'none') +
  scale_color_manual(values = col) +
  ylab("")+ xlab("Log odds") +
  guides(color = guide_legend(title  = "FDR < significance threshold (0.025)")) +
  ggtitle("") 
ggsave(file = paste0(path,"/output/7_NK_clusters/NK_cell_composition_sccomp.png"), 
       width=3.5, height=1.8, units="in",bg = 'white')

# *** end of NK clusters ***