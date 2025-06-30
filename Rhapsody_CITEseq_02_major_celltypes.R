#### major cell type clustering 

library(Seurat)
library(ggplot2)
library(dplyr)
library(Signac)
library(tidyr)
library(SingleR)
library(celldex)
library(data.table)

rm(list=ls())

#load the data
pbmc <- readRDS(file = "output/1_integration/pbmc_integrated.rds")

## use 18 ADTs for major cell type clustering 
markers <- c("CCR7","CD11b","CD11c","CD16","CD19","CD20","CD2","CD36","CD38",
             "CD45RA","CD45RO","CD56","CD9","HLA.DR","CD3", "CD4", "CD8", "CD14")
DefaultAssay(pbmc) <- "ADT"
pbmc.subset <- pbmc[rownames(pbmc) %in% markers, ]
adt_assay <- GetAssay(pbmc.subset, assay = "integratedADT")
pbmc[["ADT18"]] <- adt_assay

DefaultAssay(pbmc) <- "ADT18"
pbmc <- pbmc %>% FindVariableFeatures()  %>% ScaleData() %>% RunPCA() 
#ElbowPlot(pbmc, ndims = 25)
pbmc <- pbmc %>% RunUMAP(dims = 1:13) %>% FindNeighbors() %>% FindClusters(resolution=c(0.2, 0.3, 0.4),
                                                                           cluster.name = c("cluster_0.2", "cluster_0.3", "cluster_0.4"), random.seed=42) %>% identity()

pbmc$cluster_0.2 <- as.factor(as.numeric(as.character(pbmc$cluster_0.2)) + 1)
pbmc$cluster_0.3 <- as.factor(as.numeric(as.character(pbmc$cluster_0.3)) + 1)
pbmc$cluster_0.4 <- as.factor(as.numeric(as.character(pbmc$cluster_0.4)) + 1)

cluster.name = c("cluster_0.2", "cluster_0.3","cluster_0.4")

for (i in cluster.name){
  DimPlot(pbmc, group.by = i, label = T, repel = T, pt.size = 1, label.size = 7, )
  ggsave(file = paste0("output/2_major_clusters/umap_comparison1/DimPlot_pbmc_", i, ".png"),
         height = 7, width = 8, units = "in", bg = 'white')
}

## *Figure S2B*
DimPlot(pbmc, group.by ="cluster_0.3",label = T, repel=TRUE, label.size = 5) + 
  ggsci::scale_color_d3(palette="category20") +ggtitle("")
ggsave(file = "output/2_major_clusters/UMAP_major_clusters_initial.png", 
       width=10, height=8, units = "in", bg = 'white')

saveRDS(pbmc, file = "output/2_major_clusters/pbmc_major_celltype_initial.rds")

# Feature plots
for(i in markers){
  FeaturePlot(pbmc, features =i, cols =  c(alpha("lightgrey",0.4), "blue"), min.cutoff =0 ,max.cutoff =5)
  ggsave(file = paste0("output/2_major_clusters/FeaturePlot_1/FeaturePlot_", i, ".png"), 
         width=6, height=6, units = "in", bg = 'white')
}

pbmc <- subset(pbmc, subset = cluster_0.3 %in% c("9", "12", "14"), invert = T) 

DefaultAssay(pbmc) <- "ADT18"
pbmc <- pbmc %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
#ElbowPlot(pbmc)
pbmc <- pbmc %>%
  RunUMAP(dims = 1:13) %>% FindNeighbors() %>% FindClusters(resolution=c(0.2, 0.25, 0.3),
                                                            cluster.name = c("cluster_0.2","cluster_0.25", "cluster_0.3"), random.seed=42) %>% identity()

pbmc$cluster_0.2 <- as.factor(as.numeric(as.character(pbmc$cluster_0.2)) + 1)
pbmc$cluster_0.25 <- as.factor(as.numeric(as.character(pbmc$cluster_0.25)) + 1)
pbmc$cluster_0.3 <- as.factor(as.numeric(as.character(pbmc$cluster_0.3)) + 1)

cluster.name = c("cluster_0.2", "cluster_0.25","cluster_0.3")

for (i in cluster.name){
  DimPlot(pbmc, group.by = i, label = T, repel = T, pt.size = 1, label.size = 7)
  ggsave(file = paste0("output/2_major_clusters/umap_comparison2/DimPlot_pbmc_", i, ".png"),
         height = 7, width = 8, units = "in", bg = 'white')
}

DimPlot(pbmc, group.by ="cluster_0.2",label = T, repel=TRUE, label.size = 5) + 
  ggsci::scale_color_d3(palette="category20") +ggtitle("")
ggsave(file = "output/2_major_clusters/UMAP_major_clusters_cleaned.png", 
       width=10, height=8, units = "in", bg = 'white')

# *Figure S2B* 
for(i in markers){
  FeaturePlot(pbmc, features =i, cols =  c(alpha("lightgrey",0.2), "blue"),  
                 min.cutoff =0 ,max.cutoff =5) +
    theme(plot.title = element_text(size = 24)) +
    NoAxes()
  ggsave(file = paste0("output/2_major_clusters/FeaturePlot_2/FeaturePlot_", i, ".png"), 
         width=4.5, height=3.7, units = "in", bg = 'white')
}

# save RDS file 
saveRDS(pbmc, file = "output/2_major_clusters/pbmc_major_celltype_cleaned.rds")

## Single R annotation
monaco.ref <- MonacoImmuneData()

# convert Seurat object to SCE object
sce <- as.SingleCellExperiment(DietSeurat(pbmc), assay = 'RNA')
pred.monaco.main <- SingleR(sce, 
                            ref = monaco.ref, labels = monaco.ref$label.main, de.method = 'classic')
#table(pred.monaco.main$pruned.labels)
df_monaco <- data.frame(monaco.main =pred.monaco.main$pruned.labels)
rownames(df_monaco) <- pred.monaco.main@rownames
cell_list <- c("Monocytes","Dendritic cells", "B cells", "CD8+ T cells", "NK cells","T cells","CD4+ T cells" )
df_monaco$monaco.main <- ifelse(!df_monaco$monaco.main %in% cell_list, "Others", df_monaco$monaco.main)

# add cell type as metadata 
pbmc <- AddMetaData(pbmc, metadata = df_monaco)

# *Figure S2C*
DimPlot(pbmc, group.by ="monaco.main", label = F, pt.size =1) + 
  scale_color_brewer(palette = "Set2") +
  theme(legend.text = element_text(size = 20),
        axis.title = element_text(size = 20))+
  guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("")
ggsave(file = "output/2_major_clusters/UMAP_monaco_SingleR.png", 
       width=10, height=8, units = "in", bg = 'white')

# Add major cell type label 
pbmc$major_celltype <- ""
pbmc$major_celltype[pbmc$cluster_0.2=="1"] <- "Classical mono"
pbmc$major_celltype[pbmc$cluster_0.2=="2"] <- "Memory CD8T"
pbmc$major_celltype[pbmc$cluster_0.2=="3"] <- "Memory CD4T"
pbmc$major_celltype[pbmc$cluster_0.2=="4"] <- "Natural killer"
pbmc$major_celltype[pbmc$cluster_0.2=="5"] <- "Nonclassical mono" 
pbmc$major_celltype[pbmc$cluster_0.2=="6"] <- "Naive CD4T"
pbmc$major_celltype[pbmc$cluster_0.2=="7"] <- "B"
pbmc$major_celltype[pbmc$cluster_0.2=="8"] <- "Naive CD8T"
pbmc$major_celltype[pbmc$cluster_0.2=="9"] <- "Dendritic"

# *Figure S2B*
pbmc$major_celltype <-factor(pbmc$major_celltype, 
                             levels = c("Classical mono","Memory CD8T","Memory CD4T","Natural killer","Nonclassical mono","Naive CD4T","B","Naive CD8T","Dendritic"))

DimPlot(pbmc, group.by ="major_celltype") + ggsci::scale_color_d3(palette="category20") + 
  theme(legend.text = element_text(size = 20),
        axis.title = element_text(size = 20))+
  guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("")
ggsave(file = "output/2_major_clusters/UMAP_major_celltype.png", 
       width=10, height=8, units = "in", bg = 'white')

# *Figure 1B*
DimPlot(pbmc, group.by ="major_celltype", label = T, label.size = 5, repel = T) + 
  NoAxes()+ ggsci::scale_color_d3(palette="category20") +
  theme(legend.text = element_text(size = 20),
        plot.title = element_text(hjust = 1, size = 16))+
  guides(color = guide_legend(override.aes = list(size=5))) + 
  ggtitle("Major cell types (66,995 cells)")
ggsave(file = "output/2_major_clusters/UMAP_major_celltype2.png",  
       width=7.5, height=5, units = "in", bg = 'white')

## save RDS file 
saveRDS(pbmc, file = "output/2_major_clusters/pbmc_major_celltype_annotated.rds")

## DGE *Excel file 3* 
dge <- FindAllMarkers(pbmc, assay = "RNA", slot = "data", 
                      group.by = "major_celltype", only.pos = T, test.use = "MAST", latent.vars = "SampleName") 
write.csv(dge, file = "output/2_major_clusters/major_celltype_dge_MAST.csv")

# Heatmap (Fig S3)
dge$p_val_adj[which(dge$p_val_adj==0)] <- 1e-300
dge$rank <- dge$avg_log2FC*(-1)*log10(dge$p_val_adj)

top.genes <- dge %>% group_by(cluster) %>% top_n(5, rank)
top.genes <- top.genes$gene %>% unique()
col <- ggsci::pal_d3("category10")(9)
names(col) <- unique(sort(pbmc$major_celltype))

DoHeatmap(subset(pbmc, downsample = min(table(pbmc$major_celltype))),
          group.by = 'major_celltype', 
          features = top.genes, assay = 'RNA', slot = 'scale.data',
          group.colors = col, label = T) +
  scale_fill_gradientn(colors = c("navy","white","red")) +
  theme(legend.text = element_text(size = 14, face = 'bold'),
        axis.text.y = element_text(size = 14, face = 'bold'))
ggsave(file = "output/2_major_clusters/DoHM_pbmc_cluster_gene_markers.png",
       height = 10, width = 10, units = 'in', bg = 'white')

## Heatmap for ADT expression 
pbmc_df <- pbmc@assays$ADT18@data %>% t() %>% as.data.frame() 
pbmc_df$cluster <- paste0("C", pbmc$cluster_0.2,"_", pbmc$major_celltype)

avg <- pbmc_df %>% group_by(cluster) %>% summarise_all(funs(mean)) %>% data.table() # average expression 
cluster <- avg$cluster
avg$cluster <- NULL
avg <- avg %>% t() %>% as.data.frame()
colnames(avg) <- cluster
avg <- avg %>% as.matrix()

avg <- avg %>% t() 
newnames_row <- lapply(
  rownames(avg) ,
  function(x) bquote(bold(.(x))))
newnames_cols <- lapply(
  colnames(avg),
  function(x) bquote(bold(.(x))))

pal_cols <- colorRampPalette(c("navy","white","red"))

## *Figure S2A* 
png(filename = "output/2_major_clusters/heatmap_major_celltype_ADT.png", height =5, width =12, res=350, units = "in")
pheatmap::pheatmap(avg, 
                   color = pal_cols(100), 
                   fontsize = 18, fontsize_col = 18, fontsize_row = 18, fontsize_number=18, 
                   labels_row = as.expression(newnames_row),
                   labels_col = as.expression(newnames_cols),
                   cluster_cols = FALSE, 
                   border_color = "white",
                   ) 
dev.off()

## Cell counts for major cell types per patient (*Excel File 2*) 
cell_count <- pbmc@meta.data %>% group_by(major_celltype, SampleName) %>% summarise(num=n())
cell_count <- cell_count %>% pivot_wider(names_from = "major_celltype", values_from = "num", values_fill = 0)
write.csv(cell_count, file = "output/2_major_clusters/cell_count_major_celltype_per_patient.csv",
          row.names = F)

## cell composition analysis (*Fig 2C*)
library(sccomp)

pbmc$SampleType <- factor(pbmc$SampleType, levels = c("Normal HTx", "CAV"))
set.seed(47)
sccomp_result <- pbmc %>% 
  sccomp_estimate(formula_composition = ~SampleType,
                  formula_variability = ~SampleType,
                  .sample = SampleName, 
                  .cell_group = major_celltype,
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
  ylab("")+ xlab("Log (odds ratio)") +
  ggtitle("") 
ggsave(file = "/output/2_major_clusters/pbmc_cell_composition_sccomp.png", 
       width=3, height=2.2, units="in",bg = 'white') 


# *** end of the major cell type clustering ***
