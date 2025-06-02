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

## *Figure S1B*
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

# *Figure S1C*
DimPlot(pbmc, group.by ="monaco.main", label = F, pt.size =1) + 
  scale_color_brewer(palette = "Set2") +
  theme(legend.text = element_text(size = 20),
        axis.title = element_text(size = 20))+
  guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("")
ggsave(file = "output/2_major_clusters/UMAP_monaco_SingleR.png", 
       width=10, height=8, units = "in", bg = 'white')

# Add major cell type label 
pbmc$major_celltype <- ""
pbmc$major_celltype[pbmc$cluster_0.2=="1"] <- "Classical_mono"
pbmc$major_celltype[pbmc$cluster_0.2=="2"] <- "Memory_CD8T"
pbmc$major_celltype[pbmc$cluster_0.2=="3"] <- "Memory_CD4T"
pbmc$major_celltype[pbmc$cluster_0.2=="4"] <- "Natural_killer"
pbmc$major_celltype[pbmc$cluster_0.2=="5"] <- "Nonclassical_mono" 
pbmc$major_celltype[pbmc$cluster_0.2=="6"] <- "Naive_CD4T"
pbmc$major_celltype[pbmc$cluster_0.2=="7"] <- "B"
pbmc$major_celltype[pbmc$cluster_0.2=="8"] <- "Naive_CD8T"
pbmc$major_celltype[pbmc$cluster_0.2=="9"] <- "Dendritic"

# *Figure S1B*
pbmc$major_celltype <-factor(pbmc$major_celltype, 
                             levels = c("Classical_mono","Memory_CD8T","Memory_CD4T","Natural_killer","Nonclassical_mono","Naive_CD4T","B","Naive_CD8T","Dendritic"))

DimPlot(pbmc, group.by ="major_celltype") + ggsci::scale_color_d3(palette="category20") + 
  theme(legend.text = element_text(size = 20),
        axis.title = element_text(size = 20))+
  guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("")
ggsave(file = "output/2_major_clusters/UMAP_major_celltype.png", 
       width=10, height=8, units = "in", bg = 'white')

# *Figure 1B*
DimPlot(pbmc, group.by ="major_celltype", label = T, label.size = 10, repel = T) + 
  NoLegend()+ NoAxes()+ ggsci::scale_color_d3(palette="category20") + ggtitle("") 
ggsave(file = "output/2_major_clusters/UMAP_major_celltype_Fig1B.png", 
       width=9, height=7, units = "in", bg = 'white')

## save RDS file 
saveRDS(pbmc, file = "output/2_major_clusters/pbmc_major_celltype_annotated.rds")

## DGE *Excel file 3* 
dge <- FindAllMarkers(pbmc, assay = "RNA", slot = "data", 
                      group.by = "major_celltype", only.pos = T, test.use = "MAST", latent.vars = "SampleName") 
write.csv(dge, file = "output/2_major_clusters/major_celltype_dge_MAST.csv")

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

pal_cols <- colorRampPalette(c("blue","white","red"))

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

# *** end of the major cell type clustering ***
