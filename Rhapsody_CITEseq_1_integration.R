#### 1. Normalization & batch effect corrected by Seurat integration (CCA) 
# data includes 3 disease groups (CONTROL, CAV, NGD) 
# remove NGD data in the beginning  

library(Seurat) #v5.2.0
library(ggplot2) #v3.5.1

rm(list=ls())

path = "~/Documents/UCSD_CAV/CITEseq"

#Read the data
pbmc.rna <- as.sparse(read.csv(paste0(path, "/raw_data/Genes_raw.csv"), sep=",", header = TRUE, row.names = 1))
pbmc.rna <- CollapseSpeciesExpressionMatrix(pbmc.rna)
pbmc.rna <- t(as.matrix(pbmc.rna))

#Seurat Object Setup 
pbmc <- CreateSeuratObject(counts=pbmc.rna)

#Adding metadata to the seurat object 
meta <- read.csv(paste0(path, "/raw_data/Metadata.csv"), row.names = 1)
pbmc <- AddMetaData(pbmc, metadata = meta)

# remove NGD 
pbmc <- subset(pbmc, subset = SampleType == "NGD", invert = T)

pbmc.list <- SplitObject(pbmc, split.by ="Library")
# LogNormalization of RNA data
pbmc.list <- lapply(X=pbmc.list, FUN = function(x){
  x <- NormalizeData(x) 
  x <- FindVariableFeatures(x, selection.method = "vst")
})

# RNA Integration
features <- SelectIntegrationFeatures(object.list = pbmc.list)
immune.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features) # ~30min
immune.combined <- IntegrateData(anchorset = immune.anchors)
immune.combined <- JoinLayers(immune.combined, assay = "RNA")

# umap: visualize integration
immune.combined <- ScaleData(immune.combined)
immune.combined <- RunPCA(immune.combined, npcs = 30)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)

# save DimPlot for integrated libraries for RNA data (*Figure S1A*)
DimPlot(immune.combined, reduction = "umap", group.by = "Library") + ggtitle("Integrated:RNA")
ggsave(file = paste0(path, "/output/1_integration/UMAP_Library_Integrated_RNA.png"), width=10, height=8, units = "in", bg = 'white')

# save RDS file 
saveRDS(immune.combined, file = paste0(path, "/output/1_integration/immune.combined_rna.rds"))

## ADT integration
pbmc.adt <- as.sparse(read.csv(paste0(path, "/raw_data/Ab_raw.csv"), sep = ",", header = TRUE, row.names = 1))
pbmc.adt <- CollapseSpeciesExpressionMatrix(pbmc.adt)
pbmc.adt <- t(as.matrix(pbmc.adt))
#all.equal(colnames(pbmc.rna), colnames(pbmc.adt))  

# change ADT label
adtname <- rownames(pbmc.adt)
adtname <- strsplit(adtname, "[.]")
adtname <- do.call(rbind, adtname)
adtname[,1] <- gsub("HLA","HLA.DR", adtname[,1])
rownames(pbmc.adt) <- adtname[,1]

# create ADT assay 
objectadt<- CreateSeuratObject(counts=pbmc.adt)
objectadt <- AddMetaData(objectadt, metadata = meta)
# remove NGD
objectadt <- subset(objectadt, subset = SampleType == "NGD", invert = T)

objectadt.list <- SplitObject(objectadt, split.by ="Library")

# CLR normalization for ADT data
objectadt.list <- lapply(X=objectadt.list, FUN = function(x){
  x <- NormalizeData(x,normalization.method="CLR", margin=2) 
  x <- FindVariableFeatures(x, selection.method = "vst")
})

# integration
features <- SelectIntegrationFeatures(object.list = objectadt.list)
immune.anchors <- FindIntegrationAnchors(object.list = objectadt.list, anchor.features = features) # ~15min
immune.combined_adt <- IntegrateData(anchorset = immune.anchors)
immune.combined_adt <- JoinLayers(immune.combined_adt, assay = "RNA")

# umap: visualize integration
immune.combined_adt <-ScaleData(immune.combined_adt)
immune.combined_adt <- RunPCA(immune.combined_adt, npcs = 30)
immune.combined_adt <- RunUMAP(immune.combined_adt, reduction = "pca", dims = 1:30)

# save DimPlot for integrated libraries for ADT data (*Figure S1A*)
DimPlot(immune.combined_adt, reduction = "umap", group.by = "Library") + ggtitle("Integrated:ADT")
ggsave(file = paste0(path, "/output/1_integration/UMAP_Library_Integrated_ADT.png"), width=10, height=8, units = "in", bg = 'white')

# save RDS file 
saveRDS(immune.combined_adt, file = paste0(path, "/output/1_integration/immune.combined_adt.rds"))

# read previously created objects 
immune.combined <- readRDS(file = paste0(path, "/output/1_integration/immune.combined_rna.rds"))
immune.combined_adt <- readRDS(file = paste0(path, "/output/1_integration/immune.combined_adt.rds"))

# add ADT
adt_assay <- GetAssay(immune.combined_adt, assay = "RNA")
immune.combined[["ADT"]] <- adt_assay
integratedadt_assay <- GetAssay(immune.combined_adt, assay = "integrated")
immune.combined[["integratedADT"]] <- integratedadt_assay

# modify metadata 
immune.combined$SampleName <- gsub("CONTROL", "Normal_HTx", immune.combined$SampleName)
immune.combined$SampleType <- gsub("CONTROL", "Normal_HTx", immune.combined$SampleType)

saveRDS(immune.combined, file = paste0(path, "/output/1_integration/pbmc_integrated.rds"))


# *** end of the integration protocol***


