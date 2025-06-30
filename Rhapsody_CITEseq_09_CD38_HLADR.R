#### CD38 and HLA.DR expression analysis 

library(Seurat)
library(ggplot2)
library(dplyr)

rm(list=ls())

#load data
pbmc <- readRDS(file = "output/2_major_clusters/pbmc_major_celltype_annotated.rds")
CD8T <- readRDS(file = "output/8_diffusion_map/CD8T_diffusion_map.rds")

# Biaxial plot to show thresholds for CD38 and HLA.DR: major cell types 
color <- ggsci::pal_d3("category10")(9)
names(color) <- unique(pbmc$major_celltype)[order(unique(pbmc$major_celltype))]

Lym <- subset(pbmc, major_celltype %in% c("Naive CD4T", "Memory CD4T", "Naive CD8T", "Memory CD8T", "Natural killer")) 
DefaultAssay(Lym) <- "integratedADT"
# *Figure S9A* 
FeatureScatter(Lym, feature1 ="CD3",feature2 ="CD38", slot = 'data',
                     group.by = "major_celltype", cols = color)+   
  theme(text = element_text(face = "bold", size = 22), axis.text.x=element_text(size=15), 
        axis.title = element_text(size=20,face="bold"), title = element_blank(), 
        axis.title.y.right = element_text(size =15),
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        axis.line = element_line(linewidth =1))+
  geom_density2d(col = 'black')+
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_hline(yintercept =1.3)+
  xlab("CD3")+ ylab("CD38")+
  ggtitle("") 
ggsave(file = "output/9_CD38_HLADR/CD3_vs_CD38_Biaxial_Plot_y1.3_pbmc.png",
       width=8, height=6, units = "in", bg = 'white')


Myl <- subset(pbmc, major_celltype %in% c("Classical mono", "Nonclassical mono", "B", "Dendritic")) 
DefaultAssay(Myl) <- "integratedADT"
# *Figure S9A* 
FeatureScatter(Myl, feature1 ="CD14",feature2 ="HLA.DR", slot = 'data',
                     group.by = "major_celltype", cols = color)+   
  theme(text = element_text(face = "bold", size = 22), axis.text.x=element_text(size=15), 
        axis.title = element_text(size=20,face="bold"), title = element_blank(), 
        axis.title.y.right = element_text(size =15),
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        axis.line = element_line(linewidth =1))+
  geom_density2d(col = 'black')+
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_hline(yintercept =2.4)+
  xlab("CD14")+ ylab("HLA.DR")+
  ggtitle("")
ggsave(file = "output/9_CD38_HLADR/CD14_vs_HLADR_Biaxial_Plot_y2.4_pbmc.png",
       width=8, height=6, units = "in", bg = 'white')

# Add metadata 
dat <- CD8T@assays$integratedADT$data %>% t() %>% as.data.frame()
dat$CD38HLADR <- ""
dat$CD38HLADR[dat$CD38 >1.3 & dat$HLA.DR >2.4] <- "CD38+HLADR+"
dat$CD38HLADR[dat$CD38 >1.3 & dat$HLA.DR <2.4] <- "CD38+HLADR-"
dat$CD38HLADR[dat$CD38 <1.3 & dat$HLA.DR >2.4] <- "CD38-HLADR+"
dat$CD38HLADR[dat$CD38 <1.3 & dat$HLA.DR <2.4] <- "CD38-HLADR-"
CD8T$CD38HLADR <- dat$CD38HLADR

CD8T$CD38HLADR <- factor(CD8T$CD38HLADR, levels = c("CD38+HLADR-","CD38-HLADR-","CD38-HLADR+","CD38+HLADR+"))

# DimPlot for CD38 & HLA.DR expression (*Figure 2G*) 
DimPlot(CD8T, group.by ="CD38HLADR", reduction = "wnn.umap", label = F, pt.size =1) + 
  scale_color_manual(values = c("#F8EEC1","#85827A","#E8A761","#1B6393"))+
  theme(legend.position = "top",
        legend.text = element_text(size = 10, face = "bold"),
        legend.key.spacing = unit(0, "in")) +
  NoAxes() +
  guides(color = guide_legend(override.aes = list(size=3), nrow = 2, title.hjust=0.5)) +
  ggtitle("")
ggsave(file = "output/9_CD38_HLADR/DimPlot_CD8T_CD38HLADR.png",
       width=4, height=4, units = "in", bg = 'white')


# save RDS file
saveRDS(CD8T, file = "output/9_CD38_HLADR/CD8T_diffusion_map_CD38HLADR.rds")


# *** end of CD38_HLADR analysis ***