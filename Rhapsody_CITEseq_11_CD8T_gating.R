#### CD8 T cell gating simulation 

library(Seurat)
library(ggplot2)

rm(list=ls())


#load the data
CD8T <- readRDS(file = "output/4_CD8T_clusters/CD8T_clusters.rds")

## CD161 vs CD56 
color_d3 <- c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#8C564BFF",
              "#D62728FF","#E377C2FF","#9467BDFF")

color1 <- c("gray","gray","gray","gray",
            "gray","#E377C2FF","gray")
names(color1)<- sort(unique(CD8T$cluster_name0.3))
DefaultAssay(CD8T) <- "integratedADT"
# *Figure 6A*
FeatureScatter(CD8T, feature1 ="CD56",feature2 ="CD161", group.by = "cluster_name0.3",pt.size = 1.7)+   
  theme(axis.text.x=element_text(size=15), 
        axis.title = element_text(size=20,face="bold"), title = element_blank(), 
        axis.title.y.right = element_text(size =20),
        legend.position = 'right', legend.text = element_text(size = 16, face = 'bold'),
        axis.line = element_line(size=1))+
  scale_color_manual(values = color1)+
  geom_hline(yintercept =2.5)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  ggtitle("")
ggsave(file = "output/11_CD8T_gating/FeatureScatter_CD8T_CD161_vs_CD56.png",
       height = 5, width = 7, units = 'in', bg = 'white')


## CCR7 vs CD27 
# remove C6
CD161n_CD8T <- subset(CD8T, cluster_0.3 !="6")
CD161n_CD8T$CCR7CD27 <- ifelse(CD161n_CD8T$cluster_0.3=="3","3:naive",
                               ifelse(CD161n_CD8T$cluster_0.3 %in% c("2","4","7"), "2+4+7:EM","1+5:EMRA"))

color2 <- c("#BCBD22FF","#17BECFFF","#2CA02CFF")
names(color2) <- sort(unique(CD161n_CD8T$CCR7CD27))

# *Figure 6A*
FeatureScatter(CD161n_CD8T, feature1 ="CCR7",feature2 ="CD27", group.by = "CCR7CD27", pt.size = 1.7)+   
  theme(axis.text.x=element_text(size=15), 
        axis.title = element_text(size=20,face="bold"), title = element_blank(), 
        axis.title.y.right = element_text(size =20),
        legend.position = 'right', legend.text = element_text(size = 16, face = 'bold'),
        axis.line = element_line(size=1))+
  scale_color_manual(values = color2)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_vline(xintercept =1.1)+
  geom_hline(yintercept =1.9)+
  ggtitle("")
ggsave(file = "output/11_CD8T_gating/FeatureScatter_CD8T_CD27_vs_CCR7.png",
       height = 5, width = 7, units = 'in', bg = 'white')

## CD38 vs HLA.DR
# TEM
Tem_CD8T <- subset(CD161n_CD8T, cluster_0.3 %in% c("2","4","7"))

color3 <- c("#FF7F0EFF","#8C564BFF","#9467BDFF")
names(color3) <- sort(unique(Tem_CD8T$cluster_name0.3))

# *Figure 6A*
FeatureScatter(Tem_CD8T, feature1 ="CD38",feature2 ="HLA.DR", group.by = "cluster_name0.3", pt.size = 1.7)+   
  theme(axis.text.x=element_text(size=15), 
       axis.title = element_text(size=20,face="bold"), title = element_blank(), 
       axis.title.y.right = element_text(size =20),
       legend.position = 'right', legend.text = element_text(size = 16, face = 'bold'),
       axis.line = element_line(size=1))+
  scale_color_manual(values = color3)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_vline(xintercept =1.3)+
  geom_hline(yintercept =2.4)+
  ggtitle("")
ggsave(file = "output/11_CD8T_gating/FeatureScatter_CD8T_CD38_vs_HLADR.png",
       height = 5, width = 7, units = 'in', bg = 'white')

## CD56 vs CD38 
# TEMRA
Temra_CD8T <- subset(CD161n_CD8T, cluster_0.3 %in% c("1","5"))

color4 <- c("#1F77B4FF","#D62728FF")
names(color4) <- sort(unique(Temra_CD8T$cluster_name0.3))

# *Figure 6A* 
FeatureScatter(Temra_CD8T, feature1 ="CD38",feature2 ="CD56", group.by = "cluster_name0.3", pt.size = 1.7)+   
  theme(axis.text.x=element_text(size=15), 
        axis.title = element_text(size=20,face="bold"), title = element_blank(), 
        axis.title.y.right = element_text(size =20),
        legend.position = 'right', legend.text = element_text(size = 16, face = 'bold'),
        axis.line = element_line(size=1))+
  scale_color_manual(values = color4)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_hline(yintercept =2.7)+
  ggtitle("") 
ggsave(file = "output/11_CD8T_gating/FeatureScatter_CD8T_CD38_vs_CD56.png",
       height = 5, width = 7, units = 'in', bg = 'white')

# batch effect (*Figure S9B)
DefaultAssay(CD161n_CD8T) <- "integratedADT"
FeatureScatter(CD161n_CD8T, feature1 ="CCR7",feature2 ="CD27", group.by = "Library")+   
  theme(axis.text.x=element_text(size=15), 
        axis.title = element_text(size=20,face="bold"), title = element_blank(), 
        axis.title.y.right = element_text(size =20),
        legend.position = 'right', legend.text = element_text(size = 16, face = 'bold'),
        axis.line = element_line(size=1))+
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_vline(xintercept =1.1)+
  geom_hline(yintercept =1.9)+
  ggtitle("")
ggsave(file = "output/11_CD8T_gating/FeatureScatter_CD8T_CD27_vs_CCR7_Library_integrated.png",
       height = 5, width = 7, units = 'in', bg = 'white')

DefaultAssay(CD161n_CD8T) <- "ADT"
FeatureScatter(CD161n_CD8T, feature1 ="CCR7",feature2 ="CD27", group.by = "Library")+   
  theme(axis.text.x=element_text(size=15), 
        axis.title = element_text(size=20,face="bold"), title = element_blank(), 
        axis.title.y.right = element_text(size =20),
        legend.position = 'right', legend.text = element_text(size = 16, face = 'bold'),
        axis.line = element_line(size=1))+
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_vline(xintercept =1.1)+
  geom_hline(yintercept =1.9)+
  ggtitle("")
ggsave(file = "output/11_CD8T_gating/FeatureScatter_CD8T_CD27_vs_CCR7_Library_raw.png",
       height = 5, width = 7, units = 'in', bg = 'white')


# *** end of CD8T gating ***