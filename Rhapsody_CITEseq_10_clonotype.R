#### clonotype analysis  

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(tidyr)

rm(list=ls())

#load the data
CD8T <- readRDS(file =  "output/9_CD38_HLADR/CD8T_diffusion_map_CD38HLADR.rds")
vdj_df <- read.csv(file = "raw_data/vdj_combined.csv", row.names = 1)
vdj_df <- vdj_df %>% select(-contains("BCR"))

meta <- CD8T@meta.data
vdj_df2 <- meta %>% merge(vdj_df, by = 0)

# use only paired TCR chain data 
vdj_df2 <- vdj_df2[vdj_df2$TCR_Paired_Chains == "TRUE",]
vdj_df2$CDR3AB <- paste0(vdj_df2$TCR_Alpha_Gamma_CDR3_Translation_Dominant,"_",vdj_df2$TCR_Beta_Delta_CDR3_Translation_Dominant)

CDR3_clonotype <- vdj_df2 %>% 
  group_by(CDR3AB, SampleName)%>%
  summarise(cell_count=n())

CDR3_clonotype <- CDR3_clonotype[order(CDR3_clonotype$cell_count, decreasing = T),]
CDR3_clonotype$clonality <- ifelse(CDR3_clonotype$cell_count >=2, "Expanded", "Non_expanded")
CDR3_clonotype$clone <- factor(1:length(rownames(CDR3_clonotype)))
CDR3_clonotype$CDR3AB_SampleName <- paste0(CDR3_clonotype$CDR3AB,"_",CDR3_clonotype$SampleName)

vdj_df2$CDR3AB_SampleName <- paste0(vdj_df2$CDR3AB,"_",vdj_df2$SampleName)
vdj_df2 <- vdj_df2 %>% merge(CDR3_clonotype, by.x = c("CDR3AB_SampleName","CDR3AB","SampleName"), 
                             by.y = c("CDR3AB_SampleName","CDR3AB","SampleName")) #, all = TRUE

df_add <- vdj_df2[,c("CDR3AB","cell_count","clonality","clone","cell_id","CDR3AB_SampleName","TCR_Alpha_Gamma_V_gene_Dominant")]
rownames(df_add) <- df_add$cell_id

CD8T <- AddMetaData(CD8T, metadata = df_add)

## TRAV
prop.table(table(CD8T$cluster_name0.3, CD8T$TCR_Alpha_Gamma_V_gene_Dominant), margin = 1)
#TRAV1-2*01: C1-0.09, C2-0.05, C3-0.03, C4-0.01, C5-0.03, C6-0.44, C7-0.07 

# save RDS file 
saveRDS(CD8T, file = "output/10_clonotype/CD8T_clonotype.rds")

## Subset data with paired TCR data
CD8T <- subset(CD8T, clonality != "NA") 
vdj_df <- CD8T@meta.data

## TCR diversity comparison between the disease types 
library(iNEXT)

clone_df  <- vdj_df %>%
  group_by(SampleName, CDR3AB) %>%
  summarise(cell_count = n())

ls_clone <- split(clone_df$cell_count, clone_df$SampleName)
out <- iNEXT(ls_clone, q=1, datatype = "abundance") #q = Hill number 

# Rarefaction/extrapolation curve (*Figure S8A*) #type = 3: coverage based R/E curve 
ggiNEXT(out, type=3, color.var = "Assemblage", se=FALSE) +
  scale_color_manual(values = c(rep("maroon",6),rep("darkgray",12))) +
  theme_classic()+  
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(face = "bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        legend.position = 'right'
  ) +  
  scale_shape_manual(values = rep(20,18)); 
ggsave(file = "output/10_clonotype/R_E_curve_q1_CD8T_type3.png",
       width = 7, height=5, units = "in", bg = 'white')

# Richness, Shannon index, and Simpson index 
out_Richness <- ChaoRichness(ls_clone, datatype = "abundance")
out_Shannon <- ChaoShannon(ls_clone, datatype = "abundance", B=3000)
out_Simpson <- ChaoSimpson(ls_clone, datatype = "abundance", B=3000)
out_Shannon$SampleType <- c(rep("CAV",6),rep("Normal_HTx",12))
out_Richness$SampleType <- c(rep("CAV",6),rep("Normal_HTx",12))
out_Simpson$SampleType <- c(rep("CAV",6),rep("Normal_HTx",12))

library(ggpubr)
library(rstatix)

shannon_stat <- out_Shannon %>% group_by(SampleType) %>%
  summarise(med=median(Estimator),
            iqr0.25 = quantile(Estimator,0.25),
            iqr0.75 = quantile(Estimator, 0.75))
# CAV:          4.37    [3.92,    4.87]
# Normal_HTx:   6.08    [5.67,    6.61]

# Shannon
stat_test <- out_Shannon %>%
  wilcox_test(Estimator ~SampleType) %>%
  add_significance() %>% add_xy_position()

gg1<- ggplot(out_Shannon, aes(x=SampleType, y=Estimator)) +
  geom_boxplot(aes(fill=SampleType), width = 0.3) +
  theme_classic()+
  geom_point(aes(group = SampleType), position = position_dodge(0.75), size=0.5)+
  theme(axis.text.x = element_text(face = "bold", size=7, angle = 20, vjust = 0.5), 
        axis.text.y = element_text(face = "bold", size=7),
        axis.title.y = element_text(face = "bold", size=10),
        legend.position = 'none'
        )+
  scale_fill_manual(values = c("maroon","gray"))+
  stat_pvalue_manual(stat_test, label = 'p={round(p,3)}', hide.ns = FALSE, label.size = 3)+
  ylim(c(0,8))+
  xlab("") + ylab("Shannon Index");gg1 

# Richness
stat_test <- out_Richness %>%
  wilcox_test(Estimator ~SampleType) %>%
  add_significance() %>% add_xy_position()

gg2 <- ggplot(out_Richness, aes(x=SampleType, y=Estimator)) +
  geom_boxplot(aes(fill=SampleType), width = 0.3) +
  theme_classic()+
  geom_point(aes(group = SampleType), position = position_dodge(0.75), size=0.5)+
  theme(axis.text.x = element_text(face = "bold", size=7, angle = 20, vjust = 0.5), 
        axis.text.y = element_text(face = "bold", size=7),
        axis.title.y = element_text(face = "bold", size=10),
        legend.position = 'none'
  )+
  scale_fill_manual(values = c("maroon","gray"))+
  stat_pvalue_manual(stat_test, label = 'p={round(p,3)}', hide.ns = FALSE, label.size = 3)+
  ylim(c(0,5500))+
  xlab("") + ylab("Species Richness") ;gg2

# Simpson
stat_test <- out_Simpson %>%
  wilcox_test(Estimator ~SampleType) %>%
  add_significance() %>% add_xy_position()
stat_test$y.position <- 1.02

gg3 <- ggplot(out_Simpson, aes(x=SampleType, y=Estimator)) +
  geom_boxplot(aes(fill=SampleType), width = 0.3) +
  theme_classic()+
  geom_point(aes(group = SampleType), position = position_dodge(0.75), size=0.5)+
  theme(axis.text.x = element_text(face = "bold", size=7, angle = 20, vjust = 0.5), 
        axis.text.y = element_text(face = "bold", size=7),
        axis.title.y = element_text(face = "bold", size=10),
        legend.position = 'none')+
  scale_fill_manual(values = c("maroon","gray"))+
  stat_pvalue_manual(stat_test, label = 'p={round(p,3)}', hide.ns = FALSE, label.size = 3)+
  ylim(c(0,1.1))+
  xlab("") + ylab("Simpson Index") ;gg3

# *Figure 4A* 
ggarrange(gg2,gg1,gg3, ncol = 3)
ggsave(file = "output/10_clonotype/CD8T_diversity_SampleType.png",
       width = 4, height=2.3, units = "in", bg = 'white') 

## TCR diversity comparison between the CD38HLADR groups
clone_df <- vdj_df %>%
  group_by(CD38HLADR, CDR3AB) %>%
  summarise(cell_count = n())

ls_clone <- split(clone_df$cell_count, clone_df$CD38HLADR)
out_Shannon <- ChaoShannon(ls_clone, datatype = "abundance", B=3000)
out_Shannon$cluster <- rownames(out_Shannon)

out_Shannon$cluster <- factor(out_Shannon$cluster, levels = c("CD38+HLADR-","CD38-HLADR-","CD38-HLADR+","CD38+HLADR+"))
# *Figure 5F* 
ggplot(out_Shannon,aes(x = cluster, y = Estimator, fill = cluster))+
  geom_col(width = 0.5) +
  scale_fill_manual(values = c("#F8EEC1","#85827A","#E8A761","#1B6393")) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = 'none',
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(angle=-90, hjust = 0, vjust = 0.5, size = 16, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        plot.title = element_text(hjust = 0.5))+
  ylab("Shannon index") + xlab("") + ggtitle("")
ggsave(file = "output/10_clonotype/DPT_CD8T_Shannon_CD38HLADR.png",
       width=4, height=4.5, units = "in", bg = 'white')


## log10(clonotype size) 
CD8T$log_cell_count <- log10(CD8T$cell_count)

# *Figure 4B*
FeaturePlot(CD8T, features = 'log_cell_count', reduction = 'wnn.umap', order = T,
                  pt.size = 1.5,
                  cols = c("blue","yellow")) +
  theme(text = element_text(face = "bold", size = 20))+
  NoAxes()+
  labs(color = "Size of\nclonotypes\nlog10 value") +
  ggtitle("") 
ggsave(file = "output/10_clonotype/FeaturePlot_CD8T_clone_log10.png",
       width = 7, height=5, units = "in", bg = 'white')

## Comparison of size of expanded clonotype between clusters 
vdj_df <- CD8T@meta.data 
vdj_df <- vdj_df %>%
  filter(clonality=="Expanded") 

stat <- vdj_df %>%
  group_by(cluster_name0.3) %>%
  summarise(avg_log=mean(log_cell_count))

col <- ggsci::pal_d3("category20")(20)
fillset <- c(col[11:13], col[16],col[14],col[17],col[15])
names(fillset) <- sort(unique(vdj_df$cluster_name0.3))

# *Figure 4C* 
ggplot() +
  geom_quasirandom(data=vdj_df, aes(x=cluster_name0.3, y=log_cell_count, color = cluster_name0.3)) +
  scale_color_manual(values = fillset) +
  geom_point(data= stat, aes(x=cluster_name0.3, y=avg_log), size=2)+
  theme_classic() +
  theme(legend.position = 'right',
        legend.text = element_text(size = 11, face = 'bold'),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title.y = element_text(size = 11, face = 'bold'),
        legend.key.height = unit(0.2, "in"),
        legend.key.width = unit(0.2, "in"),
        legend.key.spacing = unit(0, "in")
        ) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  xlab("") +ylab("Size of clonotypes (log10)")
ggsave(file = "output/10_clonotype/CD8T_clone_log10count_cluster.png",
       width = 4.5, height=3, units = "in", bg = 'white')

## Clonotype overlaps: C4 vs others  
vdj_df <- CD8T@meta.data
C4_clone <- vdj_df[vdj_df$cluster_0.3 == "4", ]$clone %>% unique()
vdj_df$C4_clone <- ifelse(vdj_df$clone %in% C4_clone, "True", "False")

df_summary <- vdj_df %>% 
  filter(C4_clone =="True") %>%
  group_by(clone, cluster_0.3) %>%
  summarise(clone_num = n())

df_summary$cluster_0.3 <- paste0("C",df_summary$cluster_0.3)
df3 <- vdj_df %>%
  group_by(cluster_0.3) %>%
  summarise(cluster_cell = n())
df3$cluster_0.3 <- paste0("C", df3$cluster_0.3)

df_summary <- df_summary %>%
  merge(df3, by.x=c("cluster_0.3", by.y=c("cluster_0.3")))

df_summary <- df_summary %>%  mutate(prop=clone_num*100/cluster_cell)
df_summary <- df_summary %>% group_by(clone) %>% mutate(summed = sum(clone_num))

df_summary <- df_summary %>% group_by(clone) %>%
  mutate(C4_FC = prop/.data$prop[grepl("C4",.data$cluster_0.3)],
         C4_FC_inv = 1/C4_FC, 
         C4_log2FC = log2(C4_FC_inv))
df_summary$comparison <- ifelse(df_summary$C4_log2FC == 0, "Equally expanded", 
                                 ifelse(df_summary$C4_log2FC < 0, "Less expanded in C4", 
                                        "More expanded in C4" ))
df_summary2 <- df_summary %>%
  filter(cluster_0.3 != "C4")
df_summary2$cluster <- paste0(df_summary2$cluster_0.3," vs C4")

color_pal <- c("black","red","blue")
names(color_pal) <- c("Equally distributed", "More expanded in C4","Less expanded in C4")

# *Figure 4D*
ggplot(df_summary2, aes(x=cluster, y=C4_log2FC, colour = comparison)) +
  geom_quasirandom(aes(size=summed)) +
  scale_color_manual(values = color_pal, name = "") +
  scale_radius(range = c(1,5), name = "clonotype size") +
  theme_minimal() + 
  theme(axis.text.y = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.text.x = element_text(size = 14, face = 'bold',angle = 270, vjust = 0.5, hjust = 0),
        legend.text = element_text(size = 14, face = 'bold'),
        legend.title=element_text(size = 14, face = 'bold'), 
        axis.line = element_line(color = 'black'),
        #axis.ticks = element_line(color = 'black', linewidth = 0.2),
        legend.spacing = unit(0, "in"),
        legend.key.width = unit(0, "in"),
        legend.key.height = unit(0.1, "in"),
        legend.position = "right")+ 
  ylab("Fold Change (log2)") +xlab("") +ggtitle("")
ggsave(file = "output/10_clonotype/CD8T_clone_overlap_log2FC_cluster.png",
       width = 5.5, height=4.5, units = "in", bg = 'white')

## Repertoire overlap among CD8T clusters (Jaccard index) 
vdj_df <- CD8T@meta.data
vdj_df$cluster_0.3 <- paste0("C", vdj_df$cluster_0.3)

df_clone <- vdj_df %>%
  filter(clonality=="Expanded") %>%
  group_by(cluster_0.3, CDR3AB) %>%
  summarise(clone_num = n())

df_clone <- df_clone %>%
  pivot_wider(names_from = "cluster_0.3", values_from = clone_num,
              values_fill = 0)

library(vegan)
library(reshape2)
df_mat <- as.matrix(df_clone[,2:8]) %>% t()
dissimilarity <- vegdist(df_mat, method = 'jaccard', na.rm = T, upper = T, diag = T)
similarity <- rep(1, length(dissimilarity)) - dissimilarity
df_out <- melt(as.matrix(similarity), na.rm = T) 
mat_result <- as.matrix(df_out[,c(-1,-2)]) 
hm_mat <- matrix(mat_result, ncol = 7, nrow = 7)
colnames(hm_mat) <- sort(unique(vdj_df$cluster_0.3))
rownames(hm_mat) <- sort(unique(vdj_df$cluster_0.3))
hm_mat[upper.tri(hm_mat, diag = T)] <- NA 
hm_mat <- t(hm_mat)
df_out <- melt(as.matrix(hm_mat), na.rm = T)
df_out$value <- round(df_out$value, 3)

# Jacaard index heatmap (*Figure 4E*)
ggplot(df_out, aes(Var2, Var1, , fill=value))+
  geom_tile(color='black', linewidth = 0.5)+
  geom_text(aes(label = value))+
  scale_fill_gradient(high = 'red', low = 'white',
                      name = 'Jaccard\nindex') +
  theme_minimal() +
  theme(axis.text = element_text(size = 16, face = 'bold'),
        legend.title = element_text(size=16))+
  xlab("") + ylab("")+
  coord_fixed()
ggsave(file = "output/10_clonotype/CD8T_heatmap_jaccard_clusters.png",
       width = 5, height=5, units = "in", bg = 'white')


## Jaccard index calculation per patient 
df_clone <- vdj_df %>%
  filter(clonality=="Expanded") %>%
  group_by(cluster_0.3, CDR3AB, SampleName) %>%
  dplyr::summarise(clone_num = n())

df_clone <- df_clone %>%
  pivot_wider(names_from = "cluster_0.3", values_from = clone_num,
              values_fill = 0)

Var1 <- c("C1","C2","C3","C4","C5","C6","C7")
Var2 <- c("C1","C2","C3","C4","C5","C6","C7")
output <- expand.grid("Var1" = Var1, "Var2" = Var2)

for (name in sort(unique(vdj_df$SampleName))){
  df <- df_clone %>% filter(SampleName == name)
  df <- as.matrix(df[,3:9]) %>% t()
  dissimilarity <- vegdist(df, method = 'jaccard', na.rm = F, upper = T, diag = T)
  similarity <- rep(1, length(dissimilarity)) - dissimilarity
  df_out <- melt(as.matrix(similarity), na.rm = F) 
  colnames(df_out)[3] <- name
  output <- output %>% merge(df_out, by.x = c("Var1","Var2"), by.y = c("Var1","Var2"))
}
output[is.na(output)] <- 0
output2 <- output[output$Var1 != "C4" & output$Var2 == "C4" ,]

rownames(output2)<- paste0(output2$Var1, " vs ", output2$Var2)
output2$Var1 <- NULL
output2$Var2 <- NULL
output2 <- output2 %>% t() %>% as.data.frame()
output2$SampleType <- c(rep("CAV",6), rep("Normal_HTx",12))
output2 <- output2 %>% pivot_longer(cols = 1:6, names_to = "cluster", values_to = "jaccard")

# Jaccard index: C4 vs the other (per patient) (*Figure 4F*)
ggplot(output2, aes(x=cluster, y=jaccard, color=SampleType)) +
  geom_quasirandom(aes(shape=SampleType),size=3) +
  scale_color_manual(values = c("maroon","gray50"))+
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        axis.text.x = element_text(size = 15, face = 'bold',angle = 270, vjust = 0.5, hjust = 0),
        legend.text = element_text(size = 15, face = 'bold'),
        legend.title=element_blank(), 
        axis.line = element_line(color = 'black'),
        legend.spacing = unit(0, "in"),
        legend.key.width = unit(0, "in"),
        legend.key.height = unit(0.1, "in"),
        legend.position = "right")+ 
  guides(colour = guide_legend(override.aes = list(size=4)))+
  xlab("") +ylab("Jaccard index")
ggsave(file = "output/10_clonotype/CD8T_jaccard_C4_vs_other.png",
       width = 5.5, height=4.5, units = "in", bg = 'white')

## Overlapping clonotypes 
df_vdj <- CD8T@meta.data
df_vdj <- df_vdj %>% filter(!is.na(CDR3AB)) #5777 cells
unique_clone <- df_vdj %>% group_by(CDR3AB, SampleName) %>%
  summarise(clone_num=n()) #2926 clonotypes 
unique_clone2 <- unique_clone %>% group_by(CDR3AB) %>%
  summarise(clone_num=n())
shared_clone <- unique_clone2 %>% filter(clone_num==2) #11 clonotypes 
df_vdj_shared <- df_vdj %>% filter(CDR3AB %in% shared_clone$CDR3AB) #262 cells
df_vdj_shared <- df_vdj_shared %>%
  group_by(CDR3AB, SampleName) %>%
  summarise(clone_num = n())

## TCR diversity comparison: CAV vs Normal HTx (CD4T) 
rm(list = ls())

#load the data
CD4T <- readRDS(file = "output/5_CD4T_clusters/CD4T_clusters.rds")
vdj_df <- read.csv(file = "raw_data/vdj_combined.csv", row.names = 1)
vdj_df <- vdj_df %>% select(-contains("BCR"))

meta <- CD4T@meta.data
vdj_df2 <- meta %>% merge(vdj_df, by = 0)

# use only paired TCR chain data 
vdj_df2 <- vdj_df2[vdj_df2$TCR_Paired_Chains == "TRUE",]
vdj_df2$CDR3AB <- paste0(vdj_df2$TCR_Alpha_Gamma_CDR3_Translation_Dominant,"_",vdj_df2$TCR_Beta_Delta_CDR3_Translation_Dominant)

CDR3_clonotype <- vdj_df2 %>% 
  group_by(CDR3AB, SampleName)%>%
  summarise(cell_count=n())

CDR3_clonotype <- CDR3_clonotype[order(CDR3_clonotype$cell_count, decreasing = T),]
CDR3_clonotype$clonality <- ifelse(CDR3_clonotype$cell_count >=2, "Expanded", "Non_expanded")
CDR3_clonotype$clone <- factor(1:length(rownames(CDR3_clonotype)))
CDR3_clonotype$CDR3AB_SampleName <- paste0(CDR3_clonotype$CDR3AB,"_",CDR3_clonotype$SampleName)

vdj_df2$CDR3AB_SampleName <- paste0(vdj_df2$CDR3AB,"_",vdj_df2$SampleName)
vdj_df2 <- vdj_df2 %>% merge(CDR3_clonotype, by.x = c("CDR3AB_SampleName","CDR3AB","SampleName"), 
                             by.y = c("CDR3AB_SampleName","CDR3AB","SampleName")) 

df_add <- vdj_df2[,c("CDR3AB","cell_count","clonality","clone","cell_id","CDR3AB_SampleName")]
rownames(df_add) <- df_add$cell_id

CD4T <- AddMetaData(CD4T, metadata = df_add)

# save RDS file 
saveRDS(CD4T, file = "output/10_clonotype/CD4T_clonotype.rds")

CD4T <- subset(CD4T, clonality != "NA") 
vdj_df <- CD4T@meta.data

clone_df  <- vdj_df %>%
  group_by(SampleName, CDR3AB) %>%
  summarise(cell_count = n())

ls_clone <- split(clone_df$cell_count, clone_df$SampleName)

# Richness, Shannon index, and Simpson index 
out_Richness <- ChaoRichness(ls_clone, datatype = "abundance")
out_Shannon <- ChaoShannon(ls_clone, datatype = "abundance", B=3000)
out_Simpson <- ChaoSimpson(ls_clone, datatype = "abundance", B=3000)
out_Shannon$SampleType <- c(rep("CAV",6),rep("Normal_HTx",12))
out_Richness$SampleType <- c(rep("CAV",6),rep("Normal_HTx",12))
out_Simpson$SampleType <- c(rep("CAV",6),rep("Normal_HTx",12))

shannon_stat <- out_Shannon %>% group_by(SampleType) %>%
  summarise(med=median(Estimator),
            iqr0.25 = quantile(Estimator,0.25),
            iqr0.75 = quantile(Estimator, 0.75))
# CAV:          8.60    [7.99,    9.27]
# Normal_HTx:   9.81    [9.46,   10.4]

# Shannon
stat_test <- out_Shannon %>%
  wilcox_test(Estimator ~SampleType) %>%
  add_significance() %>% add_xy_position()

gg1<- ggplot(out_Shannon, aes(x=SampleType, y=Estimator)) +
  geom_boxplot(aes(fill=SampleType), width = 0.5) +
  theme_classic()+
  geom_point(aes(group = SampleType), position = position_dodge(0.75), size=0.5)+
  theme(axis.text.x = element_text(face = "bold", size=7, angle = 20, vjust = 0.5), 
        axis.text.y = element_text(face = "bold", size=7),
        axis.title.y = element_text(face = "bold", size=10),
        legend.position = 'none')+
  scale_fill_manual(values = c("maroon","gray"))+
  stat_pvalue_manual(stat_test, label = 'p={round(p,3)}', hide.ns = FALSE, label.size = 3)+
  ylim(c(0,11.5))+
  xlab("") + ylab("Shannon Index");gg1 

# Richness
stat_test <- out_Richness %>%
  wilcox_test(Estimator ~SampleType) %>%
  add_significance() %>% add_xy_position()

gg2 <- ggplot(out_Richness, aes(x=SampleType, y=Estimator)) +
  geom_boxplot(aes(fill=SampleType), width = 0.5) +
  theme_classic()+
  geom_point(aes(group = SampleType), position = position_dodge(0.75), size=0.5)+
  theme(axis.text.x = element_text(face = "bold", size=7, angle = 20, vjust = 0.5), 
        axis.text.y = element_text(face = "bold", size=7),
        axis.title.y = element_text(face = "bold", size=10),
        legend.position = 'none')+
  scale_fill_manual(values = c("maroon","gray"))+
  stat_pvalue_manual(stat_test, label = 'p={round(p,3)}', hide.ns = FALSE, label.size = 3)+
  ylim(c(0,65000))+
  xlab("") + ylab("Richness"); gg2 

# Simpson
stat_test <- out_Simpson %>%
  wilcox_test(Estimator ~SampleType) %>%
  add_significance() %>% add_xy_position()
stat_test$y.position <- 1.02

gg3 <- ggplot(out_Simpson, aes(x=SampleType, y=Estimator)) +
  geom_boxplot(aes(fill=SampleType), width = 0.5) +
  theme_classic()+
  geom_point(aes(group = SampleType), position = position_dodge(0.75), size=0.5)+
  theme(axis.text.x = element_text(face = "bold", size=7, angle = 20, vjust = 0.5), 
        axis.text.y = element_text(face = "bold", size=7),
        axis.title.y = element_text(face = "bold", size=10),
        legend.position = 'none')+
  scale_fill_manual(values = c("maroon","gray"))+
  stat_pvalue_manual(stat_test, label = 'p={round(p,3)}', hide.ns = FALSE, label.size = 3)+
  ylim(c(0,1.1))+
  xlab("") + ylab("Simpson Index") ;gg3

# *Figure S8B*
fig <- ggarrange(gg2,gg1,gg3, ncol = 3) 
annotate_figure(fig, top = text_grob("CD4T cell", size = 11, face = 'italic'))
ggsave(file = "output/10_clonotype/CD4T_diversity_SampleType.png",
       width = 4.2, height=2.5, units = "in", bg = 'white')


## BCR diversity comparison: CAV vs Normal HTx 
rm(list = ls())

# load data
Bcell <- readRDS(file = "output/6_Bcell_clusters/Bcell_clusters.rds")
vdj_df <- read.csv("raw_data/vdj_combined.csv", row.names = 1)

vdj_df <- vdj_df %>% select(-contains("TCR"))

meta <- Bcell@meta.data
vdj_df2 <- meta %>% merge(vdj_df, by = 0)

# use only BCR Heavy chain 
CDR3_clonotype <- vdj_df2 %>%
  group_by(BCR_Heavy_CDR3_Translation_Dominant, SampleName)%>%
  summarise(cell_count=n())

CDR3_clonotype <- CDR3_clonotype[order(CDR3_clonotype$cell_count, decreasing = T),]
CDR3_clonotype$clonality <- ifelse(CDR3_clonotype$cell_count >=2, "Expanded", "Non_expanded")
CDR3_clonotype$clone <- factor(1:length(rownames(CDR3_clonotype)))

CDR3_clonotype$BCR_SampleName <- paste0(CDR3_clonotype$BCR_Heavy_CDR3_Translation_Dominant,"_",CDR3_clonotype$SampleName)

vdj_df2$BCR_SampleName <- paste0(vdj_df2$BCR_Heavy_CDR3_Translation_Dominant,"_",vdj_df2$SampleName)
vdj_df2 <- vdj_df2 %>% merge(CDR3_clonotype, by.x = c("BCR_SampleName","BCR_Heavy_CDR3_Translation_Dominant","SampleName"), 
                             by.y = c("BCR_SampleName","BCR_Heavy_CDR3_Translation_Dominant","SampleName")) #, all = TRUE

df_add <- vdj_df2[,c("BCR_Heavy_CDR3_Translation_Dominant","cell_count","clonality","clone","cell_id","BCR_SampleName")]
rownames(df_add) <- df_add$cell_id

Bcell <- AddMetaData(Bcell, metadata = df_add)

# save RDS file 
saveRDS(Bcell, file = "output/10_clonotype/Bcell_clonotype.rds")

Bcell <- subset(Bcell, clonality != "NA") 
vdj_df <- Bcell@meta.data

clone_df  <- vdj_df %>%
  group_by(SampleName, BCR_Heavy_CDR3_Translation_Dominant) %>%
  summarise(cell_count = n())

ls_clone <- split(clone_df$cell_count, clone_df$SampleName)

# Richness, Shannon index, and Simpson index 
out_Richness <- ChaoRichness(ls_clone, datatype = "abundance")
out_Shannon <- ChaoShannon(ls_clone, datatype = "abundance", B=3000)
out_Simpson <- ChaoSimpson(ls_clone, datatype = "abundance", B=3000)
out_Shannon$SampleType <- c(rep("CAV",6),rep("Normal_HTx",12))
out_Richness$SampleType <- c(rep("CAV",6),rep("Normal_HTx",12))
out_Simpson$SampleType <- c(rep("CAV",6),rep("Normal_HTx",12))

shannon_stat <- out_Shannon %>% group_by(SampleType) %>%
  summarise(med=median(Estimator),
            iqr0.25 = quantile(Estimator,0.25),
            iqr0.75 = quantile(Estimator, 0.75))
# CAV:          5.96    [4.74,    6.53]
# Normal_HTx:   6.49    [5.92    8.68]

# Shannon
stat_test <- out_Shannon %>%
  wilcox_test(Estimator ~SampleType) %>%
  add_significance() %>% add_xy_position()

gg1 <- ggplot(out_Shannon, aes(x=SampleType, y=Estimator)) +
  geom_boxplot(aes(fill=SampleType), width = 0.5) +
  theme_classic()+
  geom_point(aes(group = SampleType), position = position_dodge(0.75), size=0.5)+
  theme(axis.text.x = element_text(face = "bold", size=7, angle = 20, vjust = 0.5), 
        axis.text.y = element_text(face = "bold", size=7),
        axis.title.y = element_text(face = "bold", size=10),
        legend.position = 'none')+
  scale_fill_manual(values = c("maroon","gray"))+
  stat_pvalue_manual(stat_test, label = 'p={round(p,3)}', hide.ns = FALSE, label.size = 3)+
  ylim(c(0,10))+
  xlab("") + ylab("Shannon Index");gg1 

# Richness
stat_test <- out_Richness %>%
  wilcox_test(Estimator ~SampleType) %>%
  add_significance() %>% add_xy_position()

gg2 <- ggplot(out_Richness, aes(x=SampleType, y=Estimator)) +
  geom_boxplot(aes(fill=SampleType), width = 0.5) +
  theme_classic()+
  geom_point(aes(group = SampleType), position = position_dodge(0.75), size=0.5)+
  theme(axis.text.x = element_text(face = "bold", size=7, angle = 20, vjust = 0.5), 
        axis.text.y = element_text(face = "bold", size=7),
        axis.title.y = element_text(face = "bold", size=10),
        legend.position = 'none')+
  scale_fill_manual(values = c("maroon","gray"))+
  stat_pvalue_manual(stat_test, label = 'p={round(p,3)}', hide.ns = FALSE, label.size = 3)+
  ylim(c(0,24000))+
  xlab("") + ylab("Richness"); gg2 

# Simpson
stat_test <- out_Simpson %>%
  wilcox_test(Estimator ~SampleType) %>%
  add_significance() %>% add_xy_position()
stat_test$y.position <- 1.02

gg3 <- ggplot(out_Simpson, aes(x=SampleType, y=Estimator)) +
  geom_boxplot(aes(fill=SampleType), width = 0.5) +
  theme_classic()+
  geom_point(aes(group = SampleType), position = position_dodge(0.75), size=0.5)+
  theme(axis.text.x = element_text(face = "bold", size=7, angle = 20, vjust = 0.5), 
        axis.text.y = element_text(face = "bold", size=7),
        axis.title.y = element_text(face = "bold", size=10),
        legend.position = 'none')+
  scale_fill_manual(values = c("maroon","gray"))+
  stat_pvalue_manual(stat_test, label = 'p={round(p,3)}', hide.ns = FALSE, label.size = 3)+
  ylim(c(0,1.1))+
  xlab("") + ylab("Simpson Index") ;gg3

# *Figure S8C*
fig <- ggarrange(gg2,gg1,gg3, ncol = 3)
annotate_figure(fig, top = text_grob("B cell", size = 11, face = "italic"))
ggsave(file = "output/10_clonotype/Bcell_diversity_SampleType.png",
       width = 4.2, height=2.5, units = "in", bg = 'white')


# ***end of clonotype analysis ***

