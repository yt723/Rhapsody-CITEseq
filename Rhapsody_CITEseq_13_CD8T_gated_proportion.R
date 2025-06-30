#### CD8 T cell proportion comparison determined by gates 

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list=ls())

#load the data
CD8T <- readRDS(file = "output/4_CD8T_clusters/CD8T_clusters.rds")
#CD161: 2.5
#CCR7: 1.1
#CD27: 1.9
#CD38: 1.3
#HLA.DR: 2.4
#CD56: 2.7 

df <- CD8T@assays$integratedADT@data %>% t() %>% as.data.frame()
df <- df[,c("CD161","CCR7","CD27","CD38","HLA.DR","CD56")]
df$cluster <- "others"
df$cluster[df$CD161 >2.5] <- "C6"
df$cluster[df$CD161 <=2.5 & df$CCR7 >1.1 & df$CD27 >1.9] <- "C3"
df$cluster[df$CD161 <=2.5 & df$CCR7 <=1.1 & df$CD27 >1.9 & df$CD38 <=1.3 & df$HLA.DR <=2.4] <- "C7"
df$cluster[df$CD161 <=2.5 & df$CCR7 <=1.1 & df$CD27 >1.9 & df$CD38 <=1.3 & df$HLA.DR >2.4] <- "C2"
df$cluster[df$CD161 <=2.5 & df$CCR7 <=1.1 & df$CD27 >1.9 & df$CD38 >1.3 & df$HLA.DR >2.4] <- "C4"
df$cluster[df$CD161 <=2.5 & df$CCR7 <=1.1 & df$CD27 <=1.9 & df$CD56 <=2.7] <- "C1"
df$cluster[df$CD161 <=2.5 & df$CCR7 <=1.1 & df$CD27 <=1.9 & df$CD56 >2.7] <- "C5"

table(df$cluster)
all(rownames(CD8T@meta.data)==rownames(df)) 
CD8T$gate_cluster <- df$cluster

# gated proportion comparison
df_prop <- CD8T@meta.data %>% group_by(gate_cluster, SampleName, SampleType) %>% summarise(num=n())
df_prop <- df_prop %>% pivot_wider(names_from = 'gate_cluster', values_from = 'num', values_fill = 0)

mat_prop <- df_prop[, 3:ncol(df_prop)] %>% as.matrix()
mat_prop <- round(prop.table(mat_prop, margin = 1) *100, 2)
df_prop[,3:ncol(df_prop)] <- mat_prop 

write.csv(df_prop, file = "output/13_gated_proportion/CD8T_gate_cluster_proportion_patient.csv", row.names = F) 

## similarity between the unsupervised and gated clusters 
df <- CD8T@meta.data %>% select(cluster_name0.3, gate_cluster)
df_sum <- df %>% group_by(cluster_name0.3, gate_cluster) %>% summarise(num=n())
df_sum <- df_sum %>% group_by(gate_cluster) %>% mutate(total = sum(num), prop = round(num/total,2))
df_mat <- df_sum[,c("cluster_name0.3","gate_cluster","prop")]
df_mat <- df_mat %>% pivot_wider(names_from = cluster_name0.3, values_from = prop, values_fill = 0)
df_mat <- df_mat %>% pivot_longer(cols = 2:8, names_to = "cluster_name0.3", values_to = "prop")
df_mat$gate_cluster <- paste0("gate ", df_mat$gate_cluster)

ggplot(df_mat, aes(y=cluster_name0.3, x=gate_cluster, fill = prop)) +
  geom_tile(color='black', linewidth = 0.5) +
  geom_text(aes(label=prop)) +
  scale_fill_gradient(high = 'red', low = 'white',
                      name = 'Cluster\noverlap') +
  theme_minimal() +
  theme(axis.text = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 0.75),
        legend.title = element_text(size=16))+
  xlab("") + ylab("")+
  coord_fixed()  
ggsave(file= "output/13_gated_proportion/heatmap_CD8T_gate_cluster_overlap.png",
       height = 6, width = 7, units = 'in', bg = 'white')


# *** end of CD8T gated clusters ***