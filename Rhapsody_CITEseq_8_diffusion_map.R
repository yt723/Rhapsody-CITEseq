#### Diffusion map 

library(Seurat)
library(ggplot2)
library(destiny)
library(dplyr)
library(tidyr)

rm(list=ls())

path = "~/Documents/UCSD_CAV/CITEseq"

#load data
CD8T <- readRDS(file = paste0(path, "/output/4_CD8T_clusters/CD8T_clusters.rds"))

naive <- subset(CD8T, cluster_0.3 =="3")

# determine the start point (in the naive cluster)
p <- DimPlot(naive)
cell.located <- CellSelector(plot=p, ident = "centered_cell")
# "L6_4925
cell.located <- "L6_4925"
# run diffusion map
counts <- CD8T@assays$RNA$data %>% as.matrix()
colnames(counts) <- CD8T$cluster_0.3

set.seed(123)
dm <- DiffusionMap(t(counts))
# use the first diffusion component (DC1) as a measure of pseudotime 
DCrank <- data.frame(pseudotime_diffusionmap = rank(eigenvectors(dm)[,1])) #rank cells by their dpt 
rownames(DCrank) <- rownames(CD8T@meta.data)
CD8T <- AddMetaData(CD8T, metadata = DCrank)

#  Diffusion pseudotime (DPT)
position <- grep(cell.located, rownames(DCrank)) #10992L 
dpt <- DPT(dm, tips =position) 

DPTrank <- data.frame(pseudotime_dpt = rank(dpt$dpt)) #rank cells by their dpt 
rownames(DPTrank) <- rownames(CD8T@meta.data)
CD8T <- AddMetaData(CD8T, metadata = DPTrank)

# Display DPT on UMAP (*Fig 5A*)
FeaturePlot(CD8T, features = 'pseudotime_dpt', reduction = 'wnn.umap',
            pt.size = 1) +
  theme(legend.text = element_text(face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.32),
        legend.key.height = unit(0.15, "in")) +
  scale_color_gradient(low = "blue", high = "yellow", name = "pseudotime") +
  NoAxes() + ggtitle("")
ggsave(file = paste0(path, "/output/8_diffusion_map/FeaturePlot_CD8T_DPT.png"), 
       width=4, height=4, units = "in", bg = 'white')


# Violin plot to compare clusters (*Fig 3C*)
col <- ggsci::pal_d3("category20")(20)
fillset <- c(col[11:13], col[16],col[14],col[17],col[15])

# *Figure 5B* 
ggplot()+
  geom_violin (data=CD8T@meta.data,
               aes(x = pseudotime_dpt, y = cluster_name0.3, 
                   colour = cluster_name0.3, fill = cluster_name0.3), linewidth=1) +
  geom_point(data=CD8T@meta.data,
             aes(x = pseudotime_dpt, y = cluster_name0.3), 
             stat = "summary", size=3) +
  scale_color_manual(values = fillset) +
  scale_fill_manual(values = fillset) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11, face = 'bold'),
        axis.title.x = element_text(size = 12, face = 'bold')) +
  xlab("Diffusion pseudotime") + ylab("") + ggtitle("")
ggsave(file = paste0(path, "/output/8_diffusion_map/DPT_CD8T_diffusion_cluster_name0.3.png"), 
       width=4.5, height=3, units = "in", bg = 'white')



# save RDS
saveRDS(CD8T, file = paste0(path, "/output/8_diffusion_map/CD8T_diffusion_map.rds"))

## DPT vs expression of CD8T activation markers
ADT <- CD8T@assays$ADT$scale.data %>% as.matrix() %>% t()
add <- CD8T$pseudotime_dpt %>% as.matrix()
colnames(add) <- "pseudotime_dpt"
df_ADT <- cbind(ADT, add) %>% as.data.frame()

RNA <- CD8T@assays$RNA$scale.data %>% as.matrix() %>% t()
add <- CD8T$pseudotime_dpt %>% as.matrix()
colnames(add) <- "pseudotime_dpt"
df_RNA <- cbind(RNA, add) %>% as.data.frame()

df_ADT2 <- df_ADT[,c("pseudotime_dpt","CD38","HLA.DR")]
df_RNA2 <- df_RNA[,c("IFNG","PRF1","GZMB")]

df <- cbind(df_ADT2, df_RNA2)
df <- df %>% pivot_longer(cols = c("CD38","HLA.DR","IFNG","GZMB","PRF1"), values_to = "expression",names_to = "marker")
df$marker <- factor(df$marker, levels = c("CD38","HLA.DR","IFNG","GZMB","PRF1"))


dpt <- CD8T@meta.data %>% filter(cluster_0.3=="4")
q1 <- quantile(dpt$pseudotime_dpt, 0.25)
q3 <- quantile(dpt$pseudotime_dpt, 0.75)

ggplot(df, aes(x=pseudotime_dpt, y=expression, colour = marker)) +
  geom_smooth(aes(fill = marker)) + 
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  theme_classic() +
  geom_vline(xintercept = c(q1,q3), linetype = "dashed", color = "grey", linewidth = 1) +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.15,0.8),
        legend.key.height = unit(0.15, "in"),
        legend.text = element_text(size = 12, face = 'bold'),
        axis.title.x = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        axis.text  = element_text(size = 11)) +
  xlab("Diffusion pseudotime") + ylab("Expression")
ggsave(file = paste0(path, "/output/8_diffusion_map/DPT_CD8T_marker_vs_DPT.png"), 
       width=4.5, height=3, units = "in", bg = 'white')

# *** end of diffusion map analysis ***
