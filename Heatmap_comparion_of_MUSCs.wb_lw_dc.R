### script from 'MuSCs_similarity.R'

# MuSCs homology
library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(dendextend)
library(ape)
library(ggdendro)
# ------------------------------------------------------
## all sub MuSCs hclust
## (1) laiwu MuSCs
ph.MuSCs.integrated <- readRDS('lw_MUSCs.rds')
## (2) duroc MuSCs
pt.MuSCs.integrated <- readRDS('dc_MUSCs.rds')
## (3) wild MuSCs
pw.MuSCs.integrated <- readRDS('wb_MUSCs.rds')


ph.MuSCs.integrated$cell_name_1128 <- paste0('lw ',ph.MuSCs.integrated$anno_subtype_1128)
pt.MuSCs.integrated$cell_name_1128 <- paste0('dc ',pt.MuSCs.integrated$anno_subtype_1128)
pw.MuSCs.integrated$cell_name_1128 <- paste0('wb ',pw.MuSCs.integrated$anno_sub_1027)


DefaultAssay(ph.MuSCs.integrated) <- 'integrated'
DefaultAssay(pt.MuSCs.integrated) <- 'integrated'
DefaultAssay(pw.MuSCs.integrated) <- 'integrated'

var_ph <- VariableFeatures(ph.MuSCs.integrated)
var_pt <- VariableFeatures(pt.MuSCs.integrated)
var_pw <- VariableFeatures(pw.MuSCs.integrated)
# (1) 
ph_MuSCs_markers <- read.table('~/test/0407/20210412/ph_MuSCs_cca/ph.MuSCs.markers-2.txt',header=T,row.names=1,sep='\t')
ph_MuSCs_markers <- ph_MuSCs_markers[ph_MuSCs_markers$p_val_adj < 0.05,]

pt_MuSCs_markers <- read.table('~/test/0407/20210412/pt_MuSCs_cca/pt.MuSCs.markers-2.txt',header=T,row.names=1,sep='\t')
pt_MuSCs_markers <- pt_MuSCs_markers[pt_MuSCs_markers$p_val_adj < 0.05,]

pw_MuSCs_markers <- read.table('~/test/0407/20210412/pw_MuSCs_cca/pw.MuSCs.markers-2.txt',header=T,row.names=1,sep='\t')
pw_MuSCs_markers <- pw_MuSCs_markers[pw_MuSCs_markers$p_val_adj<0.05,]

ph_gene <- intersect(var_ph,ph_MuSCs_markers$gene)
pt_gene <- intersect(var_pt,pt_MuSCs_markers$gene)
pw_gene <- intersect(var_pw,pw_MuSCs_markers$gene)

hvg <- intersect(intersect(ph_gene,pt_gene),pw_gene)
# (2)
#hvg <- intersect(intersect(var_ph,var_pt),var_pw)
out.dir <- '~/test/0407/20221027_fig_plot/'
mylist <- list()
for (c in unique(ph.MuSCs.integrated$cell_name_1128)){
  temp <- names(ph.MuSCs.integrated$cell_name_1128[which(ph.MuSCs.integrated$cell_name_1128==c)])
  temp <- ph.MuSCs.integrated@assays$integrated@scale.data[hvg,temp]
  temp <- apply(temp,1,mean)
  mylist[[c]] <- temp
}
for (c in unique(pt.MuSCs.integrated$cell_name_1128)){
  temp <- names(pt.MuSCs.integrated$cell_name_1128[which(pt.MuSCs.integrated$cell_name_1128==c)])
  temp <- pt.MuSCs.integrated@assays$integrated@scale.data[hvg,temp]
  temp <- apply(temp,1,mean)
  mylist[[c]] <- temp
}
for (c in unique(pw.MuSCs.integrated$cell_name_1128)){
  temp <- names(pw.MuSCs.integrated$cell_name_1128[which(pw.MuSCs.integrated$cell_name_1128==c)])
  temp <- pw.MuSCs.integrated@assays$integrated@scale.data[hvg,temp]
  temp <- apply(temp,1,mean)
  mylist[[c]] <- temp
}
data_f_h <- as.data.frame(mylist)
data_f_h[data_f_h > 2] <- 2
data_f_h[data_f_h < -2] <- -2

# pheatmap(data_f_h)
p_clustered <- 
  pheatmap(t(data_f_h),clustering_distance_cols='correlation',
           color = colorRampPalette(colors = c("blue","white","red"))(100),
           fontsize = 7,
  )

data_f_h.trans <- t(data_f_h)
data_f_h.trans <- data_f_h.trans[p_clustered$tree_row$order,p_clustered$tree_col$order]

rNames <- sub('\\.\\.fast\\.','\\(fast\\)',rownames(data_f_h.trans))
rNames <- sub('\\.\\.slow\\.','\\(slow\\)',rNames)
rNames <- sub('\\.\\.','\\+ ',rNames)
rNames <- gsub('\\.',' ',rNames)
rownames(data_f_h.trans) <- rNames
pheatmap(data_f_h.trans,cluster_rows = F,cluster_cols = F,fontsize_row = 10,show_colnames = F,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         fontsize = 7,filename = paste0(out.dir,'Fig5B.heatmap.pdf'),width = 8.3,height = 7.2)


d <- dist(data_f_h.trans/2)
hc <- hclust(d)
dend <- as.dendrogram(hc)
dend <- dend %>% set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 6) %>%  # node point size
  # set("leaves_col",colorRampPalette(brewer.pal(n=9,name='YlOrRd'))(36)) %>%
  set("leaves_col",rainbow(36)) %>%
  set("labels_cex",0.6)
ggd1 <- as.ggdend(dend)
# ggplot(ggd1) 
ggplot(ggd1,horiz=T,offset_labels = -15) + scale_x_reverse()

