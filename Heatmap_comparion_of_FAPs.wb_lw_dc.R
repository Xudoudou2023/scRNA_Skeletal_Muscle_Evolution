############################################## Fig3 BC###################
library(pheatmap)
library(RColorBrewer)
out.dir <- "~/fig_plot/"
########### 
# wild FAPs
wb_FAPs <- readRDS('wb_FAPs.rds')
# laiwu FAPs
lw_FAPs <- readRDS('lw_FAPs.rds'))
# duroc FAPs
dc_FAPs <- readRDS('dc_FAPs.rds'))
######## hclust
markers_dc_FAPs <- read.table(paste0(out.dir,'pt.FAPs.markers-2.txt'),header = T,row.names = 1,sep = '\t')
markers_dc_FAPs <- markers_dc_FAPs[markers_dc_FAPs$p_val_adj < 0.05,]

markers_lw_FAPs <- read.table(paste0(out.dir,'ph.FAPs.markers-2.txt'),header = T,row.names = 1,sep = '\t')
markers_lw_FAPs <- markers_lw_FAPs[markers_lw_FAPs$p_val_adj < 0.05,]

markers_wb_FAPs <- read.table(paste0(out.dir,'pw.FAPs.markers-2.txt'),header = T,row.names = 1,sep = '\t')
markers_wb_FAPs <- markers_wb_FAPs[markers_wb_FAPs$p_val_adj < 0.05,]

DefaultAssay(lw_FAPs.filtered) <- 'integrated'
DefaultAssay(dc_FAPs.filtered) <- 'integrated'
DefaultAssay(wb_FAPs.filtered) <- 'integrated'
var_ph <- VariableFeatures(lw_FAPs.filtered)
var_pt <- VariableFeatures(dc_FAPs.filtered)
var_pw <- VariableFeatures(wb_FAPs.filtered)

ph_gene <- intersect(var_ph,markers_lw_FAPs$gene)
pt_gene <- intersect(var_pt,markers_dc_FAPs$gene)
pw_gene <- intersect(var_pw,markers_wb_FAPs$gene)

hvg <- intersect(intersect(ph_gene,pt_gene),pw_gene)

mylist <- list()
for (c in unique(lw_FAPs.filtered$anno_20221031)){
  temp <- names(lw_FAPs.filtered$anno_20221031[which(lw_FAPs.filtered$anno_20221031==c)])
  temp <- lw_FAPs.filtered@assays$integrated@scale.data[hvg,temp]
  temp <- apply(temp,1,mean)
  mylist[[c]] <- temp
}
for (c in unique(dc_FAPs.filtered$anno_20221031)){
  temp <- names(dc_FAPs.filtered$anno_20221031[which(dc_FAPs.filtered$anno_20221031==c)])
  temp <- dc_FAPs.filtered@assays$integrated@scale.data[hvg,temp]
  temp <- apply(temp,1,mean)
  mylist[[c]] <- temp
}
for (c in unique(wb_FAPs.filtered$anno_20221031)){
  temp <- names(wb_FAPs.filtered$anno_20221031[which(wb_FAPs.filtered$anno_20221031==c)])
  temp <- wb_FAPs.filtered@assays$integrated@scale.data[hvg,temp]
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
           fontsize = 7)

data_f_h.trans <- t(data_f_h)
data_f_h.trans <- data_f_h.trans[p_clustered$tree_row$order,p_clustered$tree_col$order]
rNames <- sub('\\.',' ',rownames(data_f_h.trans))
rNames <- sub('\\.','-',rNames)
rownames(data_f_h.trans) <- rNames
pheatmap(data_f_h.trans,cluster_rows = F,cluster_cols = F,fontsize_row = 10,show_colnames = F,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         fontsize = 7,filename = paste0(out.dir,'Fig3C.heatmap.pdf'),width = 8.3,height = 7.2)

all.FAPs.anno <- read.table('~/test/0407/20221027_fig_plot/all.FAPs.annotation.txt',header = T,row.names = 1,sep = '\t')
all.FAPs.anno.ordered <- all.FAPs.anno[rownames(data_f_h.trans),,drop=F]
data_f_h.trans_2 <- data_f_h.trans
rownames(data_f_h.trans_2) <- all.FAPs.anno.ordered$annotation # change subtype names

library(dendextend)
d <- dist(data_f_h.trans_2/2)
hc <- hclust(d)
dend <- as.dendrogram(hc)
dend <- dend %>% set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 5) %>%  # node point size
  # set("leaves_col",colorRampPalette(brewer.pal(n=9,name='YlOrRd'))(36)) %>%
  set("leaves_col",rainbow(36)) %>% 
  set("labels_cex",0.8)
ggd1 <- as.ggdend(dend)
# ggplot(ggd1) 
ggplot(ggd1,horiz=T,offset_labels = -15) + scale_x_reverse()

