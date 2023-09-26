
####### Fig 9D
### search specific/shared genes
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(pheatmap)
### read species data
sub_pig <- readRDS('./sub.pig.rds')
sub_hsa <- readRDS('./sub.hsa.rds')
sub_mus <- readRDS('./sub.mus.rds')
## homo gene
homo_gene <- read.table('./homo.gene_hsa_mmu_ssu.txt',header = T,row.names = 1)
homo_gene$homo <- paste(homo_gene$hsa,homo_gene$mmu,sep = '/')

all.markers.pig <- read.table('.//pig/all.marers.txt',
                              sep = '\t',header = TRUE,row.names = 1)
all.markers.hsa <- read.table('.//hsa/all.markers.hsa.txt',
                              sep = '\t',header = T,row.names = 1)
all.markers.mus <- read.table('.//mus/all.markers.mus.txt',
                              sep = '\t',header = T, row.names = 1)
all.markers.pig <- all.markers.pig[all.markers.pig$gene %in% homo_gene$ssu,]
all.markers.hsa <- all.markers.hsa[all.markers.hsa$gene %in% homo_gene$hsa,]
all.markers.mus <- all.markers.mus[all.markers.mus$gene %in% homo_gene$mmu,]
top50.pig <- all.markers.pig %>% group_by(cluster) %>% top_n(50,wt = avg_logFC)
top50.hsa <- all.markers.hsa %>% group_by(cluster) %>% top_n(50,wt = avg_logFC)
top50.mus <- all.markers.mus %>% group_by(cluster) %>% top_n(50,wt = avg_logFC)

top50.pig$gene <- as.vector(top50.pig$gene)
top50.hsa$gene <- as.vector(top50.hsa$gene)
top50.mus$gene <- as.vector(top50.mus$gene)
########################## 
## overlap celltype
## pig <- FAPs,MuSCs,PDGFRb+ cells,Endothelial cells,Macrophages+Neutrophils,T+B cells
## hsa <- LUM+ FAP cells/PRG4+ FAP cells,Satellite cells,Smooth Muscle cells,Endothelial cells/PCV Endothelial cells,Myeloid cells,T/B cells/NK cells
## mus <- FAP/Scx+,MuSCs,ITGA7+ VCAM1-,Endothelial,Macrophages/Neutrophils,T/B
# (1) FAPs
FAPs.pig <- top50.pig[grep('FAP',top50.pig$cluster),]$gene
FAPs.pig.h <- homo_gene[homo_gene$ssu %in% FAPs.pig,]$homo

FAPs.hsa <- top50.hsa[grep('FAP',top50.hsa$cluster),]$gene
FAPs.hsa.h <- homo_gene[homo_gene$hsa %in% FAPs.hsa,]$homo

FAPs.mus <- top50.mus[grep('FAP|Scx+',top50.mus$cluster),]$gene
FAPs.mus.h <- homo_gene[homo_gene$mmu %in% FAPs.mus,]$homo

FAP.shared <- intersect(intersect(FAPs.pig.h,FAPs.hsa.h),FAPs.mus.h)

FAP.specific.pig <- FAPs.pig.h[!(FAPs.pig.h %in% FAPs.hsa.h | FAPs.pig.h %in% FAPs.mus.h)]
FAP.specific.hsa <- FAPs.hsa.h[!(FAPs.hsa.h %in% FAPs.pig.h | FAPs.hsa.h %in% FAPs.mus.h)]
FAP.specific.mus <- FAPs.mus.h[!(FAPs.mus.h %in% FAPs.pig.h | FAPs.mus.h %in% FAPs.hsa.h)]

# (2) MuSCs
MuSCs.pig <- top50.pig[grep('MuSCs',top50.pig$cluster),]$gene
MuSCs.pig.h <- homo_gene[homo_gene$ssu %in% MuSCs.pig,]$homo

MuSCs.hsa <- top50.hsa[grep('Satellite Cells',top50.hsa$cluster),]$gene
MuSCs.hsa.h <- homo_gene[homo_gene$hsa %in% MuSCs.hsa,]$homo

MuSCs.mus <- top50.mus[grep('MuSCs',top50.mus$cluster),]$gene
MuSCs.mus.h <- homo_gene[homo_gene$mmu %in% MuSCs.mus,]$homo

MuSCs.shared <- intersect(intersect(MuSCs.pig.h,MuSCs.hsa.h),MuSCs.mus.h)

MuSCs.specific.pig <- MuSCs.pig.h[!(MuSCs.pig.h %in% MuSCs.hsa.h | MuSCs.pig.h %in% MuSCs.mus.h)]
MuSCs.specific.hsa <- MuSCs.hsa.h[!(MuSCs.hsa.h %in% MuSCs.pig.h | MuSCs.hsa.h %in% MuSCs.mus.h)]
MuSCs.specific.mus <- MuSCs.mus.h[!(MuSCs.mus.h %in% MuSCs.pig.h | MuSCs.mus.h %in% MuSCs.hsa.h)]

# (3) muscle
muscle.pig <- top50.pig[grep('PDGFRb\\+cells',top50.pig$cluster),]$gene
muscle.pig.h <- homo_gene[homo_gene$ssu %in% muscle.pig,]$homo

muscle.hsa <- top50.hsa[grep('Smooth Muscle cells',top50.hsa$cluster),]$gene
muscle.hsa.h <- homo_gene[homo_gene$hsa %in% muscle.hsa,]$homo

muscle.mus <- top50.mus[grep('ITGA7\\+VCAM1-',top50.mus$cluster),]$gene
muscle.mus.h <- homo_gene[homo_gene$mmu %in% muscle.mus,]$homo

muscle.shared <- intersect(intersect(muscle.pig.h,muscle.hsa.h),muscle.mus.h)

muscle.specific.pig <- muscle.pig.h[!(muscle.pig.h %in% muscle.hsa.h | muscle.pig.h %in% muscle.mus.h)]
muscle.specific.hsa <- muscle.hsa.h[!(muscle.hsa.h %in% muscle.pig.h | muscle.hsa.h %in% muscle.mus.h)]
muscle.specific.mus <- muscle.mus.h[!(muscle.mus.h %in% muscle.pig.h | muscle.mus.h %in% muscle.hsa.h)]

# (4) Endo
Endo.pig <- top50.pig[grep('Endothelial',top50.pig$cluster),]$gene
Endo.pig.h <- homo_gene[homo_gene$ssu %in% Endo.pig,]$homo

Endo.hsa <- top50.hsa[grep('Endothelial',top50.hsa$cluster),]$gene
Endo.hsa.h <- homo_gene[homo_gene$hsa %in% Endo.hsa,]$homo

Endo.mus <- top50.mus[grep('Endothelial',top50.mus$cluster),]$gene
Endo.mus.h <- homo_gene[homo_gene$mmu %in% Endo.mus,]$homo

Endo.shared <- intersect(intersect(Endo.pig.h,Endo.hsa.h),Endo.mus.h)

Endo.specific.pig <- Endo.pig.h[!(Endo.pig.h %in% Endo.hsa.h | Endo.pig.h %in% Endo.mus.h)]
Endo.specific.hsa <- Endo.hsa.h[!(Endo.hsa.h %in% Endo.pig.h | Endo.hsa.h %in% Endo.mus.h)]
Endo.specific.mus <- Endo.mus.h[!(Endo.mus.h %in% Endo.pig.h | Endo.mus.h %in% Endo.hsa.h)]

# (5) Myeloid
Myeloid.pig <- top50.pig[grep('Macrophage|Neutrophils',top50.pig$cluster),]$gene
Myeloid.pig.h <- homo_gene[homo_gene$ssu %in% Myeloid.pig,]$homo

Myeloid.hsa <- top50.hsa[grep('Myeloid',top50.hsa$cluster),]$gene
Myeloid.hsa.h <- homo_gene[homo_gene$hsa %in% Myeloid.hsa,]$homo

Myeloid.mus <- top50.mus[grep('Macrophage|Neutrophils',top50.mus$cluster),]$gene
Myeloid.mus.h <- homo_gene[homo_gene$mmu %in% Myeloid.mus,]$homo

Myeloid.shared <- intersect(intersect(Myeloid.pig.h,Myeloid.hsa.h),Myeloid.mus.h)

Myeloid.specific.pig <- Myeloid.pig.h[!(Myeloid.pig.h %in% Myeloid.hsa.h | Myeloid.pig.h %in% Myeloid.mus.h)]
Myeloid.specific.hsa <- Myeloid.hsa.h[!(Myeloid.hsa.h %in% Myeloid.pig.h | Myeloid.hsa.h %in% Myeloid.mus.h)]
Myeloid.specific.mus <- Myeloid.mus.h[!(Myeloid.mus.h %in% Myeloid.pig.h | Myeloid.mus.h %in% Myeloid.hsa.h)]

# (6) lymphocyte
lymphocyte.pig <- top50.pig[grep('T\\+B cells',top50.pig$cluster),]$gene
lymphocyte.pig.h <- homo_gene[homo_gene$ssu %in% lymphocyte.pig,]$homo

lymphocyte.hsa <- top50.hsa[grep('T/B|NK cells',top50.hsa$cluster),]$gene
lymphocyte.hsa.h <- homo_gene[homo_gene$hsa %in% lymphocyte.hsa,]$homo

lymphocyte.mus <- top50.mus[top50.mus$cluster=='T' | top50.mus$cluster=='B',]$gene
lymphocyte.mus.h <- homo_gene[homo_gene$mmu %in% lymphocyte.mus,]$homo

lymphocyte.shared <- intersect(intersect(lymphocyte.pig.h,lymphocyte.hsa.h),lymphocyte.mus.h)

lymphocyte.specific.pig <- lymphocyte.pig.h[!(lymphocyte.pig.h %in% lymphocyte.hsa.h | lymphocyte.pig.h %in% lymphocyte.mus.h)]
lymphocyte.specific.hsa <- lymphocyte.hsa.h[!(lymphocyte.hsa.h %in% lymphocyte.pig.h | lymphocyte.hsa.h %in% lymphocyte.mus.h)]
lymphocyte.specific.mus <- lymphocyte.mus.h[!(lymphocyte.mus.h %in% lymphocyte.pig.h | lymphocyte.mus.h %in% lymphocyte.hsa.h)]

#### integrate gene
shared.genes <- c(FAP.shared,MuSCs.shared,muscle.shared,
                  Endo.shared,Myeloid.shared,lymphocyte.shared)
specific.gene.pig <- c(FAP.specific.pig,MuSCs.specific.pig,muscle.specific.pig,
                       Endo.specific.pig,Myeloid.specific.pig,lymphocyte.specific.pig)
specific.gene.hsa <- c(FAP.specific.hsa,MuSCs.specific.hsa,muscle.specific.hsa,
                       Endo.specific.hsa,Myeloid.specific.hsa,lymphocyte.specific.hsa)
specific.gene.mus <- c(FAP.specific.mus,MuSCs.specific.mus,muscle.specific.mus,
                       Endo.specific.mus,Myeloid.specific.mus,lymphocyte.specific.mus)

color <- colorRampPalette(rev(c("yellow", "gray", "blue")))(100)
################################ pig ####################################
ann_cell_pig <- Idents(sub_pig) # cell annotation
ordered_pig <- c() # same cell type put together
for(L in levels(Idents(sub_pig))){
  ordered_pig <- c(ordered_pig,which(ann_cell_pig== L))
}
# ordered <- c(which(anno_cell=='FAPs'),which(anno_cell=='MuSCs'),
#             which(anno_cell=='Endothelial cells'),which(anno_cell=='PDGFRb+cells'),
#             which(anno_cell=='Macrophages+Neutrophils'),which(anno_cell=='T+B cells'))
pig.data <- as.matrix(sub_pig@assays$RNA@scale.data)[,ordered_pig]
### cal gaps index
sum <- 0
gaps_row <- c()
for(L in levels(Idents(sub_pig))){
  sum <- sum+length(which(ann_cell_pig==L))
  gaps_row <- c(gaps_row,sum)
}
gaps_col <- c()

#### pig heatmap plot
s.g <- gsub('/.*$','',shared.genes)
sp.g <- gsub('/.*$','',specific.gene.pig)
sh.g <- gsub('/.*$','',specific.gene.hsa)
sm.g <- gsub('/.*$','',specific.gene.mus)
select.gene <- c(s.g,sp.g,sh.g,sm.g)
p1 <- pheatmap(t(pig.data[unique(select.gene),]),
               color = color,
               cluster_rows = F,cluster_cols = F,
               show_rownames = F,
               breaks = seq(-1,2,length.out = 100),
               annotation_row = as.data.frame(ann_cell_pig),
               gaps_row = gaps_row)
################################ hsa #############################
ann_cell_hsa <- Idents(sub_hsa)
ordered_hsa <- c()
for(L in levels(Idents(sub_hsa))){
  ordered_hsa <- c(ordered_hsa,which(ann_cell_hsa== L))
}
hsa.data <- as.matrix(sub_hsa@assays$RNA@scale.data)[,ordered_hsa]
### cal gaps index
sum <- 0
gaps_row <- c()
for(L in levels(Idents(sub_hsa))){
  sum <- sum+length(which(ann_cell_hsa==L))
  gaps_row <- c(gaps_row,sum)
}
gaps_row
gaps_row <- c(584,942,1055,2315,2352)
gaps_col <- c()
p2 <- pheatmap(t(hsa.data[unique(select.gene),]),
               color = color,
               cluster_rows = F,cluster_cols = F,
               show_rownames = F,
               breaks = seq(-1,2,length.out = 100),
               annotation_row = as.data.frame(ann_cell_hsa),
               gaps_row = gaps_row)
###################################### mus #################
ann_cell_mus <- Idents(sub_mus)
ordered_mus <- c()
for(L in levels(Idents(sub_mus))){
  ordered_mus <- c(ordered_mus,which(ann_cell_mus== L))
}
mus.data <- as.matrix(sub_mus@assays$RNA@scale.data)[,ordered_mus]
### cal gaps index
sum <- 0
gaps_row <- c()
for(L in levels(Idents(sub_mus))){
  sum <- sum+length(which(ann_cell_mus==L))
  gaps_row <- c(gaps_row,sum)
}
gaps_row
gaps_row <- c(1818,2794,3334,4334,4868)
s.g <- gsub('^.*/','',shared.genes)
sp.g <- gsub('^.*/','',specific.gene.pig)
sh.g <- gsub('^.*/','',specific.gene.hsa)
sm.g <- gsub('^.*/','',specific.gene.mus)
select.gene <- c(s.g,sp.g,sh.g,sm.g)
select.gene <- select.gene[select.gene %in% rownames(sub_mus)]
p3 <- pheatmap(t(mus.data[unique(select.gene),]),
               color = color,
               cluster_rows = F,cluster_cols = F,
               show_rownames = F,
               breaks = seq(-1,2,length.out = 100),
               annotation_row = as.data.frame(ann_cell_mus),
               gaps_row = gaps_row)

