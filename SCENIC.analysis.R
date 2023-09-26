#### Fig 9 
common.TFs <- read.csv('hsa_mus_ssu.shared.regulons.csv',header = T)
common.TFs <- as.character(gsub('\\(.*','',common.TFs$human.pig.mouse.shared))

ssu <- readRDS('./p.merged.rds')
Idents(ssu) <- 'cell_type'
levels(ssu) <- c('FAPs','MuSCs','Myocytes','Myoblasts','PDGFRb+cells','Glial cells',
                 'Endothelial cells','Macrophages+Neutrophils','T+B cells')
ssu.counts <- ssu[['RNA']]@data
ssu.counts <- ssu.counts[common.TFs,]
ssu.meta <- ssu@meta.data
ssu.meta$anno <- Idents(ssu)
ssu.TFs.data <- apply(ssu.counts,1,function(x) tapply(x,ssu.meta$anno,mean))
ssu.TFs.data <- t(ssu.TFs.data)
ssu.TFs.data <- as.data.frame(ssu.TFs.data)

tmp <- t(apply(ssu.TFs.data,1,scale))
colnames(tmp) <- colnames(ssu.TFs.data)
tmp <- as.data.frame(tmp)
g.tmp <- c()
for(i in levels(ssu)){
  temp <- tmp[which(tmp[,i]>1),i,drop=F]
  temp <- temp[order(-temp[,1]),,drop=F]
  g.tmp <- c(g.tmp,rownames(temp))
  
}
g.tmp <- unique(g.tmp)
tmp <- tmp[g.tmp,]
tmp_2 <- ssu.TFs.data[g.tmp,]
col <- brewer.pal(n = 11,name = 'RdYlBu')

pheatmap(tmp_2,scale = 'row',cluster_rows = F,cluster_cols = F,breaks = seq(-2,3,length.out = 100),
         fontsize_row = 4,
         color = colorRampPalette(colors = c("#74ADD1","white","red"))(100))

hsa <- readRDS('./hsa.anno.rds')
levels(hsa) <- c('LUM+ FAP cells','PRG4+ FAP cells','Satellite Cells','Smooth Muscle cells','Pericytes',
                 'Endothelial cells','PCV Endothelial cells','Myeloid cells','T/B cells','NK cells')
hsa.counts <- hsa[['RNA']]@data
hsa.counts <- hsa.counts[common.TFs,]
hsa.meta <- hsa@meta.data
hsa.meta$anno <- Idents(hsa)
hsa.TFs.data <- apply(hsa.counts,1,function(x) tapply(x,hsa.meta$anno,mean))
hsa.TFs.data <- t(hsa.TFs.data)
write.csv(hsa.TFs.data,'20201103-all.cell.pyscenic/hsa.regulons.exp.csv')
pheatmap(hsa.TFs.data,scale = 'row',cluster_cols = F,cluster_rows = F)

tmp_hsa <- hsa.TFs.data[g.tmp,]
pheatmap(tmp_hsa,scale = 'row',cluster_rows = F,cluster_cols = F,breaks = seq(-2,3,length.out = 100),
         fontsize_row = 4, color = colorRampPalette(colors = c("#74ADD1","white","red"))(100))


homo.gene <- read.table('homo.gene_hsa_mmu_ssu.txt',header = T,row.names = 1)
homo.gene <- data.frame(homo.gene,row.names = 1)
mus.common.TFs <- as.character(homo.gene[common.TFs,]$mmu)
mus <- readRDS('mus.cca.rds')
mus <- subset(mus,idents = '14',invert = T)
levels(mus) <- c('FAP','Scx+','MuSCs','Glial','Endothelial','Macrophage','Neutrophils','T','B','ITGA7+VCAM1-')
mus.counts <- mus[['RNA']]@data
mus.counts <- mus.counts[mus.common.TFs,]
mus.meta <- mus@meta.data
mus.meta$anno <- Idents(mus)
mus.TFs.data <- apply(mus.counts,1,function(x) tapply(x,mus.meta$anno,mean))
mus.TFs.data <- t(mus.TFs.data)


g.tmp_mus <- as.character(homo.gene[g.tmp,]$mmu)
tmp_mus <- mus.TFs.data[g.tmp_mus,]
pheatmap(tmp_mus,scale = 'row',cluster_rows = F,cluster_cols = F,breaks = seq(-2,3,length.out = 100),
         fontsize_row = 4, color = colorRampPalette(colors = c("#74ADD1","white","red"))(100))
