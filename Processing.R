################ All cells 
library(Seurat)
library(cowplot)
library(dplyr)

cellranger.dir <- '~/cr_output/'
samples <- list.dirs(cellranger.dir)
obj.list <- list()
for(s in samples){
  obj.list[[s]] <- readRDS(paste(cellranger.dir,s,'outs/filtered_gene_bc_matrices',sep = '/'))
}
merged <- merge(x=obj.list[[1]],y=obj.list[2:length(obj.list)])
merged[['percent.mt']] <- PercentageFeatureSet(merged,pattern = '^MT-')
VlnPlot(merged,features = c('nCount_RNA','nFeature_RNA','percent.mt'))
merged <- subset(merged, subset = nFeature_RNA>200 & nFeature_RNA <6000 & percent.mt <10)
# SCT
merged <- SCTransform(merged, vars.to.regress = c("percent.mt",'nCount_RNA'),variable.features.n = 2000)
## 4. linear dimensional reduction
merged <- RunPCA(merged)
merged <- RunTSNE(merged, dims = 1:20)

merged <- FindNeighbors(merged, dims = 1:20)
merged <- FindClusters(merged, resolution = 0.5)

p4 <- DimPlot(merged, reduction = "tsne",label = T)
p5 <- DimPlot(merged, reduction = "tsne",group.by = 'orig.ident')
plot_grid(p4,p5)
##### cell annotation
cell_type_all <- c('FAPs','Satellite cells', 'Myoblasts','Myocytes', 'PDGFRb+ cells', 
                   'Endothelial cells','Glial cells', 'Myeloid', 'Lymphoid')
celltype_gene <- c('DCN','COL1A1','PAX7','MYOD1','ACTC1','TPM1','RGS5','PECAM1','CDH5','PLP1', 
                   'MBP', 'CSF1R', 'S100A8', 'PTPRC', 'CD3G')
marker_merged <- FindAllMarkers(merged,only.pos = T,logfc.threshold = 0.4)

top30 <- marker_merged %>% group_by(cluster) %>% top_n(30,wt = avg_logFC)
merged <- RenameIdents(merged,`1`='FAPs',`2`='FAPs',`3`='FAPs',`4`='FAPs',`10`='FAPs',`18`='FAPs',`21`='FAPs',
                       `0`='Satellite cells',`7`='Satellite cells',`9`='Satellite cells',`5`='Myocytes',`6`='Myocytes',`8`='Myocytes',
                       `12`='Myoblasts',`14`='Myoblasts',`11`='PDGFRb+ cells',`17`='Endothelial cells',
                       `16`='Lymphoid',`13`='Myeloid',`19`='Myeloid',
                       `15`='Glial cells',`20`='Glial cells')
merged$cell_type <- Idents(merged)
saveRDS(merged,'merged.rds')



