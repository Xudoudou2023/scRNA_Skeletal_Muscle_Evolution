#### add doublets predicted by Scrublet to Seurat object
p.merged <- readRDS('./p.merged.rds')

### doublet.meta
doublets.dir <- list.files('./rawdata/cr_outs/')
doublets.meta <- c()
for(i in doublets.dir){
  temp.dou <- read.csv(paste0('./rawdata/cr_outs/',i,'/doublet.txt'),header = TRUE,row.names = 1)
  doublets.meta <- rbind(doublets.meta,temp.dou)
}

### Wild FAPs
wild.FAPs <- readRDS('./p/FAPs/pw.FAPs.integrated.rds')
DimPlot(wild.FAPs,label = TRUE)
wild.FAPs$barcode <- colnames(wild.FAPs)

wild.doublets.meta <- doublets.meta[colnames(wild.FAPs),]
all(rownames(wild.doublets.meta)==colnames(wild.FAPs))
wild.FAPs@meta.data <- cbind(wild.FAPs@meta.data,wild.doublets.meta)

p1 <- DimPlot(wild.FAPs,label = TRUE)
p2 <- DimPlot(wild.FAPs,group.by = 'predicted_doublets') + labs(title = 'Predicted doublets')
Cairo::CairoPDF('./Doublets_detection-Wild.FAPs.pdf',width = 12,height = 5)
plot_grid(p1,p2,ncol = 2)
dev.off()

### Laiwu FAPs
laiwu.FAPs <- readRDS('./p/FAPs/ph.FAPs.integrated.rds')
DimPlot(laiwu.FAPs,label = TRUE)
laiwu.FAPs$barcode <- colnames(laiwu.FAPs)

laiwu.doublets.meta <- doublets.meta[colnames(laiwu.FAPs),]
all(rownames(laiwu.doublets.meta)==colnames(laiwu.FAPs))
laiwu.FAPs@meta.data <- cbind(laiwu.FAPs@meta.data,laiwu.doublets.meta)

p1 <- DimPlot(laiwu.FAPs,label = TRUE)
p2 <- DimPlot(laiwu.FAPs,group.by = 'predicted_doublets') + labs(title = 'Predicted doublets')
Cairo::CairoPDF('./Doublets_detection-laiwu.FAPs.pdf',width = 12,height = 5)
plot_grid(p1,p2,ncol = 2)
dev.off()

### duroc FAPs
duroc.FAPs <- readRDS('./p/FAPs/pt.FAPs.integrated.rds')
DimPlot(duroc.FAPs,label = TRUE)
duroc.FAPs$barcode <- colnames(duroc.FAPs)

duroc.doublets.meta <- doublets.meta[colnames(duroc.FAPs),]
all(rownames(duroc.doublets.meta)==colnames(duroc.FAPs))
duroc.FAPs@meta.data <- cbind(duroc.FAPs@meta.data,duroc.doublets.meta)

p1 <- DimPlot(duroc.FAPs,label = TRUE)
p2 <- DimPlot(duroc.FAPs,group.by = 'predicted_doublets') + labs(title = 'Predicted doublets')
Cairo::CairoPDF('./Doublets_detection-duroc.FAPs.pdf',width = 12,height = 5)
plot_grid(p1,p2,ncol = 2)
dev.off()