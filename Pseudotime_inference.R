library(pheatmap)
library(monocle)

######## Fig 3A
FAPs <- readRDS('FAPs.rds')
Idents(FAPs) <- 'cell_anno'
all.markers <- FindAllMarkers(FAPs,only.pos = TRUE)
# t <- 20
# top <- all.markers%>% group_by(cluster) %>% top_n(t,wt = avg_logFC)
FAPs <- FindVariableFeatures(FAPs,nfeatures = 200)
top <- VariableFeatures(FAPs)

expr_matrix <- FAPs@assays$SCT@counts
sample.info <- FAPs@meta.data
gene_anno <- data.frame(gene_short_name=rownames(FAPs))
rownames(gene_anno) <- gene_anno$gene_short_name
pd <- new("AnnotatedDataFrame", data = sample.info)
fd <- new("AnnotatedDataFrame", data = gene_anno)
cds <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
# cds_ordering_genes <- unique(top$gene)
cds_ordering_genes <- unique(top)
cds <-setOrderingFilter(cds,ordering_genes = cds_ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds,color_by = 'cell_anno',cell_size = 0.5)

# plot_cell_trajectory(cds,color_by = 'cell_anno',cell_size = 0.5) + theme(legend.position="right")
for(i in unique(plot.data$cell_anno)){
  
  plot.data$gg <- ifelse(plot.data$cell_anno==i & plot.data$Type=='Laiwu',1,
                         ifelse(plot.data$cell_anno == i & plot.data$Type=='Duroc',2,
                                ifelse(plot.data$cell_anno ==i & plot.data$Type=='Wild',3,4)))
  # ggplot(plot.data,aes(x=Component1,y=Component2)) + geom_point()
  plot.data$gg <- factor(plot.data$gg,levels = c(1,2,3,4),ordered = TRUE)
  p_tmp <- plot.data[plot.data$gg != 4,]
  p_tmp$Type <- factor(p_tmp$Type,levels = c('Wild','Laiwu','Duroc'),ordered = TRUE)
  
  p <- ggscatter(plot.data,x='Component1',y='Component2',color = 'lightgrey',size = 0.5,legend = 'right') +
    geom_point(data = p_tmp,aes(x=Component1,y=Component2,color=Type),size=0.2) +
    scale_color_manual(values = c('#387FB9','#E5352F','#5CAD4B')) + labs(title = i)
  print(p)
}

saveRDS(cds,'FAPs.monocle.rds')

########### heatmap 
FAPs_cds <- readRDS('FAPs.monocle.rds')
plot_cell_trajectory(FAPs_cds)

BEAM_res <- BEAM(FAPs_cds, branch_point = 4, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
hm_gene <- row.names(subset(BEAM_res,qval < 1e-4))
pheat <- plot_genes_branched_heatmap(FAPs_cds[hm_gene,],
                                     branch_point = 4,
                                     num_clusters = 5,
                                     cores = 4,
                                     use_gene_short_name = T,
                                     show_rownames = F,
                                     return_heatmap = TRUE
)

ph_orderd <- pheat$heatmap_matrix # extract orderd matrix

# heatmap color function
table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
  x <- seq(0, 1, length.out = n)
  y <- rep(0, length(x))
  sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
  sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
  y[sill.min:sill.max] <- 1
  base.min <- round((n - 1) * (mid - base / 2)) + 1
  base.max <- round((n - 1) * (mid + base / 2)) + 1
  xi <- base.min:sill.min
  yi <- seq(0, 1, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  xi <- sill.max:base.max
  yi <- seq(1, 0, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  height * y
}

rgb.tables <- function(n,
                       red = c(0.75, 0.25, 1),
                       green = c(0.5, 0.25, 1),
                       blue = c(0.25, 0.25, 1))
{
  rr <- do.call("table.ramp", as.list(c(n, red)))
  gr <- do.call("table.ramp", as.list(c(n, green)))
  br <- do.call("table.ramp", as.list(c(n, blue)))
  rgb(rr, gr, br)
}

matlab.like <- function(n) rgb.tables(n)

matlab.like2 <- function(n){
  rgb.tables(n,
             red = c(0.8, 0.2, 1),
             green = c(0.5, 0.4, 0.8),
             blue = c(0.2, 0.2, 1))
}

blue2green2red <- matlab.like2
exp_rng <- range(ph_orderd) #bks is based on the expression range
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
hmcols <- blue2green2red(length(bks) - 1) # final color

row_dist <- as.dist((1 - cor(Matrix::t(ph_orderd)))/2) # cal row distance
row_dist[is.na(row_dist)] <- 1

anno_row <- pheat$annotation_row
anno_col <- data.frame(r=colnames(ph_orderd),
                       Celltype=rep(c('cellfate1','cellfate2','cellfate2','cellfate3'),each=50),
                       row.names = 1)
p_tmp <- pheatmap(ph_orderd,cluster_cols = F,show_rownames = F,show_colnames = F,
                  cutree_rows = 5,gaps_col = c(50,150),
                  annotation_row = anno_row,annotation_col = anno_col,
                  breaks=bks,color=hmcols,
                  clustering_method = 'ward.D2',clustering_distance_rows = row_dist)
row_cluster <- cutree(p_tmp$tree_row,k=5)
gene_heat <- p_tmp$tree_row$labels[p_tmp$tree_row$order]
gene_heat2 <- data.frame(gene_heat,gene=gene_heat,row.names = 1)
gene_heat2$cluster <- row_cluster[rownames(gene_heat2)]
