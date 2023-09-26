#### Fig 9G
# q-q plot
# tmp <- read.csv('Becker_Muscular_Dystrophy.csv',header = TRUE)
files <- list.files('.',pattern = 'csv')
genes_disease.df <- c()
for (f in files){
  sample <- sub('.csv','',f)
  tmp <- read.csv(f,header = TRUE,stringsAsFactors = F)
  tmp <- tmp[,c('Disease','Gene')]
  genes_disease.df <- rbind(genes_disease.df,tmp)
}
homo_gene <- read.table('./homo.gene_hsa_mmu_ssu.txt',sep = '\t',header = T)
genes_disease.df <- genes_disease.df[genes_disease.df$Gene %in% homo_gene$hsa,]


write.table(genes_disease.df,'genes_disease.df.txt',sep = '\t',quote = F)

homo_gene <- data.frame(row=homo_gene$hsa,homo_gene,row.names = 1,stringsAsFactors = F)

hsa <- readRDS('./hsa.anno.rds')
hsa_counts <- hsa@assays$RNA@counts
mus <- readRDS('./mus.cca.rds')
mus_counts <- mus@assays$RNA@counts
ssu <- readRDS('./p.merged.rds')
ssu_counts <- ssu@assays$RNA@counts

mus_gene <- rownames(mus_counts)
genes_disease <- genes_disease.df$Gene
inter_gene <- intersect(mus_gene,homo_gene$mmu[homo_gene$hsa %in% genes_disease])
genes_disease.df <- genes_disease.df[genes_disease.df$Gene%in% homo_gene$hsa[homo_gene$mmu%in%inter_gene],]
## test 
n <- 1
plot.list <- list()
for(disease in unique(genes_disease.df$Disease)){
  # disease <- unique(genes_disease.df$Disease)[1]
  gene <- genes_disease.df$Gene[genes_disease.df$Disease==disease]
  gene_mus <- as.character(homo_gene[gene,]$mmu)
  
  all(gene %in% rownames(hsa_counts))
  all(gene %in% rownames(ssu_counts))
  all(gene_mus %in% rownames(mus_counts))
  
  hsa_dis <- hsa_counts[gene,,drop=F]
  ssu_dis <- ssu_counts[gene,,drop=F]
  mus_dis <- mus_counts[gene_mus,,drop=F]
  
  ssu_hsa <- data.frame(ssu=log2(rowMeans(ssu_dis)+1),hsa = log2(rowMeans(hsa_dis)+1))
  p1 <- ggscatter(ssu_hsa,x='ssu',y='hsa',add = 'reg.line',add.params = list(color = 'red'),size = 1,
                  xlim = c(0,ceiling(max(ssu_hsa$ssu))),ylim = c(0,max(ceiling(ssu_hsa$hsa))))
  t.model <- lm(ssu~hsa,data = ssu_hsa)
  r2 <- round(summary(t.model)$r.squared,digits = 2)
  p1 <- p1 + annotate("text",label=bquote(r^2 == .(r2)),x=1,y=3.6,color='red',size=3)+
    theme(axis.title = element_blank()) 
  plot.list[[paste0('sh_',n)]] <- p1
  
  
  ssu_mus <- data.frame(ssu=log2(rowMeans(ssu_dis)+1),mus = log2(rowMeans(mus_dis)+1))
  p1 <- ggscatter(ssu_mus,x='ssu',y='mus',add = 'reg.line',add.params = list(color = 'red'),size = 1,
                  xlim = c(0,ceiling(max(ssu_mus$ssu))),ylim = c(0,max(ceiling(ssu_mus$mus)))) 
  t.model <- lm(ssu~mus,data = ssu_mus)
  r2 <- round(summary(t.model)$r.squared,digits = 2)
  p1 <- p1 + annotate("text",label=bquote(r^2 == .(r2)),x=1,y=3.6,color='red',size=3)+
    theme(axis.title = element_blank())
  plot.list[[paste0('sm_',n)]] <- p1
  
  hsa_mus <- data.frame(hsa=log2(rowMeans(hsa_dis)+1),mus = log2(rowMeans(mus_dis)+1))
  p1 <- ggscatter(hsa_mus,x='hsa',y='mus',add = 'reg.line',add.params = list(color = 'red'),size = 1,
                  xlim = c(0,ceiling(max(hsa_mus$hsa))),ylim = c(0,max(ceiling(hsa_mus$mus))))
  t.model <- lm(hsa~mus,data = hsa_mus)
  r2 <- round(summary(t.model)$r.squared,digits = 2)
  p1 <- p1 + annotate("text",label=bquote(r^2 == .(r2)),x=1,y=3.6,color='red',size=3)+
    theme(axis.title = element_blank()) 
  plot.list[[paste0('hm_',n)]] <- p1
  n <- n+1
}
plot_grid(plotlist = plot.list,ncol = 3,nrow = 6)

