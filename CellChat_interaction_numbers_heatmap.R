## Fig 6
## heatmap
library(Cairo)
library(pheatmap)
# laiwu
setwd('~/test/0407/cellChat3/Laiwu/')
cellchat_Laiwu <- readRDS('cellchat3.rds')
count_Laiwu <- cellchat_Laiwu@net$count

Laiwu_count <- count_Laiwu
for(row in 1:nrow(Laiwu_count)){
  for(col in 1:ncol(Laiwu_count)){
    if(row == col){
      next
    }
    else if(row < col){Laiwu_count[row,col] <- sum(Laiwu_count[row,col],Laiwu_count[col,row])}
    else{Laiwu_count[row,col] <- Laiwu_count[col,row]}
  }
}
ordered_Laiwu <- c(3,1,5,7,4,9,11,8,2,6,10,13,14,12,15,19,16,18,17,23,22,20,21,27,28,29,25,26,24)
Laiwu_count <- Laiwu_count[ordered_Laiwu,ordered_Laiwu]
anno_col <- data.frame(colnames(Laiwu_count),
                       type = c(rep('FAPs',11),rep('Myogenic linage',8),rep('Immune',10)),row.names = 1)

CairoPNG(filename = 'Laiwu.interaction.heatmap-2.png',width = 1200,height = 900)
print(pheatmap(Laiwu_count,cluster_rows = F,cluster_cols = F,annotation_col = anno_col,main = 'Laiwu',
               breaks = seq(0,200,length.out = 100)))
dev.off()

# duroc
setwd('~/test/0407/cellChat3/Duroc/')
cellchat_Duroc <- readRDS('cellchat3.rds')
count_Duroc <- cellchat_Duroc@net$count

Duroc_count <- count_Duroc
for(row in 1:nrow(Duroc_count)){
  for(col in 1:ncol(Duroc_count)){
    if(row == col){
      next
    }
    else if(row < col){Duroc_count[row,col] <- sum(Duroc_count[row,col],Duroc_count[col,row])}
    else{Duroc_count[row,col] <- Duroc_count[col,row]}
  }
}
ordered_Duroc <- c(1,2,8,4,6,10,11,5,7,9,3,15,17,13,16,18,14,12,20,22,21,19,27,29,26,25,28,23,24)
Duroc_count <- Duroc_count[ordered_Duroc,ordered_Duroc]

anno_col <- data.frame(colnames(Duroc_count),
                       type = c(rep('FAPs',11),rep('Myogenic linage',7),rep('Immune',11)),row.names = 1)
CairoPNG(filename = 'Duroc.interaction.heatmap-2.png',width = 1200,height = 900)
print(pheatmap(Duroc_count,cluster_rows = F,cluster_cols = F,annotation_col = anno_col,main = 'Duroc',
               breaks = seq(0,200,length.out = 100)))
dev.off()

# wild
setwd('~/test/0407/cellChat3/Wild/')
cellchat_Wild <- readRDS('cellchat3.rds')
count_Wild <- cellchat_Wild@net$count

Wild_count <- count_Wild
for(row in 1:nrow(Wild_count)){
  for(col in 1:ncol(Wild_count)){
    if(row == col){
      next
    }
    else if(row < col){Wild_count[row,col] <- sum(Wild_count[row,col],Wild_count[col,row])}
    else{Wild_count[row,col] <- Wild_count[col,row]}
  }
}
ordered_Wild <- c(2,1,7,6,4,5,8,3,13,15,9,12,16,18,11,14,10,17,22,21,19,20,25,26,27,28,29,23,24)
Wild_count <- Wild_count[ordered_Wild,ordered_Wild]

anno_col <- data.frame(colnames(Wild_count),
                       type = c(rep('FAPs',8),rep('Myogenic linage',10),rep('Immune',11)),row.names = 1)
CairoPNG(filename = 'Wild.interaction.heatmap-2.png',width = 1200,height = 900)
print(pheatmap(Wild_count,cluster_rows = F,cluster_cols = F,annotation_col = anno_col,main = 'Wild',
               breaks = seq(0,200,length.out = 100)))
dev.off()
