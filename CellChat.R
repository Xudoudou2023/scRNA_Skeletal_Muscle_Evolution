## Fig 6
# cell-cell interaction (CellChat)
library(Seurat)
library(dplyr)
library(CellChat)
library(Cairo)
library(ggplot2)
library(ggalluvial)

setwd('~/test/0407/cellChat3')
##############################
##  cellChat
options(stringsAsFactors = FALSE)
setwd('~/test/0407/cellChat3')
seurat_object <- readRDS('Duroc.FAPs_Myogenic_Immune.rds')
# Idents(seurat_object) <- seurat_object$cell_type
Idents(seurat_object) <- seurat_object$cell_anno
##### Part1 Data input & processing and initialization of CellChat object
data.input <- GetAssayData(seurat_object, assay = "SCT", slot = "data") # normalized data matrix
labels <- Idents(seurat_object)
identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(data = data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
# showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 10) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#### Part2 Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
### Part3, Visualization and systems analysis of cell-cell communication network
# Create a directory to save figures
data.dir <- '~/test/0407/cellChat'
dir.create(data.dir)
setwd(data.dir)
pathways <- cellchat@netP$pathways
for(p in pathways){ 
  # pathways.show <- c("TGFb") 
  
  pathways.show <- p
  levels(cellchat@idents)
  vertex.receiver = seq(1,10) # a numeric vector
  
  CairoPDF(paste0(p,'_pahtways.pdf'),width = 20)
  # Hierarchy plot
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
  # Circle plot
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
  
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  print(netAnalysis_contribution(cellchat, signaling = pathways.show))
  
  # Identify signaling roles of cell groups
  cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the signaling roles of cell groups
  netVisual_signalingRole(cellchat, signaling = pathways.show)
  dev.off()
}
#### 4,Identify global communication patterns and major signals for specific cell groups
# Identify and visualize outgoing communication pattern of secreting cells
CairoPDF(file = 'outgoing.pdf')
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

### Identify and visualize incoming communication pattern of target cells
CairoPDF(file = 'incoming.pdf')
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

##### Manifold and classification learning analysis of signaling networks

CairoPDF(file = 'Embedding.pdf')
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional")
netVisual_embeddingZoomIn(cellchat, type = "functional")

#### Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural")
netVisual_embeddingZoomIn(cellchat, type = "structural")
dev.off()
saveRDS(cellchat, 'cellchat3.rds')


