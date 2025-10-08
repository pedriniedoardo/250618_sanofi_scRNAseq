# Load the required libraries ---------------------------------------------
library(patchwork)
library(tidyverse)
library(CellChat)
library(Matrix)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

# cellchat_CTRL <- readRDS("out/cellchat_CTRL_full.rds")

# set up the objects ------------------------------------------------------
# read in the seurat object
seurat_CTRL <- readRDS("../scRNAseq_seurat/data/subset_ctrl.rds")
# fix the P/F naming to avoid downstream errors
ref_meta <- seurat_CTRL@meta.data %>%
  rownames_to_column("barcode")

# ref_meta_fix <- ref_meta %>%
#   mutate(class_fix2 = case_when(class_fix=="P/F"~"PF",
#                                 T ~class_fix))
# 
# # chage the label of teh cells type based on the paper
# # table(ref_meta$seurat_clusters)
# seurat_CTRL$class_fix2 <- ref_meta_fix$class_fix2
# Create a CellChat object ------------------------------------------------
# cellchat_AK <- createCellChat(object = mat_AK, meta = meta_AK, group.by = "class")

cellchat_CTRL <- createCellChat(object = seurat_CTRL, group.by = "cell_type",assay = "RNA")
# number of cells in each cell group
groupSize <- as.numeric(table(cellchat_CTRL@idents))

# Set the ligand-receptor interaction database ----------------------------
# use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
ggsave("images/CellChatDB.pdf",width = 5,height = 5)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB$interaction %>%
  group_by(annotation) %>%
  summarise()

# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat_CTRL@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication an --------
# To infer the cell state-specific communications, we identify over-expressed ligands or receptors in one cell group and then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.

cellchat_CTRL <- subsetData(cellchat_CTRL)

cellchat_CTRL <- identifyOverExpressedGenes(cellchat_CTRL)
cellchat_CTRL <- identifyOverExpressedInteractions(cellchat_CTRL)
# project gene expression data onto PPI network (optional)
cellchat_CTRL <- projectData(cellchat_CTRL, PPI.mouse)

# Part II: Inference of cell-cell communication network -------------------
# CellChat infers the biologically significant cell-cell communication by assigning each interaction with a probability value and peforming a permutation test. CellChat models the probability of cell-cell communication by integrating gene expression with prior known knowledge of the interactions between signaling ligands, receptors and their cofactors using the law of mass action.

# Compute the communication probability and infer cellular communi --------
# if there is an error in this step try to check:
# 1) Use the correct CellChatDB (mouse or human)
# 2) the input data matrix cellchat@data is present and non 0 values
# 3) the data matrix cellchat@data.signaling is present and non 0 values
# 4) the cell group information is correct: unique(cellchat@idents): 1 to 11.
# 5) try to change the name of the idents
cellchat_CTRL <- computeCommunProb(cellchat_CTRL)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
# to include AMA I reduced the min.cells from 10 to 5
cellchat_CTRL <- filterCommunication(cellchat_CTRL, min.cells = 10)

# Extract the inferred cellular communication network as a data fr --------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,

# Infer the cell-cell communication at a signaling pathway level ----------
# CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.
cellchat_CTRL <- computeCommunProbPathway(cellchat_CTRL)
cellchat_CTRL@net
cellchat_CTRL@netP

# Calculate the aggregated cell-cell communication network ----------------
# We can calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.

cellchat_CTRL <- aggregateNet(cellchat_CTRL)
# save the object up to this poitn
# saveRDS(cellchat,file = "cellchat.rds")

# try to plot some communications
groupSize <- as.numeric(table(cellchat_CTRL@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("images/01_number_of_interactions_cellchat_CTRL.pdf",width = 5,height = 5)
netVisual_circle(cellchat_CTRL@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("images/01_strength_of_interactions_cellchat_CTRL.pdf",width = 5,height = 5)
netVisual_circle(cellchat_CTRL@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

# focus on the IMMUNE cells communication
mat <- cellchat_CTRL@net$weight
# pdf("images/01_grid_interactions_cellchat_CTRL.pdf",width = 10,height = 10)
# par(mfrow = c(3,3), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
# dev.off()

pdf("images/01_grid_interactions_cellchat_CTRL.pdf",width = 10,height = 10)
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = F, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# partIII: visualization of ccc networks ----------------------------------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
cellchat_CTRL <- netAnalysis_computeCentrality(cellchat_CTRL, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_CTRL, pattern = "outgoing",height = 16)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_CTRL, pattern = "incoming",height = 16)

pdf("images/02_contribution_cellchat_CTRL.pdf",width = 11,height = 16)
ht1 + ht2
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_CTRL)
#> Signaling role analysis on the cell-cell communication network from user's input
gg1
ggsave("images/02_contribution_all_cellchat_CTRL.pdf",width = 5,height = 5)

# -------------------------------------------------------------------------
# #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# # Signaling role analysis on the cell-cell communication networks of interest
# gg2_IMM <- netAnalysis_signalingRole_scatter(cellchat_CTRL, signaling = c("VISFATIN", "SPP1","PSAP"))
# # gg2 <- netAnalysis_signalingRole_scatter(cellchat_CTRL, signaling = c("VISFATIN"))
# # gg3 <- netAnalysis_signalingRole_scatter(cellchat_CTRL, signaling = c("SPP1"))
# # gg4 <- netAnalysis_signalingRole_scatter(cellchat_CTRL, signaling = c("PSAP"))
# 
# gg2_IMM
# ggsave("images/02_contribution_subset_IMM_OUT_cellchat_CTRL.pdf",width = 5,height = 5)
# # gg3
# # ggsave("images/02_contribution_SPP1_cellchat_CTRL.pdf",width = 5,height = 5)
# # gg4
# # ggsave("images/02_contribution_PSAP_cellchat_CTRL.pdf",width = 5,height = 5)
# 
# #
# pathways.show_IMM <- c("SPP1","PSAP")
# # pathways.show <- c("VISFATIN")
# # Hierarchy plot
# # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
# vertex.receiver = seq(1,3) # a numeric vector. 
# pdf("images/02_contribution_hierarchy_IMM_OUT_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pathways.show_IMM, function(x){
#   netVisual_aggregate(cellchat_CTRL, signaling = x,  vertex.receiver = vertex.receiver,layout = "hierarchy")
# })
# dev.off()
# # also other approach in terms of aesthtuchs
# vertex.receiver = seq(1,3) # a numeric vector. 
# pdf("images/02_contribution_chord_IMM_OUT_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pathways.show_IMM, function(x){
#   netVisual_aggregate(cellchat_CTRL, signaling = x,  vertex.receiver = vertex.receiver,layout = "chord")
# })
# dev.off()
# #
# vertex.receiver = seq(1,3) # a numeric vector. 
# pdf("images/02_contribution_circle_IMM_OUT_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pathways.show_IMM, function(x){
#   netVisual_aggregate(cellchat_CTRL, signaling = x,  vertex.receiver = vertex.receiver,layout = "circle")
# })
# dev.off()
# 
# # Heatmap
# pdf("images/02_contribution_heatmap_IMM_OUT_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pathways.show_IMM, function(x){
#   netVisual_heatmap(cellchat_CTRL, signaling = x, color.heatmap = "Reds")
# })
# dev.off()
# # par(mfrow=c(1,1))
# #> Do heatmap based on a single object
# 
# # Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
# netAnalysis_contribution(cellchat_CTRL, signaling = pathways.show_IMM)
# ggsave("images/02_contribution_braplotLR_IMM_OUT_cellchat_CTRL.pdf",width = 5,height = 5)
# 
# #
# pairLR <- extractEnrichedLR(cellchat_CTRL, signaling = pathways.show_IMM, geneLR.return = FALSE)
# # LR.show <- pairLR[1,] # show one ligand-receptor pair
# # Hierarchy plot
# vertex.receiver = seq(1,3) # a numeric vector
# pdf("images/02_contribution_hierarcyLR_IMM_OUT_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pairLR$interaction_name, function(x){
#   netVisual_individual(cellchat_CTRL, signaling = pathways.show_IMM,  pairLR.use = x, vertex.receiver = vertex.receiver,layout = "hierarchy")
# })
# dev.off()
# #
# pdf("images/02_contribution_circleLR_IMM_OUT_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pairLR$interaction_name, function(x){
#   netVisual_individual(cellchat_CTRL, signaling = pathways.show_IMM,  pairLR.use = x, vertex.receiver = vertex.receiver,layout = "circle")
# })
# dev.off()
# #
# pdf("images/02_contribution_chordLR_IMM_OUT_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pairLR$interaction_name, function(x){
#   netVisual_individual(cellchat_CTRL, signaling = pathways.show_IMM,  pairLR.use = x, vertex.receiver = vertex.receiver,layout = "chord")
# })
# dev.off()
# # -------------------------------------------------------------------------
# # do the same for the one incoming into AST
# # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# # gg1 <- netAnalysis_signalingRole_scatter(cellchat_CTRL)
# #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# # Signaling role analysis on the cell-cell communication networks of interest
# gg2_AST <- netAnalysis_signalingRole_scatter(cellchat_CTRL, signaling = c("CD99", "SPP1","COLLAGEN","LAMININ","FN1","CDH","BMP"))
# #> Signaling role analysis on the cell-cell communication network from user's input
# gg2_AST
# ggsave("images/02_contribution_subset_AST_IN_cellchat_CTRL.pdf",width = 5,height = 5)
# 
# #
# pathways.show_AST <- c("CD99", "SPP1","COLLAGEN","LAMININ","FN1","CDH","BMP")
# # pathways.show <- c("VISFATIN")
# # Hierarchy plot
# # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
# vertex.receiver = seq(1,3) # a numeric vector. 
# pdf("images/02_contribution_hierarchy_AST_IN_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pathways.show_AST, function(x){
#   netVisual_aggregate(cellchat_CTRL, signaling = x,  vertex.receiver = vertex.receiver,layout = "hierarchy")
# })
# dev.off()
# # also other approach in terms of aesthtuchs
# vertex.receiver = seq(1,3) # a numeric vector. 
# pdf("images/02_contribution_chord_AST_IN_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pathways.show_AST, function(x){
#   netVisual_aggregate(cellchat_CTRL, signaling = x,  vertex.receiver = vertex.receiver,layout = "chord")
# })
# dev.off()
# #
# vertex.receiver = seq(1,3) # a numeric vector. 
# pdf("images/02_contribution_circle_AST_IN_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pathways.show_AST, function(x){
#   netVisual_aggregate(cellchat_CTRL, signaling = x,  vertex.receiver = vertex.receiver,layout = "circle")
# })
# dev.off()
# 
# # Heatmap
# pdf("images/02_contribution_heatmap_AST_IN_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pathways.show_AST, function(x){
#   netVisual_heatmap(cellchat_CTRL, signaling = x, color.heatmap = "Reds")
# })
# dev.off()
# # par(mfrow=c(1,1))
# #> Do heatmap based on a single object
# 
# # Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
# netAnalysis_contribution(cellchat_CTRL, signaling = pathways.show_AST)
# ggsave("images/02_contribution_braplotLR_AST_IN_cellchat_CTRL.pdf",width = 5,height = 10)
# 
# #
# pairLR <- extractEnrichedLR(cellchat_CTRL, signaling = pathways.show_AST, geneLR.return = FALSE)
# # LR.show <- pairLR[1,] # show one ligand-receptor pair
# # Hierarchy plot
# vertex.receiver = seq(1,3) # a numeric vector
# pdf("images/02_contribution_hierarcyLR_AST_IN_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pairLR$interaction_name, function(x){
#   netVisual_individual(cellchat_CTRL, signaling = pathways.show_AST,  pairLR.use = x, vertex.receiver = vertex.receiver,layout = "hierarchy")
# })
# dev.off()
# #
# pdf("images/02_contribution_circleLR_AST_IN_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pairLR$interaction_name, function(x){
#   netVisual_individual(cellchat_CTRL, signaling = pathways.show_AST,  pairLR.use = x, vertex.receiver = vertex.receiver,layout = "circle")
# })
# dev.off()
# #
# pdf("images/02_contribution_chordLR_AST_IN_cellchat_CTRL.pdf",width = 10,height = 7)
# lapply(pairLR$interaction_name, function(x){
#   netVisual_individual(cellchat_CTRL, signaling = pathways.show_AST,  pairLR.use = x, vertex.receiver = vertex.receiver,layout = "chord")
# })
# dev.off()

# -------------------------------------------------------------------------
# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# Bubble plot
# We can also show all the significant interactions (L-R pairs) from some cell groups to other cell groups using netVisual_bubble.

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
levels(cellchat_CTRL@idents)
netVisual_bubble(cellchat_CTRL, sources.use = 13, targets.use = c(1:18), remove.isolate = FALSE)
ggsave("images/03_contribution_bubble_from_ENDO_OUT_LR_cellchat_CTRL.pdf",width = 5,height = 17)
# look at what comes into ASTRO
levels(cellchat_CTRL@idents)
netVisual_bubble(cellchat_CTRL, sources.use = c(1:18), targets.use = 13, remove.isolate = FALSE)
ggsave("images/03_contribution_bubble_to_ENDO_IN_LR_cellchat_CTRL.pdf",width = 5,height = 17)

# # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# # show all the interactions sending from Inflam.FIB
# pdf("images/03_contribution_circle_from_IMM_OUT_LR_cellchat_CTRL.pdf",width = 7,height = 7)
# netVisual_chord_gene(cellchat_CTRL, sources.use = 2, targets.use = c(1:7), lab.cex = 0.5,legend.pos.y = 30)
# dev.off()
# 
# pdf("images/03_contribution_circle_to_AST_IN_LR_cellchat_CTRL.pdf",width = 10,height = 10)
# netVisual_chord_gene(cellchat_CTRL, sources.use = c(1:7), targets.use = 1, lab.cex = 0.5,legend.pos.y = 30)
# dev.off()

# Part IV: Systems analysis of cell-cell communication network ------------
# To facilitate the interpretation of the complex intercellular communication networks, CellChat quantitively measures networks through methods abstracted from graph theory, pattern recognition and manifold learning.

# Compute and visualize the network centrality scores
# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# cellchat_CTRL <- netAnalysis_computeCentrality(cellchat_CTRL, slot.name = "netP")

# # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# pdf("images/04_contribution_role_IMM_OUT_cellchat_CTRL.pdf",width = 5,height = 5)
# netAnalysis_signalingRole_network(cellchat_CTRL, signaling = pathways.show_IMM, width = 8, height = 2.5, font.size = 10)
# dev.off()
# 
# pdf("images/04_contribution_role_AST_IN_cellchat_CTRL.pdf",width = 5,height = 5)
# netAnalysis_signalingRole_network(cellchat_CTRL, signaling = pathways.show_AST, width = 8, height = 2.5, font.size = 10)
# dev.off()

# Identify signals contributing most to outgoing or incoming signa --------
# We can also answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups.

# Identify global communication patterns to explore how multiple c --------
# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
# In addition to exploring detailed communications for individual pathways, an important question is how multiple cell groups and signaling pathways coordinate to function. CellChat employs a pattern recognition method to identify the global communication patterns.

# Identify and visualize outgoing communication pattern of secreti --------
# Outgoing patterns reveal how the sender cells (i.e. cells as signal source) coordinate with each other as well as how they coordinate with certain signaling pathways to drive communication.

# library(NMF)
# library(ggalluvial)
# Here we run selectK to infer the number of patterns.
selectK(cellchat_CTRL, pattern = "outgoing")
ggsave("images/selectK_outgoing_cellchat_CTRL.pdf",width = 8,height = 4)
dev.off()
# Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 6.
nPatterns <- 6
pdf("images/patterns_outgoing_cellchat_CTRL.pdf",width = 15,height = 17)
cellchat_CTRL <- identifyCommunicationPatterns(cellchat_CTRL, pattern = "outgoing", k = nPatterns,height = 17,width = 10)
dev.off()
# river plot
netAnalysis_river(cellchat_CTRL, pattern = "outgoing")
ggsave("images/river_outgoing_cellchat_CTRL.pdf",width = 8,height = 10)
# Identify and visualize incoming communication pattern of target  --------
# Incoming patterns show how the target cells (i.e. cells as signal receivers) coordinate with each other as well as how they coordinate with certain signaling pathways to respond to incoming signals.
selectK(cellchat_CTRL, pattern = "incoming")
ggsave("images/selectK_incoming_cellchat_CTRL.pdf",width = 8,height = 4)
dev.off()
# siluette values begin to drop when the number of incoming patterns is 4.
nPatterns <- 5
pdf("images/patterns_incoming_cellchat_CTRL.pdf",width = 15,height = 17)
cellchat_CTRL <- identifyCommunicationPatterns(cellchat_CTRL, pattern = "incoming", k = nPatterns,height = 17,width = 10)
dev.off()
# river plot
netAnalysis_river(cellchat_CTRL, pattern = "incoming")
ggsave("images/river_incoming_cellchat_CTRL.pdf",width = 8,height = 10)
# Manifold and classification learning analysis of signaling netwo --------
# Manifold and classification learning analysis of signaling networks

# Identify signaling groups based on structure similarity
cellchat_CTRL <- computeNetSimilarity(cellchat_CTRL, type = "functional")
cellchat_CTRL <- netEmbedding(cellchat_CTRL, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat_CTRL <- netClustering(cellchat_CTRL, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat_CTRL, type = "structural", label.size = 3.5)
ggsave("images/netVisual_structural_cellchat_CTRL.pdf",width = 5,height = 5)

#
cellchat_CTRL <- computeNetSimilarity(cellchat_CTRL, type = "structural")
cellchat_CTRL <- netEmbedding(cellchat_CTRL, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat_CTRL <- netClustering(cellchat_CTRL, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat_CTRL, type = "structural", label.size = 3.5)
ggsave("images/netVisual_functional_cellchat_CTRL.pdf",width = 5,height = 5)

# Part V: Save the CellChat object ----------------------------------------
saveRDS(cellchat_CTRL, file = "out/cellchat_CTRL_full.rds")