# downstream analysis



# GSEA --------------------------------------------------------------------


work.path <- "/share/genomics/cwx/wang";setwd(work.path)
raw.path <- file.path(work.path, "RawData")
out.path <- file.path(work.path, "OutPut")

geneset <- readRDS("/share/genomics/cwx/public/msigdb.rds")
geneset <- subset(geneset, gs_cat %in% c("H", "C2"))
geneset <- split(geneset$gene_symbol, geneset$gs_name)

DEgenes <- lapply(unique(seu$Cell_type), function(celltype){
  FindMarkers(object = seu, 
              group.by = "Group", 
              ident.1 = "remote", ident.2 = "proximal", subset.ident = celltype, 
              logfc.threshold = 0)
})
names(DEgenes) <- unique(seu$Cell_type)
DEpathway <- lapply(DEgenes, function(DEgene){
  DEgene$feature = rownames(DEgene)
  genes.list = DEgene[, c("feature","avg_log2FC")]
  genes.list = dplyr::arrange(genes.list,-genes.list[,2])
  rank = tibble::deframe(genes.list)
  Res = fgsea::fgseaMultilevel(pathways = geneset, stats = rank)
  Res = subset(Res, padj < 0.05) 
  Res = dplyr::arrange(Res, -Res$padj)
})


# GSVA --------------------------------------------------------------------

mtx <- AverageExpression(seu, group.by = "Cell_type")[[1]]
GSVA.celltype <- GSVA::gsva(expr = mtx, gset.idx.list = geneset)

# Cellchat ----------------------------------------------------------------

library(CellChat)
cellchat = createCellChat(object = seu, group.by = "Cell_type")
CellChatDB <-CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)

CellChat.group <- SplitObject(seu, "Group")
CellChat.group <- lapply(CellChat.group, function(x) RunCellChat(object = x, feature = "Cell_type"))
mergeCellChat.group <- mergeCellChat(CellChat.group)
netVisual_heatmap(mergeCellChat.group)
mergeCellChat.group <- computeNetSimilarityPairwise(mergeCellChat.group, type = "functional")
mergeCellChat.group <- netEmbedding(mergeCellChat.group, type = "functional")
mergeCellChat.group <- netClustering(mergeCellChat.group, type = "functional")
netVisual_embeddingPairwise(mergeCellChat.group, type = "functional", label.size = 3.5)
rankNet(mergeCellChat.group, mode = "comparison", stacked = T, do.stat = TRUE)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

