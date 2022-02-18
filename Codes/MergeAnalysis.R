library(Seurat)
library(magrittr)

# data <- t(as.matrix(seu@assays$RNA@counts))
# a <- SeuratObject::as.Neighbor(seu@graphs$RNA_nn)
# knn <- a@nn.idx
# batch.estimate <- kBET(df = data, batch = seu$Batch, knn = knn, k0 = 20)


# Overview of the Atlas ---------------------------------------------------


seu <- readRDS("/share/genomics/cwx/wang/RawData/seurat.rds")
cell.col <- c("Epithelial cells" = "#E41A1C", "T lymphocytes" = "#377EB8",
              "NK cells" = "#4DAF4A", "Fibroblasts" = "#45AADD", 
              "Myeloid cells" = "#984EA3", "B lymphocytes" = "#FF7F00",
              "Endothlial cells" = "#A65628", "MAST cells" = "#F781BF", "Mixed" = "#999999")

seu <- RunUMAP(seu, reduction = "mnn", dims = 1:20)
seu <- RunTSNE(seu, reduction = "mnn", dims = 1:20)
DimPlot(seu, reduction = "tsne", group.by = "seurat_clusters", label = T)
DimPlot(seu, reduction = "tsne", group.by = "new_celltype", label = T)
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = T)
DimPlot(seu, reduction = "umap", group.by = "new_celltype", label = T)

seu$sampleid <- plyr::mapvalues(x = seu$sampleid,
                                from = c("N1", "N2", "N3", "N4", "T1", "T2", "T3", "T4"),
                                to = c("001-PRO", "002-PRO", "003-PRO", "004-PRO",
                                       "001-RMT", "002-RMT", "003-RMT", "004-RMT"))
seu$group <- plyr::mapvalues(x = seu$group, from = c("N", "T"), to = c("PRO", "RMT"))
seu$celltypes <- seu$seurat_clusters
seu$celltypes <- plyr::mapvalues(x = seu$celltypes,
                                 from = c(7, 9, 19, 15, 17, 20,
                                        10, 11, 6, 23,
                                        3, 8, 13, 14, 16, 21, 22,
                                        5, 12,
                                        0, 1, 2, 4, 18),
                                 to = c(rep("Epithelial cells", 6),
                                        "Endothlial cells", "Fibroblasts", "MAST cells", "Mixed",
                                        rep("Myeloid cells", 7),
                                        rep("B lymphocytes", 2),
                                        rep("T lymphocytes", 4),
                                        "NK cells"
                                 ))
DimPlot(seu, reduction = "umap", group.by = "celltypes", cols = cell.col) 

data = seu@meta.data[, c("sampleid", "group", "celltypes")]
plot.data <- data[, c("sampleid")]


# Epithelial --------------------------------------------------------------

Epi.seu <- subset(seu, celltypes %in% "Epithelial cells")
Epi.seu <- FindVariableFeatures(Epi.seu) %>% ScaleData()
Epi.seu <- RunPCA(Epi.seu)
ElbowPlot(Epi.seu)
Epi.seu <- RunUMAP(Epi.seu, dims = 1:15)
Epi.seu <- FindNeighbors(Epi.seu, dims = 1:15)
Epi.seu <- FindClusters(Epi.seu, resolution = 0.4)
DimPlot(Epi.seu, label = T)
FeaturePlot(Epi.seu, features = c("SFTPC", "LAMP3", "AGER", "FOXJ1", "RFX2", "SCGB1A1", "MUC16", "KRT19", "PECAM1"))
Epi.seu$sub_celltype <- Epi.seu$seurat_clusters
Epi.seu$sub_celltype <- plyr::mapvalues(x = Epi.seu$sub_celltype,
                                        from = c(0, 2, 4, 5, 8, 1, 3,
                                                 11, 4, 7),
                                        to = c(rep("AT2", 5), rep("AT1", 2),
                                               "Ciliated", "Club"))

# Stromal -----------------------------------------------------------------


Fibro.seu <- subset(seu, celltypes %in% c("Fibroblasts"))
Fibro.seu <- FindVariableFeatures(Fibro.seu) %>% ScaleData()
Fibro.seu <- RunPCA(Fibro.seu)
ElbowPlot(Fibro.seu)
Fibro.seu <- RunUMAP(Fibro.seu, dims = 1:15)
Fibro.seu <- FindNeighbors(Fibro.seu, dims = 1:15)
Fibro.seu <- FindClusters(Fibro.seu, resolution = 0.4)
DimPlot(Fibro.seu, label = T)
# "COL14A1", "GSN", "PI16", "CYGB", "PRRX1",
# "COL13A1", "TCF21", "ITGA8", "CXCL14", "NPNT",
# "TAGLN", "ACTA2", "ACTG2", "MYH1", "MYLK",
# "SYNPO2", "CRYAB", "CNN1", "DES",
# "UPK3B", "MSLN", "CALB2", 
# "CYP1B1", "APOD",
# "RGS5", "CSPG4", "ABCC9", "KCNJ8"
FeaturePlot(Fibro.seu, features = c("COL14A1", "COL13A1", "ACTA2", "PDGFRA", "RGS5"))
FeaturePlot(Fibro.seu, features = c("CXCL12", "SOD2", "COL1A2", "TAGLN"))

EC.seu <- subset(seu, celltypes %in% c("Endothlial cells"))
EC.seu <- FindVariableFeatures(EC.seu) %>% ScaleData()
EC.seu <- RunPCA(EC.seu)
ElbowPlot(EC.seu)
EC.seu <- RunUMAP(EC.seu, dims = 1:15)
EC.seu <- FindNeighbors(EC.seu, dims = 1:15)
EC.seu <- FindClusters(EC.seu, resolution = 0.4)
DimPlot(EC.seu, label = T)
FeaturePlot(EC.seu, features = c("CD3D", "CD4", "CD8A", "FOXP3", "NKG7"))

# Myeloid -----------------------------------------------------------------


Mye.seu <- subset(seu, celltypes %in% "Myeloid cells")
Mye.seu <- FindVariableFeatures(Mye.seu) %>% ScaleData()
Mye.seu <- RunPCA(Mye.seu)
ElbowPlot(Mye.seu)
Mye.seu <- RunUMAP(Mye.seu, dims = 1:15)
Mye.seu <- FindNeighbors(Mye.seu, dims = 1:15)
Mye.seu <- FindClusters(Mye.seu, resolution = 0.4)
DimPlot(Mye.seu, label = T)
FeaturePlot(T.seu, features = c("CD3D", "CD4", "CD8A", "FOXP3", "NKG7"))

# immune cells ------------------------------------------------------------

T.seu <- subset(seu, celltypes %in% c("NK cells", "T lymphocytes"))
T.seu <- FindVariableFeatures(T.seu) %>% ScaleData()
T.seu <- RunPCA(T.seu)
ElbowPlot(T.seu)
T.seu <- RunUMAP(T.seu, dims = 1:15)
T.seu <- FindNeighbors(T.seu, dims = 1:15)
T.seu <- FindClusters(T.seu, resolution = 0.8)
DimPlot(T.seu, label = T)
FeaturePlot(T.seu, features = c("CD3D", "CD4", "CD8A", "FOXP3", "NKG7", "CCR7", "PDCD1", "NCR3"))
VlnPlot(T.seu, features = c("CD3D", "CD4", "CD8A", "FOXP3", "NKG7", "CCR7", "PDCD1", "NCR3"))
plot.marker <- c("CD3D", "CD3E", "CD3G",
                 "XCL1", "FCGR3A", "KLRD1", "KLRF1",
                 "CD8A", "CD8B",
                 "CD4", "IL7R",
                 "TCF7", "SELL", "LEF1", "CCR7",
                 "LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4",
                 "IL2", "GZMA", "GNLY", "PRF1", "GZMB", "GZMK", "IFNG", "NKG7",
                 "IL2RA", "FOXP3", "IKZF2", "TGFB1", "TGFB3", "TGFBI", "TGFBR1",
                 "MAF", "CXCL13", "CXCR5", "PDCD1",
                 "IRF4", "CREM", "NR4A2",
                 "STAT4", "IL12RB2", "IFNG",
                 "GATA3", "STAT6", "IL4",
                 "TRDC", "TRGC2", "TRGC1")
mtx <- standarize.fun(indata = AverageExpression(T.seu, features = plot.marker)[[1]], halfwidth = 2)
T.seu$sub_celltype <- T.seu$seurat_clusters
T.seu$sub_celltype <- plyr::mapvalues(x = T.seu$seurat_clusters,
                                      from = c(7, 10, 11, 9, 14,
                                               2, 4, 6, 8, 16, 17, 18,
                                               1, 12, 15,
                                               0, 3, 5, 13
                                               ),
                                      to = c(rep("NK cells", 3), "Treg", "Exhausted CD8T",
                                             rep("Cytotoxic CD8+ T", 6), "Cytotoxic CD4+ T",
                                             rep("Naive T", 3),
                                             rep("NKT cells", 3), "undetermined"))

B.seu <- subset(seu, celltypes %in% c("B lymphocytes"))
B.seu <- FindVariableFeatures(B.seu) %>% ScaleData()
B.seu <- RunPCA(B.seu)
ElbowPlot(B.seu)
B.seu <- RunUMAP(B.seu, dims = 1:15)
B.seu <- FindNeighbors(B.seu, dims = 1:15)
B.seu <- FindClusters(B.seu, resolution = 0.4)
DimPlot(B.seu, label = T)
FeaturePlot(B.seu, features = c("MS4A1", "HLA-DRA","IGHA1", "JCHAIN", "IGHG1"))
B.seu$sub_celltype <- B.seu$seurat_clusters
B.seu$sub_celltype <- plyr::mapvalues(x = B.seu$sub_celltype,
                                      from = c(2, 4, 9,
                                               0, 1, 3, 5, 6, 7, 8),
                                      to = c(rep("Plasma cells", 3),
                                             rep("B cells", 7)))
