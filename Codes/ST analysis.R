# ST analysis

library(Seurat)
library(harmony)
library(magrittr)
library(SingleR)
library(clustree)
source("/share/genomics/cwx/public/utils.R")

work.path <- "/share/genomics/cwx/wang";setwd(work.path)
raw.path <- file.path(work.path, "RawData")
out.path <- file.path(work.path, "OutPut")

if (!dir.exists(out.path)) dir.create(out.path)


# Data Preparation --------------------------------------------------------


ST.seu <- Load10X_Spatial(data.dir = file.path(raw.path, "HT2021-12389_Report_2022_01_06/1.SpaceRanger/N3"))
VlnPlot(ST.seu, features = c("nCount_Spatial", "nFeature_Spatial"))
SpatialFeaturePlot(ST.seu, features = c("nCount_Spatial", "nFeature_Spatial"))

# ST.seu <- SCTransform(ST.seu, assay = "Spatial")
# ST.seu <- RunPCA(ST.seu, assay = "SCT", verbose = FALSE)
ST.seu <- NormalizeData(ST.seu) %>% FindVariableFeatures() %>% ScaleData()
ST.seu <- RunPCA(ST.seu)
ST.seu <- FindNeighbors(ST.seu, reduction = "pca", dims = 1:30)
ST.seu <- FindClusters(ST.seu, verbose = FALSE)
ST.seu <- RunUMAP(ST.seu, reduction = "pca", dims = 1:30)
SpatialFeaturePlot(ST.seu, features = marker$lineage$Common)
SpatialFeaturePlot(ST.seu, features = c("PECAM1", "CDH5"))

DefaultAssay(ST.seu) <- "Spatial"

# ST.seu <- FindSpatiallyVariableFeatures(ST.seu, assay = "SCT", features = VariableFeatures(ST.seu)[1:1000],
#                                        selection.method = "markvariogram")
# top.features <- head(SpatiallyVariableFeatures(ST.seu, selection.method = "markvariogram"), 6)

anchors <- FindTransferAnchors(reference = seu, query = ST.seu, normalization.method = "LogNormalize")
predictions.assay <- TransferData(anchorset = anchors, refdata = seu$Cell_type, prediction.assay = TRUE,
                                  weight.reduction = ST.seu[["pca"]], dims = 1:30)
DimPlot(ST.seu, group.by = "seurat_clusters")
ST.seu[["predictions"]] <- predictions.assay

VlnPlot(ST.seu, features = names(table(seu$Cell_type)))
DefaultAssay(ST.seu) <- "predictions"
FeaturePlot(ST.seu, features = c("CD4 Tcell"))
SpatialFeaturePlot(ST.seu, features = c("CD4 Tcell", "Mono-Macro", "Endothelial", "Epithelial", "Fibroblast"))

DimPlot(ST.seu, group.by = "seurat_clusters") +
  SpatialDimPlot(ST.seu, group.by = "seurat_clusters")

mtx <- lapply(levels(ST.seu$seurat_clusters), function(x){
  ST.seu$seurat_clusters == x
})
mtx <- do.call(cbind, mtx)
ST.ref <- read.table(file.path(raw.path, "HT2021-12389_Report_2022_01_06", "5.CellType", "T3-celltype.xls"))
rownames(ST.ref) <- gsub("(\\w+)(-)(\\w+)", "\\3", rownames(ST.ref))
all(rownames(ST.ref)==gsub("(\\w+)(-)(\\w+)", "\\1", colnames(ST.seu)))
mtx <- t(as.matrix(ST.ref[, c(2:12)])) %*% mtx

# mtx <- lapply(levels(ST.seu$seurat_clusters), function(x){
#   ST.seu$seurat_clusters == x
# })
# mtx <- do.call(cbind, mtx)
# mtx <- as.matrix(ST.seu@assays$predictions@data) %*% mtx
# colnames(mtx) <- as.character(levels(ST.seu$seurat_clusters))

# saveRDS(ST.seu, file.path(out.path, "T3.rds"))
# saveRDS(ST.seu, file.path(out.path, "N3.rds"))