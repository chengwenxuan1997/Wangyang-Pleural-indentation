# a exploratory analysis for sample 1 and sample 2
# single cell transcriptome
# 2022/01/13

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


# SampleToAnalyse <- c("N1", "N2", "T1", "T2")
SampleToAnalyse <- c("N1", "N2", "N3", "N4", "T1", "T2", "T3", "T4")
# SampleToAnalyse <- c("N3", "T3")

dirs <- list.dirs(file.path(raw.path, "HT2021-12392_Report/1.CellRanger"), full.names = F)
dirs <- dirs[grep("filtered_feature_bc_matrix", dirs)]
dirs <- dirs[-grep("ipynb", dirs)]
dirs <- dirs[grep(pattern = paste(SampleToAnalyse, collapse = "|"), x = dirs)]
dirs <- file.path(raw.path, "HT2021-12392_Report/1.CellRanger", dirs)

expr <- lapply(dirs, function(x){
  File <- Read10X(x)
  colnames(File) <- paste0(basename(dirname(x)), "_", 
                           gsub("([A-Z]+)(-)(\\d)", "\\1", colnames(File)))
  return(File)
})

expr <- do.call(cbind, expr)
summary(colSums(expr))
summary(colSums(expr>0))
summary(rowSums(expr>0))
dim(expr)


# QC, Clustering and visualization --------------------------------------------------


seu <- CreateSeuratObject(expr, min.cells = 10, min.features = 3)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu$Sample <- substr(colnames(seu), 1, 2)
Idents(seu) <- seu$Sample
seu$Batch <- Idents(seu)
seu$Batch <- plyr::mapvalues(x = seu$Batch,
                             from = c("N1", "N2", "N3", "N4", "T1", "T2", "T3", "T4"),
                             to = c("1", "1", "2", "2", "1", "1", "2", "2"))
seu$Group <- Idents(seu)
seu$Group <- plyr::mapvalues(x = seu$Group,
                             from = c("N1", "N2", "N3", "N4", "T1", "T2", "T3", "T4"),
                             to = c(rep("proximal", 4), rep("remote", 4)))
seu$Patient <- Idents(seu)
seu$Patient <- paste0("Patient", substr(seu$Sample, 2,2))

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
seu <- subset(seu, nFeature_RNA >= 200 & nFeature_RNA <=6000 & percent.mt < 30)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
# seu <- seu[, colnames(seu)%in%rownames(annotation)]
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu)
ElbowPlot(seu, ndims = 50)
reduction.to.use = "pca"
dim.to.use = 1:(which.min(diff(seu@reductions[[reduction.to.use]]@stdev[10:50]))+9)
seu <- RunUMAP(seu, dims = dim.to.use, reduction = reduction.to.use, reduction.name = "umap")
seu <- RunTSNE(seu, dims = dim.to.use, reduction = reduction.to.use)
# seu <- RunTSNE(seu, dims = dim.to.use)
FeaturePlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
DimPlot(seu, reduction = "umap", group.by = "Batch")

seu <- FindNeighbors(seu, dims = dim.to.use, reduction = reduction.to.use, graph.name = NULL)
# seu <- CalKSpearmanScore(seu = seu, new_feature = "RawKSscore")

seu <- FindClusters(seu, resolution = 0.4)
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = T)
# seu <- CellCycleScoring(seu,
#                         s.features = cc.genes$s.genes, 
#                         g2m.features = cc.genes$g2m.genes, 
#                         set.ident = TRUE)
# DimPlot(seu, group.by = "Phase")

# Annotation provided by SingleR and OE -----------------------------------

ref.annotation <- readRDS(file.path(raw.path, "annotation.rds"))
ref.annotation$fullid <- paste0(ref.annotation$sampleid, 
                                "_", 
                                do.call(rbind, strsplit(ref.annotation$orig.ident, "-"))[,1])
sum(ref.annotation$fullid %in% colnames(expr))
seu$ref <- ref.annotation$new_celltype[match(colnames(seu), ref.annotation$fullid)]

hpca.se <- reference$HumanPrimaryCellAtlasData
blueprint <- reference$BlueprintEncodeData
typea <- as.data.frame(SingleR(test = seu@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main))
typea <- subset(typea, pruned.labels %in% names(table(typea$pruned.labels))[table(typea$pruned.labels) > 100])
typeb <- as.data.frame(SingleR(test = seu@assays$RNA@data, ref = blueprint, labels = blueprint$label.main))
typeb <- subset(typeb, pruned.labels %in% names(table(typeb$pruned.labels))[table(typeb$pruned.labels) > 100])
seu$hpca.ref <- typea$pruned.labels[match(colnames(seu), rownames(typea))]
seu$blue.ref <- typeb$pruned.labels[match(colnames(seu), rownames(typeb))]
DimPlot(seu, reduction = "umap", group.by = "ref", label = T)
DimPlot(seu, reduction = "umap", group.by = "hpca.ref", label = T)
DimPlot(seu, reduction = "umap", group.by = "blue.ref", label = T)

DimPlot(seu, group.by = "seurat_clusters", label = T)
DimPlot(seu, group.by = "Cell_type", label = T)
FeaturePlot(seu, features = marker$lineage$Common)
FeaturePlot(seu, features = marker$lineage$`Thymus lymphocyte`)
FeaturePlot(seu, features = marker$`Thymus lymphocyte`$NKT)
FeaturePlot(seu, features = marker$lineage$`Bone lymphocyte`)
FeaturePlot(seu, features = marker$lineage$`Monocyte-Macrophage`)
FeaturePlot(seu, features = marker$Macrophage$M2)
FeaturePlot(seu, features = marker$lineage$NK)
FeaturePlot(seu, features = marker$lineage$Neutrophils)
FeaturePlot(seu, features = marker$lineage$Tumor)
FeaturePlot(seu, features = marker$lineage$Mast)

seu$Cell_type <- seu$seurat_clusters
# seu$Cell_type <- plyr::mapvalues(x = seu$Cell_type,
#                                  from = c(0, 2, 3,
#                                           1, 10,
#                                           4, 6, 15,
#                                           7, 8, 9, 12, 13,
#                                           5, 11, 14),
#                                  to = c("NK&T", "NK&T", "NK&T",
#                                         "Bcell", "Bcell", 
#                                         "Mono-Macro", "Mono-Macro", "Neutrophils",
#                                         "Endothelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial",
#                                         "Mast", "Fibroblast", "Epithelial"))
## patient3
seu$Cell_type <- plyr::mapvalues(x = seu$Cell_type,
                                 from = c(1, 2, 6, 0,
                                          10, 4, 3, 8,
                                          5, 7, 12, 9,
                                          11, 13, 15, 14),
                                 to = c("NK&T", "NK&T", "Bcell", "Neutrophils",
                                        "Mono-Macro", "Mono-Macro", "Mono-Macro", "Mono-Macro",
                                        "Epithelial", "Epithelial", "Epithelial", "Epithelial",
                                        "Mast", "Fibroblast", "Endothelial", "Mix"))
Idents(seu) <- seu$seurat_clusters
DEGs <- FindAllMarkers(seu, )
DEG <- DEGs[DEGs$cluster == 4,]

# T cell ------------------------------------------------------------------

T.seu <- subset(seu, Cell_type %in% c("CD4 Tcell", "CD8 Tcell", "NK&T"))
T.seu <- NormalizeData(T.seu) %>% FindVariableFeatures() %>% ScaleData()
T.seu <- RunPCA(T.seu)
# T.seu <- RunHarmony(T.seu, group.by.vars = "Patient")
ElbowPlot(T.seu)
reduction.to.use = "pca"
dim.to.use = 1:15
T.seu <- RunUMAP(T.seu, dims = dim.to.use, reduction = reduction.to.use)
T.seu <- FindNeighbors(T.seu, dims = dim.to.use)
T.seu <- FindClusters(T.seu, resolution = 0.8)
DimPlot(T.seu, group.by = "seurat_clusters", label = T)
FeaturePlot(T.seu, features = marker$`Thymus lymphocyte`$CD4T)
FeaturePlot(T.seu, features = marker$`Thymus lymphocyte`$CD8T)
FeaturePlot(T.seu, features = marker$`Thymus lymphocyte`$NKT)
FeaturePlot(T.seu, features = marker$`Thymus lymphocyte`$Treg)

plot.genes <- c("FOXP3", "CD8A", "CD4", "CD3D", "MKI67", "KLRB1")
# lapply(plot.genes, function(x) CheckMarkers(object = T.seu, feature = x))

T.seu$Cell_subtype <- T.seu$RNA_snn_res.0.8
T.seu$Cell_subtype <- plyr::mapvalues(x = T.seu$Cell_subtype,
                                      from = c(7, 5, 12, 14, 16,
                                               0, 1, 9, 11, 15, 17,
                                               6),
                                      to = c("Treg", rep("CD8T", 3), "CD4T", "CD4 cycling",
                                             rep("NKT", 5),
                                             "CD4"))

res <- GSVA::gsva(expr = expr,
                  gset.idx.list = geneset,
                  parallel.sz=40,
                  kcdf = "Gaussian")
res

pheatmap(mat = mtx,
         color = colorRampPalette(c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F"))(100),
         breaks = seq(0, 5, 0.05),
         angle_col = 45)
grid.lines(x = c(2,3), y = c(4,5))
