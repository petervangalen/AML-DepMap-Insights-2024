# Peter van Galen, 220502
# Visualize the latest reference

library(tidyverse)
library(Seurat)
library(limma)
library(data.table)

# Set up
setwd("~/DropboxMGB/GitHub/AML-DepMap-Insights-2024/AnalysisPeter/1_Projection")
rm(list=ls())
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load data ---------------------------------------------------------------------------------------
# This script is broken because these files no longer exist
seu.all <- readRDS("~/DropboxMGB/Projects/DrugSynergy/Single-cell/AnalysisYoke/5BM_Integrated_annotated.rds")
seu.tnk <- readRDS("~/DropboxMGB/Projects/DrugSynergy/Single-cell/AnalysisYoke/5BM_Integrated_T cells_annotated.rds")

# Quick visualization
DimPlot(seu.all) + theme(aspect.ratio = 1)
DimPlot(seu.tnk) + theme(aspect.ratio = 1)

# Save marker genes
all_markers <- FindAllMarkers(seu.all, slot = "data", only.pos = T)
all_markers.dt.ls <- lapply(split(all_markers, f = all_markers$cluster), function(x) data.table(x))
all_markers.dt.ls <- lapply(all_markers.dt.ls, function(x) setorder(x, -avg_log2FC))
all_markers.df <- do.call(cbind, lapply(all_markers.dt.ls, function(x) x$gene[1:50]))
all_markergenesOrder.df <- all_markers.df[,c("HSC", "EarlyEryth", "LateEryth", "ProMono", "Mono", "ncMono", "cDC", "pDC",
  "ProB", "PreB", "B cells", "Plasma", "CD4 Naive", "CD4 Memory", "CD8 Memory", "gamma delta-like T", "NK-like T", "NK cells")]
write.table(all_markergenesOrder.df, file = "All_markerGenes.txt", sep = "\t", col.names = T, row.names = F, quote = F)

tnk_markers <- FindAllMarkers(seu.tnk, slot = "data", only.pos = T)
tnk_markers.dt.ls <- lapply(split(tnk_markers, f = tnk_markers$cluster), function(x) data.table(x))
tnk_markers.dt.ls <- lapply(tnk_markers.dt.ls, function(x) setorder(x, -avg_log2FC))
tnk_markers.df <- do.call(cbind, lapply(tnk_markers.dt.ls, function(x) x$gene[1:50]))
tnk_markergenesOrder.df <- tnk_markers.df[,c("CD4 Naive", "CD4 Memory", "CD8 Naive", "CD8 Memory", "CD8 Term Exh", "Treg", "gamma delta-like", "NK cells", "NKT cells")]
write.table(tnk_markergenesOrder.df, file = "TNK_markerGenes.txt", sep = "\t", col.names = T, row.names = F, quote = F)


# Merge cell type annotations ---------------------------------------------------------------------
all( colnames(seu.tnk) %in% colnames(seu.all) )

# Make metadata tibbles
all_metadata_tib <- as_tibble(seu.all@meta.data, rownames = "cell")
all_metadata_tib$CellType <- seu.all@active.ident
tnk_metadata_tib <- as_tibble(seu.tnk@meta.data, rownames = "cell")
tnk_metadata_tib$CellType <- seu.tnk@active.ident

# Join
merge_metadata_tib <- all_metadata_tib %>% left_join(select(tnk_metadata_tib, cell, CellType), by = "cell")
merge_metadata_tib <- merge_metadata_tib %>% mutate(CellType = ifelse(is.na(CellType.y), as.character(CellType.x), as.character(CellType.y)))
merge_metadata_tib <- merge_metadata_tib %>% mutate(CellType = gsub("HSC", "HSPC", CellType))
merge_metadata_tib$CellType %>% table(useNA = "always")

# Overwrite with more granular clusters
seu.all$CellType <- NULL
seu.all@meta.data$CellType <- factor(merge_metadata_tib$CellType, levels = c("HSPC", "EarlyEryth", "LateEryth", "ProMono", "Mono", "ncMono", "cDC", "pDC",
  "ProB", "PreB", "B cells", "Plasma",
  "CD4 Naive", "CD4 Memory", "CD8 Naive", "CD8 Memory", "CD8 Term Exh", "Treg",
  "gamma delta-like", "NKT cells", "NK cells"))
seu.all@active.ident <- seu.all$CellType

# Checks
DefaultAssay(seu.all) <- "RNA"
GetAssayData(seu.all, slot = "count")[1:10,1:3]
GetAssayData(seu.all, slot = "data")[1:10,1:3] # This doesn't make sense. The data slot is not normalized
expr_mat <- GetAssayData(seu.all, slot = "count")
plot(colSums(expr_mat), pch = ".", col = factor(seu.all$orig.ident))
seu.all$replicate <- cutf(colnames(seu.all), d = "_", f = 2)
seu.all$replicate %>% table
cutf(colnames(seu.all), d= "-", f = 2) %>% table

# Improve cell IDs
new_names <- colnames(seu.all)
new_names <- gsub("BM191227_5", "BM191227_1",
  gsub("BM191227_6", "BM191227_2",
  gsub("SBM1064_4", "SBM1064_1",
  gsub("SBM1064_5", "SBM1064_2", new_names))))
seu.all <- RenameCells(seu.all, new.names = new_names)
# Improve metadata
seu.all$replicate <- cutf(colnames(seu.all), d= "_", f = 2)
seu.all@meta.data %>% head
seu.all@meta.data <- seu.all@meta.data[,c("orig.ident", "replicate", "nCount_RNA", "nFeature_RNA",
                                          "nCount_SCT", "nFeature_SCT", "CellType")]
seu.all@meta.data %>% head
# What is the data?
plot(colSums(GetAssayData(seu.all, slot = "count")), pch = ".",
     col = factor(paste0(seu.all$orig.ident, "_", seu.all$replicate)))
plot(colSums(GetAssayData(seu.all, slot = "data")), pch = ".",
     col = factor(paste0(seu.all$orig.ident, "_", seu.all$replicate)))


# Visualize ---------------------------------------------------------------------------------------

# From https://www.r-bloggers.com/2013/02/the-paul-tol-21-color-salute/
tol21rainbow <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

#pdf(file = "CellTypes.pdf")
DimPlot(seu.all, cols = tol21rainbow) + theme(aspect.ratio = 1)
#dev.off()

FeaturePlot(seu.all, features = "GZMK") + theme(aspect.ratio = 1)

VlnPlot(seu.all, features = c("GAPDH", "KDM1A"))

