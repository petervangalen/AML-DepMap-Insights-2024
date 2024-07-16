# Peter van Galen, 220105
# Troubleshoot addModuleScore function

library(tidyverse)
library(Seurat)

setwd("~/DropboxMGB/Projects/FLT3_inhibition/Single-cell/AnalysisYoke/")

rm(list=ls())

seu <- readRDS(file = "5HDBM_Integrated_pca_clusterID.rds")

# Specify genes
genes.ch <- c("CLEC9A", "CLNK", "BATF3", "THBD", "FLT3", "IRF8") # "XCR1" is not in the Seurat object
genes.ch %in% rownames(seu)
length(genes.ch)

# Calculate signature scores
seu <- AddModuleScore(object = seu, features = list(genes.ch), assay = "RNA", name = "cDC1_genes")
seu@meta.data %>% head

FeaturePlot(seu, features = "cDC1_genes1")