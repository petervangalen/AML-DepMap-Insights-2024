# Peter van Galen, 220111
# Integrate & cluster bone marrow data

library(tidyverse)
library(Seurat)
library(sctransform)
library(readxl)
library(gridExtra)
#library(patchwork)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/FLT3_inhibition/Single-cell/AnalysisPeter/")

# Function to split character string
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load data ---------------------------------------------------------------------------------------

# Load up RDS files (with PCA analyses) for 5 healthy BM samples
BM180221 = readRDS("../AnalysisYoke/BM180221.rds")
BM191119 = readRDS("../AnalysisYoke/BM191119.rds")
BM191227 = readRDS("../AnalysisYoke/BM191227.rds")
SBM1014 = readRDS("../AnalysisYoke/SBM1014.rds")
SBM1064 = readRDS("../AnalysisYoke/SBM1064.rds")

# Avoid duplicate cell names later
BM180221 <- RenameCells(BM180221, new.names = gsub("-1", "-BM180221", colnames(BM180221)))
BM191119 <- RenameCells(BM191119, new.names = gsub("-1", "-BM191119", colnames(BM191119)))
BM191227 <- RenameCells(BM191227, new.names = gsub("-1", "-BM191227", colnames(BM191227)))
SBM1014 <- RenameCells(SBM1014, new.names = gsub("-1", "-SBM1014", colnames(SBM1014)))
SBM1064 <- RenameCells(SBM1064, new.names = gsub("-1", "-SBM1064", colnames(SBM1064)))

bm.list <- list(BM191119, BM191227, BM180221, SBM1014, SBM1064)
rm(BM191119, BM191227, BM180221, SBM1014, SBM1064)


# Normalize, scale, integrate ---------------------------------------------------------------------

# sctransform normalization (https://satijalab.org/seurat/articles/sctransform_vignette.html)
bm.list <- lapply(X = bm.list, FUN = function(x) {
  x <- SCTransform(x)
})

# Integration
features <- SelectIntegrationFeatures(object.list = bm.list)
bm.list <- PrepSCTIntegration(object.list = bm.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = bm.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

# Initial dimensionality reduction
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)

# Score signatures --------------------------------------------------------------------------------

# Load signature genes & colors
signatures <- read_excel("~/DropboxMGB/Pipelines/AuxiliaryFiles/200727_signatures.xlsx")[-1,]
heatcol <- c("#083160", "#2668aa", "#4794c1", "#94c5dd", "#d2e5ef", "#f7f7f7", "#fcdbc8", "#f2a585", "#d46151", "#b01b2f", "#660220")

# Change the tibble to a list (simultaneously exclude NAs)
signatures <- lapply(as.list(signatures), function(x) x[!is.na(x)])

# Only keep genes that are also in the Seurat object
# Note that the current default matrix in the Seurat object only has 2000 genes..!
signatures <- lapply(signatures, function(x) x[x %in% rownames(immune.combined.sct)])
signatures <- signatures[grep("Griffin", names(signatures))]

# Add scores for each signature
for (s in seq_along(signatures)) {
  print(names(signatures[s]))
  immune.combined.sct <- AddModuleScore(object = immune.combined.sct,
                                        features = signatures[s],
                                        assay = "SCT", name = names(signatures)[s])
}

# UMAP, cluster & visualize -----------------------------------------------------------------------

# Test UMAP and cluster generation with different parameters
dims <- c(30, 40, 50)
resolutions <- c(0.5, 1)

# Test UMAP and clustering with different parameters
for (d in dims) {
  for (r in resolutions) {
    #d <- 15
    #r <- 0.5
    
    # Calculate UMAP coordinates and clusters
    immune.combined.tests <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:d)
    immune.combined.tests <- FindNeighbors(immune.combined.tests, reduction = "pca", dims = 1:d)
    immune.combined.tests <- FindClusters(immune.combined.tests, resolution = r)
    
    # Generate UMAPs colored by sample and cluster
    g1 <- DimPlot(immune.combined.tests, reduction = "umap", group.by = "orig.ident") +
      theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    g2 <- DimPlot(immune.combined.tests, reduction = "umap", label = TRUE, repel = TRUE) +
      ggtitle("clusters") +
      theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    # Generate plots colored by signature score    
    sign_plots <- lapply(names(signatures), function(x)
      FeaturePlot(immune.combined.tests, features = paste0(x, "1"), cols = heatcol, raster = T) +
        ggtitle(cutf(x, d = "_", f = 3)) +
        theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              legend.position = "none", plot.title = element_text(size = 10)))
    
    # Save plots
    pdf(paste0("220111_Test_dims_res/dim", d, "_res", r, ".pdf"))
    grid.arrange(g1, g2, ncol = 1)
    grid.arrange(grobs = sign_plots, ncol = 4)
    dev.off()
    
    rm(immune.combined.tests)
    }
}





