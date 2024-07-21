# Peter van Galen, 240718
# For eIF4A1, 2, and 3, plot expression in AML cell lines

# Set up
library(tidyverse)
library(ComplexHeatmap)
library(janitor)

setwd("~/DropboxMGB/GitHub/AML-DepMap-Insights-2024/AnalysisPeter/2_Expression")

rm(list=ls())

# My favorites
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))
heat_colors <- c("#083160", "#2668aa", "#4794c1", "#94c5dd", "#d2e5ef", "#f7f7f7",
                 "#fcdbc8", "#f2a585", "#d46151", "#b01b2f", "#660220")

# Plot expression ---------------------------------------------------------------------------------

# Load DepMap data
expr <- read_csv("../../DataDownloads/24Q2/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")
expr <- rename(expr, ModelID = `...1`)
models <- read_csv("../../DataDownloads/24Q2/Model.csv", guess_max = 2000)
models$OncotreePrimaryDisease %>% tabyl %>% arrange(desc(n)) %>% head(15)

# Inspect cell lines, subset for AML, filter for relevant info
models_subset <- models %>% filter(OncotreePrimaryDisease == "Acute Myeloid Leukemia") %>%
  select(ModelID, OncotreePrimaryDisease, CellLineName)

# Filter for cell lines and genes of interest
mygenes <- factor(c("EIF4A1", "EIF4A2", "EIF4E", "EIF4G1"))
expr2 <- expr %>% filter(ModelID %in% models_subset$ModelID) %>%
  left_join(models_subset) %>%
  rename_with(cutf, d = " ") %>%
  select(CellLineName, all_of(mygenes))

# Create matrix to plot
expr3 <- as.matrix(expr2[,2:ncol(expr2)])
rownames(expr3) <- expr2$CellLineName
head(expr3)

# Generate the heatmap
h1 <- Heatmap(t(expr3),
  name = "Expression",
  col = heat_colors,
  cluster_rows = F,
  cluster_columns = T,
  show_row_names = T,
  show_column_names = T,
  column_names_rot = 45,
  border = T)

# This is nice, but I want to reverse the columns
h1

# Extract & reverse the column dendrogram
ht_drawn <- column_dend(draw(h1))
column_dend_reversed <- rev(column_dend)

# Generate the heatmap with the reversed dendrogram
h1_reversed <- Heatmap(t(expr3),
  name = "Expression",
  col = heat_colors,
  cluster_rows = F,
  cluster_columns = column_dend_reversed,  # Use reversed dendrogram
  show_row_names = T,
  show_column_names = T,
  column_names_rot = 45,
  border = T)

# Draw the heatmap with the reversed dendrogram
h1_reversed

# Save plot
pdf(file = "240720_AML-CellLine-Expr.pdf", width = 10, height = 3)
h1_reversed
dev.off()

