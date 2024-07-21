# Peter van Galen, 240511
# Assess the effect of eIF4A1, 2, and 3 in AML cell lines

# Set up
library(tidyverse)

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

# Extract cell line IDs from models
cellLine_ids <- models %>% filter(grepl("MV4-11|HL-60|U-937|OCI-AML3|MOLM-13", CellLineName)) %>%
  select(ModelID, CellLineName)

# Filter for cell lines and genes of interest
mygenes <- c("MCL1", "BCL2", "EIF4A3", "EIF4A2", "EIF4A1", "GAPDH", "ACTB")
expr2 <- expr %>% filter(ModelID %in% cellLine_ids$ModelID) %>%
  left_join(cellLine_ids) %>%
  rename_with(cutf, d = " ") %>%
  select(CellLineName, all_of(mygenes))

# Plot heatmap
p1 <- expr2 %>%
  pivot_longer(cols = -CellLineName, names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = mygenes)) %>%
  ggplot(aes(x = CellLineName, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradientn(colours = heat_colors) +
  labs(title = "Gene Expression Heatmap", x = "Cell Line", y = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

pdf(file = "240511_Expression.pdf", width = 4, height = 4)
p1
dev.off()


