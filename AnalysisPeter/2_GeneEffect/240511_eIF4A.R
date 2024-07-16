# Peter van Galen, 240511
# Assess the effect of eIF4A1, 2, and 3 in AML cell lines

# Set up
library(tidyverse)

setwd("~/DropboxMGB/GitHub/DepMap/AnalysisPeter/2_GeneEffect")

rm(list=ls())

# My favorites
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))
heat_colors <- c("#083160", "#2668aa", "#4794c1", "#94c5dd", "#d2e5ef", "#f7f7f7",
                 "#fcdbc8", "#f2a585", "#d46151", "#b01b2f", "#660220")

# Plot expression ---------------------------------------------------------------------------------

# Load DepMap data
expr <- read_csv("../../ExpressionAndDependencyData2024/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
expr <- rename(expr, ModelID = `...1`)
models <- read_csv("../../ExpressionAndDependencyData2024/Model.csv", guess_max = 2000)

# Extract cell line IDs from models
cellLine_ids <- models %>% filter(grepl("MV4-11|HL-60|U-937|OCI-AML3|MOLM-13|HAP1", CellLineName)) %>%
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

pdf(file = "Expression.pdf", width = 4, height = 4)
p1
dev.off()


# Plot dependencies for each cell line ------------------------------------------------------------

# Load DepMap data
df <- read_csv("../../ExpressionAndDependencyData2024/CRISPRGeneEffect.csv")
df <- rename(df, ModelID = `...1`)

# Filter for cell lines of interest, transpose, and simplify gene names
df_wrangled <- df %>% filter(ModelID %in% cellLine_ids$ModelID) %>%
  pivot_longer(cols = -ModelID, names_to = "gene", values_to = "dependency") %>%
  mutate(gene = cutf(gene, d = " ")) %>%
  left_join(cellLine_ids) %>%
  select(CellLineName, gene, dependency)

# Plot dependencies
p2 <- df_wrangled %>% ggplot(aes(x = dependency)) +
  geom_histogram(bins = 1000) +
  geom_vline(data = filter(df_wrangled, grepl("EIF4A", gene)), aes(xintercept = dependency, color = gene)) +
  facet_wrap(~CellLineName, ncol = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Save plot
pdf(file = "Dependencies.pdf")
p2
dev.off()


# Plot dependencies for each gene -----------------------------------------------------------------

# AML cells are not more dependent on eIF4A than other cancer cell lines
df %>% select(ModelID, `EIF4A1 (1973)`, `EIF4A2 (1974)`, `EIF4A3 (9775)`) %>%
  pivot_longer(cols = -ModelID, names_to = "gene", values_to = "dependency") %>%
  mutate(Class = factor(ifelse(ModelID %in% filter(models, OncotreeCode == "AML")$ModelID,
                               yes = "AML", no = "Other"), levels = c("Other", "AML"))) %>%
  ggplot(aes(x = dependency, color = Class, fill = Class)) +
  geom_histogram(bins = 1100) +
  facet_wrap(~gene, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())


