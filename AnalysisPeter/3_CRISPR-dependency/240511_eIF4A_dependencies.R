# Peter van Galen, 240511
# Assess the effect of eIF4A1, 2, and 3 in AML cell lines

# Set up
library(tidyverse)

# fmt: skip
setwd("~/Library/CloudStorage/OneDrive-MassGeneralBrigham//GitHub/AML-DepMap-Insights-2024/AnalysisPeter/3_CRISPR-dependency")

# Clear environment variables
rm(list = ls())

# My favorites
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load DepMap model data & extract cell line IDs from models
models <- read_csv("../../DataDownloads/24Q2/Model.csv", guess_max = 2000)
cellLine_ids <- models %>%
  filter(grepl("MV4-11|HL-60|U-937|OCI-AML3|MOLM-13", CellLineName)) %>%
  select(ModelID, CellLineName)

# Load DepMap dependency data
df <- read_csv("../../DataDownloads/24Q2/CRISPRGeneEffect.csv")
df <- rename(df, ModelID = `...1`)

# Filter for cell lines of interest, transpose, and simplify gene names
df_wrangled <- df %>%
  filter(ModelID %in% cellLine_ids$ModelID) %>%
  pivot_longer(cols = -ModelID, names_to = "gene", values_to = "dependency") %>%
  mutate(gene = cutf(gene, d = " ")) %>%
  left_join(cellLine_ids) %>%
  select(CellLineName, gene, dependency)

# Plot dependencies
p2 <- df_wrangled %>%
  ggplot(aes(x = dependency)) +
  geom_histogram(bins = 1000) +
  geom_vline(
    data = filter(df_wrangled, grepl("EIF4A", gene)),
    aes(xintercept = dependency, color = gene)
  ) +
  facet_wrap(~CellLineName, ncol = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Save plot
pdf(file = "240511_Dependencies.pdf")
p2
dev.off()


# Plot dependencies for all AML lines -----------------------------------------

# Select all AML lines
cellLine_ids <- models %>%
  filter(OncotreeCode == "AML") %>%
  select(ModelID, CellLineName)

# Filter for cell lines of interest, transpose, and simplify gene names
df_wrangled <- df %>%
  filter(ModelID %in% cellLine_ids$ModelID) %>%
  pivot_longer(cols = -ModelID, names_to = "gene", values_to = "dependency") %>%
  mutate(gene = cutf(gene, d = " ")) %>%
  left_join(cellLine_ids) %>%
  select(CellLineName, gene, dependency)

# Plot dependencies
df_wrangled %>%
  ggplot(aes(x = dependency)) +
  geom_histogram(bins = 1000) +
  geom_vline(
    data = filter(df_wrangled, grepl("EIF4A", gene)),
    aes(xintercept = dependency, color = gene)
  ) +
  coord_cartesian(xlim = c(-3, 1)) +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text = element_text(color = "black"))

# Save plot
ggsave("250926_Dependencies.pdf", width = 6, height = 3)
