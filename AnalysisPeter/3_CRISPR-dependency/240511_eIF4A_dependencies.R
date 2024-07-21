# Peter van Galen, 240511
# Assess the effect of eIF4A1, 2, and 3 in AML cell lines

# Set up
library(tidyverse)

setwd("~/DropboxMGB/GitHub/AML-DepMap-Insights-2024/AnalysisPeter/3_CRISPR-dependency")

rm(list=ls())

# My favorites
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load DepMap model data & extract cell line IDs from models
models <- read_csv("../../DataDownloads/24Q2/Model.csv", guess_max = 2000)
cellLine_ids <- models %>% filter(grepl("MV4-11|HL-60|U-937|OCI-AML3|MOLM-13", CellLineName)) %>%
  select(ModelID, CellLineName)

# Load DepMap dependency data
df <- read_csv("../../DataDownloads/24Q2/CRISPRGeneEffect.csv")
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
pdf(file = "240511_Dependencies.pdf")
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


