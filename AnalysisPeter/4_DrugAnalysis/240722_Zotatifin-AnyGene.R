# Peter van Galen, 240719
# Plot correlation between Zotatifin efficacy and eIF4A1 expression

# Set up
library(tidyverse)
library(ggrepel)

setwd("~/DropboxMGB/GitHub/AML-DepMap-Insights-2024/AnalysisPeter/4_DrugAnalysis")

rm(list=ls())

# My favorites
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Plot expression ---------------------------------------------------------------------------------

# Load DepMap data
expr <- read_csv("../../DataDownloads/24Q2/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")
expr <- rename(expr, ModelID = `...1`)
models <- read_csv("../../DataDownloads/24Q2/Model.csv", guess_max = 2000)
drugs <- read_csv("../../DataDownloads/24Q2/Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv")
drugs <- rename(drugs, DrugID = `...1`)
anno <- read_csv("../../DataDownloads/24Q2/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")

# Optional: subset for AML cell lines
expr <- expr %>% subset(ModelID %in% filter(models, OncotreePrimaryDisease == "Acute Myeloid Leukemia")$ModelID)

# Extract zotatifin drug ID from table
anno %>% filter(Drug.Name == "ZOTATIFIN")
zotatifin.id <- filter(anno, Drug.Name == "ZOTATIFIN")$IDs

# Extract zotatifin data from drug screen
zot_tib <- drugs %>% filter(DrugID == zotatifin.id) %>%
  pivot_longer(cols = -DrugID, names_to = "ModelID", values_to = "ZotatifinLFC") %>%
  select(-DrugID)
zot_df <- data.frame(zot_tib, row.names = 1)

# Create expression matrix from expression tibble
expr_mat <- expr %>% rename_with(cutf, d = " ") %>% select(-ModelID) %>% as.matrix
rownames(expr_mat) <- expr$ModelID
expr_mat[1:3,1:3]

# Put the expression matrix in the same order as the tibble
shared_models <- intersect(rownames(na.omit(zot_df)), rownames(expr_mat))
expr_mat_order <- expr_mat[shared_models,]
zot_df_order <- zot_df[shared_models,,drop=F]
expr_mat_order[1:3,1:3]
head(zot_df_order)

# Correlate sensitivity with expression
all_cors <- apply(expr_mat_order, 2, function(x) cor(x, zot_df_order$ZotatifinLFC))
all_cors_tib <- tibble(Gene = colnames(expr_mat_order),
                       Correlation = all_cors) %>%
  arrange(desc(Correlation))

# Genes to highlight
genes_of_interest <- c(all_cors_tib$Gene[1:5],
                       all_cors_tib$Gene[(nrow(all_cors_tib)-5):nrow(all_cors_tib)],
                       "EIF4A1", "EIF4A2", "EIF4A3", "AKT1", "AKT2", "AKT3", "STAT5A", "STAT5B")

# Visualize ranked sensitivities
all_cors_tib %>% mutate(Rank = 1:n()) %>%
  mutate(Label = ifelse(Gene %in% genes_of_interest, yes = Gene, no = "")) %>%
  ggplot(aes(x = Rank, y = Correlation, label = Label)) +
  geom_point() +
  geom_text_repel(aes(label = Label), max.overlaps = 10000) +
  theme_bw()

# Look at TP53 correlation more closely
merge(zot_df_order, expr_mat_order[,"TP53",drop=F], by = "row.names") %>%
  rename(ModelID = "Row.names") %>%
  left_join(select(models, ModelID, CellLineName)) %>%
  ggplot(aes(x = ZotatifinLFC, y = TP53, label = CellLineName)) +
  geom_point() +
  geom_text_repel() +
  theme_bw() +
  theme(aspect.ratio = 1)
  
# Save
write_tsv(all_cors_tib, file = "240722_Zotatifin_sensitivity_vs_expression.txt")

