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

# Subset cell lines for AML, select relevant info
models_subset <- models %>% filter(OncotreePrimaryDisease == "Acute Myeloid Leukemia") %>%
  select(ModelID, CellLineName)

# Wrangle expression
mygene <- "EIF4A1"
ylabel <- paste(mygene, "Expression")
expr_subset <- expr %>% rename_with(cutf, d = " ") %>%
  select(ModelID, all_of(mygene))

# Extract zotatifin drug ID from table
anno %>% filter(Drug.Name == "ZOTATIFIN")
zotatifin.id <- filter(anno, Drug.Name == "ZOTATIFIN")$IDs

# Wrangle drugs
drug_t <- drugs %>%
  filter(DrugID == zotatifin.id) %>%
  pivot_longer(cols = -DrugID, names_to = "ModelID", values_to = "Zotatifin_LFC") %>%
  select(-DrugID)

# Combine all
Join <- models_subset %>% left_join(drug_t) %>% na.omit %>%
  left_join(expr_subset) %>%
  mutate(Label = ifelse(CellLineName %in% c("HL-60", "HEL", "THP-1", "NOMO-1", "OCI-AML3", "MOLM-13",
                                            "SKM-1", "TF-1", "U-937", "MONO-MAC-1"),
                        yes = CellLineName, no = NA))

# Plot
p1 <- Join %>%
  ggplot(aes(x = Zotatifin_LFC, y = EIF4A1, label = Label)) +
  geom_point(color = "#696969") +
  geom_smooth(method = "lm", col = "#4682b4", alpha = 0.2) +  # Optional: Add a regression line
  labs(title = paste0("Correlation: ", round(cor(Join$Zotatifin_LFC, Join$EIF4A1), 2),
       ", p = ", round(cor.test(Join$Zotatifin_LFC, Join$EIF4A1)$p.value, 4))) +
  ylab(ylabel) +
  geom_text_repel() +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))

pdf("240719_Zotatifin-eIF4A.pdf", width = 3.5, height = 3.5)
p1
dev.off()

