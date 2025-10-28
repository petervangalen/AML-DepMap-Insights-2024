# Peter van Galen, 240719
# Plot correlation between Zotatifin efficacy and eIF4A1 expression

# Set up
library(tidyverse)
library(ggrepel)

# fmt: skip
setwd("~/Library/CloudStorage/OneDrive-MassGeneralBrigham//GitHub/AML-DepMap-Insights-2024/AnalysisPeter/4_DrugAnalysis")

rm(list = ls())

# My favorites
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load DepMap data
expr <- read_csv(
  "../../DataDownloads/24Q2/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv"
)
expr <- rename(expr, ModelID = `...1`)
models <- read_csv("../../DataDownloads/24Q2/Model.csv", guess_max = 2000)
drugs <- read_csv(
  "../../DataDownloads/24Q2/Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv"
)
drugs <- rename(drugs, DrugID = `...1`)
anno <- read_csv(
  "../../DataDownloads/24Q2/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"
)

# Subset cell lines for AML, select relevant info
models_subset <- models %>%
  filter(OncotreePrimaryDisease == "Acute Myeloid Leukemia") %>% # optional
  select(ModelID, CellLineName)

# Wrangle expression
mygene <- "EIF4A1"
mygene <- "JUN"
mygene <- "JUNB"
mygene <- "BCL2"
mygene <- "ERBB2"
expr_subset <- expr %>%
  rename_with(cutf, d = " ") %>%
  select(ModelID, all_of(mygene))

# Extract ID for drug of interest from the table
anno %>% filter(Drug.Name %in% c("ZOTATIFIN", "VENETOCLAX"))
mydrug <- "ZOTATIFIN"
mydrug <- "VENETOCLAX"
mydrug <- "LAPATINIB" # Positive control with gene ERBB2
drug.id <- filter(anno, Drug.Name == mydrug)$IDs

# Wrangle drugs
drug_t <- drugs %>%
  filter(DrugID == drug.id) %>%
  pivot_longer(
    cols = -DrugID,
    names_to = "ModelID",
    values_to = "Drug_LFC"
  ) %>%
  select(-DrugID)

# Combincke all
Join <- models_subset %>%
  left_join(drug_t) %>%
  na.omit %>%
  left_join(expr_subset) %>%
  mutate(
    Label = ifelse(
      CellLineName %in%
        c(
          "HL-60",
          "HEL",
          "THP-1",
          "NOMO-1",
          "OCI-AML3",
          "MOLM-13",
          "SKM-1",
          "TF-1",
          "U-937",
          "MONO-MAC-1"
        ),
      yes = CellLineName,
      no = NA
    )
  )

# Plot
p1 <- Join %>%
  ggplot(aes(x = Drug_LFC, y = .data[[mygene]], label = Label)) +
  geom_point(color = "#696969") +
  geom_smooth(method = "lm", col = "#4682b4", alpha = 0.2) + # Optional: Add a regression line
  labs(
    title = paste0(
      "Correlation: ",
      round(cor(Join$Drug_LFC, Join[[mygene]], use = "complete.obs"), 2),
      ", p = ",
      round(cor.test(Join$Drug_LFC, Join[[mygene]])$p.value, 4)
    )
  ) +
  ylab(paste(mygene, "Expression")) +
  xlab(paste(str_to_sentence(mydrug), "LFC")) +
  geom_text_repel() +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  )

# View
p1

# Save
ggsave(
  paste0("240719_", str_to_sentence(mydrug), "-", mygene, ".pdf"),
  width = 3.5,
  height = 3.5
)
