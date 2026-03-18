library(KEGGREST)
library(dplyr)
library(readxl)
library(pheatmap)

# Load data
KO_hits <- read_excel("C:/Users/../../Github_Rscripts/r-scripts-data-analysis-visualisation/data/kegg_visualisation_bins.xlsx", sheet = "MAGs")

# Convert KO IDs to KEGG format
KO_hits$KO <- paste0("ko:", KO_hits$KO)

# Arrange bins in a manual order if needed
col_manual_order <- c(
  "bin.1",
  "bin.2",
  "bin.3",
  "bin.4",
  "bin.5")

KO_hits <- KO_hits[, c("KO", col_manual_order)]

# Retrieve KO → pathway mapping

pathway2ko <- keggLink("pathway", "ko")

pathway2ko_df <- data.frame(
  KO = names(pathway2ko),
  Pathway = pathway2ko,
  stringsAsFactors = FALSE
)

# Standardise pathway IDs

pathway2ko_df$Pathway <- gsub("^path:ko", "path:map", pathway2ko_df$Pathway)

pathway2ko_df <- distinct(pathway2ko_df)

# Merge all pathway KOs with sample data
KO_hits_full <- pathway2ko_df %>%
  left_join(KO_hits, by = "KO")  # keeps all KOs in KEGG pathways

# Replace NAs (KOs not detected in any sample) with 0
sample_cols <- col_manual_order  # your bin columns
KO_hits_full[sample_cols][is.na(KO_hits_full[sample_cols])] <- 0

# Define pathways of interest
KO_hits_in_relevant_pathways <- KO_hits_full%>% 
  filter(Pathway %in% c("path:map00010", "path:map00030")) %>% 
  arrange(Pathway)

# Calculate coverage to check pathway completeness
# Coverage = (# detected KOs) / (total KOs in KEGG pathway)
KO_coverage_per_pathway <- KO_hits_in_relevant_pathways %>%
  group_by(Pathway) %>%
  summarise(
    n_KOs_total = n_distinct(KO),   # all KOs in KEGG for that pathway
    across(all_of(sample_cols), ~ sum(.x == 1) / n_KOs_total)
  ) %>%
  ungroup()

# Convert coverage to presence/absence using a threshold (e.g., 15%)
KO_pathway_presence <- KO_coverage_per_pathway %>%
  mutate(across(
    bin.1:bin.5,
    ~ as.numeric(.x >= 0.03)   # ≥3% of KEGG KOs must be detected, in the paper we chose 15%
  ))

# Define order of pathways for plotting
pathway_manual_order <- c("path:map00030", "path:map00010")

# Reorder rows based on manual pathway order
KO_pathway_presence<- KO_pathway_presence[match(pathway_manual_order, KO_pathway_presence$Pathway), ]

# Clean pathway labels (remove "path:" prefix)
KO_pathway_presence$Pathway <- sub("^path:", "", KO_pathway_presence$Pathway)

# Convert data to matrix for heatmap (exclude metadata columns)
KO_matrix_pathway <- as.matrix(
  KO_pathway_presence[, 3:(ncol(KO_pathway_presence))])

# Ensure columns are in the correct manual order
KO_matrix_pathway <- KO_matrix_pathway[, col_manual_order]

# Define color palette for heatmap (white = 0, grey = 1)
heat_colors <- colorRampPalette(c("white", "#ADADAD"))(100)

# Plot pathway presence/absence heatmap
Heatmap_pathways <- pheatmap(
  KO_matrix_pathway,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = heat_colors,
  border_color = "black",
  fontsize = 30,
  legend = FALSE,
  cellwidth = 26.6,
  cellheight = 50,
  labels_row = KO_pathway_presence$Pathway,
  filename = "Heatmap_pathways.png"
)
