library(KEGGREST)
library(dplyr)
library(readxl)
library(pheatmap)

# Load data 
KO_hits <- read_excel("C:/Users/../Github_Rscripts/r-scripts-data-analysis-visualisation/data/kegg_visualisation_bins.xlsx", sheet = "MAGs")

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

# Retrieve KO descriptions from KEGG
ko2name <- keggList("ko")

ko_df <- data.frame(
  KO = names(ko2name),
  KO_description = ko2name,
  stringsAsFactors = FALSE
)

# Add "ko:" prefix to match KO format
ko_df$KO <- paste0("ko:", ko_df$KO)

# Retrieve KO → pathway mapping from KEGG
pathway2ko <- keggLink("pathway", "ko")

pathway2ko_df <- data.frame(
  KO = names(pathway2ko),
  Pathway = pathway2ko,
  stringsAsFactors = FALSE
)

# Standardise pathway IDs (koXXXXX → mapXXXXX)
pathway2ko_df$Pathway <- gsub("^path:ko", "path:map", pathway2ko_df$Pathway)

pathway2ko_df <- distinct(pathway2ko_df)

# Retrieve pathway descriptions
pathway2name <- keggList("pathway")

pathway2name_df <- data.frame(
  Pathway = names(pathway2name),
  Pathway_description = unname(pathway2name),
  stringsAsFactors = FALSE
)

# Add "path:" prefix for consistency
pathway2name_df$Pathway <- paste0("path:", pathway2name_df$Pathway)

# Combine KO descriptions to create a full KEGG annotation table
ko2descriptions2pathway <- ko_df %>%
  full_join(pathway2ko_df, by = "KO", relationship = "many-to-many") %>% 
  full_join(pathway2name_df, by = "Pathway", relationship = "many-to-many")

# Merge KEGG annotations with your KO hit data
KO_hits <- KO_hits %>%
  left_join(ko2descriptions2pathway, by = "KO", relationship = "many-to-many")

# Filter to KOs of interest (subset for plotting)
KO_hits_relevant <- KO_hits %>% 
  filter(KO %in% c("ko:K01866",
                   "ko:K01883",
                   "ko:K02996",
                   "ko:K00287",
                   "ko:K01875",
                   "ko:K01724"))

# Create KO-level presence/absence matrix
# Converts any non-zero value → 1 (presence)
KO_heatmap <- KO_hits_relevant %>%
  group_by(KO, KO_description) %>%
  summarise(across(bin.1:bin.5, ~ as.numeric(any(!is.na(.x) & .x != 0)))) %>%
  ungroup()

# Define KO order for plotting
KO_manual_order <- c("ko:K01724",
                     "ko:K00287",
                     "ko:K01875",
                     "ko:K02996",
                     "ko:K01883",
                     "ko:K01866")

# Reorder KO rows
KO_heatmap <- KO_heatmap[match(KO_manual_order, KO_heatmap$KO), ]

# Remove rows without valid descriptions
KO_heatmap <- KO_heatmap %>% 
  filter(!is.na(KO_description), KO_description != "NA", KO_description != "")

# Clean KO IDs for display (remove "ko:")
KO_heatmap$KO_clean <- sub("^ko:", "", KO_heatmap$KO)

# Extract short enzyme name (before ";")
KO_heatmap$Short_name <- ifelse(
  grepl(";", KO_heatmap$KO_description),
  sub(";.*", "", KO_heatmap$KO_description),
  NA
)

# Extract full description (remove EC numbers)
KO_heatmap$Description_clean <- ifelse(
  grepl(";", KO_heatmap$KO_description),
  sub(" \\[EC:.*\\]", "",
      sub("^[^;]*;\\s*", "", KO_heatmap$KO_description)),
  sub(" \\[EC:.*\\]", "", KO_heatmap$KO_description)
)

# Create final row labels
KO_heatmap$Row_label <- paste(
  KO_heatmap$KO_clean, ";",
  ifelse(is.na(KO_heatmap$Short_name), " ;", paste0(KO_heatmap$Short_name, " ;")),
  KO_heatmap$Description_clean
)

# Convert KO table to matrix for heatmap
n_cols <- ncol(KO_heatmap)
KO_matrix <- as.matrix(KO_heatmap[, 3:(n_cols-4)])

# Define grayscale color palette
heat_colors <- colorRampPalette(c("white", "#ADADAD"))(100)

KO_heatmap <- pheatmap(
  KO_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = heat_colors,
  border_color = "black",
  fontsize = 30,
  legend = FALSE,
  cellwidth = 26.6,
  cellheight = 50,
  labels_row = KO_heatmap$Row_label,
  filename = "heatmap_KO_level.png"
)
