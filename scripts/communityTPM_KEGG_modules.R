#Load libraries

library(KEGGREST)
library(dplyr)
library(stringr)
library(readxl)
library(tidyr)
library(pheatmap)
library(purrr)


# Load data and add subscripts to make it match KEGGLink data
KO_hits <- read_excel("/../Documents/R/GitHub_Rscripts/communityTPM_kegg_visualisation.xlsx", sheet = "kofam_tpm")

KO_hits$KO <- paste0("ko:", KO_hits$KO)

# Load KeggLinks which extract the up to date information from the KEGG website

# List all modules and their associated KOs
ko2module <- keggLink("module", "ko")

ko2module_df <- data.frame(
  KO = names(ko2module),
  Module = ko2module,
  stringsAsFactors = FALSE
)

# List all modules and their descriptions 
module_desc <- keggList("module")

module_desc_df <- data.frame(
  Module = names(module_desc),
  Description = as.character(module_desc),
  stringsAsFactors = FALSE
)

# Add md: to all entries for it to match other data sets
module_desc_df$Module <- paste0("md:", module_desc_df$Module)

# List all pathways and their associated KOs
pathway2ko <- keggLink("pathway", "ko")

pathway2ko_df <- data.frame(
  KO = names(pathway2ko),
  Pathway = pathway2ko,
  stringsAsFactors = FALSE
)

# Substitute path:ko with path:map, to match other data sets
pathway2ko_df$Pathway <- gsub("^path:ko", "path:map", pathway2ko_df$Pathway)

# Remove duplicates
pathway2ko_df <- distinct(pathway2ko_df)

# Merge data sets
ko2module2pathway <- ko2module_df %>%
  left_join(module_desc_df, by = "Module")

# Merge data sets to get one final data set with all KOs associated to pathways, modules along with module descriptions
ko2module2pathway <- ko2module2pathway %>%
  full_join(pathway2ko_df, by = "KO", relationship = "many-to-many")

# Merge KO hits for your samples with all KEGG database curated above
KO_hits <- KO_hits %>%
  left_join(ko2module2pathway, by = "KO", relationship = "many-to-many")


# To visualize the data on a module-level, group & summarise the TPM values for all KO hits that belong to a particular module
# We intentionally do NOT normalize by the number of KOs per module because KEGG modules are relatively small 
# and similar in size, so total KO abundance provides a meaningful proxy for module-level functional potential.
# Duplicate KOs are retained (i.e., counted in multiple modules) since many enzymes participate in more than one 
# biochemical route. This redundancy reflects true biological multifunctionality rather than technical duplication
KO_hits_summarised <- KO_hits %>% 
  group_by(Module, Description, Pathway) %>% 
  summarise(across(sample_1:sample_3, sum, na.rm = TRUE)) %>% 
  ungroup()

# Filter out for pathways of interest
KO_hits_summarised_relevant_pathways <- KO_hits_summarised %>% 
  filter(Pathway %in% c("path:map00040", "path:map00053")) %>% 
  arrange(Pathway)

# Filter out for modules of interest
KO_hits_summarised_relevant_modules <- KO_hits_summarised_relevant_pathways %>% 
  filter(Module %in% c("md:M00014", "md:M00129", "md:M00761", "md:M00997", "md:M00999"))

# To plot a module-level heat map, group by module and their descriptions
KO_modules_tpm <- KO_hits_summarised_relevant_modules %>%
  group_by(Module, Description) %>%
  summarise(across(sample_1:sample_3, sum, na.rm = TRUE)) %>%
  ungroup()

# Define a manual order in which you want the modules to be plotted
manual_order <- c("md:M00014",
                  "md:M00129",
                  "md:M00761",
                  "md:M00997",
                  "md:M00999")

# Adjust data set to match the manual order
KO_modules_tpm <- KO_modules_tpm[match(manual_order, KO_modules_tpm$Module), ]

# If you have multiple samples, you can define the manual order in which they are plotted
col_manual_order <- c("sample_3", "sample_2", "sample_1")

# The data set needs to be converted to a matrix in order to be plotted as a heat map
KO_matrix_module <- as.matrix(KO_modules_tpm[, -1:-2])
KO_matrix_module <- KO_matrix_module[, col_manual_order]

# Define colors of the heat map
heat_colors <- colorRampPalette(c("#F5F5DC","#b8d8be", "#c2d6f6", "#9466B0"))(100)

# Final plot
Sample_Heatmap <- pheatmap(
  KO_matrix_module,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = heat_colors,
  fontsize = 15,
  gaps_col = 1:2,
  cellwidth = 20,
  cellheight = 20,
  angle_col = 90,
  labels_row = KO_modules_tpm$Description,
  labels_col = c("Sample 3", "Sample 2", "Sample 1"),
  filename = "Heatmap_modulelevel.png")
