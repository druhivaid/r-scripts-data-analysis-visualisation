# Load packages

library(tidyverse)
library(dplyr)
library(ggplot2)
library(randomcoloR)
library(patchwork)
library(cowplot)

# Load & clean data file
genus_data_fp <- "//../R/GitHub_Rscripts/genus_rel_abundance.xlsx"

genus_data <- readxl::read_xlsx(genus_data_fp, "genera")

genus_data <- janitor::clean_names(genus_data)

# Use factor so taxa plot is always in the order presented in the excel
# E.g. allows you to keep the Unknowns or any other taxa of interest on top of the bar plot etc. 
genus_data$genus <- factor(genus_data$genus, levels = unique(genus_data$genus))

# Manually set colors for taxa of interest
# Need NOT be done for all, those not defined will be colored using the randomcoloR package
genus_manual_colors <- c(
  'genus_1' = "#1B1212",
  'genus_2' = "#e0d6ff")

# Identify remaininy taxa not associated with a color, assign color to them then make a final color list of all taxa
genus_remaining_taxa <- setdiff(genus_data$genus, names(genus_manual_colors))

genus_random_colors <- setNames(distinctColorPalette(length(genus_remaining_taxa)), genus_remaining_taxa)

genus_final_colors <- c(genus_manual_colors, genus_random_colors)

# Calculate sum of relative abundance of a particular taxa across all samples

genus_data <- genus_data %>%
  group_by(genus) %>%
  mutate(sum_percent = rowSums(across(starts_with("sample_"))))

# Filter out taxa that have a relative abundace of > 5 % across samples to be put in the Legend
# This is helpful when there are several taxa across samples
get_genus_labels <- function(df) {
  df %>%
    filter(sum_percent > 5) %>%
    pull(genus)
}

genus_legend <- get_genus_labels(genus_data)

# Pivot data set to be able to plot it, then rename the columns
genus_data  <- genus_data %>% 
  pivot_longer(cols = sample_1:sample_4)

colnames(genus_data)[3] <- 'sample_name'
colnames(genus_data)[4] <- 'rel_abundance'

# First plot
genus_plot_1 <- ggplot(data = genus_data)+
  geom_col(mapping = aes( x = sample_name, y = rel_abundance, fill = genus),
           position = "stack",
           width = 0.4)+
  
  scale_fill_manual(values = genus_final_colors)+
  
  scale_x_discrete(limits=c("sample_1", "sample_2"), 
                   labels=c(sample_1 = "Sample 1", sample_2 = "Sample 2"))+
  
  theme_bw()+
  
  labs(
    subtitle = 'Site 1',
    x = '',
    y = 'Relative Abundance')+
  
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.position = "none",
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
  )

# Second plot
genus_plot_2 <- ggplot(data = genus_data)+
  geom_col(mapping = aes( x = sample_name, y = rel_abundance, fill = genus),
           position = "stack",
           width = 0.4)+
  
  scale_fill_manual(values = genus_final_colors)+
  
  scale_x_discrete(limits=c("sample_3", "sample_4"), 
                   labels=c(sample_3 = "Sample 3", sample_4 = "Sample 4"))+
  
  theme_bw()+
  
  labs(
    subtitle = 'Site 2',
    x = '',
    y = 'Relative Abundance')+
  
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.position = "none",
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 18)
  )

# Only legend plot with the relavent i.e. > 5% relative abundance taxa across samples displayed in the Legend
only_legend_plot <- ggplot(data = genus_data)+
  geom_col(mapping = aes( x = sample_name, y = rel_abundance, fill = genus))+
  
  scale_fill_manual(values = genus_final_colors,
                    breaks = genus_legend)+
  theme_bw()+
  
  theme(legend.position = "right",
        legend.text = element_text(size = 18),
        legend.title = element_text(size =18))+
  
  guides(fill = guide_legend(ncol = 1, title = "Genus"))

# Extract the legend from the only_legend_plot
legend <- get_legend(only_legend_plot)

legend_plot <- plot_grid(legend)

# Merge all plots together 
merged_plot <- (genus_plot_1 + genus_plot_2 + legend_plot) + plot_spacer() + plot_layout(nrow = 1)

# Save plot as a high resolution png
ggsave("genera_plot.png", merged_plot, device="png", width = 30, height = 5, dpi = 600, bg = "transparent")


