#!/usr/bin/Rscript

# Investigate convergence for the Bayesian chains for the Yeast datasets. 
# Visualise the distributions of the continuous parameters, posterior similarity
# matrices and time series for the predicted clustering in the timecourse data.

# For similarity matrix
library(mdiHelpR)

# Pipe
library(magrittr)

# Data manipulation
library(dplyr)

# Tibbles can hold list columns
library(tibble)

# Pivot
library(tidyr)

# String manipulation
library(stringr)

# Visualisation
library(ggplot2)
library(pheatmap)
library(patchwork)
library(scales)

set.seed(1)
setMyTheme()

# === Read in data =============================================================


# Original data modelled
data_dir <- "./Data/Yeast/Original_data/"
datasets <- c("Timecourse", "ChIP-chip", "PPI")
data_files <- list.files(data_dir, full.names = T)[c(2, 1, 3)]
orig_data <- data_files %>%
  lapply(read.csv, row.names = 1) %>%
  set_names(datasets)

# Save plots in this directory
save_dir <- "./SupplementaryMaterial/Images/Yeast/"

timecourse_col_pal <- dataColPal()
bin_col_pal <- simColPal()

timecourse_breaks <- defineDataBreaks(orig_data$Timecourse, timecourse_col_pal)
bin_breaks <- defineBreaks(bin_col_pal, lb = 0, ub = 1)

col_pal <- list(timecourse_col_pal, bin_col_pal, bin_col_pal)
breaks <- list(timecourse_breaks, bin_breaks, bin_breaks)

compareSimilarityMatricesAnnotated(matrices = orig_data, col_pal = col_pal, breaks = breaks)

# for (dataset in datasets){
N <- nrow(orig_data$Timecourse)
item_names <- row.names(orig_data$Timecourse)

# Set up the X and Y co-ordinates for the entries on the heatmap
row_order <- findOrder(orig_data$Timecourse)
item_order <- item_names[row_order]

timecourse_plt_data <- orig_data$Timecourse %>%
  rownames_to_column(var = "Item_i") %>%
  pivot_longer(-Item_i, values_to = "Value", names_to = "Timepoint") %>%
  mutate(Dataset = "Timecourse",
         Y = match(Item_i, item_order),
         X = as.numeric(str_remove_all(Timepoint, "X"))
         )

if(prod(dim(orig_data$Timecourse)) != nrow(timecourse_plt_data)){
  stop("Timecourse: Mismatch in size. Data has lost information in transformation.\n")
}

p_timecourse <- timecourse_plt_data %>% 
  ggplot(aes(x = X, y = Y, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#146EB4", mid = "white", high = "#FF9900") +
  labs(
    title = "Timecourse data",
    # subtitle = paste0("Posterior similarity matrices (simulation ", j, ")"),
    x = "Timepoint",
    y = "Gene",
    fill = "Gene expression"
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_text(size = 10.5),
    axis.title.x = element_text(size = 10.5),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 10.5),
    legend.text = element_text(size = 10.5),
    panel.border =element_blank()
  )

# Set up the X and Y co-ordinates for the entries on the heatmap
# row_order <- findOrder(orig_data$PPI)
col_order <- findOrder(t(orig_data$PPI))
item_order <- item_names[row_order]
item_col_order <- colnames(orig_data$PPI)[col_order]

ppi_plt_data <- orig_data$PPI %>%
  rownames_to_column(var = "Item_i") %>%
  pivot_longer(-Item_i, values_to = "Interaction", names_to = "Item_j") %>%
  mutate(Dataset = "PPI", Y = match(Item_i, item_order), X = match(Item_j, item_col_order))

# ppi_plt_data$Item_i %>% unique() %>% length()

if(prod(dim(orig_data$PPI)) != nrow(ppi_plt_data)){
  stop("PPI: Mismatch in size. Data has lost information in transformation.\n")
}

p_ppi <- ppi_plt_data %>% 
  ggplot(aes(x = X, y = Y, fill = Interaction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  labs(
    title = "PPI data",
    # subtitle = paste0("Posterior similarity matrices (simulation ", j, ")"),
    x = "Protein",
    y = "",
    fill = "Protein interaction"
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_text(size = 10.5),
    axis.title.x = element_text(size = 10.5),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 10.5),
    legend.text = element_text(size = 10.5),
    panel.border =element_blank()
  )


# Set up the X and Y co-ordinates for the entries on the heatmap
# row_order <- findOrder(orig_data$`ChIP-chip`)
col_order <- findOrder(t(orig_data$`ChIP-chip`))
item_order <- item_names[row_order]
item_col_order <- colnames(orig_data$`ChIP-chip`)[col_order]

chip_chip_plt_data <- orig_data$`ChIP-chip` %>%
  rownames_to_column(var = "Item_i") %>%
  pivot_longer(-Item_i, values_to = "Interaction", names_to = "Item_j") %>%
  mutate(Dataset = "ChIP-chip", Y = match(Item_i, item_order), X = match(Item_j, item_col_order))

if(prod(dim(orig_data$`ChIP-chip`)) != nrow(chip_chip_plt_data)){
  stop("ChIP-chip: Mismatch in size. Data has lost information in transformation.\n")
}
# chip_chip_plt_data$Item_i %>% unique() %>% length()

p_chipchip <- chip_chip_plt_data %>% 
  ggplot(aes(x = X, y = Y, fill = Interaction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  labs(
    title = "ChIP-chip data",
    # subtitle = paste0("Posterior similarity matrices (simulation ", j, ")"),
    x = "Protein",
    y = "",
    fill = "Protein interaction"
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_text(size = 10.5),
    axis.title.x = element_text(size = 10.5),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 10.5),
    legend.text = element_text(size = 10.5),
    panel.border =element_blank()
  )

p_timecourse + p_ppi + p_chipchip +
  plot_layout(guides = "collect")
  
