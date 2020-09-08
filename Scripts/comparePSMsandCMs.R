
library(Matrix)
library(mdiHelpR)
library(magrittr)
library(pheatmap)
library(patchwork)

sim_col <- simColPal()
breaks <- defineBreaks(sim_col, lb = 0, ub = 1)

# Read in the files
# Data being clustered
data <- read.csv("./Data/small_n_large_p_base/dataset_1.csv")
truth <- readRDS("./Data/small_n_large_p_base/cluster_IDs_1.Rds")

# PSMs
psm_dir <- "./Data/small_n_large_p_base/PSMs/"
psm_files <- list.files(psm_dir, full.names = T, pattern = ".txt$") %>%
  stringr::str_sort(numeric = T)

# CMs
cm_dir <- "./Data/small_n_large_p_base/CMs/"
cm_files <- list.files(cm_dir, full.names = T, pattern = ".txt$") %>%
  stringr::str_sort(numeric = T)

# Convert from sparse format
psms <- lapply(psm_files, function(x) {
  as.matrix(readMM(x))
})


cms <- lapply(cm_files, function(x) {
  as.matrix(readMM(x)) * 1
})

# Subset and add the pooled samples
cms_used <- cms[c(1, 3, 5, 6, 8, 10, 11, 12, 14)]
pooled_psm <- Reduce(`+`, psms[1:7]) / 7

psms_used <- psms
psms_used[[8]] <- pooled_psm
#
# compareSimilarityMatricesAnnotated(matrices = psms_used, cluster_IDs = truth, matrix_imposing_order = 8)
# compareSimilarityMatricesAnnotated(matrices = cms_used,
#   matrix_imposing_order = length(cms_used),
#   n_col = 3,
#   cluster_IDs = truth
#
# )


# The order for the CMs (rows adn columns)
cm_order <- findOrder(cms_used[[9]])

ph_list <- cms_used %>%
  lapply(function(x) {
    pheatmap(x[cm_order, cm_order],
      silent = T,
      color = sim_col,
      legend = F,
      cluster_rows = F,
      cluster_cols = F
    )
  }) %>%
  lapply(extract2, "gtable")

ph_patch <- patchwork::wrap_plots(ph_list)

# Order for the PSMs
psm_order <- mdiHelpR::findOrder(pooled_psm)
psm_ph_list <- psms %>%
  lapply(function(x) {
    pheatmap(x[psm_order, psm_order],
      silent = T,
      color = sim_col,
      legend = F,
      cluster_rows = F,
      cluster_cols = F
    )
  }) %>%
  lapply(extract2, "gtable")

psm_patch <- patchwork::wrap_plots(psm_ph_list[1:6])

# Pooled PSM
pooled_gg <- pooled_psm[psm_order, psm_order] %>%
  pheatmap(
    silent = T,
    color = sim_col,
    legend = F,
    cluster_rows = F,
    cluster_cols = F
  ) %>%
  extract2("gtable") %>%
  patchwork::wrap_plots()

# Example CM for CC(100, 50)
cm_gg <- cms_used[[9]][psm_order, psm_order] %>%
  pheatmap(
    silent = T,
    color = sim_col,
    legend = F,
    cluster_rows = F,
    cluster_cols = F
  ) %>%
  extract2("gtable") %>%
  patchwork::wrap_plots()

# Final patchwork of plots
p1 <- (pooled_gg + cm_gg) / psm_patch +
  plot_annotation(title = "Small N, large P", subtitle = "Simulation 1", theme = theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14)
  ))

p1
ggsave("./Images/Simulations/small_n_large_p_base/comp_psms_cm.png", p1,
  height = 6,
  width = 6
)
