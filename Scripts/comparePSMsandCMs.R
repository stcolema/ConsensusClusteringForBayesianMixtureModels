
library(Matrix)
library(mdiHelpR)
library(magrittr)
library(pheatmap)
library(patchwork)

library(stringr)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)

setMyTheme()
sim_col <- simColPal()
breaks <- defineBreaks(sim_col, lb = 0, ub = 1)

# Read in the files
# Data being clustered
data_dir <- "./Data/Simulations/"
scenarios <- list.dirs(data_dir, recursive = F, full.names = F)
scn_dirs <- list.dirs(data_dir, recursive = F, full.names = T)

# scenarios <- scenarios[-10]
# scn_dirs <- scn_dirs[-10]

scenario_str <- scenarios %>%
  str_replace_all("_", " ") %>%
  str_to_sentence()

title_str <- c(
  "Base case",
  expression(paste("Irrelevant features (", P[n], " = 10)")),
  expression(paste("Irrelevant features (", P[n], " = 100)")),
  expression(paste("Irrelevant features (", P[n], " = 20)")),
  expression(paste("Large standard deviation (", sigma^{2}, " = 9)")),
  expression(paste("Large standard deviation (", sigma^{2}, " = 25)")),
  "No structure",
  "2D",
  expression(paste("Small ", N, " large ", P, " (", Delta, mu, " = 1.0)")),
  expression(paste("Small ", N, " large ", P, " (", Delta, mu, " = 0.2)")),
  expression(paste("Varying proportions (", Delta, mu, " = 1.0)")),
  expression(paste("Varying proportions (", Delta, mu, " = 0.4)"))
)

n_scn <- length(scenarios)

sims <- c(1, 10, 100)


# CMs
cm_dir <- paste0(data_dir, scenarios[1], "/CMs/")
cm_files <- list.files(cm_dir, full.names = T, pattern = ".txt$") %>%
  stringr::str_sort(numeric = T)

cms <- lapply(cm_files, function(x) {
  as.matrix(readMM(x)) * 1
})

psm_df_full <- NULL
for (i in 1:n_scn) {
  for (j in sims) {
    data <- read.csv(paste0(data_dir, scenarios[i], "/dataset_", j, ".csv"), row.names = 1)
    truth <- readRDS(paste0(data_dir, scenarios[i], "/cluster_IDs_", j, ".Rds"))

    N <- nrow(data)
    item_names <- row.names(data)

    # PSMs
    psm_dir <- paste0(data_dir, scenarios[i], "/PSMs/")
    psm_files <- list.files(psm_dir, full.names = T, pattern = paste0("Simulation", j, "PSMN")) %>%
      stringr::str_sort(numeric = T)

    chains <- psm_files %>%
      str_match("S([:digit:]+)") %>%
      magrittr::extract(, 2) %>%
      as.numeric()
    n_chains <- length(chains)

    # Convert from sparse format
    psms <- lapply(psm_files, function(x) {
      as.matrix(readMM(x))
    })

    # Create a long data.frame of the PSM entries ready for ggplot2 heatmapping

    for (s in 1:n_chains) {
      .dat <- psms[[s]] %>%
        set_rownames(item_names) %>%
        set_colnames(item_names) %>%
        as.data.frame() %>%
        rownames_to_column(var = "Item_i") %>%
        pivot_longer(-Item_i, values_to = "Prop", names_to = "Item_j") %>%
        mutate(Chain = chains[s])

      if (s == 1) {
        psm_long_df <- .data
      } else {
        psm_long_df <- rbind(psm_long_df, .dat)
      }
    }

    psm_long_df$Scenario <- scenario_str[i]
    psm_long_df$Simulation <- j

    # The method and model ID
    psm_long_df$Method <- "Bayesian"
    psm_long_df$Model <- paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain)

    chain_labels <- paste0("Chain ", chains) %>%
      set_names(chains)

    # factor(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain),
    #        levels = unique(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain))
    # )

    # Set up the X and Y co-ordinates for the entries on the heatmap
    row_order <- findOrder(psms[[1]])
    item_order <- item_names[row_order]
    match(psm_long_df$Item_i, item_order)

    psm_long_df <- psm_long_df %>%
      mutate(X = match(Item_i, item_order), Y = N - match(Item_j, item_order))

    if (is.null(psm_df_full)) {
      psm_df_full <- psm_long_df
    } else {
      psm_df_full <- rbind(psm_df_full, psm_long_df)
    }

    psm_plt <- psm_long_df %>%
      ggplot(aes(x = X, y = Y, fill = Prop)) +
      geom_tile() +
      facet_wrap(~Chain, labeller = labeller(Chain = chain_labels)) +
      scale_fill_gradient(low = "white", high = "#146EB4") +
      labs(
        title = title_str[i],
        subtitle = paste0("Posterior similarity matrices (simulation ", j, ")"),
        x = "Item",
        y = "Item",
        fill = "Coclustering\nproportion"
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
        legend.text = element_text(size = 10.5)
      )

    ggsave(paste0("./SupplementaryMaterial/Images/Simulations/PSMs/", scenarios[i], "Sim", j, ".png"),
      plot = psm_plt,
      height = 6, width = 8
    )
  }
}


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
