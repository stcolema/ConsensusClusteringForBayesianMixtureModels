
# Investigate and compare the CMs and PSMs for the simulations

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
  expression(paste("Large standard deviation (", sigma^{
    2
  }, " = 9)")),
  expression(paste("Large standard deviation (", sigma^{
    2
  }, " = 25)")),
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

# === PSMs =====================================================================

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
      as.matrix(readMM(x)) * 1
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
        psm_long_df <- .dat
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
    # match(psm_long_df$Item_i, item_order)

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
      height = 5, width = 6
    )
  }
}

# === CMs ======================================================================

s_plotted <- c(1, 10, 30, 50, 100)
r_plotted <- c(1, 10, 100, 1000, 10000)

cm_df_full <- NULL
for (i in 1:n_scn) {
  for (j in sims) {
    data <- read.csv(paste0(data_dir, scenarios[i], "/dataset_", j, ".csv"), row.names = 1)
    truth <- readRDS(paste0(data_dir, scenarios[i], "/cluster_IDs_", j, ".Rds"))

    N <- nrow(data)
    item_names <- row.names(data)

    # CMs
    cm_dir <- paste0(data_dir, scenarios[i], "/CMs/")
    cm_files <- list.files(cm_dir, full.names = T, pattern = paste0("Simulation", j, "ConsensusMatrix")) %>%
      stringr::str_sort(numeric = T)

    # Find the ensembles present
    S_used <- cm_files %>%
      str_match("S([:digit:]+)") %>%
      magrittr::extract(, 2) %>%
      as.numeric()

    R_used <- cm_files %>%
      str_match("N([:digit:]+)") %>%
      magrittr::extract(, 2) %>%
      as.numeric()

    # The ensemble details (expect this to be the same in each scenario)
    ensembles <- data.frame("S" = S_used, "R" = R_used)
    n_ensembles <- nrow(ensembles)

    # Convert from sparse format
    cms <- lapply(cm_files, function(x) {
      as.matrix(readMM(x)) * 1
    })

    if (length(cms) != n_ensembles) {
      stop("Number of ensembles and number of CMs not matching.\n")
    }

    # Create a long data.frame of the PSM entries ready for ggplot2 heatmapping

    for (k in 1:n_ensembles) {
      if (ensembles$S[k] %in% s_plotted & ensembles$R[k] %in% r_plotted) {
        .dat <- cms[[k]] %>%
          set_rownames(item_names) %>%
          set_colnames(item_names) %>%
          as.data.frame() %>%
          rownames_to_column(var = "Item_i") %>%
          pivot_longer(-Item_i, values_to = "Prop", names_to = "Item_j") %>%
          mutate(S = ensembles$S[k], R = ensembles$R[k])

        if (k == 1) {
          cm_long_df <- .dat
        } else {
          cm_long_df <- rbind(cm_long_df, .dat)
        }
      }
    }

    cm_long_df$Scenario <- scenario_str[i]
    cm_long_df$Simulation <- j

    # The method and model ID
    cm_long_df$Method <- "Consensus clustering"
    cm_long_df$Model <- paste0("CC(", cm_long_df$R, ", ", cm_long_df$S, ")")

    R_labels <- paste0("D = ", unique(R_used)) %>%
    # R_labels <- c("1st iteration used", paste0(unique(R_used)[-1], "th iteration")) %>%
      set_names(unique(R_used))

    # S_labels <- c("1 chain used", paste0(unique(S_used)[-1], " chains used")) %>%
    #   set_names(unique(S_used))
    
    S_labels <- paste0("W = ", unique(S_used)) %>%
      set_names(unique(S_used))

    # factor(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain),
    #        levels = unique(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain))
    # )

    # Set up the X and Y co-ordinates for the entries on the heatmap
    row_order <- findOrder(cms[[n_ensembles]])
    item_order <- item_names[row_order]
    # match(cm_long_df$Item_i, item_order)

    cm_long_df <- cm_long_df %>%
      mutate(X = match(Item_i, item_order), Y = N - match(Item_j, item_order))

    if (is.null(cm_long_df)) {
      cm_df_full <- cm_long_df
    } else {
      cm_df_full <- rbind(cm_df_full, cm_long_df)
    }

    cm_plt <- cm_long_df %>%
      ggplot(aes(x = X, y = Y, fill = Prop)) +
      geom_tile() +
      facet_grid(R ~ S, labeller = labeller(R = R_labels, S = S_labels)) +
      scale_fill_gradient(low = "white", high = "#146EB4") +
      labs(
        title = title_str[i],
        subtitle = paste0("Consensus matrices (simulation ", j, ")"),
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

    ggsave(paste0("./SupplementaryMaterial/Images/Simulations/CMs/", scenarios[i], "Sim", j, ".png"),
      plot = cm_plt,
      height = 6, width = 8
    )
  }
}

# === Old comparison ===========================================================

# write.csv(psm_df_full, "./Data/Simulations/PSMLongData.csv")

curr_scn <- "small_n_large_p_base"
scn_ind <- which(scenarios == curr_scn)
sim_num <- 1

data <- read.csv(paste0(data_dir, scenarios[scn_ind], "/dataset_", sim_num, ".csv"), row.names = 1)
truth <- readRDS(paste0(data_dir, scenarios[scn_ind], "/cluster_IDs_", sim_num, ".Rds"))

N <- nrow(data)
item_names <- row.names(data)

# PSMs
psm_dir <- paste0(data_dir, scenarios[scn_ind], "/PSMs/")
psm_files <- list.files(psm_dir, full.names = T, pattern = paste0("Simulation", sim_num, "PSMN")) %>%
  stringr::str_sort(numeric = T)

chains <- psm_files %>%
  str_match("S([:digit:]+)") %>%
  magrittr::extract(, 2) %>%
  as.numeric()
n_chains <- length(chains)

# Convert from sparse format
psms <- lapply(psm_files, function(x) {
  as.matrix(readMM(x)) * 1
})

# CMs
cm_dir <- paste0(data_dir, scenarios[scn_ind], "/CMs/")
cm_files <- list.files(cm_dir, full.names = T, pattern = paste0("Simulation", sim_num, "ConsensusMatrix")) %>%
  stringr::str_sort(numeric = T)

# Find the ensembles present
S_used <- cm_files %>%
  str_match("S([:digit:]+)") %>%
  magrittr::extract(, 2) %>%
  as.numeric()

R_used <- cm_files %>%
  str_match("N([:digit:]+)") %>%
  magrittr::extract(, 2) %>%
  as.numeric()

# The ensemble details (expect this to be the same in each scenario)
ensembles <- data.frame("S" = S_used, "R" = R_used)
n_ensembles <- nrow(ensembles)

# Convert from sparse format
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
      cluster_cols = F,
      border_color = NA
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
    cluster_cols = F,
    border_color = NA
  ) %>%
  extract2("gtable")
# %>%
  # patchwork::wrap_plots()

# Example CM for CC(100, 50)
cm_gg <- cms_used[[9]][psm_order, psm_order] %>%
  pheatmap(
    silent = T,
    color = sim_col,
    legend = F,
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  ) %>%
  extract2("gtable") 
# %>%
  # patchwork::wrap_plots()
#

plt_lst <- vector("list", 8)
plt_lst[3:8] <- psm_ph_list[1:6]
plt_lst[[1]] <- cm_gg
plt_lst[[2]] <- pooled_gg

layout <- "
AAABBB
AAABBB
CCDDEE
FFGGHH
"

patchwork::wrap_plots(cm_gg, 
                      pooled_gg,
                      psm_ph_list[[1]],
                      psm_ph_list[[2]],
                      psm_ph_list[[3]],
                      psm_ph_list[[4]],
                      psm_ph_list[[5]],
                      psm_ph_list[[6]],
                      design = layout)  + plot_annotation(tag_levels = 'A')

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


# 9 corresponds to CC(100, 50)
.dat <- cms_used[[9]] %>%
      set_rownames(item_names) %>%
      set_colnames(item_names) %>%
      as.data.frame() %>%
      rownames_to_column(var = "Item_i") %>%
      pivot_longer(-Item_i, values_to = "Prop", names_to = "Item_j") %>%
      mutate(S = ensembles$S[k], R = ensembles$R[k])
    
    if (k == 1) {
      cm_long_df <- .dat
    } else {
      cm_long_df <- rbind(cm_long_df, .dat)
    }
  }
}

cm_long_df$Scenario <- scenario_str[i]
cm_long_df$Simulation <- j

# The method and model ID
cm_long_df$Method <- "Consensus clustering"
cm_long_df$Model <- paste0("CC(", cm_long_df$R, ",", cm_long_df$S, ")")

R_labels <- paste0("R = ", unique(R_used)) %>%
  set_names(unique(R_used))

S_labels <- paste0("S = ", unique(S_used)) %>%
  set_names(unique(S_used))

# factor(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain),
#        levels = unique(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain))
# )

# Set up the X and Y co-ordinates for the entries on the heatmap
row_order <- findOrder(cms[[n_ensembles]])
item_order <- item_names[row_order]
# match(cm_long_df$Item_i, item_order)

cm_long_df <- cm_long_df %>%
  mutate(X = match(Item_i, item_order), Y = N - match(Item_j, item_order))

if (is.null(cm_long_df)) {
  cm_df_full <- cm_long_df
} else {
  cm_df_full <- rbind(cm_df_full, cm_long_df)
}

# === New comparison ===========================================================

curr_scn <- "small_n_large_p_base"
scn_ind <- which(scenarios == curr_scn)
sim_num <- 1
cm_comp_df <- NULL


data <- read.csv(paste0(data_dir, scenarios[scn_ind], "/dataset_", sim_num, ".csv"), row.names = 1)
truth <- readRDS(paste0(data_dir, scenarios[scn_ind], "/cluster_IDs_", sim_num, ".Rds"))

N <- nrow(data)
item_names <- row.names(data)

# CMs
cm_dir <- paste0(data_dir, scenarios[scn_ind], "/CMs/")
cm_files <- list.files(cm_dir, full.names = T, pattern = paste0("Simulation", sim_num, "ConsensusMatrix")) %>%
  stringr::str_sort(numeric = T)

# Find the ensembles present
S_used <- cm_files %>%
  str_match("S([:digit:]+)") %>%
  magrittr::extract(, 2) %>%
  as.numeric()

R_used <- cm_files %>%
  str_match("N([:digit:]+)") %>%
  magrittr::extract(, 2) %>%
  as.numeric()

# The ensemble details (expect this to be the same in each scenario)
ensembles <- data.frame("S" = S_used, "R" = R_used)
n_ensembles <- nrow(ensembles)

# Convert from sparse format
cms <- lapply(cm_files, function(x) {
  as.matrix(readMM(x)) * 1
})

if (length(cms) != n_ensembles) {
  stop("Number of ensembles and number of CMs not matching.\n")
}

# Create a long data.frame of the PSM entries ready for ggplot2 heatmapping

for (k in 1:n_ensembles) {
  if (ensembles$S[k] %in% s_plotted & ensembles$R[k] %in% r_plotted) {
    .dat <- cms[[k]] %>%
      set_rownames(item_names) %>%
      set_colnames(item_names) %>%
      as.data.frame() %>%
      rownames_to_column(var = "Item_i") %>%
      pivot_longer(-Item_i, values_to = "Prop", names_to = "Item_j") %>%
      mutate(S = ensembles$S[k], R = ensembles$R[k])
    
    if (k == 1) {
      cm_long_df <- .dat
    } else {
      cm_long_df <- rbind(cm_long_df, .dat)
    }
  }
}

cm_long_df$Scenario <- scenario_str[scn_ind]
cm_long_df$Simulation <- sim_num
cm_long_df$Chain <- NA

# The method and model ID
cm_long_df$Method <- "Consensus clustering"
cm_long_df$Model <- paste0("CC(", cm_long_df$R, ", ", cm_long_df$S, ")")

R_labels <- paste0("R = ", unique(R_used)) %>%
  set_names(unique(R_used))

S_labels <- paste0("S = ", unique(S_used)) %>%
  set_names(unique(S_used))

# factor(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain),
#        levels = unique(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain))
# )

# Set up the X and Y co-ordinates for the entries on the heatmap
# row_order <- findOrder(cms[[n_ensembles]])
# item_order <- item_names[row_order]
# # match(cm_long_df$Item_i, item_order)
# 
# cm_long_df <- cm_long_df %>%
#   mutate(X = match(Item_i, item_order), Y = N - match(Item_j, item_order))

# PSMs
psm_dir <- paste0(data_dir, scenarios[scn_ind], "/PSMs/")
psm_files <- list.files(psm_dir, full.names = T, pattern = paste0("Simulation", sim_num, "PSMN")) %>%
  stringr::str_sort(numeric = T)

chains <- psm_files %>%
  str_match("S([:digit:]+)") %>%
  magrittr::extract(, 2) %>%
  as.numeric()

n_chains <- length(chains)

# Convert from sparse format
psms <- lapply(psm_files, function(x) {
  as.matrix(readMM(x)) * 1
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
    psm_long_df <- .dat
  } else {
    psm_long_df <- rbind(psm_long_df, .dat)
  }
}

psm_long_df$Scenario <- scenario_str[scn_ind]
psm_long_df$Simulation <- sim_num

# The method and model ID
psm_long_df$Method <- "Bayesian"
psm_long_df$Model <- paste0("Chain ", psm_long_df$Chain)

psm_long_df$S <- 1
psm_long_df$R <- 1e6


chain_labels <- paste0("Chain ", chains) %>%
  set_names(chains)


psm_pooled <- Reduce(`+`, psms) / n_chains

.dat <- psm_pooled %>%
  set_rownames(item_names) %>%
  set_colnames(item_names) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Item_i") %>%
  pivot_longer(-Item_i, values_to = "Prop", names_to = "Item_j") %>%
  mutate(Scenario = scenario_str[scn_ind],
         Simulation = sim_num,
         Chain = NA,
         Method = "Bayesian",
         Model = "Bayesian (Pooled)",
         S = n_chains,
         R = 1e6)

psm_long_df <- rbind(psm_long_df, .dat)

plt_df <- rbind(psm_long_df, cm_long_df)

# factor(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain),
#        levels = unique(paste0(psm_long_df$Method, ": chain ", psm_long_df$Chain))
# )

# Set up the X and Y co-ordinates for the entries on the heatmap
row_order <- findOrder(psm_pooled)
item_order <- item_names[row_order]
# match(psm_long_df$Item_i, item_order)

plt_df <- plt_df %>%
  mutate(X = match(Item_i, item_order), Y = N - match(Item_j, item_order))

p1 <- plt_df %>% 
  filter(Model %in% c("Bayesian (Pooled)", "CC(100, 50)")) %>% 
  ggplot(aes(x = X, y = Y, fill = Prop)) +
  geom_tile() +
  facet_wrap(~Model, nrow = 1) +
  # facet_grid(R ~ S, labeller = labeller(R = R_labels, S = S_labels)) +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  labs(
    # title = title_str[i],
    # subtitle = paste0("Consensus matrices (simulation ", j, ")"),
    x = "Item",
    y = "Item",
    fill = "Coclustering\nproportion"
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_text(size = 10.5),
    axis.title.x = element_blank(), #element_text(size = 10.5),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 10.5),
    legend.text = element_text(size = 10.5)
  )

p2 <- plt_df %>% 
  filter(Model %in% c(paste0("Chain ", 1:8))) %>% 
  ggplot(aes(x = X, y = Y, fill = Prop)) +
  geom_tile() +
  facet_wrap(~Model) +
  # facet_grid(R ~ S, labeller = labeller(R = R_labels, S = S_labels)) +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  labs(
    # title = title_str[i],
    # subtitle = paste0("Consensus matrices (simulation ", j, ")"),
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

p1 / p2 +
  plot_layout(guides = 'collect') +
  plot_annotation(title = "Small N, large P",
                  subtitle = "Comparison of similarity matrices")

ggsave(paste0("./Images/Simulations/small_n_large_p_base/comp_psms_cm.png"),
       # plot = cm_plt,
       height = 6, width = 6
)  


