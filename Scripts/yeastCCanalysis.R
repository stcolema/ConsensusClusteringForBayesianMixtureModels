#!/usr/bin/Rscript

# Analysis of the stability of the consensus clustering ensembles, visualisation
# of CMs and predicted clusters in the timecourse data.

# For similarity matrix
library(mdiHelpR)

# Pipe
library(magrittr)

# Data manipulation
library(dplyr)

# Tibbles can hold list columns
library(tibble)

# Convergence
library(stableGR)

# Pivot
library(tidyr)

# String manipulation
library(stringr)

# Visualisation
library(ggplot2)
library(pheatmap)

# Predicted clustering
library(mcclust)

# MCMC diagnostic
library(coda)

# Reading in data
library(data.table)

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
save_dir <- "./Images/Yeast/ConsensusClustering/"

# Number of items
N <- orig_data[[1]] %>%
  nrow()

# Number of datasets
L <- orig_data %>%
  length()

# Number of phi parameters from MDI
n_phis <- L * (L - 1) / 2

# Total number of continuous parameters recorded for MDI with L datasets
n_cont_parameters <- L + L * (L - 1) / 2

# Directory for MDI output for Consensus clustering
orig_output_dir <- c(
  "./Data/Yeast/MDIoutput/ConsensusClustering/",
  "./Data/Yeast/MDIoutput/ConsensusClusteringR10001/"
)

out_files <- orig_output_dir %>%
  lapply(list.files) %>%
  lapply(str_sort, numeric = T)

n_analyses <- length(out_files)

n_files <- out_files %>%
  lapply(length)

# The thinning factor (will be needed for correctly extracting the rth sample
# from the output .csv)
thin <- out_files %>%
  lapply(function(x) {
    thin <- x %>%
      magrittr::extract(1) %>%
      str_match("T([:digit:]+)") %>%
      magrittr::extract(1, 2) %>%
      as.numeric()
    thin
  })
thin

# The ensemble parameters of interest
# R <- list(c(0, 1, 5, 10) * 100 + 1, c(2:10) * 1000 + 1)
# S <- list(c(1, c(1, 5, 10) * 100), c(1, c(1, 5, 10) * 100))

R <- list(c(0, 1, 5) * 100 + 1, c(1, 5:10) * 1000 + 1)
S <- list(c(1, c(1, 5, 10) * 100), c(1, c(1, 5, 10) * 100))

# The ensembles available
ensembles <- list()
for (i in 1:n_analyses) {
  ensembles[[i]] <- expand.grid(S[[i]], R[[i]]) %>%
    set_colnames(c("S", "R"))
}

n_ensembles <- ensembles %>%
  lapply(nrow)

# Put the samples from MDI into a matrix for each ensemble of interest
cc_samples <- list()
ensemble_ind <- 0

# Iterate over chain length
for (i in 1:n_analyses) {
  for (r in R[[i]]) {

    # The row index in the output corresponding to the chain depth of interest
    row_ind <- floor(r / thin[[i]] + 1)
    sample_matrix <- data.table(matrix(0, max(S[[i]]), n_cont_parameters + L * N))

    for (s in 1:min(n_files[[i]], max(S[[i]]))) {
      .x <- fread(paste0(orig_output_dir[i], out_files[[i]][s]), nrows = row_ind)

      if (s == 1) {
        colnames(sample_matrix) <- colnames(.x)
      }

      sample_matrix[s, ] <- .x[row_ind, ]

      if (s %in% S[[i]]) {
        ensemble_ind <- ensemble_ind + 1
        cc_samples[[ensemble_ind]] <- sample_matrix[1:s, , drop = F]
      }
    }
  }
}

# Check that the number of samples in each level is correct
s_used <- lapply(cc_samples, nrow) %>%
  unlist()

if (any(s_used != ensembles$S)) {
  stop("Something is wrong with the ensemble list.")
}

if (ensemble_ind != do.call(sum, n_ensembles)) {
  stop("Something is mismatched between the number of ensembles observed and expected.")
}

# Set the names of the entries in the list to be easier to follow
names(cc_samples) <- paste0("R", ensembles$R, "S", ensembles$S)

# Index of continuous variables in csv files
continuous_params_col_index <- 1:(L + n_phis)

# Index of the dataset specific allocations in the sample matrix
dataset_samples <- list()
for (l in 1:L) {
  dataset_samples[[l]] <- grep(paste0("Dataset", l), colnames(sample_matrix))
}

# The names of each item that was clustered
item_names <- colnames(sample_matrix)[dataset_samples[[1]]] %>%
  str_remove("Dataset1_")

ensembles <- do.call(rbind, ensembles)
n_ensembles <- nrow(ensembles)

# Split out the samples into the continuous variables and the dataset specific
# allocations
continuous_data <- .cont_data <- tibble(
  Samples = list(),
  R = integer(),
  S = integer()
)

alloc_data <- .alloc_tib <- tibble(
  Samples = list(),
  R = integer(),
  S = integer(),
  Dataset = character(),
  CM = list(),
  Cl = list()
)

for (i in 1:n_ensembles) {

  # MDI output
  .curr_data <- cc_samples[[i]]

  # The phis and mass parameters
  continuous_data <- continuous_data %>%
    rbind(tibble(
      Samples = list(.curr_data[, ..continuous_params_col_index]),
      R = ensembles$R[i],
      S = ensembles$S[i]
    ))

  # Dataset specific allocations
  for (l in 1:L) {

    # Allocation samples
    alloc_cols <- colnames(.curr_data)[dataset_samples[[l]]]
    alloc <- .curr_data[, ..alloc_cols]

    # Consensus matrix
    cm <- createSimilarityMat(as.matrix(alloc)) %>%
      set_rownames(item_names) %>%
      set_colnames(item_names)

    # Predicted clustering (k is set to an arbitrarily large number)
    cl <- mcclust::maxpear(cm, max.k = 275)$cl

    # Save to a tibble
    .alloc_tib <- tibble(
      Samples = list(alloc),
      R = ensembles$R[i],
      S = ensembles$S[i],
      Dataset = datasets[[l]],
      CM = list(cm),
      Cl = list(cl)
    )

    alloc_data <- rbind(alloc_data, .alloc_tib)
  }
}

# # Check out the behaviour of the merged samples for the CC
# merged_tib <- tibble(
#   Samples = list(),
#   R = character(),
#   S = integer(),
#   Dataset = character(),
#   CM = list(),
#   Cl = list()
# )
# s_final <- 1000
# r_merged <- c(5:10)*1000 + 1
# 
# # The phis and mass parameters
# merged_cont_data <- continuous_data$Samples[which(continuous_data$R %in% r_merged & 
#                                 continuous_data$S == s_final)] %>%
#   unlist() %>% 
#   matrix(byrow = T, ncol = n_cont_parameters) %>% 
#   set_colnames(colnames(continuous_data$Samples[[1]]))
# 
# merged_cont_tib <- tibble(Samples = list(merged_cont_data), S = 1000* 6, R= "Merged")
# 
# for(dataset in datasets){
# ensembles_merged <- which(alloc_data$R %in% r_merged & 
#                             alloc_data$S == s_final &
#                             alloc_data$Dataset == dataset)
# 
# merged_samples <- alloc_data[ensembles_merged, ]$Samples %>% 
#   do.call(rbind, .)
# 
# merged_cm <- createSimilarityMat(as.matrix(merged_samples)) %>%
#   set_rownames(item_names) %>%
#   set_colnames(item_names)
# 
# # Predicted clustering (k is set to an arbitrarily large number)
# cl <- mcclust::maxpear(merged_cm, max.k = 275)$cl
# 
# .alloc_tib <- tibble(
#   Samples = list(merged_samples),
#   R = "Merged",
#   S = s_final * length(ensembles_merged),
#   Dataset = datasets[[l]],
#   CM = list(merged_cm),
#   Cl = list(cl)
# )
# 
# merged_tib <- rbind(merged_tib, .alloc_tib)
# }

# merged_tib$CM[[3]] %>% pheatmap()
# library(mcclust)
# comp_tib <- alloc_data[alloc_data$R == 10001 & alloc_data$S == 1000, ]
# arandi(comp_tib$Cl[[1]], merged_tib$Cl[[1]])
# arandi(comp_tib$Cl[[2]], merged_tib$Cl[[2]])
# arandi(comp_tib$Cl[[3]], merged_tib$Cl[[3]])
# 
# for(l in 1:L){
#   cl_dat <- data.frame(Cl = comp_tib$Cl[[l]],
#                        Dataset = datasets[l],
#                        R = comp_tib$R[l], 
#                        S = comp_tib$S[l], 
#                        Item = 1:N, 
#                        Model = paste0("CC(", comp_tib$R[l], ",", comp_tib$S[l], ")")
#                        )
#   
#   cl_dat_2 <- data.frame(Cl = merged_tib$Cl[[l]],
#                          Dataset = datasets[l],
#                          R = merged_tib$R[l],
#                          S = merged_tib$S[l], 
#                          Item = 1:N,
#                          Model = "Merged")
#   
#   if(l == 1){
#     cl_plt_data <- rbind(cl_dat, cl_dat_2)
#   } else{
#     cl_plt_data <- rbind(cl_plt_data, cl_dat, cl_dat_2)
#   }
#   
# }
# 
# cl_plt_data %>% 
#   ggplot(aes(x = Item, y = Cl, colour = Model)) +
#   geom_point() +
#   geom_jitter() +
#   facet_wrap(~Dataset)
# 
# table(comp_tib$Cl[[1]])
# table(merged_tib$Cl[[1]])
# arandi(comp_tib$Cl[[2]], merged_tib$Cl[[2]])
# arandi(comp_tib$Cl[[3]], merged_tib$Cl[[3]])
# 
# # Compare the consensus matrix for the 
# p1 <- compareSimilarityMatricesAnnotated(comp_tib$CM[[1]], merged_tib$CM[[1]])
# p1_alt <- compareSimilarityMatricesAnnotated(comp_tib$CM[[1]], merged_tib$CM[[1]], matrix_imposing_order = 2)
# p2 <-compareSimilarityMatricesAnnotated(comp_tib$CM[[2]], merged_tib$CM[[2]])
# p2_alt <-compareSimilarityMatricesAnnotated(comp_tib$CM[[2]], merged_tib$CM[[2]], matrix_imposing_order = 2)
# p3 <-compareSimilarityMatricesAnnotated(comp_tib$CM[[3]], merged_tib$CM[[3]])
# p3_alt <-compareSimilarityMatricesAnnotated(comp_tib$CM[[3]], merged_tib$CM[[3]], matrix_imposing_order = 2)
# 
# library(patchwork)
# p1 / p1_alt
# p2 / p2_alt
# p3 / p3_alt
#
# 
# rbind(alloc_data, merged_tib)
# rbind(continuous_data, merged_cont_tib)

# Save the tibbles
saveRDS(alloc_data, file = "./Data/Yeast/CCAllocTibble.rds")
saveRDS(continuous_data, file = "./Data/Yeast/CCContParamsTibble.rds")

# === Visualisation ============================================================

plt_params <- list(
  R = c(1001, 5001, 10001),
  S = c(100, 500, 1000)
)

plt_ensembles <- expand.grid(plt_params$R, plt_params$S) %>%
  set_colnames(c("R", "S"))

# Facet labels for continuous variables
param_labels <- c(paste0("alpha[", 1:3, "]"), paste0("phi[", 1:3, "]"))
names(param_labels) <- colnames(continuous_data$Samples[[1]])

# === Time series ==============================================================

# Time series of the time course data separated out into predicted clusters

# for(dataset in dataset_names){
dataset <- "Timecourse"
curr_inds <- which(alloc_data$Dataset == dataset)

# The data to be plotted
plt_data <- orig_data[[dataset]]
plt_data$Gene <- row.names(plt_data)

for (i in curr_inds) {
  r <- alloc_data$R[i]
  s <- alloc_data$S[i]

  plt_data_transformed <- plt_data %>%
    mutate(Cluster = alloc_data$Cl[[i]]) %>%
    add_count(Cluster) %>%
    pivot_longer(cols = -c(Gene, Cluster, n), values_to = "Expression") %>%
    mutate(Time = stringr::str_extract(name, "[:digit:]+"))

  # Plot time series by cluster
  p_time_series <- plt_data_transformed %>%
    ggplot(aes(x = as.numeric(Time), y = Expression, group = Gene)) +
    geom_line(alpha = 0.3) +
    facet_wrap(~Cluster) +
    labs(
      title = dataset,
      subtitle = paste0("Clustering predicted by CC(", r, ",", s, ")"),
      x = "Time"
    ) +
    theme(
      axis.text.y = element_text(size = 10.5),
      axis.text.x = element_text(size = 10.5),
      axis.title.y = element_text(size = 10.5),
      axis.title.x = element_text(size = 10.5),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14),
      strip.text.x = element_text(size = 10.5),
      legend.text = element_text(size = 10.5)
    )

  ggsave(paste0(save_dir, dataset, "TimeSeriesClusterR", r, "S", s, ".png"),
    plot = p_time_series,
    height = 12,
    width = 14
  )

  # Plot time series by cluster
  p_time_series_no_singletons <- plt_data_transformed %>%
    filter(n > 1) %>%
    ggplot(aes(x = as.numeric(Time), y = Expression, group = Gene)) +
    geom_line(alpha = 0.3) +
    facet_wrap(~Cluster) +
    labs(
      title = dataset,
      subtitle = paste0("Clustering predicted by CC(", r, ",", s, ")"),
      x = "Time"
    ) +
    theme(
      axis.text.y = element_text(size = 10.5),
      axis.text.x = element_text(size = 10.5),
      axis.title.y = element_text(size = 10.5),
      axis.title.x = element_text(size = 10.5),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14),
      strip.text.x = element_text(size = 10.5),
      legend.text = element_text(size = 10.5)
    )

  ggsave(paste0(save_dir, dataset, "TimeSeriesClusterR", r, "S", s, "NoSingletons.png"),
    plot = p_time_series_no_singletons,
    height = 12,
    width = 14
  )
}

# The CC(10001, 1000) analysis. This plot is used in the supplementary materials
r <- alloc_data$R[118]
s <- alloc_data$S[118]

plt_data_transformed <- plt_data %>%
  mutate(Cluster = alloc_data$Cl[[118]]) %>%
  add_count(Cluster) %>%
  pivot_longer(cols = -c(Gene, Cluster, n), values_to = "Expression") %>%
  mutate(Time = stringr::str_extract(name, "[:digit:]+"))

# Plot time series by cluster
p_time_series <- plt_data_transformed %>%
  ggplot(aes(x = as.numeric(Time), y = Expression, group = Gene)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~Cluster, nrow = 4) +
  labs(
    title = dataset,
    subtitle = paste0("Clustering predicted by CC(", r, ",", s, ")"),
    x = "Time"
  ) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 21, face = "bold"),
    plot.subtitle = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

ggsave(paste0("./SupplementaryMaterial/Images/Yeast/TimeSeriesClusterR", r, "S", s, ".png"),
       plot = p_time_series,
       height = 10,
       width = 14
)

# === Consensus matrices =======================================================

# Compare consensus matrices for the different ensembles

# The labels of the variables
R_labels <- paste0("R = ", unlist(R)) %>% set_names(unlist(R))
S_labels <- paste0("S = ", unlist(S)) %>% set_names(unlist(S))

# Iterate over each dataset and save a grid of heatmaps
for (dataset in datasets) {

  # The indices and CMs of the tibble corresponding to the current dataset
  curr_inds <- which(alloc_data$Dataset == dataset)
  curr_cms <- alloc_data$CM[curr_inds]

  # Find the order for the rows and oclumns of the consensus matrices based on
  # the largest & deepest ensemble (here Consensus(1001, 1000))
  row_order <- findOrder(curr_cms[[n_ensembles]])

  # Set up the re-ordered items (this will be used to align the data when in
  # long format)
  item_order <- item_names[row_order]

  # Iterate over the different ensembles converting the consensus matrics to
  # long format and prepare them for ggplot
  for (i in 1:n_ensembles) {
    .cm <- curr_cms[[i]]

    .df <- .cm[row_order, row_order] %>%

      # Convert to tibble adding the row names as a variable
      as_tibble(rownames = "Gene_i") %>%

      # Pivot to long form
      pivot_longer(cols = -Gene_i, names_to = "Gene_j", values_to = "Prop") %>%

      # Add variables to indicate the position of the grid for (gene_i, gene_j)
      # on the heatmap based on the ordering defined previously and add the
      # ensemble information
      mutate(
        X = match(Gene_i, item_order),
        Y = N - match(Gene_j, item_order),
        R = ensembles$R[i],
        S = ensembles$S[i]
      )

    # Bind the data frames together
    if (i == 1) {
      cm_plt <- .df
    } else {
      cm_plt <- rbind(cm_plt, .df)
    }
  }

  # Add a dataset label
  cm_plt$Dataset <- dataset

  if (dataset == datasets[[1]]) {
    cm_plt_data_all <- cm_plt
  } else {
    cm_plt_data_all <- rbind(cm_plt_data_all, cm_plt)
  }

  # Plot a grid of consensus matrices where column number is increasing the
  # number of samples used and row number corresponds to a value of chain depth
  p_cm <- cm_plt %>%
    ggplot(aes(X, Y, fill = Prop)) +
    geom_tile() +
    facet_grid(
      rows = vars(R),
      cols = vars(S),
      labeller = labeller(R = R_labels, S = S_labels)
    ) +
    scale_fill_gradient(low = "white", high = "#146EB4") +
    labs(
      title = dataset,
      subtitle = "Consensus matrices",
      x = "Gene",
      y = "Gene",
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

  # Save!
  ggsave(paste0(save_dir, dataset, "CMcomparison.png"),
    plot = p_cm,
    height = 14,
    width = 12
  )
}

supp_plts <- list()
for (l in 1:L) {
  dataset <- datasets[l]
  supp_plts[[l]] <- p <- cm_plt_data_all %>%
    filter(
      R %in% c(1001, 5001, 10001),
      S %in% c(100, 500, 1000),
      Dataset == dataset
    ) %>%
    ggplot(aes(X, Y, fill = Prop)) +
    geom_tile() +
    # facet_wrap(~Dataset)+
    facet_grid(
      rows = vars(R),
      cols = vars(S),
      labeller = labeller(R = R_labels, S = S_labels)
    ) +
    scale_fill_gradient(low = "white", high = "#146EB4") +
    labs(
      title = dataset,
      subtitle = "Consensus matrices",
      x = "Gene",
      y = "Gene",
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

  ggsave(paste0("./SupplementaryMaterial/Images/Yeast/", dataset, "CMcomparison.png"),
    plot = p,
    height = 8,
    width = 7
  )
}

# === Density plots ============================================================

n_plts <- nrow(plt_ensembles)

# PLace continuous parameters in a single data.frame suitable for gglot2
for (i in 1:n_plts) {
  .ind <- which(continuous_data$R == plt_ensembles$R[i] & continuous_data$S == plt_ensembles$S[i])

  if (length(.ind) > 1) {
    stop("Too many ensembles match criteria, check tibble for oddities.")
  }

  .curr_data <- continuous_data$Samples[[.ind]] %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "Parameter")

  .curr_data$R <- continuous_data$R[[.ind]]
  .curr_data$S <- continuous_data$S[[.ind]]

  if (i == 1) {
    cont_plt_data <- .curr_data
  } else {
    cont_plt_data <- rbind(cont_plt_data, .curr_data)
  }
}
cont_plt_data$Model <- paste0("CC(", cont_plt_data$R, ", ", cont_plt_data$S, ")")

# Facet labels for continuous variables
param_labels <- c(paste0("alpha[", 1:3, "]"), paste0("phi[", c(12, 13, 23), "]"))
names(param_labels) <- colnames(continuous_data$Samples[[n_ensembles]])


p_param_density <- cont_plt_data %>%
  # filter(Model == "CC(10001, 1000)") %>%
  ggplot(aes(x = value, fill = Model)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  labs(
    title = "Parameter density",
    x = "Value",
    y = "Density"
  ) +
  scale_fill_viridis_d()

p_param_density +
  theme(
    axis.text.y = element_text(size = 10.5),
    axis.text.x = element_text(size = 10.5),
    axis.title.y = element_text(size = 10.5),
    axis.title.x = element_text(size = 10.5),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 10.5),
    legend.text = element_text(size = 10.5)
  )
