#!/usr/bin/Rscript

# Redundant file for the analysis of the consensus clustering for the longer 
# chains

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
orig_output_dir <- "./Data/Yeast/MDIoutput/ConsensusClusteringR10001/"
out_files <- list.files(orig_output_dir) %>%
  str_sort(numeric = T)

n_files <- length(out_files)

# The thinning factor (will be needed for correctly extracting the rth sample
# from the output .csv)
thin <- out_files[1] %>%
  str_match("T([:digit:]+)") %>%
  magrittr::extract2(1, 2) %>%
  as.numeric()

# The ensemble parameters of interest
R <- c(2:10) * 1000 + 1
S <- c(1, c(1, 5, 10) * 100)

# The ensembles available
ensembles <- expand.grid(S, R) %>%
  set_colnames(c("S", "R"))

n_ensembles <- nrow(ensembles)

# Put the samples from MDI into a matrix for each ensemble of interest
cc_samples <- list()
ensemble_ind <- 0

# Iterate over chain length
for (r in R) {
  
  # The row index in the output corresponding to the chain depth of interest
  row_ind <- floor(r / thin + 1)
  sample_matrix <- matrix(0, max(S), n_cont_parameters + L * N)
  
  for (s in 1:n_files) {
    .x <- fread(paste0(orig_output_dir, out_files[s]), nrows = row_ind)
    
    if (s == 1) {
      colnames(sample_matrix) <- colnames(.x)
    }
    
    sample_matrix[s, ] <- round(as.matrix(.x[row_ind, ], nrow = 1))
    
    if (s %in% S) {
      ensemble_ind <- ensemble_ind + 1
      cc_samples[[ensemble_ind]] <- sample_matrix[1:s, , drop = F]
    }
  }
}

# Check that the number of samples in each level is correct
s_used <- lapply(cc_samples, nrow) %>%
  unlist()

if (any(s_used != ensembles$S)) {
  stop("Something is wrong with the ensemble list.")
}

if (ensemble_ind != n_ensembles) {
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
      Samples = list(.curr_data[, continuous_params_col_index]),
      R = ensembles$R[i],
      S = ensembles$S[i]
    ))
  
  # Dataset specific allocations
  for (l in 1:L) {
    
    # Allocation samples
    alloc <- .curr_data[, dataset_samples[[l]], drop = F]
    
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

# Save the tibbles
saveRDS(alloc_data, file = "./Data/Yeast/CCR10001AllocTibble.rds")
saveRDS(continuous_data, file = "./Data/Yeast/CCR10001ContParamsTibble.rds")

# === Visualisation ============================================================

# === Time series ==============================================================

# Time series of the time course data separated out into predicted clusters

# Facet labels for continuous variables
param_labels <- c(paste0("alpha[", 1:3, "]"), paste0("phi[", 1:3, "]"))
names(param_labels) <- colnames(continuous_data$Samples[[1]])

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
    mutate(Cluster = alloc_data$Cl[[i]])  %>% 
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

# === Consensus matrices =======================================================

# Compare consensus matrices for the different ensembles

# The labels of the variables
R_labels <- paste0("R = ", R) %>% set_names(R)
S_labels <- paste0("S = ", S) %>% set_names(S)

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
  
  if(dataset == datasets[[1]]){
    cm_plt_data_all <- cm_plt
  } else{
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
  ggsave(paste0(save_dir, dataset, "CMcomparisonR10001.png"),
         plot = p_cm,
         height = 14,
         width = 12
  )
}

cm_plt %>%
  filter(R %in% c(101, 501), S %in% c(100, 500)) %>% 
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
