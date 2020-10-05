#!/usr/bin/Rscript

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

set.seed(1)
setMyTheme()

# === Functions ================================================================

gewekeDF <- function(x,
                     frac_1 = 0.1,
                     frac_2 = 0.5,
                     n_bins = 20,
                     p_value_threshold = 0.05) {
  # The preferred object type for interacting with coda functions
  x <- coda::as.mcmc.list(x)
  
  # The vector of start iterations to calculate the Geweke statistic for
  start_iter_vec <- floor(seq(
    from = stats::start(x),
    to = (stats::start(x) + stats::end(x)) / 2,
    length = n_bins
  ))
  
  # The matrix that will hold the Geweke stat
  geweke_mat <- matrix(nrow = length(start_iter_vec), ncol = coda::nvar(x), dimnames = list(start_iter_vec, coda::varnames(x)))
  
  for (n in 1:length(start_iter_vec)) {
    curr_geweke_diag <- coda::geweke.diag(stats::window(x, start = start_iter_vec[n]),
                                          frac1 = frac_1,
                                          frac2 = frac_2
    )
    
    geweke_mat[n, ] <- curr_geweke_diag[[1]]$z
  }
  
  # The 1.96 threshold for 0.05 significance on a standard normal distribution
  c_limit <- stats::qnorm(1 - p_value_threshold / 2)
  
  # The variables to gather when moving from wide to long data (these are our
  # parameters)
  vars_to_gather <- coda::varnames(x)
  
  # The data.frame we will plot (transform to long data to use the ggplot2
  # framework)
  geweke_df <- data.frame(Start_iteration = start_iter_vec) %>%
    cbind(geweke_mat) %>%
    tidyr::gather_("Parameter", "Geweke_statistic", vars_to_gather)
}

# === Read in data =============================================================


# Original data modelled
data_dir <- "./Data/Yeast/Original_data/"
datasets <- c("Timecourse", "ChIP-chip", "PPI")
data_files <- list.files(data_dir, full.names = T)[c(2, 1, 3)]
orig_data <- data_files %>%
  lapply(read.csv, row.names = 1) %>%
  set_names(datasets)

# Save plots in this directory
save_dir <- "./Images/Yeast/Bayesian/"

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

# All subdirs
all_dirs <- list.dirs(main_dir, recursive = F)

# Dirs for Bayesian inference
bayes_dirs <- all_dirs[grep(all_dirs, pattern = "*BayesianYeastN*")]

# MDI output files
bayes_files <- bayes_dirs %>%
  lapply(list.files, pattern = "*BayesianYeastN676000T1000Seed*", full.names = T) %>%
  str_sort(numeric = T) %>%
  lapply(read.csv)

# Number of chains used
n_chains <- length(bayes_files)

dataset_samples <- list()
for (l in 1:L) {
  dataset_samples[[l]] <- grep(paste0("Dataset", l), colnames(bayes_files[[1]]))
}

# The names of each item that was clustered
item_names <- colnames(bayes_files[[1]])[dataset_samples[[1]]] %>%
  str_remove("Dataset1_")

# Split out the samples into the continuous variables and the dataset specific
# allocations
continuous_data <- .cont_data <- tibble(Samples = list(), Seed = integer())
alloc_data <- .alloc_tib <- tibble(Samples = list(), Seed = integer(), Dataset = character(), PSM = list(), Cl = list())
for (i in 1:n_chains) {
  
  # MDI output
  .curr_data <- bayes_files[[i]]
  
  # The phis and mass parameters
  continuous_data <- continuous_data %>%
    rbind(tibble(Seed = i, Samples = list(.curr_data[, continuous_params_col_index])))
  
  # Dataset specific allocations
  for (l in 1:L) {
    
    # Allocation samples
    alloc <- .curr_data[, dataset_samples[[l]]]
    
    # PSM
    psm <- createSimilarityMat(as.matrix(alloc)) %>% 
      set_rownames(item_names) %>% 
      set_colnames(item_names)
    
    # Predicted clustering (k is set to an arbitrarily large number)
    cl <- mcclust::maxpear(psm, max.k = 275)$cl
    
    # Save to a tibble
    .alloc_tib <- tibble(Samples = list(alloc), Seed = i, Dataset = datasets[[l]], PSM = list(psm), Cl = list(cl))
    alloc_data <- rbind(alloc_data, .alloc_tib)
  }
}

# The number of iterations and thinning factor
alloc_data$R <- 676000
alloc_data$Thin <- 1000

# Save the tibbles
saveRDS(alloc_data, file = "./Data/yeastBayesAllocTibble.rds")
saveRDS(continuous_data, file = "./Data/yeastBayesContParamsTibble.rds")

# === Bayes Visualisation ======================================================

# # Facet labels for continuous variables
# param_labels <- c(paste0("alpha[", 1:3, "]"), paste0("phi[", 1:3, "]"))
# names(param_labels) <- colnames(continuous_data$Samples[[1]])
# 
# # Plot the timecourse data as a time series
# plt_data <- orig_data$Timecourse
# plt_data$Gene <- row.names(plt_data)
# plt_data$Cluster <- alloc_data$Cl[[1]]
# plt_data <- plt_data %>%
#   pivot_longer(cols = -c(Gene, Cluster), values_to = "Expression") %>%
#   mutate(Time = stringr::str_extract(name, "[:digit:]+"))
# 
# plt_data %>%
#   ggplot(aes(x = Time, y = Expression, group = Gene)) +
#   geom_line()
# 
# # Plot time series by cluster
# plt_data %>%
#   ggplot(aes(x = Time, y = Expression, group = Gene)) +
#   geom_line(alpha = 0.3) +
#   facet_wrap(~Cluster)
# 
# pheatmap(orig_data$Timecourse, color = dataColPal(), cluster_cols = F)
# 
# annotatedHeatmap(orig_data$Timecourse, alloc_data$Cl[[1]], cluster_cols = F)
# 
# cl_order <- order(alloc_data$Cl[[1]])
# annotatedHeatmap(orig_data$Timecourse[cl_order, ], alloc_data$Cl[[1]][cl_order],
#                  cluster_cols = F,
#                  cluster_rows = F
# )
# 
# row_order <- findOrder(orig_data$Timecourse)
# compareMatrices(orig_data$Timecourse[row_order, ], alloc_data$PSM[[1]][row_order, row_order], order_cols = F, order_rows = F)
# 
# compareSimilarityMatricesAnnotated(
#   matrices = list(orig_data$Timecourse[row_order, ], alloc_data$PSM[[1]][row_order, row_order]),
#   order_rows = F,
#   order_cols = F,
#   # cluster_IDs =  alloc_data$Cl[[1]],
#   collect_legend = F,
#   col_pal = list(dataColPal(), simColPal()),
#   breaks = list(
#     defineDataBreaks(orig_data$Timecourse[row_order, ], col_pal = dataColPal()),
#     defineBreaks(lb = 0, col_pal = simColPal())
#   )
# )

# PLace continuous parameters in a single data.frame suitable for gglot2
for (i in 1:n_chains) {
  .curr_data <- continuous_data$Samples[[i]] %>%
    pivot_longer(cols = everything(), names_to = "Parameter")
  
  .curr_data$Seed <- continuous_data$Seed[[i]]
  
  if (i == 1) {
    cont_plt_data <- .curr_data
  } else {
    cont_plt_data <- rbind(cont_plt_data, .curr_data)
  }
}
cont_plt_data$Seed <- factor(cont_plt_data$Seed)

p_param_density <- cont_plt_data %>%
  ggplot(aes(x = value, fill = Seed)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  labs(
    title = "Parameter density",
    x = "Value",
    y = "Density"
  ) +
  scale_fill_viridis_d()

# === MCMC diagnostics =========================================================

# Save the samples from each chain as an ``coda::mcmc`` object
mcmc_lst <- list()
for (i in 1:n_chains) {
  mcmc_lst[[i]] <- as.mcmc(continuous_data$Samples[[i]])
}


geweke_lst <- mcmc_lst %>%
  lapply(gewekeDF)


for (i in 1:n_chains) {
  geweke_lst[[i]]$Seed <- i
  if (i == 1) {
    geweke_df <- geweke_lst[[i]]
  } else {
    geweke_df <- rbind(geweke_df, geweke_lst[[i]])
  }
}
geweke_df$Seed <- factor(geweke_df$Seed)

# The 1.96 threshold for 0.05 significance on a standard normal distribution
c_limit <- stats::qnorm(1 - 0.05 / 2)

# For this kind of plot I prefer unifrom axis, thus we find the y-axis limits
y_limit <- max(c(c_limit, abs(geweke_df$Geweke_statistic)))

p_geweke <- geweke_df %>%
  ggplot(aes(x = Start_iteration, y = Geweke_statistic, color = Seed)) +
  geom_line() +
  facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  ggplot2::geom_hline(yintercept = c_limit, linetype = "dashed", color = "grey") +
  ggplot2::geom_hline(yintercept = -c_limit, linetype = "dashed", color = "grey") +
  ggplot2::ylim(-y_limit, y_limit) +
  ggplot2::labs(
    x = "First iteration in segment",
    y = "Geweke's convergence dianostic",
    title = "Within chain convergence"
  ) +
  scale_color_viridis_d()
p_geweke

geweke_df %>%
  filter(Seed %in% c(2, 3, 5, 6, 7, 8, 10)) %>% 
  ggplot(aes(x = Start_iteration, y = Geweke_statistic, color = Seed)) +
  geom_line() +
  facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  ggplot2::geom_hline(yintercept = c_limit, linetype = "dashed", color = "grey") +
  ggplot2::geom_hline(yintercept = -c_limit, linetype = "dashed", color = "grey") +
  # ggplot2::ylim(-y_limit, y_limit) +
  ggplot2::labs(
    x = "First iteration in segment",
    y = "Geweke's convergence dianostic",
    title = "Within chain convergence"
  ) +
  scale_color_viridis_d()

p_geweke_reduced <- geweke_df %>%
  filter(Seed %in% c(2, 3, 5, 7, 8, 10)) %>% 
  ggplot(aes(x = Start_iteration, y = Geweke_statistic, color = Seed)) +
  geom_line() +
  facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  ggplot2::geom_hline(yintercept = c_limit, linetype = "dashed", color = "grey") +
  ggplot2::geom_hline(yintercept = -c_limit, linetype = "dashed", color = "grey") +
  # ggplot2::ylim(-y_limit, y_limit) +
  ggplot2::labs(
    x = "First iteration in segment",
    y = "Geweke's convergence dianostic",
    title = "Within chain convergence"
  ) +
  scale_color_viridis_d()

chains_to_keep <- c(2, 3, 5, 7, 8, 10)

mcmc_lst <- mcmc_lst[chains_to_keep]

# Gelman-Rubin values for each parameter
p_gelman <- gelmanPlot(mcmc_lst) +
  facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) 

# Gelman-Rubin diagnostic using stable variance estimators
stableGR::stable.GR(mcmc_lst)

# Check convergence
n_eff <- n.eff(mcmc_lst)

converged <- n_eff$converged

# Check MCMC convergence
cat(paste0("Chains are converged: ", converged))
p_geweke +
  theme( strip.text = element_text( size = 12, face = "bold" )) +
  theme(axis.text.y=element_text(size = 10.5),
        axis.text.x=element_text(size = 10.5),
        axis.title.y=element_text(size = 10.5),
        axis.title.x=element_text(size = 10.5),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 10.5),
        legend.text = element_text(size = 10.5)
  )

ggsave("./Images/Yeast/gewekePlot.png",
       height = 4, width = 6)

p_geweke_reduced +
  theme( strip.text = element_text( size = 12, face = "bold" )) +
  theme(axis.text.y=element_text(size = 10.5),
        axis.text.x=element_text(size = 10.5),
        axis.title.y=element_text(size = 10.5),
        axis.title.x=element_text(size = 10.5),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 10.5),
        legend.text = element_text(size = 10.5)
  )

ggsave("./Images/Yeast/gewekePlotReduced.png",
       height = 4, width = 6)

p_gelman +
  theme( strip.text = element_text( size = 12, face = "bold" ))+
  ggplot2::scale_color_manual(values = c("grey", "black")) +
  ggplot2::scale_linetype_manual(values = c(2, 1)) +
  theme(axis.text.y=element_text(size = 10.5),
        axis.text.x=element_text(size = 10.5),
        axis.title.y=element_text(size = 10.5),
        axis.title.x=element_text(size = 10.5),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 10.5),
        legend.text = element_text(size = 10.5)
  )

ggsave("./Images/Yeast/gelmanPlot.png",
       height = 4, width = 6)

p_param_density +
  theme(axis.text.y=element_text(size = 10.5),
        axis.text.x=element_text(size = 10.5),
        axis.title.y=element_text(size = 10.5),
        axis.title.x=element_text(size = 10.5),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 10.5),
        legend.text = element_text(size = 10.5)
  )
ggsave("./Images/Yeast/densityPlot.png")

p_cc_timecourse + 
  theme(axis.text.y=element_text(size = 10.5),
        axis.text.x=element_text(size = 10.5),
        axis.title.y=element_text(size = 10.5),
        axis.title.x=element_text(size = 10.5),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 10.5),
        legend.text = element_text(size = 10.5)
  )

ggsave("./Images/Yeast/timecourseClustering.png",
       height = 12, width = 14)

# == NEw =================

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
  
  chain <- alloc_data$Seed[i]
  
  plt_data_transformed <- plt_data %>%
    mutate(Cluster = alloc_data$Cl[[i]]) %>%
    pivot_longer(cols = -c(Gene, Cluster), values_to = "Expression") %>%
    mutate(Time = stringr::str_extract(name, "[:digit:]+"))
  
  # Plot time series by cluster
  p_time_series <- plt_data_transformed %>%
    ggplot(aes(x = as.numeric(Time), y = Expression, group = Gene)) +
    geom_line(alpha = 0.3) +
    facet_wrap(~Cluster) +
    labs(
      title = dataset,
      subtitle = paste0("Clustering predicted by chain ", chain),
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
  
  ggsave(paste0(save_dir, dataset, "TimeSeriesClusterChain", chain, ".png"),
         plot = p_time_series,
         height = 12,
         width = 14
  )
}

# Compare PSMs for the different chains
# Reduce to OK chains
alloc_data_reduced <- alloc_data[alloc_data$Seed %in% chains_to_keep, ]
n_kept <- length(chains_to_keep)

# The labels of the variables
chain_labels <- paste0("Chain ", 1:n_chains) %>% set_names(1:n_chains)

# Iterate over each dataset and save a grid of heatmaps
for (dataset in datasets) {
  
  # The indices and CMs of the tibble corresponding to the current dataset
  curr_inds <- which(alloc_data_reduced$Dataset == dataset)
  curr_psms <- alloc_data_reduced$PSM[curr_inds]
  chains <- alloc_data_reduced$Seed[curr_inds]
  
  # Find the order for the rows and oclumns of the consensus matrices based on
  # the largest & deepest ensemble (here Consensus(1001, 1000))
  row_order <- findOrder(curr_psms[[n_kept]])
  
  # Set up the re-ordered items (this will be used to align the data when in
  # long format)
  item_order <- item_names[row_order]
  
  # Iterate over the different ensembles converting the consensus matrics to
  # long format and prepare them for ggplot
  for (i in 1:n_kept) {
    .psm <- curr_psms[[i]]
    
    .df <- .psm[row_order, row_order] %>%
      
      # Convert to tibble adding the row names as a variable
      as_tibble(rownames = "Gene_i") %>%
      
      # Pivot to long form
      pivot_longer(cols = -Gene_i, names_to = "Gene_j", values_to = "Prop") %>%
      
      # Add variables to indicate the position of the grid for (gene_i, gene_j)
      # on the heatmap based on the ordering defined previously and add the
      # ensemble information
      mutate(
        X = match(Gene_i, item_order),
        Y = match(Gene_j, item_order),
        Chain = chains[i]
      )
    
    # Bind the data frames together
    if (i == 1) {
      psm_plt <- .df
    } else {
      psm_plt <- rbind(psm_plt, .df)
    }
  }
  
  # Add a dataset label
  psm_plt$Dataset <- dataset
  
  if(dataset == datasets[[1]]){
    psm_plt_data_all <- psm_plt
  } else{
    psm_plt_data_all <- rbind(psm_plt_data_all, psm_plt)
  }
  
  # Plot a grid of consensus matrices where column number is increasing the
  # number of samples used and row number corresponds to a value of chain depth
  p_psm <- psm_plt %>%
    ggplot(aes(X, Y, fill = Prop)) +
    geom_tile() +
    facet_wrap(~Chain, labeller = labeller(Chain = chain_labels)
    ) +
    scale_fill_gradient(low = "white", high = "#146EB4") +
    labs(
      title = dataset,
      subtitle = "Posterior similarity matrices",
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
  ggsave(paste0(save_dir, dataset, "PSMcomparisonReduced.png"),
         plot = p_psm,
         height = 14,
         width = 12
  )
}

psm_plt %>%
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
