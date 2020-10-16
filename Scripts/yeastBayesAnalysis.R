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

# Convergence
library(stableGR)

# Pivot
library(tidyr)

# String manipulation
library(stringr)

# Visualisation
library(ggplot2)
library(pheatmap)
library(patchwork)
library(scales)

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
save_dir <- "./SupplementaryMaterial/Images/Yeast/"

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

# Directory containing MCMC samples
main_dir <- "./Data/Yeast/MDIoutput/Bayesian/"

# MDI output files
bayes_files <- main_dir %>%
  list.files(pattern = "*BayesianYeastN676000T1000Seed*", full.names = T) %>%
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
  .curr_data <- bayes_files[[i]][-(1:10), ]

  # The phis and mass parameters
  continuous_data <- continuous_data %>%
    rbind(tibble(Seed = i, Samples = list(.curr_data[, 1:n_cont_parameters])))

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

# === Bayes Visualisation ======================================================

# Facet labels for continuous variables
param_labels <- c(paste0("alpha[", 1:3, "]"), paste0("phi[", c(12, 13, 23), "]"))
names(param_labels) <- colnames(continuous_data$Samples[[1]])

# The labels of the variables
chain_labels <- paste0("Chain ", 1:n_chains) %>% set_names(1:n_chains)

# === MCMC diagnostics =========================================================

# Save the samples from each chain as an ``coda::mcmc`` object
mcmc_lst <- list()
for (i in 1:n_chains) {
  mcmc_lst[[i]] <- as.mcmc(continuous_data$Samples[[i]])
}


mdiHelpR::findESS(mcmc_lst)

# setMyTheme()
# library(fitR) # devtools::install_github("sbfnk/fitR")
# ess_plts <- mdiHelpR::makeESSPlots(mcmc_lst)

# ess_plts[[1]]
# ess_plts[[2]]
# ess_plts[[3]]
# ess_plts[[4]]
# ess_plts[[5]]
# ess_plts[[6]]
# ess_plts[[7]]
# ess_plts[[8]]
# ess_plts[[9]]
# ess_plts[[10]]
# 
# burn_in <- 200
# 
# mcmc_lst <- mcmc_lst %>% 
#   lapply(magrittr::extract, -(1:burn_in), ) %>% 
#   lapply(as.mcmc)

geweke_lst <- mcmc_lst %>%
  lapply(gewekeDF)


for (i in 1:n_chains) {
  geweke_lst[[i]]$Seed <- i
  if (i == 1) {
    geweke_df <- geweke_lst[[i]]
  } else {
    geweke_df <- rbind(geweke_df, geweke_lst[[i]]
    )
  }
}
geweke_df$Seed <- factor(geweke_df$Seed)
# geweke_df$Parameter <- factor(geweke_df$Parameter)

# The 1.96 threshold for 0.05 significance on a standard normal distribution
c_limit <- stats::qnorm(1 - 0.05 / 2)

# For this kind of plot I prefer unifrom axis, thus we find the y-axis limits
y_limit <- max(c(c_limit, abs(geweke_df$Geweke_statistic)))

p_geweke <- geweke_df %>%
  mutate(Parameter_expr = recode_factor(Parameter, 
                                   MassParameter_1 = "alpha[1]",
                                   MassParameter_2 = "alpha[2]",
                                   MassParameter_3 = "alpha[3]",
                                   Phi_12 = "phi[12]",
                                   Phi_13 = "phi[13]",
                                   Phi_23 = "phi[23]")) %>% 
  ggplot(aes(x = Start_iteration, y = Geweke_statistic, color = Parameter_expr)) +
  geom_line() +
  # facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  facet_wrap(~Seed, labeller = as_labeller(chain_labels)) +
  ggplot2::geom_hline(yintercept = c_limit, linetype = "dashed", color = "grey") +
  ggplot2::geom_hline(yintercept = -c_limit, linetype = "dashed", color = "grey") +
  ggplot2::ylim(-y_limit, y_limit) +
  ggplot2::labs(
    x = "First iteration in segment",
    y = "Geweke's convergence dianostic",
    title = "Within chain convergence",
    color = "Parameter"
  ) +
  # scale_colour_discrete(labels = parse_format())
  scale_color_viridis_d(labels = parse_format())

p_geweke

p_gweke_reduced_partially <- geweke_df %>%
  filter(Seed %in% c(1, 2, 3, 5, 6, 7, 8, 10)) %>%
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

p_gweke_reduced_partially

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

p_geweke_reduced

p_phi_geweke <- geweke_df %>%
  filter(Seed %in% c(1, 2, 3, 4, 5, 6, 7, 8, 10), Parameter %in% c("Phi_12", "Phi_13", "Phi_23")) %>%
  ggplot(aes(x = Start_iteration, y = Geweke_statistic)) +
  geom_line() +
  facet_grid(rows = vars(Seed),
  cols = vars(Parameter),
  labeller = labeller(Parameter = as_labeller(param_labels, label_parsed), Seed = chain_labels)
  ) + 
  # facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  ggplot2::geom_hline(yintercept = c_limit, linetype = "dashed", color = "grey") +
  ggplot2::geom_hline(yintercept = -c_limit, linetype = "dashed", color = "grey") +
  # ggplot2::ylim(-y_limit, y_limit) +
  ggplot2::labs(
    x = "First iteration in segment",
    y = "Geweke's convergence dianostic",
    title = "Within chain convergence"
  ) +
  scale_color_viridis_d()

p_phi_geweke

geweke_df %>%
  filter(Seed %in% c(3, 5, 7, 8, 10), Parameter %in% c("Phi_12", "Phi_13", "Phi_23")) %>%
  ggplot(aes(x = Start_iteration, y = Geweke_statistic)) +
  geom_line() +
  facet_grid(rows = vars(Seed),
             cols = vars(Parameter),
             labeller = labeller(Parameter = as_labeller(param_labels, label_parsed), Seed = chain_labels)
  ) + 
  # facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  ggplot2::geom_hline(yintercept = c_limit, linetype = "dashed", color = "grey") +
  ggplot2::geom_hline(yintercept = -c_limit, linetype = "dashed", color = "grey") +
  # ggplot2::ylim(-y_limit, y_limit) +
  ggplot2::labs(
    x = "First iteration in segment",
    y = "Geweke's convergence dianostic",
    title = "Within chain convergence"
  ) +
  scale_color_viridis_d()

geweke_df %>%
  filter(Seed %in% c(3, 7, 8, 10), Parameter %in% c("Phi_12", "Phi_13", "Phi_23")) %>%
  ggplot(aes(x = Start_iteration, y = Geweke_statistic)) +
  geom_line() +
  facet_grid(rows = vars(Seed),
             cols = vars(Parameter),
             labeller = labeller(Parameter = as_labeller(param_labels, label_parsed), Seed = chain_labels)
  ) + 
  # facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  ggplot2::geom_hline(yintercept = c_limit, linetype = "dashed", color = "grey") +
  ggplot2::geom_hline(yintercept = -c_limit, linetype = "dashed", color = "grey") +
  # ggplot2::ylim(-y_limit, y_limit) +
  ggplot2::labs(
    x = "First iteration in segment",
    y = "Geweke's convergence dianostic",
    title = "Within chain convergence"
  ) +
  scale_color_viridis_d()

# chains_to_keep <- c(2, 3, 5, 7, 8, 10)
chains_to_keep <- c(3, 5, 7, 8, 10)

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
  theme(strip.text = element_text(size = 12, face = "bold")) +
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

ggsave("./SupplementaryMaterial/Images/Yeast/Convergence/gewekePlot.png",
  height = 6, width = 6
)

p_geweke_reduced +
  theme(strip.text = element_text(size = 12, face = "bold")) +
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

ggsave("./SupplementaryMaterial/Images/Yeast/Convergence/gewekePlotReduced.png",
  height = 4, width = 6
)

p_gelman +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  ggplot2::scale_color_manual(values = c("grey", "black")) +
  ggplot2::scale_linetype_manual(values = c(2, 1)) +
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

ggsave("./SupplementaryMaterial/Images/Yeast/Convergence/gelmanPlot.png",
  height = 4, width = 6
)

p_gelman$data %>% 
  filter(Parameter %in% c("Phi_12", "Phi_13", "Phi_23")) %>% 
  group_by(Parameter) %>% 
  summarise(quantiles = quantile(Shrinkage_factor))

p_gelman$data %>% 
  ggplot(aes(x = Parameter, y = Shrinkage_factor)) +
  geom_boxplot() +
  geom_hline(yintercept = 1.25)

p_phi_geweke  +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  ggplot2::scale_color_manual(values = c("grey", "black")) +
  ggplot2::scale_linetype_manual(values = c(2, 1)) +
  theme(
    axis.text.y = element_text(size = 10.5),
    axis.text.x = element_text(size = 10.5),
    axis.title.y = element_text(size = 10.5),
    axis.title.x = element_text(size = 10.5),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 10.5),
    strip.text = element_text(size = 10.5),
    legend.text = element_text(size = 10.5)
  )
 
ggsave("./SupplementaryMaterial/Images/Yeast/Convergence/gewekePhiChain.png",
       height = 10, width = 6
)

# Add column indicating if the chain should be used in the analysis
alloc_data$Use_chain <- alloc_data$Seed %in% chains_to_keep
continuous_data$Use_chain <- continuous_data$Seed %in% chains_to_keep

# Save the tibbles
saveRDS(alloc_data, file = "./Data/Yeast/BayesAllocTibble.rds")
saveRDS(continuous_data, file = "./Data/Yeast/BayesContParamsTibble.rds")

# == Parameter densities =======================================================

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
ggsave("./Images/Yeast/densityPlot.png")

p_param_density_reduced <- cont_plt_data %>%
  filter(Seed %in% chains_to_keep) %>%
  ggplot(aes(x = value, fill = Seed)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  labs(
    title = "Parameter density",
    x = "Value",
    y = "Density"
  ) +
  scale_fill_viridis_d() +
  theme(
    axis.text.y = element_text(size = 10.5),
    axis.text.x = element_text(size = 10.5),
    axis.title.y = element_text(size = 10.5),
    axis.title.x = element_text(size = 10.5),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 10.5),
    legend.text = element_text(size = 10.5)
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

ggsave("./SupplementaryMaterial/Images/Yeast/densityPlotReduced.png",
  plot = p_param_density_reduced,
  height = 4, width = 6
)

# === Time series ==============================================================

# Time series of the time course data separated out into predicted clusters

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

# === PSMs =====================================================================

# Compare PSMs for the different chains
# Reduce to OK chains
alloc_data_reduced <- alloc_data[alloc_data$Seed %in% chains_to_keep, ]
n_kept <- length(chains_to_keep)

# Iterate over each dataset and save a grid of heatmaps
for (dataset in datasets) {

  # The indices and CMs of the tibble corresponding to the current dataset
  curr_inds <- which(alloc_data_reduced$Dataset == dataset)
  curr_psms <- alloc_data_reduced$PSM[curr_inds]
  chains <- alloc_data_reduced$Seed[curr_inds]

  # Find the order for the rows and oclumns of the consensus matrices based on
  # the largest & deepest ensemble (here Consensus(1001, 1000))
  row_order <- findOrder(curr_psms[[1]])

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
        Y = N - match(Gene_j, item_order),
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

  if (dataset == datasets[[1]]) {
    psm_plt_data_all <- psm_plt
  } else {
    psm_plt_data_all <- rbind(psm_plt_data_all, psm_plt)
  }

  # Plot a grid of consensus matrices where column number is increasing the
  # number of samples used and row number corresponds to a value of chain depth
  p_psm <- psm_plt %>%
    ggplot(aes(X, Y, fill = Prop)) +
    geom_tile() +
    facet_wrap(~Chain, labeller = labeller(Chain = chain_labels)) +
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
  # ggsave(paste0(save_dir, dataset, "PSMcomparisonReduced.png"),
  #   plot = p_psm,
  #   height = 6,
  #   width = 4
  # )
}

dataset_labels <- datasets %>% 
  set_names(datasets)

psm_plt_data_all$Dataset <- factor(psm_plt_data_all$Dataset, levels = datasets)
psm_all <- psm_plt_data_all %>%
  ggplot(aes(X, Y, fill = Prop)) +
  geom_tile() +
  facet_grid(Dataset ~ Chain, labeller = labeller(Dataset = dataset_labels,
                                                  Chain = chain_labels)
             ) +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  labs(
    title = "Yeast dataset",
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

ggsave(paste0(save_dir, "YeastPSMcomparisonReduced.png"),
  plot = psm_all,
  height = 6,
  width = 6
)

# Iterate over each dataset and save a grid of heatmaps
  

#   # The indices and CMs of the tibble corresponding to the current dataset
#   curr_inds <- which(alloc_data_reduced$Dataset == "PPI")
#   curr_psms <- alloc_data_reduced$PSM[curr_inds]
#   chains <- alloc_data_reduced$Seed[curr_inds]
#   
#   # Find the order for the rows and oclumns of the consensus matrices based on
#   # the largest & deepest ensemble (here Consensus(1001, 1000))
#   row_order <- findOrder(curr_psms[[1]])
#   
#   # Set up the re-ordered items (this will be used to align the data when in
#   # long format)
#   item_order <- item_names[row_order]
#   
#   # Iterate over the different ensembles converting the consensus matrics to
#   # long format and prepare them for ggplot
#   for (i in 1:n_kept) {
#     .psm <- curr_psms[[i]]
#     row_order <- findOrder(.psm)
#     item_order <- item_names[row_order]
#     
#     
#     .df <- .psm[row_order, row_order] %>%
#       
#       # Convert to tibble adding the row names as a variable
#       as_tibble(rownames = "Gene_i") %>%
#       
#       # Pivot to long form
#       pivot_longer(cols = -Gene_i, names_to = "Gene_j", values_to = "Prop") %>%
#       
#       # Add variables to indicate the position of the grid for (gene_i, gene_j)
#       # on the heatmap based on the ordering defined previously and add the
#       # ensemble information
#       mutate(
#         X = match(Gene_i, item_order),
#         Y = N - match(Gene_j, item_order),
#         Chain = chains[i]
#       )
#     
#     # Bind the data frames together
#     if (i == 1) {
#       psm_plt <- .df
#     } else {
#       psm_plt <- rbind(psm_plt, .df)
#     }
#   }
#   
#   # Add a dataset label
#   psm_plt$Dataset <- "PPI"
#   
#   if (dataset == datasets[[1]]) {
#     psm_plt_data_all <- psm_plt
#   } else {
#     psm_plt_data_all <- rbind(psm_plt_data_all, psm_plt)
#   }
#   
#   # Plot a grid of consensus matrices where column number is increasing the
#   # number of samples used and row number corresponds to a value of chain depth
#   p_psm <- psm_plt %>%
#     ggplot(aes(X, Y, fill = Prop)) +
#     geom_tile() +
#     facet_wrap(~Chain, labeller = labeller(Chain = chain_labels)) +
#     scale_fill_gradient(low = "white", high = "#146EB4") +
#     labs(
#       title = dataset,
#       subtitle = "Posterior similarity matrices",
#       x = "Gene",
#       y = "Gene",
#       fill = "Coclustering\nproportion"
#     ) +
#     theme(
#       axis.text = element_blank(),
#       axis.ticks = element_blank(),
#       panel.grid = element_blank(),
#       axis.title.y = element_text(size = 10.5),
#       axis.title.x = element_text(size = 10.5),
#       plot.title = element_text(size = 18, face = "bold"),
#       plot.subtitle = element_text(size = 14),
#       strip.text.x = element_text(size = 10.5),
#       legend.text = element_text(size = 10.5)
#     )
#   
#   # Save!
#   ggsave(paste0("./SupplementaryMaterial/Images/Yeast/PPIPSMcomparisonReducedOwnOrder.png"),
#          plot = p_psm,
#          height = 14,
#          width = 12
#   )
# 
# my_order <- findOrder(orig_data$PPI)
# compareSimilarityMatricesAnnotated(matrices = list(orig_data$PPI[my_order, ], .psm[my_order, my_order]),
#                                    order_cols = F, order_rows = F)
