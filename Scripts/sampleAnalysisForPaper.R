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

# Directory where MDI output files live
main_dir <- "C:/Users/stephen/Documents/PhD/Year_1/Consensus_inference/Consensus_inference_gen/Analysis/MDI_yeast_dataset/CMDLineMDI"
data_dir <- "C:/Users/stephen/Documents/PhD/Year_1/Consensus_inference/Consensus_inference_gen/Data/YeastData"
datasets <- c("Timecourse", "ChIP-ChIP", "PPI")
data_files <- list.files(data_dir, full.names = T)[c(2, 1, 3)]
orig_data <- data_files %>%
  lapply(read.csv, row.names = 1) %>%
  set_names(datasets)

# All subdirs
all_dirs <- list.dirs(main_dir, recursive = F)

# Dirs for Bayesian inference
bayes_dirs <- all_dirs[grep(all_dirs, pattern = "*BayesianYeastN*")]
cc_dir <- paste0(main_dir, "/ConsensusClustering/R500S100/")

# MDI output files
bayes_files <- bayes_dirs %>%
  lapply(list.files, pattern = "*BayesianYeastN676000T1000Seed*", full.names = T) %>%
  str_sort(numeric = T) %>%
  lapply(read.csv)

cc_file <- cc_dir %>% 
  list.files(pattern = ".csv$", full.names = T) %>% 
  read.csv()

# Number of chains used
n_chains <- length(bayes_files)

# Number of datasets
L <- colnames(bayes_files[[1]])[grep("MassParameter", colnames(bayes_files[[1]]))] %>%
  length()

# Number of phi parameters from MDI
n_phis <- L * (L - 1) / 2

# Index of continuous variables in csv files
continuous_params_col_index <- 1:(L + n_phis)

dataset_samples <- list()
for (l in 1:L) {
  dataset_samples[[l]] <- grep(paste0("Dataset", l), colnames(bayes_files[[1]]))
}

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
    psm <- createSimilarityMat(as.matrix(alloc))

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
saveRDS(alloc_data, file = "./Data/yeastAllocTibble.rds")
saveRDS(continuous_data, file = "./Data/yeastContParamsTibble.rds")

# The consensus clustering results
cc_alloc <- .alloc_tib <- tibble(Samples = list(), Dataset = character(), CM = list(), Cl = list())
# MDI output
.curr_data <- cc_file

# Continuous parameters
cc_continuous <- .curr_data[, continuous_params_col_index]

# Dataset specific allocations
for (l in 1:L) {
  
  # Allocation samples
  .alloc <- .curr_data[, dataset_samples[[l]]]
  
  # CM
  cm <- createSimilarityMat(as.matrix(.alloc))
  
  # Predicted clustering (k is set to an arbitrarily large number)
  cl <- mcclust::maxpear(cm, max.k = 275)$cl
  
  # Save to a tibble
  .alloc_tib <- tibble(Samples = list(.alloc), Dataset = datasets[[l]], CM = list(cm), Cl = list(cl))
  cc_alloc <- rbind(cc_alloc, .alloc_tib)
}

# === Bayes Visualisation ======================================================

# Facet labels for continuous variables
param_labels <- c(paste0("alpha[", 1:3, "]"), paste0("phi[", 1:3, "]"))
names(param_labels) <- colnames(continuous_data$Samples[[1]])

# Plot the timecourse data as a time series
plt_data <- orig_data$Timecourse
plt_data$Gene <- row.names(plt_data)
plt_data$Cluster <- alloc_data$Cl[[1]]
plt_data <- plt_data %>%
  pivot_longer(cols = -c(Gene, Cluster), values_to = "Expression") %>%
  mutate(Time = stringr::str_extract(name, "[:digit:]+"))

plt_data %>%
  ggplot(aes(x = Time, y = Expression, group = Gene)) +
  geom_line()

# Plot time series by cluster
plt_data %>%
  ggplot(aes(x = Time, y = Expression, group = Gene)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~Cluster)

pheatmap(orig_data$Timecourse, color = dataColPal(), cluster_cols = F)

annotatedHeatmap(orig_data$Timecourse, alloc_data$Cl[[1]], cluster_cols = F)

cl_order <- order(alloc_data$Cl[[1]])
annotatedHeatmap(orig_data$Timecourse[cl_order, ], alloc_data$Cl[[1]][cl_order],
  cluster_cols = F,
  cluster_rows = F
)

row_order <- findOrder(orig_data$Timecourse)
compareMatrices(orig_data$Timecourse[row_order, ], alloc_data$PSM[[1]][row_order, row_order], order_cols = F, order_rows = F)

compareSimilarityMatricesAnnotated(
  matrices = list(orig_data$Timecourse[row_order, ], alloc_data$PSM[[1]][row_order, row_order]),
  order_rows = F,
  order_cols = F,
  # cluster_IDs =  alloc_data$Cl[[1]],
  collect_legend = F,
  col_pal = list(dataColPal(), simColPal()),
  breaks = list(
    defineDataBreaks(orig_data$Timecourse[row_order, ], col_pal = dataColPal()),
    defineBreaks(lb = 0, col_pal = simColPal())
  )
)

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

# === CC Visualisation =========================================================

# Plot the timecourse data as a time series
plt_data <- orig_data$Timecourse
plt_data$Time <- as.numeric(plt_data$Time)
plt_data$Gene <- row.names(plt_data)
plt_data$Cluster <- cc_alloc$Cl[[1]]
plt_data <- plt_data %>%
  pivot_longer(cols = -c(Gene, Cluster), values_to = "Expression") %>%
  mutate(Time = stringr::str_extract(name, "[:digit:]+"))

plt_data %>%
  ggplot(aes(x = Time, y = Expression, group = Gene)) +
  geom_line()

# Plot time series by cluster
p_cc_timecourse <- plt_data %>%
  ggplot(aes(x = Time, y = Expression, group = Gene)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~Cluster) +
  labs(
    title = "Timecourse data",
    subtitle= "Facetted by predicted clustering from Consensus clustering"
  )

# PLace continuous parameters in a single data.frame suitable for gglot2
cc_cont_long <- cc_continuous %>%
    pivot_longer(cols = everything(), names_to = "Parameter")
  
p_param_density_cc <- cc_cont_long %>%
  ggplot(aes(x = value)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~Parameter, labeller = as_labeller(param_labels, label_parsed)) +
  labs(
    title = "Parameter density",
    x = "Value",
    y = "Density"
  ) +
  scale_fill_viridis_d()

p_param_density_cc
p_param_density

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

# === Compare merged PSM and CM ================================================
datasets
dataset_curr <- "PPI"
bayes_ind <- which(alloc_data$Dataset == dataset_curr)
merge_psm <- Reduce(`+`, alloc_data$PSM[bayes_ind]) / n_chains
pheatmap(merge_psm, color = simColPal())

pred_cl <- maxpear(merge_psm, max.k = 150)$cl
unique(pred_cl)
annotatedHeatmap(orig_data[[dataset_curr]], pred_cl, cluster_cols = F)

cc_ind <- which(cc_alloc$Dataset == dataset_curr)
compareSimilarityMatricesAnnotated(merge_psm, cc_alloc$CM[[cc_ind]])
compareSimilarityMatricesAnnotated(cc_alloc$CM[[cc_ind]], merge_psm)
annotatedHeatmap(orig_data[[dataset_curr]], cc_alloc$Cl[[cc_ind]], cluster_cols = F)

ari_mat <- matrix(0, nrow = n_chains, n_chains)
diag(ari_mat) <- 1
for(i in 1:(n_chains - 1)){
  for(j in (i+1):n_chains){
    ari_mat[i, j] <- ari_mat[j, i] <- arandi(alloc_data$Cl[bayes_ind][[i]],alloc_data$Cl[bayes_ind][[j]])
  }
  # ari_mat[i, n_chains + 1] <- ari_mat[n_chains + 1, i] <-  arandi(alloc_data$Cl[bayes_ind][[i]], pred_cl)
}
# ari_mat[n_chains, n_chains + 1] <- ari_mat[n_chains + 1, n_chains] <-  arandi(alloc_data$Cl[bayes_ind][[n_chains]], pred_cl)

ari_mat
sim_col <- simColPal()
breaks <- defineBreaks(sim_col, lb = 0)
pheatmap(ari_mat, color = sim_col, breaks = breaks)

boxplot(ari_mat)
hist(ari_mat)

alloc_data$Cl[bayes_ind] %>% 
  lapply(arandi, cc_alloc$Cl[[cc_ind]]) %>% 
  unlist() %>% 
  hist()

arandi(pred_cl, cc_alloc$Cl[[cc_ind]])

alloc_data$Cl[bayes_ind] %>% 
  lapply(arandi, pred_cl) %>% 
  unlist() %>% 
  hist()

# === Compare similarity matrices ==============================================

for(dataset_curr in datasets){
bayes_ind <- which(alloc_data$Dataset == dataset_curr)
similarity_matrices <- alloc_data$PSM[bayes_ind]

cc_ind <- which(cc_alloc$Dataset == dataset_curr)
similarity_matrices[[length(similarity_matrices) + 1]] <- cc_alloc$CM[[cc_ind]]

col_pal <- breaks <- NULL
chain_to_lead <- 1
if(dataset_curr == "Timecourse"){
  col_pal <- list(dataColPal(), simColPal(), simColPal())
  breaks = list(defineDataBreaks(x = orig_data[[dataset_curr]], col_pal = dataColPal()),
                defineBreaks(simColPal(), lb = 0),
                defineBreaks(simColPal(), lb = 0)
  )
  
  # Based on GO over representation this chain appears to be more representative of the most popular result
  chain_to_lead <- 6
}
if (dataset_curr == "ChIP-ChIP"){
  chain_to_lead <- 5
}

# if (dataset_curr == "PPI"){
#   chain_to_lead <- 5
# }

compareSimilarityMatricesAnnotated(matrices = similarity_matrices,
  title = paste0(dataset_curr, "PSMs and CM ordered by chain ", chain_to_lead),
  matrix_imposing_order = chain_to_lead
)

ggsave(paste0("./Images/Yeast/", dataset_curr, "/", dataset_curr, "SimilarityMatCompChain", chain_to_lead, ".png"),
       height = 10, 
       width = 12)

compareSimilarityMatricesAnnotated(matrices = similarity_matrices, 
                                   title = paste0(dataset_curr, "PSMs and CM ordered by CM"),
                                   matrix_imposing_order = length(similarity_matrices)
                                   )

ggsave(paste0("./Images/Yeast/", dataset_curr, "/", dataset_curr, "SimilarityMatCompCM.png"),
       height = 10, 
       width = 12)

row_order <- findOrder(similarity_matrices[[chain_to_lead]])
col_order <- findOrder(t(orig_data[[dataset_curr]]))


compareSimilarityMatricesAnnotated(orig_data[[dataset_curr]][row_order, col_order], 
  similarity_matrices[[chain_to_lead]][row_order, row_order], 
  cc_alloc$CM[[cc_ind]][row_order, row_order],
  title = paste0(dataset_curr, " data, PSM (chain ", chain_to_lead, "), CM ordered by PSM"),
  order_rows = F, 
  order_cols = F,
  col_pal = col_pal,
  breaks = breaks
)

ggsave(paste0("./Images/Yeast/", dataset_curr, "/", dataset_curr, "dataPSMandCM.png"),
       height = 6, 
       width = 6)

}

