#!/usr/env/bin/Rscript

################################################################################
#                                                                              #
# Title: Bayesian inference performance                                        #
#                                                                              #
# Aim: Analyse saved clustering output of Mason's implementation of MDI for    #
# Bayesian inference simulations run April 2020.                               #
# Compare models to the ground truth using ARI and a Frobenius norm.           #
#                                                                              #
#                                                                              #
# Author: Stephen Coleman                                                      #
# Date: 04/05/2020                                                             #
#                                                                              #
# Example: Rscript BayesLargeNModelEval.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/large_n_small_p_base/simulation_1/Bayesian/ -s ./ --geweke ~/rds/hpc-work/Analysis/Simulations/Convergence/large_n_small_p_base/Within_chain/ --truth ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/large_n_small_p_base/ #
#                                                                              #
#                                                                              #
#                                                                              #
################################################################################

# For my ggplot theme settings
library(mdiHelpR)

# File name interactions
library(stringr)

# Plotting
library(ggplot2)

# I don't think this is used
library(tibble)

# To read in data
library(data.table)

# Visualising heatmaps (no longer used)
library(pheatmap)

# For frobenius.prod
library(matrixcalc)

# For arandi and maxpear for model evaluation and predicted clustering
library(mcclust)

# For facet_wrap_paginate
library(ggforce)

# For the saving of Sparse matrices
library(Matrix)

# For command line arguments
library(optparse)

# For pipe and related functions
library(magrittr)

# Combining plots
library(patchwork)

# For filter and group_by
library(dplyr)

# === Functions ================================================================

# The Frobenius norm is used to measure model uncertainty
#' @title Normalised Frobenius similarity
#' @description The Frobenius norm of the difference between matrices A and B is
#' calculated. It is then normalised by dividing by the sum of B less it diagonal
#' entries. As B is positive definite there is no need to take the absolute value
#' first.
#' @param A An $n \times n$ coclustering matrix for the true clustering structure
#' for the dataset being clustered.
#' @param B An $n \times n$ posterior similarity or consensus matrix for some
#' clustering method. All values should be in $[0,1]$.
#' @param n The number of rows/columns in A and B; default is NULL in which case
#' ``nrow(B)`` is used.
#' @return A (non-symmetric) measure of similarity between A and B.
normalisedFrobeniusSimilarity <- function(A, B) {

  # Check matrices dimensions match
  if (any(dim(A) != dim(B))) {
    stop("Matrix dimensions do not match.")
  }


  # Our distance measure
  distance <- sum(sqrt((A - B)**2))

  distance
}

inputArguments <- function() {
  option_list <- list(

    # Data to cluster
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = "./",
      help = "File path to model output [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-s", "--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save model and predicted clustering to [default= %default]",
      metavar = "character"
    ),

    # Expected chain length
    optparse::make_option(c("-l", "--length"),
      type = "character",
      default = "1000001",
      help = "Expected chain length; will call an error if any chains are less than this [default= %default]",
      metavar = "character"
    ),

    # Thinning factor applied within each chain
    optparse::make_option(c("-t", "--thin"),
      type = "character",
      default = "1000",
      help = "Thinning factor applied within each chain [default= %default]",
      metavar = "character"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("--n_seeds"),
      type = "integer",
      default = 10,
      help = "Number of chains in analysis [default= %default]",
      metavar = "integer"
    ),

    # Simulation number
    optparse::make_option(c("--sim_num"),
      type = "integer",
      default = 1,
      help = "Current simulation [default= %default]",
      metavar = "integer"
    ),

    # Directory containing truth
    optparse::make_option(c("--truth_dir"),
      type = "character",
      default = NA,
      help = "Directory containing truth [default= %default]",
      metavar = "character"
    ),

    # Directory containing Geweke statistic data for each chain
    optparse::make_option(c("--geweke_dir"),
      type = "character",
      default = NA,
      help = "Directory containing Geweke statistic data for each chain [default= %default]",
      metavar = "character"
    ),

    # Columns of interest for convergence
    optparse::make_option(c("-c", "--cols"),
      type = "character",
      default = "1",
      help = "Columns to subset by in ``fread`` [default= %default]",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# === Initialisation ===========================================================

args <- inputArguments()

# Directory to read from and to save Sparse matrices to
data_dir <- args$dir
save_dir <- args$save_dir
truth_dir <- args$truth_dir
geweke_dir <- args$geweke_dir

dir.create(save_dir)

# data_dir <- "/Users/stephen/Desktop/Testing_pipeline/simple_2d/MDI_output/simulation_1/Bayesian/"
# truth_dir <- "/Users/stephen/Desktop/Testing_pipeline/simple_2d/Input_data/"
# geweke_dir <- "/Users/stephen/Desktop/Convergence/large_n_small_p_large_k_small_dm/"


# Details of analysis
expected_length <- args$length
thin <- args$thin %>%
  str_split(" ") %>%
  unlist() %>%
  as.numeric()

n_seeds <- args$n_seeds
sim_num <- args$sim_num

# File to save results of model performance to
results_file <- paste0(save_dir, "BayesianSimulation", sim_num, "ResultsDF.csv")

# Columns to exclude as irrelevant to clustering
subset_cols <- args$cols %>%
  str_split(" ") %>%
  unlist() %>%
  as.numeric()

# For reproducibility
set.seed(1)

# Use a subset of people in the tests
# Remember to set subset cols to -1 to exclude Mass_Parameter
# subset_cols <- -1
# if (interactive()) {
#   subset_cols <- c(2:600)
#   n_people <- length(subset_cols)
# }

# Files containing simulation data and structure
cluster_id_file <- paste0(truth_dir, "cluster_IDs_", sim_num, ".Rds")
orig_data_file <- paste0(truth_dir, "dataset_", sim_num, ".csv")

# True labelling and data being clustered
truth <- readRDS(cluster_id_file)
orig_data <- read.csv(orig_data_file)

# True coclustering matrix
true_cc <- createSimilarityMat(matrix(truth, nrow = 1))

# Read in the data containinng Geweke statistics
geweke_data <- paste0(geweke_dir, "simulation_", sim_num, "GewekeData.csv") %>%
  read.csv()

# From the Gewke data, find which chains have converged
test_data <- geweke_data %>%
  group_by(Chain) %>%
  summarise(Shapiro_p_value = shapiro.test(Geweke_statistic)$p.value) %>%
  mutate(Normal = Shapiro_p_value > 0.05)

# Extract the chain number to remove these chains
chain_number <- test_data$Chain %>%
  str_extract_all("[1-9][0-9]*") %>%
  unlist() %>%
  as.numeric()

# List input files for reading in and for anlaying
files_full <- list.files(data_dir, full.names = T) %>%
  str_sort(numeric = T)

files_short <- list.files(data_dir, full.names = F) %>%
  str_sort(numeric = T)

# Find the details of each file based on a naming convention of:
# "${Model type}N{Number of iterations}T{thinning factor}S{random seed}.csv"
file_details <- files_short %>%
  str_extract_all("[1-9][0-9]*") %>%
  unlist() %>%
  as.numeric() %>%
  matrix(ncol = 4, byrow = T) %>%
  set_colnames(c("n", "thin", "seed", "k")) %>%
  as.data.frame()

# The chains to use ordered to match files
chains_to_use <- test_data$Normal[match(file_details$seed, chain_number)]

# Remove any non-converged chains
files_short <- files_short[chains_to_use]
files_full <- files_full[chains_to_use]
file_details <- file_details[which(chains_to_use), ]

# Inspect the file details
# file_details %>%
#   head()
#
# file_details %>% nrow()

# Seed order (NOT USED)
seed_order <- file_details$seed %>%
  order()

# Files ordered by number of iterations
iter_order <- file_details$n %>%
  order()

# Inspect
# file_details[seed_order, ] %>%
# head()

# Number of seeds and iterations present
n_seeds <- length(file_details$seed)
n_iter <- max(file_details$n)

# Iterations saved across all consensus inference models in simulations run April 2020
iter_used <- c(
  seq(1e3, 1e4, 1e3),
  seq(2e4, 1e5, 1e4),
  seq(2e5, 1e6, 1e5)
)

# Used for book-keeping
n_model <- length(iter_used)


# Inspect
# files_full[seed_order] %>%
# head()

# Find dimensions of datasets by reading in a 0-row data.table
s1 <- fread(files_full[1], drop = subset_cols, nrows = 0)
n_people <- ncol(s1)

# Old object to hold psms
# n_seed_df <- n_iter_df <- matrix(0, nrow = n_people, ncol = n_people)

# seed_psm_array <- array(0, c(n_people, n_people, n_seeds, n_model))

# Number of iterations and seeds to record for model evaluation
iters_to_record <- iter_used # c(1, 10, 100, 1000, 10000)
seeds_to_record <- file_details$seed

# The number of models recorded under each measure
n_iters_recorded <- length(iters_to_record)
n_seeds_recorded <- length(seeds_to_record)

# The data.frame that holds the score of each model
results_df <- data.frame(
  N_iter = rep(iters_to_record, n_seeds_recorded),
  Seed = rep(seeds_to_record, each = n_iters_recorded),
  ARI = 0,
  Frobenius_norm = 0
)

# NOT USED
# recorded_values <- c("N_iter", "Seed", "ARI", "Frobenius_norm")


for (i in 1:n_seeds) {
  curr_seed <- seeds_to_record[i]

  # Which consensus output matches current requirements of seed and chain
  # length
  file_used <- which(file_details$seed == curr_seed)[1]

  # Read in the relevant samples
  curr_samples <- fread(files_full[[file_used]], drop = subset_cols) %>%
    as.matrix(ncol = n_people)
  
  #print("Create small PSM.")
  # curr_psm <- createSimilarityMat(curr_samples)
  #print("Small PSM made.")
  if( i == 1 ) {

     # all_seeds_psm <- curr_psm
     all_samples <- curr_samples
  } else {
     #all_seeds_psm <- all_seeds_psm + curr_psm
     all_samples <- rbind(all_samples, curr_samples)

  }
}

print("Making PSM")
all_seeds_psm <- createSimilarityMat(all_samples)
print("PSM made")

sparse_psm <-  Matrix(all_seeds_psm, sparse = TRUE)
psm_file <- paste0(
  save_dir,
  "Simulation",
  sim_num,
  "PSM_all",
  ".txt"
)

writeMM(sparse_psm, file = psm_file)
#print("Calculate ARI.")
#all_ari <- apply(all_samples, 1, arandi, truth )
#print("ARI calculated")




# Predicted clustering from these
final_pred_cl <- maxpear(all_seeds_psm, max.k = 100)$cl

# Summary accuracy
final_ARI <-  arandi(final_pred_cl, truth[1:n_people])
# final_ARI <-  mean(all_ari) # 

# Uncertainty
final_Frobenius_norm <- normalisedFrobeniusSimilarity(true_cc, all_seeds_psm)

# Put these in a data.frame
final_results <- data.frame(n_iter, "Pooled", final_ARI, final_Frobenius_norm) %>% 
  set_colnames(colnames(results_df))

# Bind them onto the other results
results_df <- final_results # rbind(results_df, final_results)

# Add the simulation number as a variable in the data.frame
results_df$Simulation <- sim_num

# Save the model performance results to a file
write.csv(results_df, file = results_file)

psm_plot_file <- paste0(
  save_dir,
  "BayesianSimulation",
  sim_num,
  "PSMs",
  ".png"
)

psm_main <- paste0("Simulation ", sim_num, ": Pooled PSM")

annotatedHeatmap(all_seeds_psm, truth, 
  col_pal = simColPal(),
  my_breaks = defineBreaks(col_pal = col_pal, lb = 0),
  main = psm_main,
  filename = psm_plot_file,  
  show_rownames = F,
  show_colnames = F
)


# === Plotting =================================================================

if (interactive()) {

  # Set labels for facet wrapping
  iter_labels <- c(paste0("Number of iterations: ", results_df$N_iter))
  names(iter_labels) <- results_df$N_iter

  seed_labels <- c(paste0("Chain ", results_df$Seed))
  names(seed_labels) <- results_df$Seed

  # Plots
  # p_ari_iter <- results_df %>%
  #   filter(N_iter %in% c(1e3, 1e4, 1e5, 1e6)) %>%
  #   ggplot(aes(x = Seed, y = ARI)) +
  #   geom_point() +
  #   facet_wrap(~N_iter, labeller = labeller(N_iter = iter_labels), ncol = 1) # +
  # labs(
  #   title = "Performance stabilitses by 100 iterations",
  #   subtitle = "CI for large n, large k, small distance between means",
  #   x = "Number of chains used"
  # )

  p_ari_seed <- results_df %>%
    ggplot(aes(x = N_iter, y = ARI)) +
    geom_line() +
    facet_wrap(~Seed, labeller = labeller(Seed = seed_labels), ncol = 1) +
    labs(
      # title = "Predictive performance",
      subtitle = "ARI between summary clustering from PSM and true labelling",
      x = "Number of iterations used"
    )


  # p_unc_iter <- results_df %>%
  #   ggplot(aes(x = Seed, y = Frobenius_norm)) +
  #   geom_line() +
  #   facet_wrap(~N_iter, labeller = labeller(N_iter = iter_labels), ncol = 1) # +
  # labs(
  #   title = "Performance stabilitses after 40 iterations",
  #   subtitle = "CI for large n, large k, small distance between means",
  #   x = "Number of chains used"
  # )

  p_unc_seed <- results_df %>%
    ggplot(aes(x = N_iter, y = Frobenius_norm)) +
    geom_line() +
    facet_wrap(~Seed, labeller = labeller(Seed = seed_labels), ncol = 1) +
    labs(
      # title = "Frobenius norm",
      # subtitle = "CI for large n, large k, small distance between means",
      x = "Number of iterations used",
      y = "Frobenius norm",
      subtitle = "Frobenius norm of difference between PSM and true coclustering matrix"
    )

  # p_ari_iter + p_ari_seed
  # p_unc_iter + p_unc_seed

  # p_ari_iter + p_unc_iter
  patch <- p_ari_seed + p_unc_seed
  patch + plot_annotation(
    title = paste0(scenario, ": simulation ", sim_num),
    subtitle = "Bayesian inference"
  )
}


#
# results_df %>%
#   ggplot(aes(x = N_iter, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~Seed, ncol = 3, nrow = 3, page =3)
#
# results_df %>%
#   ggplot(aes(x = Seed, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 4, nrow = 4, page =1) +
# labs(
#   title = "Performance stabilitses after 40 iterations",
#   subtitle = "CI for large n, large k, small distance between means",
#   x = "Number of chains used"
# )
#
# results_df %>%
#   ggplot(aes(x = N_iter, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~Seed, ncol = 4, nrow = 4, page =1) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI for large n, large k, small distance between means",
#     x = "Number of chains used"
#   )
#
# results_df %>%
#   ggplot(aes(x = Seed, y = ARI)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 4, nrow = 4, page =1,
#                       labeller = labeller(N_iter = facet_labels)) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI performance for large n, large k, small distance between means (ARI)",
#     x = "Number of chains used"
#   )
#
#
# results_df %>%
#   ggplot(aes(x = Seed, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 4, nrow = 4, page =3,
#                       labeller = labeller(N_iter = facet_labels)) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI performance for large n, large k, small distance between means (Uncertainty)",
#     x = "Number of chains used"
#   )
#
#
#
# results_df %>%
#   ggplot(aes(x = N_iter, y = ARI)) +
#   geom_point() +
#   facet_wrap_paginate(~Seed, ncol = 4, nrow = 4, page =1) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI for large n, large k, small distance between means",
#     x = "Number of chains used"
#   )
#
#
#
# results_df %>%
#   ggplot(aes(x = Seed, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 3, nrow = 3, page =2)
#
# results_df %>%
#   ggplot(aes(x = Seed, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 3, nrow = 3, page = 3)
#
#
#
#
# sim_col <- simColPal()
# my_breaks <- defineBreaks(sim_col, lb = 0)
# seed_psm_array[, , 37] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
#
# seed_psm_array[, , 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
#
# seed_psm_array[, , 10] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
#
# seed_psm_array[, , 19] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
#
#
# seed_psm_array[, , 28] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
#
#
# my_c <- maxpear(seed_psm_array[, , 6])$cl
# arandi(my_c, truth[1:n_people])

# annotatedHeatmap(orig_data[1:n_people, 2:5], my_c)
#
# my_c <- maxpear(seed_psm_array[, , 20])$cl
# annotatedHeatmap(orig_data[1:n_people, 2:5], my_c)
#
# seed_psm_array[, , 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
#
# seed_psm_array[, , 1, 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
#
# seed_psm_array[, , 10, 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
#
#
# seed_psm_array[, , 10, 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)

