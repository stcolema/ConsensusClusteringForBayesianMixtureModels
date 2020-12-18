#!/usr/env/bin/Rscript

################################################################################
#                                                                              #
# Title: Consensus inference performance                                       #
#                                                                              #
# Aim: Analyse saved clustering output of Mason's implementation of MDI for    #
# Consensus inference simulations run April 2020.                              #
# Compare models to the ground truth using ARI and a Frobenius norm.           #
#                                                                              #
#                                                                              #
# Author: Stephen Coleman                                                      #
# Date: 21/04/2020                                                             #
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
suppressMessages(library(mcclust))

# For facet_wrap_paginate
library(ggforce)

# For the saving of Sparse matrices
library(Matrix)

# For command line arguments
library(optparse)

# For pipe and related functions
library(magrittr)

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
  if(any(dim(A) != dim(B))){
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
      default = "10 100 1000 10001",
      help = "Expected chain length; will call an error if any chains are less than this [default= %default]",
      metavar = "character"
    ),

    # Thinning factor applied within each chain
    optparse::make_option(c("-t", "--thin"),
      type = "character",
      default = "1 10 100 1000",
      help = "Thinning factor applied within each chain [default= %default]",
      metavar = "character"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("--n_seeds"),
      type = "integer",
      default = 100,
      help = "Number of chains in analysis [default= %default]",
      metavar = "integer"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("--sim_num"),
      type = "integer",
      default = 1,
      help = "Current simulation [default= %default]",
      metavar = "integer"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("--truth_dir"),
      type = "character",
      default = NA,
      help = "Directory containing truth [default= %default]",
      metavar = "character"
    ),

    # Columns of interest for convergence
    optparse::make_option(c("-c", "--cols"),
      type = "character",
      default = "1",
      help = "Columns to subset by in ``fread`` [default= %default]",
      metavar = "character"
    ),


    # Scenario being analysed
    optparse::make_option(c("--scn"),
      type = "character",
      default = NA,
      help = "Scenarion being analysed [default= %default]",
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

# Create save directory if not existing
dir.create(save_dir, showWarnings = F)


# Details of analysis
expected_length <- args$length
thin <- args$thin %>%
  str_split(" ") %>%
  unlist() %>%
  as.numeric()

n_seeds <- args$n_seeds
sim_num <- args$sim_num


# scenario in question
scn <- args$scn
scn_pretty <- scn %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::str_to_sentence()


# File to save results of model performance to
results_file <- paste0(save_dir, "ConsensusSimulation", sim_num, "ResultsDF.csv")
ph_file <- paste0(save_dir, "ConsensusSimulation", sim_num, "ConsensusMatrixGrid.png")

annotated_ph_file <- paste0(save_dir, "ConsensusAnnotatedHeatmap")

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
orig_data <- read.csv(orig_data_file, row.names = 1)

# True coclustering matrix
true_cc <- createSimilarityMat(matrix(truth, nrow = 1))
N <- nrow(true_cc)

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

#if(50 %in% file_details){
#print("K might be present in some file names")
#stop("Look at output files.")
#}

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
n_seeds <- max(file_details$seed)
n_iter <- max(file_details$n)

# Iterations saved across all consensus inference models in simulations run April 2020
iter_used <- c(
  seq(1, 10, 1),
  seq(20, 100, 10),
  seq(200, 1000, 100),
  seq(2000, 10000, 1000)
)

# Used for book-keeping
n_model <- length(iter_used)

# Hold matrices to compare heatmaps
grid_matrices <- data.frame(Iter = c(1, 10, 100, 1000, 10000), N_chain = c(1, 10, 30, 50, 100))
m_list <- vector("list", nrow(grid_matrices)**2)

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
iters_to_record <- c(1, 10, 100, 500, 1000, 5000, 10000)
seeds_to_record <- c(1, 10, 30, 50, 100)

# The number of models recorded under each measure
n_iters_recorded <- length(iters_to_record)
n_seeds_recorded <- length(seeds_to_record)

# The data.frame that holds the score of each model
results_df <- data.frame(
  N_iter = rep(iters_to_record, n_seeds_recorded),
  N_seeds = rep(seeds_to_record, each = n_iters_recorded),
  ARI = 0,
  Frobenius_norm = 0
)

# NOT USED
# recorded_values <- c("N_iter", "N_seeds", "ARI", "Frobenius_norm")

# The array that hold the consensus matrix for each model
# seed_psm_array <- array(0, c(n_people, n_people, n_model))

# === Analysis =================================================================

for (j in 1:n_model) {
  # Current number of iterations used
  curr_iter <- iter_used[j]

  samples_j <- NULL # matrix(NA, nrow = n_seeds, ncol = N)

  for (i in 1:n_seeds) {

    .curr_cm <- NULL

    # Which consensus output matches current requirements of seed and chain
    # length
    file_used <- which(file_details$n > curr_iter & file_details$seed == i)[1]

    # What thinning is present in the current file (influences which sample is
    # kept, if in a file where we want the thousandth sample and there is
    # thinning of 100, we want the (1000/100 + 1)=11th sample, with the +1
    # correcting for Mason always saving the 0th iteration)
    thin_present <- file_details$thin[file_used]

    # Read in the relevant samples
    curr_samples <- fread(files_full[[file_used]], drop = subset_cols, nrows = 1 + curr_iter / thin_present) %>%
      as.matrix(ncol = n_people)

    print(curr_iter)
    print(dim(curr_samples))
    print(curr_samples[1 + curr_iter / thin_present, 1:6])

    if(is.null(samples_j)){
      samples_j <- curr_samples[1 + curr_iter / thin_present,]
    } else {
      samples_j <- rbind(samples_j, curr_samples[ 1 + curr_iter / thin_present, ])
    }

    # Create a coclustering matrix based on the current sample
    #curr_psm <- createSimilarityMat(matrix(curr_samples[1 + curr_iter / thin_present, ], nrow = 1))

    # Update the consensus matrix
    #seed_psm_array[, , j] <- .curr_cm <- seed_psm_array[, , j] * (i - 1) / i + curr_psm * (1 / i)
    
    #if(curr_iter %in% grid_matrices$Iter && i %in% grid_matrices$N_chain){
    #  # grid_matrices
    #  grid_row <- which(curr_iter == grid_matrices$Iter)
    #  grid_col <- which(i == grid_matrices$N_chain)
    #  m_entry <- (grid_row-1)*nrow(grid_matrices) + grid_col

    #  .curr_cm <- createSimilarityMat(samples_j[1:i, ])
    #  m_list[[m_entry]] <- .curr_cm
    #}

    # If current model is to to be saved, perform inference
    if (curr_iter %in% iters_to_record & i %in% seeds_to_record) {
      print(str(samples_j))
      print("Making CM.")
      
      #if(is.null(.curr_cm)){
      if(i == 1){
        .curr_cm <- createSimilarityMat(matrix(samples_j, nrow = 1))
      } else {
        .curr_cm <- createSimilarityMat(samples_j)
      }

      # Save a sparse representation of the similarity matrix
      sparse_cm <- Matrix(.curr_cm, sparse = TRUE)
      cm_file <- paste0(
        save_dir,
        "Simulation",
        sim_num,
        "ConsensusMatrixN",
        curr_iter,
        "S",
        i,
        ".txt"
      )
      writeMM(sparse_cm, file = cm_file)

      # Predicted clustering
      curr_pred_cl <- maxpear(.curr_cm)$cl

      # Record frobenius product and ARI
      record_ind <- which(results_df$N_iter == curr_iter & results_df$N_seeds == i)
      results_df$ARI[record_ind] <- arandi(curr_pred_cl, truth[1:n_people])
      results_df$Frobenius_norm[record_ind] <- normalisedFrobeniusSimilarity(true_cc, .curr_cm)

      annotatedHeatmap(orig_data, curr_pred_cl, 
        main = paste0(scn_pretty, " simulation ", sim_num, ": Consensus (", curr_iter, ", ", i, ")" ),
        filename = paste0(annotated_ph_file, "N", curr_iter, "S", i, ".png")
      )
      
    }
  }
}

# Add the simulation number as a variable in the data.frame
results_df$Simulation <- sim_num

# Save the model performance results to a file
write.csv(results_df, file = results_file)

# Save the heatmap grid to a file
#ph_grid <- suppressWarnings(compareSimilarityMatricesAnnotated(matrices = m_list,
#  title = paste0("Simulation ", sim_num, ": comparison of Consensus Matrices for different settings"),
#  cluster_IDs = truth,
#  matrix_imposing_order = 25,
#  fill = "#D3D3D3"
#)
#)
#
#ggsave(ph_file,
#       plot = ph_grid, 
#       width = 24,
#       height = 16)
#
# === Plotting =================================================================

if (interactive()) {

  # Set labels for facet wrapping
  iter_labels <- c(paste0("Number of iterations: ", results_df$N_iter))
  names(iter_labels) <- results_df$N_iter

  seed_labels <- c(paste0("Number of chains: ", results_df$N_seeds))
  names(seed_labels) <- results_df$N_seeds

  # Plots
  p_ari_iter <- results_df %>%
    ggplot(aes(x = N_seeds, y = ARI)) +
    geom_line() +
    facet_wrap(~N_iter, labeller = labeller(N_iter = iter_labels), ncol = 1) # +
  # labs(
  #   title = "Performance stabilitses by 100 iterations",
  #   subtitle = "CI for large n, large k, small distance between means",
  #   x = "Number of chains used"
  # )

  p_ari_seed <- results_df %>%
    ggplot(aes(x = N_iter, y = ARI)) +
    geom_line() +
    facet_wrap(~N_seeds, labeller = labeller(N_seeds = seed_labels), ncol = 1) # +
  # labs(
  # title = "Performance stabilitses by 100 iterations",
  # subtitle = "CI for large n, large k, small distance between means",
  # x = "Number of iterations used"
  # )


  p_unc_iter <- results_df %>%
    ggplot(aes(x = N_seeds, y = Frobenius_norm)) +
    geom_line() +
    facet_wrap(~N_iter, labeller = labeller(N_iter = iter_labels), ncol = 1) # +
  # labs(
  #   title = "Performance stabilitses after 40 iterations",
  #   subtitle = "CI for large n, large k, small distance between means",
  #   x = "Number of chains used"
  # )

  p_unc_seed <- results_df %>%
    ggplot(aes(x = N_iter, y = Frobenius_norm)) +
    geom_line() +
    facet_wrap(~N_seeds, labeller = labeller(N_seeds = seed_labels), ncol = 1)

  p_ari_iter + p_ari_seed
  p_unc_iter + p_unc_seed

  p_ari_iter + p_unc_iter
}


#
# results_df %>%
#   ggplot(aes(x = N_iter, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_seeds, ncol = 3, nrow = 3, page =3)
#
# results_df %>%
#   ggplot(aes(x = N_seeds, y = Uncertainty)) +
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
#   facet_wrap_paginate(~N_seeds, ncol = 4, nrow = 4, page =1) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI for large n, large k, small distance between means",
#     x = "Number of chains used"
#   )
#
# results_df %>%
#   ggplot(aes(x = N_seeds, y = ARI)) +
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
#   ggplot(aes(x = N_seeds, y = Uncertainty)) +
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
#   facet_wrap_paginate(~N_seeds, ncol = 4, nrow = 4, page =1) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI for large n, large k, small distance between means",
#     x = "Number of chains used"
#   )
#
#
#
# results_df %>%
#   ggplot(aes(x = N_seeds, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 3, nrow = 3, page =2)
#
# results_df %>%
#   ggplot(aes(x = N_seeds, y = Uncertainty)) +
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


