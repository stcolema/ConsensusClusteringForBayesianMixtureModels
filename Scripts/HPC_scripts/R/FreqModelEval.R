#!/usr/env/bin/Rscript

################################################################################
#                                                                              #
# Title: Frequentists inference performance                                    #
#                                                                              #
# Aim: Analyse saved clustering output of Mason's implementation of MDI for    #
# Frequentists inference simulations run April 2020.                           #
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
library(mcclust)

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

dir.create(save_dir)

# data_dir <- "/Users/stephen/Desktop/Testing_pipeline/simple_2d/MDI_output/simulation_1/Consensus/"
# truth_dir <- "/Users/stephen/Desktop/Testing_pipeline/simple_2d/Input_data/"

# Details of analysis
sim_num <- args$sim_num

# File to save results of model performance to
results_file <- paste0(save_dir, "FrequentistsSimulation", sim_num, "ResultsDF.csv")

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
K <- length(unique(truth))
orig_data <- read.csv(orig_data_file)

# True coclustering matrix
true_cc <- createSimilarityMat(matrix(truth, nrow = 1))

# List input files for reading in and for anlaying
files_full <- list.files(data_dir, full.names = T) %>%
  str_sort(numeric = T)

# Interpret MClust model results
f_model <- readRDS(files_full[1])
f_cl <- f_model$classification
K_hat <- length(unique(f_cl))

f_ari <- arandi(f_cl, truth)
f_cc <- createSimilarityMat(matrix(f_cl, nrow = 1))

f_frob <- normalisedFrobeniusSimilarity(true_cc, f_cc)

# The data.frame that holds the score of each model
results_df <- data.frame(
  N_iter = NA,
  N_seeds = NA,
  ARI = f_ari,
  Frobenius_norm = f_frob,
  K_diff = K - K_hat
)

# Add the simulation number as a variable in the data.frame
results_df$Simulation <- sim_num

# Save the model performance results to a file
write.csv(results_df, file = results_file)

