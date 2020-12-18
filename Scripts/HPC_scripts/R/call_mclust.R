#!/usr/en/bin/Rscript

##############################################################
#
# Frequentist mixture model
# Call Mclust from the command line with arguments
#
# Stephen Coleman
#
##############################################################

# To read in arguments from the command line
library(optparse)

# For mixture modelling
library(mclust)

# For the pipe
library(magrittr)

# User inputs from command line
inputArguments <- function() {
  option_list <- list(

    # Data to cluster
    optparse::make_option(c("-d", "--data"),
      type = "character",
      default = NA,
      help = "File path to data to cluster [default= %default]",
      metavar = "character"
    ),

    # Number of clusters to model within the data
    optparse::make_option(c("-k", "--n_clust"),
      type = "integer",
      default = 2,
      help = "Lower bouud on number of clusters to model within the data [default= %default]",
      metavar = "integer"
    ),

    # Upper bound on number of clusters to model within the data
    optparse::make_option(c("-g", "--n_clust_upper"),
      type = "integer",
      default = 50,
      help = "Upper bound on range of number of clusters to model within the data [default= %default]",
      metavar = "integer"
    ),

    # Model type
    optparse::make_option(c("-m", "--model_names"),
      type = "character",
      default = NULL,
      help = "Vector of model types considered by Mclust [default= %default]",
      metavar = "character"
    ),

    # File path to save model under 
    optparse::make_option(c("-f", "--model_file"),
      type = "character",
      default = "./mclustModel.Rds",
      help = "Filepath and name to save final model under [default= %default]",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


# Read in arguments from the command line
args <- inputArguments()

# Transfer these to variables
my_data <- args$data %>%
  read.csv(row.names = 1)

# The range of cluster to model
K1 <- args$n_clust
K2 <- args$n_clust_upper

# The input to Mclust
G <- seq(K1, K2)

# Model types to try in Mclust
model_names <- args$model_names

# If the user has given an input of NULL convert back from a string
if(model_names == "NULL"){
  model_names <- NULL
}

# Save details
save_dir <- args$dir
model_file <- args$model_file
cluster_file <- args$cluster_file

# Call Mclust
my_mod <- Mclust(my_data, G = G, modelNames = model_names)

# Save the model
saveRDS(my_mod, file = paste0(save_dir, model_file))

