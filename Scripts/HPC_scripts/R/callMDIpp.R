#!/usr/bin/Rscript

# Call MDI++ from within R
# Example 
# Rscript callMDIpp.R --datasets "~/rds/hpc-work/Input_data/Yeast_data/Reduced/Granovskaia_timecourse_normalised_reduced.csv ~/rds/hpc-work/Input_data/Yeast_data/Reduced/harbison_marina.csv ~/rds/hpc-work/Input_data/Yeast_data/Reduced/ppi.csv" --datatypes "GP M M" --n_iter 1000 --mdi_dir "/home/sdc56/hpc_stuff/scripts/mdipp-1.0.1/" --thin 1 -k 50 -s 1
#
# Stephen Coleman
# 15/07/2020

# String manipulation
library(stringr)

# For command line arguments
library(optparse, quietly = T) # install.packages("optparse")

# Pipe
library(magrittr)

inputArguments <- function() {
  option_list <- list(

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--datasets"),
      type = "character",
      default = NA,
      help = "names of dataset files [default= %default]",
      metavar = "character"
    ),

    # Instruction to time programme
    optparse::make_option(c("--datatypes"),
      type = "character",
      default = NULL,
      help = "Datatypes of given datasets [default= %default]",
      metavar = "character"
    ),

    # Number of iterations
    optparse::make_option(c("-r", "--n_iter"),
      type = "integer",
      default = NULL,
      help = "Number of iterations to run MDI for [default= %default]",
      metavar = "integer"
    ),

 
    # mdipp home directory
    optparse::make_option(c("--mdi_dir"),
      type = "character",
      default = NULL,
      help = "Directory where mdipp lives [default= %default]",
      metavar = "character"
    ),

    # Maximum number of clusters
    optparse::make_option(c("-k", "--n_clust"),
      type = "integer",
      default = 50,
      help = "Upper limit on number of clusters in MDI [default= %default]",
      metavar = "integer"
    ),

    # Thinning factor
    optparse::make_option(c("-t", "--thin"),
      type = "integer",
      default = 1,
      help = "Thinning factor for MCMC samples generated [default= %default]",
      metavar = "integer"
    ),

    # Random seed
    optparse::make_option(c("-s", "--seed"),
      type = "integer",
      default = 1,
      help = "Random seed controlling initialisation and sampling [default= %default]",
      metavar = "integer"
    ),

    # Output file
    optparse::make_option(c("-o", "--output"),
      type = "character",
      default = NULL,
      help = "File name for output [default= %default]",
      metavar = "character"
    ),

    # Call MDI or print call command
    optparse::make_option(c("-c", "--call"),
      type = "logical",
      default = TRUE,
      help = "Flag to call MDI; if FALSE will print command for the call [default= %default]",
      metavar = "logical"
    )


  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


runMDI <- function(datasets, dataTypes, R,
                   mdipp_dir = NULL,
                   K = 50,
                   thin = 1,
                   seed = 1,
                   output = NULL,
                   call = TRUE
                   ){
  
  L <- length(datasets)
  
  if(L != length(dataTypes)) {
    stop("Number of datasets and number of data types not equal.")
  }
  
  # options <- sprintf("-n ${num_iter} -t ${thin} -s ${seed} -c ${n_clust} > ${output}")
  options <- sprintf(" -n %s -t %s -s %s -c %s", R, thin, seed, K)

  
  cmd <- paste0(mdipp_dir, "mdipp")
  
  data_str <- ""
  for(l in 1:L){
    curr_type <- str_to_lower(dataTypes[l])
    mdippType <- NULL
    if(curr_type %in% c("gaussian", "normal", "g", "n")){
      mdippType <- "N"
    }
    if(curr_type %in% c("categorical", "multinomial", "c", "m")){
      mdippType <- "M"
    }
    if(curr_type %in% c("gaussianprocess", "gp")){
      mdippType <- "GP"
    }
    if(curr_type %in% c("bagofwords", "bw")){
      mdippType <- "BW"
    }
    if(is.null(mdippType)){
      stop(cat("Data type ", l, " is not acceptable. Please check inputs.\n"))
    }
    
    data_str <- paste(data_str, mdippType, datasets[l])
    
  }
  
  cmd <- paste0(cmd, data_str, options)
  
  if(is.null(output)){
    if(L < 5){
      
      filenames <- datasets %>% 
        strsplit("/") %>% 
        lapply(function(x) x[length(x)]) %>% 
        str_remove(".csv") %>% 
        unlist()
      
      output <- paste0(filenames, collapse = "")
    } else {
      today <- strsplit(date(), " ")
      output <- paste0(today[[1]][3], today[[1]][2], today[[1]][5])
      
    }
    
    output <- paste0(output, "R", R, "T", thin, "K", K, "S", seed, ".csv")
  }
  
  cmd <- paste0(cmd, " > ", output)

  if(call) system(cmd)

  cmd
}

# Read in arguments from the command line
args <- inputArguments()

# Take the dataset string and split
datasets <- args$datasets %>%
  strsplit(" ") %>%
  unlist()

# take data types and split
dataTypes <- args$datatypes %>%
  strsplit(" ") %>%
  unlist()

# Directory mdi lives in
mdipp_dir <- args$mdi_dir

# Parameters of MDI
R <-args$n_iter
thin <- args$thin
K <- args$n_clust
seed <- args$seed
call <- args$call

# Output file
output <- args$output

x1 <- runMDI(datasets, dataTypes, R, 
  mdipp_dir = mdipp_dir,
  K = K,
  thin = thin,
  seed = seed,
  output = output,
  call = call
)

print(x1)

