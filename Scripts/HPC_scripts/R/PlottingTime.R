#!/bin/Rscript

# Author: Stephen Coleman
# Date: 03/06/2020
# Purpose: plot time taken for each mixture model method for simulations for
# First Year Report

# Call libraries
library(magrittr)
library(ggplot2)
library(dplyr)
library(optparse)

# Set plot theme
mdiHelpR::setMyTheme()

inputArguments <- function() {
  option_list <- list(

    # Data to cluster
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = "~/rds/hpc-work/MDI_output/Simulations/Single_dataset/Time/",
      help = "File path to model output [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-s", "--save_dir"),
      type = "character",
      default = "~/rds/hpc-work/Analysis/Simulations/Time/",
      help = "Directory to save model and predicted clustering to [default= %default]",
      metavar = "character"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("--scn"),
      type = "character",
      default = "base_case",
      help = "Current scenario to analyse [default= %default]",
      metavar = "character"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("--models"),
      type = "character",
      default = "Bayesian Consensus Frequentist",
      help = "Current scenario to analyse [default= %default]",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# Input arguments from the command line
args <- inputArguments()

# Pass these to objects within R
gen_dir <- args$dir # "/Users/stephen/Desktop/Testing_pipeline/Time/"
save_dir <- args$save_dir
scn <- args$scn # "large_n_small_p_base"
inferences <- args$models %>% # c("Bayesian", "Frequentist", "Consensus")
  strsplit(" ") %>%
  unlist()

N <- list("Bayesian" = 1000001, "Consensus" = c(1, 10, 100, 1000, 10001))
# thin <- 1000
K <- 50
sims <- 1:100
seeds <- list("Bayesian" = 1:10, "Consensus" = 1:100)

# Expand grid of the consensus inference sims and seeds
all_comb <- expand.grid(sims, seeds)

# files <- paste(main_dir, scn, "sim", all_comb[, 1], "model", inf, "N", N, "seed", all_comb[, 2], "K", K, ".txt", sep = "_")

time_df <- tibble(
  "Simulation" = as.numeric(),
  "Real" = as.numeric(),
  "User" = as.numeric(),
  "Inference" = as.character()
)

# Iterate over the types of inference
for (inf in inferences) {
  skip <- 0
  if (inf == "Frequentist") {
    # The main directory for the given scneario adn inference
    main_dir <- paste0(gen_dir, scn, "/", inf, "/")
    
    skip <- 3

    for (i in sims) {
    # If inference is Frequentist (actually MLE) treat differently as naming
    # convention broke, e.g.
    # time_scn_small_n_large_p_base_sim_52_model_Frequentist{model_name}
    files <- paste("time_scn", scn, "sim", i, "model", inf, sep = "_") %>%
      paste0(main_dir, .) %>%
      paste0("{model_name}.txt")


    # Read in the files, skipping the first 3 lines in the case of MLE as Mclust
    # printed to the output file
    for (f in files) {
      .x <- tryCatch(read.table(f, skip = skip),
        error = function(s) NULL
      )

      # Split out the relevant information into a usable data shape; originally
      # in the format:
      #
      #   V1    V2
      # 1 real  4m47.247s
      # 2 user  4m46.536s
      # 3 sys   0m0.003s
      if (!is.null(.x)) {
        time_split <- .x[, 2] %>%
          as.character() %>%
          strsplit("m")

        t <- list()

        for (j in 1:2) {
          time_split[[j]][2] <- time_split[[j]][2] %>%
            stringr::str_remove("s")

          time_split[[j]] <- time_split[[j]] %>% as.numeric()

          time_split[[j]][1] <- time_split[[j]][1] * 60

          t[[j]] <- time_split[[j]] %>%
            sum()
        }

        # Add the new entry to the full df
        # Use a tibble to avoid weird behaviour around characters that
        # data.frames default to
        new_entry <- tibble(
          "Simulation" = i,
          "Real" = t[[1]],
          "User" = t[[2]],
          "Inference" = inf,
          "R" = NA
        )

        time_df <- rbind(time_df, new_entry)
      }
    }
    } 
    } else {
    # The main directory for the given scneario adn inference
    main_dir <- paste0(gen_dir, scn, "/", inf, "/")

    # Iterate over the 100 simulations within the scenario
    for (i in sims) {
      for (n in N[[inf]]) {
        
        files <- paste(scn, "sim", i, "model", inf, "N", n, "seed", seeds[[inf]], "K", K,
          sep = "_"
        ) %>%
          paste0(main_dir, .) %>%
          paste0(".txt")

        # # If inference is Frequentist (actually MLE) treat differently as naming
        # # convention broke
        # if (inf == "Frequentist") {
        #   # time_scn_small_n_large_p_base_sim_52_model_Frequentist{model_name}
        #   files <- paste("time_scn", scn, "sim", i, "model", inf, sep = "_") %>%
        #     paste0(main_dir, .) %>%
        #     paste0("{model_name}.txt")
        # }

        # Read in the files, skipping the first 3 lines in the case of MLE as Mclust
        # printed to the output file
        for (f in files) {
          .x <- tryCatch(read.table(f, skip = skip),
            error = function(s) NULL
          )

          # Split out the relevant information into a usable data shape; originally
          # in the format:
          #
          #   V1    V2
          # 1 real  4m47.247s
          # 2 user  4m46.536s
          # 3 sys   0m0.003s
          if (!is.null(.x)) {
            time_split <- .x[, 2] %>%
              as.character() %>%
              strsplit("m")

            t <- list()

            for (j in 1:2) {
              time_split[[j]][2] <- time_split[[j]][2] %>%
                stringr::str_remove("s")

              time_split[[j]] <- time_split[[j]] %>% as.numeric()

              time_split[[j]][1] <- time_split[[j]][1] * 60

              t[[j]] <- time_split[[j]] %>%
                sum()
            }

            # Add the new entry to the full df
            # Use a tibble to avoid weird behaviour around characters that
            # data.frames default to
            new_entry <- tibble(
              "Simulation" = i,
              "Real" = t[[1]],
              "User" = t[[2]],
              "Inference" = inf,
              "R" = n
            )

            time_df <- rbind(time_df, new_entry)
          }
        }
      }
    }
  }
}

# Correct "Frequentist" to "Mclust"
time_df$Inference[time_df$Inference == "Frequentist"] <- "Mclust"

# Add scenario column
time_df$scenario <- scn

# Make a nice string for the plot title from the scenario
scn_title <- scn %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::str_to_sentence()

# Plot the data
p1 <- time_df %>%
  group_by(Inference) %>%
  ggplot(aes(x = Inference, y = log(User, base = 10))) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip() +
  labs(
    title = scn_title,
    subtitle = "User time for each chain of each method",
    y = expression(log[10](t)) # ,
    # caption = "Times for different methods across simulations. \nConsensus is running a Gibbs sampler for 10,001 iterations, Bayesian for 1,000,001 (and is a factor of 10^2 slower)."
  )

# Save the plot
save_file <- paste0(save_dir, scn, "TimeComparison.png")
ggsave(save_file, plot = p1)

# Save the data
data_file <- paste0(save_dir, scn, "TimeComparisonData.csv")
write.csv(time_df, data_file)

