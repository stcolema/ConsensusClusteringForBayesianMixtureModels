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
library(stringr)

# Set plot theme
mdiHelpR::setMyTheme()

extractTime <- function(f, skip = 0) {
  .x <- tryCatch(read.table(f, skip = skip),
    error = function(s) NULL
  )

  t <- NULL

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
  }
  t
}

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
gen_dir <- args$dir #
save_dir <- args$save_dir # "./Data/"
scn <- args$scn # "no_structure"
inferences <- args$models %>% # c("Bayesian", "Frequentist", "Consensus")
  strsplit(" ") %>%
  unlist()

gen_dir <- "C:/Users/stephen/Documents/PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Time/"
save_dir <- "./Data/"
plt_dir <- "./Images/Simulations/"
inferences <- c("Bayesian", "Frequentist", "Consensus")

N <- list("Bayesian" = 1000001, "Consensus" = c(1, 10, 100, 1000, 10001))
# thin <- 1000
K <- 50
sims <- 1:100
seeds <- list("Bayesian" = 1:10, "Consensus" = 1:100)

# Expand grid of the consensus inference sims and seeds
all_comb <- expand.grid(sims, seeds)

# files <- paste(main_dir, scn, "sim", all_comb[, 1], "model", inf, "N", N, "seed", all_comb[, 2], "K", K, ".txt", sep = "_")

time_df <- tibble(
  "Real" = as.numeric(),
  "User" = as.numeric(),
  "Simulation" = as.numeric(),
  "R" = as.numeric(),
  "Inference" = as.character(),
  "Scenario" = as.character()
)


# scenarios <- list.dirs(gen_dir, recursive = F, full.names = F)

scenarios <- c(
  "base_case",
  "irrelevant_features_10",
  "irrelevant_features_20",
  "irrelevant_features_100",
  "large_standard_deviation_3",
  "large_standard_deviation_5",
  "no_structure",
  "simple_2d",
  "small_n_large_p_base",
  "small_n_large_p_small_dm",
  "varying_proportions",
  "varying_proportions_small_dm"
)

for (scn in scenarios) { # c("no_structure", "simple_2d", "large_standard_deviation_5", "large_standard_deviation_3")){
  # Iterate over the types of inference
  for (inf in inferences) {
    skip <- 0

    # The main directory for the given scneario adn inference
    main_dir <- paste0(gen_dir, scn, "/", inf, "/")

    # Columns in file names
    col_names <- c("Simulation", "R", "Seed")

    if (inf == "Frequentist") {

      # Skip the first 3 lines of the text file (Mclust print statement)
      skip <- 3

      col_names <- "Simulation"
    }
    # for (i in sims) {
    # # If inference is MLE treat differently as naming
    # # convention broke, e.g.
    # # time_scn_small_n_large_p_base_sim_52_model_Frequentist{model_name}
    # files <- paste("time_scn", scn, "sim", i, "model", inf, sep = "_") %>%
    #   paste0(main_dir, .) %>%
    #   paste0("{model_name}.txt")
    #
    #
    # # Read in the files, skipping the first 3 lines in the case of MLE as Mclust
    # # printed to the output file
    # for (f in files) {
    #   .x <- tryCatch(read.table(f, skip = skip),
    #     error = function(s) NULL
    #   )
    #
    #   # Split out the relevant information into a usable data shape; originally
    #   # in the format:
    #   #
    #   #   V1    V2
    #   # 1 real  4m47.247s
    #   # 2 user  4m46.536s
    #   # 3 sys   0m0.003s
    #   if (!is.null(.x)) {
    #     time_split <- .x[, 2] %>%
    #       as.character() %>%
    #       strsplit("m")
    #
    #     t <- list()
    #
    #     for (j in 1:2) {
    #       time_split[[j]][2] <- time_split[[j]][2] %>%
    #         stringr::str_remove("s")
    #
    #       time_split[[j]] <- time_split[[j]] %>% as.numeric()
    #
    #       time_split[[j]][1] <- time_split[[j]][1] * 60
    #
    #       t[[j]] <- time_split[[j]] %>%
    #         sum()
    #     }
    #
    #     # Add the new entry to the full df
    #     # Use a tibble to avoid weird behaviour around characters that
    #     # data.frames default to
    #     new_entry <- tibble(
    #       "Simulation" = i,
    #       "Real" = t[[1]],
    #       "User" = t[[2]],
    #       "Inference" = inf,
    #       "R" = NA,
    #       "Scenario" = scn
    #     )
    #
    #     time_df <- rbind(time_df, new_entry)
    #   }
    # }
    # }
    # }
    # else {

    files <- list.files(main_dir, pattern = ".txt$")
    if (length(files) > 0) {
      file_details <- files %>%
        str_replace_all(scn, "") %>%
        str_extract_all("_(\\d+)")


      n_col <- length(file_details[[1]])
      if (n_col == 4) col_names <- c("Simulation", "R", "Seed", "K")

      file_details <- file_details %>%
        unlist() %>%
        str_replace_all("_", "") %>%
        as.numeric() %>%
        matrix(ncol = n_col, byrow = T) %>%
        set_colnames(col_names)

      if (inf == "Frequentist") file_details$R <- NA

      files <- paste0(main_dir, files)
      n_files <- length(files)

      time_mat <- files %>%
        lapply(extractTime) %>%
        lapply(unlist)

      missing <- time_mat %>%
        lapply(function(x) length(x) == 0) %>%
        unlist()

      if (!all(missing)) {
        time_mat <- time_mat %>%
          unlist() %>%
          matrix(ncol = 2, byrow = T) %>%
          set_colnames(c("Real", "User"))

        if (sum(!missing) == 1) {
          .curr_time <- c(time_mat, file_details[!missing, 1:2]) %>%
            t() %>%
            as_tibble() %>%
            set_colnames(c("Real", "User", "Simulation", "R"))
        } else {
          .curr_time <- cbind(time_mat, file_details[!missing, 1:2]) %>%
            as_tibble()
        }

        .curr_time$Inference <- inf
        .curr_time$Scenario <- scn

        time_df <- rbind(time_df, .curr_time)
      }
    }
    # # Iterate over the 100 simulations within the scenario
    # for (i in sims) {
    #   for (n in N[[inf]]) {
    #
    #     files <- paste("time_scn", scn, "sim", i, "model", inf, "N", n, "seed", seeds[[inf]], # "K", K,
    #       sep = "_"
    #     ) %>%
    #       paste0(main_dir, .) %>%
    #       paste0(".txt")
    #
    #     # # If inference is Frequentist (actually MLE) treat differently as naming
    #     # # convention broke
    #     # if (inf == "Frequentist") {
    #     #   # time_scn_small_n_large_p_base_sim_52_model_Frequentist{model_name}
    #     #   files <- paste("time_scn", scn, "sim", i, "model", inf, sep = "_") %>%
    #     #     paste0(main_dir, .) %>%
    #     #     paste0("{model_name}.txt")
    #     # }
    #
    #     # Read in the files, skipping the first 3 lines in the case of MLE as Mclust
    #     # printed to the output file
    #     for (f in files) {
    #       .x <- tryCatch(read.table(f, skip = skip),
    #         error = function(s) NULL
    #       )
    #
    #       # Split out the relevant information into a usable data shape; originally
    #       # in the format:
    #       #
    #       #   V1    V2
    #       # 1 real  4m47.247s
    #       # 2 user  4m46.536s
    #       # 3 sys   0m0.003s
    #       if (!is.null(.x)) {
    #         time_split <- .x[, 2] %>%
    #           as.character() %>%
    #           strsplit("m")
    #
    #         t <- list()
    #
    #         for (j in 1:2) {
    #           time_split[[j]][2] <- time_split[[j]][2] %>%
    #             stringr::str_remove("s")
    #
    #           time_split[[j]] <- time_split[[j]] %>% as.numeric()
    #
    #           time_split[[j]][1] <- time_split[[j]][1] * 60
    #
    #           t[[j]] <- time_split[[j]] %>%
    #             sum()
    #         }
    #
    #         # Add the new entry to the full df
    #         # Use a tibble to avoid weird behaviour around characters that
    #         # data.frames default to
    #         new_entry <- tibble(
    #           "Simulation" = i,
    #           "Real" = t[[1]],
    #           "User" = t[[2]],
    #           "Inference" = inf,
    #           "R" = n,
    #           "Scenario" = scn
    #         )
    #
    #         time_df <- rbind(time_df, new_entry)
    #       }
    #     }
    #   }
  }
}
# }
# }

# Correct "Frequentist" to "Mclust"
time_df$Inference[time_df$Inference == "Frequentist"] <- "Mclust"

# Add scenario column
# time_df$scenario <- scn

# # Make a nice string for the plot title from the scenario
# scn_title <- scn %>%
#   stringr::str_replace_all("_", " ") %>%
#   stringr::str_to_sentence()

# ?geom_smooth()

mdiHelpR::setMyTheme(axis.text.y=element_text(hjust=0.0, size = 10.5),
           axis.text.x=element_text(angle=30, size = 10.5),
           # axis.title.y=element_blank(),
           # axis.title.x=element_blank(),
           plot.title = element_text(size = 18, face = "bold"),
           plot.subtitle = element_text(size = 14),
           strip.text.x = element_text(size = 10.5)
)

scenarios <- c(
  "base_case",
  "irrelevant_features_10",
  "irrelevant_features_20",
  "irrelevant_features_100",
  "large_standard_deviation_3",
  "large_standard_deviation_5",
  "no_structure",
  "simple_2d",
  "small_n_large_p_base",
  "small_n_large_p_small_dm",
  "varying_proportions",
  "varying_proportions_small_dm"
)

dimensions <- c(20, 30, 40, 120, 20, 20, 2, 2, 500, 500, 20, 20)
n_sim <- length(scenarios)

time_df$Dimension <- 1

for(i in 1:n_sim){
  scn <- scenarios[i]
  time_df$Dimension[time_df$Scenario == scn] <- dimensions[i]
}

time_df$Dimension <- as.factor(time_df$Dimension)

# Plot the data
p1 <- time_df %>%
  filter(Inference != "Mclust", R > 1) %>%
  group_by(R) %>%
  ggplot(aes(x = log(R, base = 10), y = log(User, base = 10), group = Dimension, colour = Dimension)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, mapping = aes(color = Dimension)) +
  # geom_boxplot(colour = "black", fill = "#FDE725FF") +
  # coord_flip() +
  labs(
    title = "Time taken for MCMC iterations",
    # subtitle = "User time for each chain",
    x = expression(log[10](R)),
    y = expression(log[10](t)),
    colour = "P"
    # caption = "Times for different methods across simulations. \nConsensus is running a Gibbs sampler for 10,001 iterations, Bayesian for 1,000,001 (and is a factor of 10^2 slower)."
  ) +
  scale_color_viridis_d() +
  theme(axis.text.y=element_text(size = 10.5),
        axis.text.x=element_text(size = 10.5),
        axis.title.y=element_text(size = 10.5),
        axis.title.x=element_text(size = 10.5),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 10.5),
        legend.text = element_text(size = 10.5)
  )
# + theme(legend.position = "bottom")

p1

# Save the plot
save_file <- paste0(plt_dir, "TimeComparison.png")
ggsave(save_file, plot = p1, height = 5.5, width = 5)

# Save the data
data_file <- paste0(save_dir, "TimeComparisonData.csv")
write.csv(time_df, data_file)
