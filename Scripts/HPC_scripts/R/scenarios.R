#!/usr/bin/env/Rscript

library(magrittr)
library(stringr)
library(mdiHelpR)

# The directory to save datasets
main_dir <- "~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/"

# The generic dataset and cluster label names
gen_data_filename <- "dataset_"
gen_cluster_IDs_filename <- "cluster_IDs_"

# The scenarios of interest
scn_table <- data.frame(
Scenario = c("Simple 2D", "No structure", "Base Case", rep("Large N, small P", 3), rep("Large standard deviation", 3), rep("Irrelevant features", 5), rep("Small N, large P", 2), rep("Varying proportions", 2), rep("Large N, large P", 2)),
  N = c(100, 100, 2e2, 1e4, 1e4, 1e4, 2e2, 2e2, 2e2, 2e2, 2e2, 2e2, 2e2, 2e2, 50, 50, 200, 200, 3000, 3000),
  P_s = c(2, 0, 20, 4, 4, 4, 20, 20, 20, 20, 20, 20, 20, 20, 500, 500, 20, 20, 3000, 3000),
  P_n = c(0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 10, 20, 100, 200, 0, 0, 0, 0, 0 ,0),
  K = c(5, 1, 5, 5, 50, 50, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 8, 8),
  Delta_mu = c(3, 0, 1, 1, 1, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.2, 1, 0.4, 1, 0.4),
  sigma2 = c(1, 1, 1, 1, 1, 1, 3, 5, 10, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  Pi = c(rep("vec(1/K)", 16), "(0.5, 0.25, 0.125, 0.0675, 0.0675)", "(0.5, 0.25, 0.125, 0.0675, 0.0675)", "vec(1/K)", "vec(1/K)")
  )

scn_table <- scn_table[c(1:3, 7, 8, 11:13, 17:18, 4:6, 15:16, 19:20), ]

scenario_names <- scn_table$Scenario %>%
  str_replace_all(",", "") %>%
  str_replace_all(" ", "_") %>%
  str_to_lower()

scenario_names[c(4,5)] %<>%
  paste0("_", c(3, 5))

scenario_names[6:8] %<>%
  paste0("_", c(10, 20, 100))

scenario_names[11:13] %<>%
  paste0("_", c("base", "large_k", "large_k_small_dm"))

scenario_names[14:15] %<>%
  paste0("_", c("base", "small_dm"))

scenario_names[10] %<>%
  paste0("_small_dm")

scenario_names[17] %<>%
  paste0("_small_dm")

# Change the scenario names to a nice saveable format
scn_table$Scenario <- scenario_names

scn_table <- scn_table[c(10, 16,17),]


# Write the scenario data.frame to a csv for future reference
write.csv(scn_table, "./my_scenarios.csv", row.names = F)

# The number of scenarios
n_scn <- nrow(scn_table)

# The number of simualtions within each scenario
n_sim <- 100

for(l in 1:n_sim){
  # Set a seed to define the lth dataset
  set.seed(l)

  # loop over the different scenarios
  for(i in 1:n_scn){

     # Create the filenames for the output
     save_dir <- paste0(main_dir, scn_table[i, 1])
     data_filename <- paste0(save_dir, "/", gen_data_filename, l, ".csv")
     cluster_IDs_filename <- paste0(save_dir, "/", gen_cluster_IDs_filename, l, ".Rds")

     # Create the save directory if it does not already exist
     dir.create(save_dir, showWarnings = FALSE)

     # Extract the current scenario parameters
     N <- scn_table[i, 2]
     P <- scn_table[i, 3]
     P_n <- scn_table[i, 4]
     K <- scn_table[i, 5]
     delta_mu <- scn_table[i, 6]
     sigma <- scn_table[i, 7]

     # Pi is a little awkward
     pi <- scn_table[i, 8] %>%
       as.character()

     if(pi == "vec(1/K)"){
       pi <- rep(1/K, K)
     } else {
       pi <- pi %>%
         str_replace_all("[()]", "") %>%
         strsplit(", ")
       pi <- pi[[1]] %>%
         as.numeric()
     }

     # Generate data
     d1 <- generateSimulationDataset(
       K = K,
       n = N,
       p = P,
       p_n = P_n,
       pi = pi,
       delta_mu = delta_mu,
       cluster_sd = sigma
     )

     # Create named vector of labels
     labels <- d1$cluster_IDs %>%
       set_names(row.names(d1$data))

     # Save labels
     saveRDS(labels, file = cluster_IDs_filename)

     # Save data
     write.csv(scale(d1$data), file = data_filename)
  
  }

}
