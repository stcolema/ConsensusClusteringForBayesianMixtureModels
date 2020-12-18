#!/bin/Rscript

## Compare the performance of the three different models for a  given scenario


library(magrittr)
library(ggplot2)
library(Matrix)
library(dplyr)
library(mdiHelpR)
library(optparse)


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

set.seed(1)
setMyTheme()
args <- inputArguments()

scn <- args$scn # "base_case"
models <- args$models %>%
  strsplit(" ") %>%
  unlist()
home_dir <- args$dir
main_dir <- paste0(home_dir, scn)
result_dir <- paste(main_dir, models, sep = "/")
save_dir <- args$save_dir
saving <- F

# I was thinking of comparing clusterigns to the truth, but 
# doesn't really fit into this script right not
#truth_dir <- args$truth_dir
#sim_used <- args$sim_example
#orig_data_f <- paste0(truth_ir, "dataset_", sim_used, ".csv")
#cluster_id_f <- paste0(truth_ir, "cluster_IDs_", sim_used, ".Rds")
#
#orig_data <- read.csv(orig_data_f, row.names= 1)
#truth <- readRDS(cluster_id_f)

scn_str <- scn %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::str_to_sentence()

gen_save_str <- paste0(save_dir, scn)


models <- c("Bayesian", "Consensus", "Frequentist")
n_mod <- length(models)
# main_dir <- paste0("/Users/stephen/Desktop/Work_update_2020_05_11/Model_performance/Model_performance/", scn, "/")
#result_dir <- paste0(main_dir, models)



files <- list()
my_data <- list()
for (i in 1:n_mod) {
  my_data[[i]] <- list()
  d <- result_dir[i]
  files[[i]] <- list.files(d, pattern = "ResultsDF.csv", full.names = T)
  n_files <- length(files[[i]])
  for (j in 1:n_files) {
    f <- files[[i]][j]
    my_data[[i]][[j]] <- read.csv(f, row.names = 1)
    my_data[[i]][[j]]$Inference <- models[i]
  }
  my_data[[i]] <- do.call(rbind, my_data[[i]])
}

bayes_pooled <- my_data[[1]][my_data[[1]]$Seed == "Pooled", ]
consensus_data <- my_data[[2]]
freq_data <- my_data[[3]]


iter_label <- consensus_data$N_iter %>%
  unique() %>%
  paste("Iteration sampled:", .) %>%
  set_names(unique(consensus_data$N_iter))

seed_label <- consensus_data$N_seeds %>%
  unique() %>%
  paste("Total number of samples:", .) %>%
  set_names(unique(consensus_data$N_seeds))

p1 <- consensus_data %>%
  ggplot(aes(x = N_iter, y = ARI, group = Simulation)) +
  facet_wrap(~N_seeds, labeller = labeller(N_seeds = seed_label)) +
  geom_line(colour = "black", alpha = 0.3) +
  labs(
    x = "Iteration used from chain",
    title = paste0(scn_str, ": Model predictive performance across simulations")
  ) +
  ylim(0, 1)


save_file <- paste0(gen_save_str, "_ci_model_performance_iterations.png")
ggsave(save_file, plot = p1)

# if (saving) ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_ci_model_performance_iteration.png") # , width = 5, height = 3.5)


p2 <- consensus_data %>%
  ggplot(aes(x = N_seeds, y = ARI, group = Simulation)) +
  facet_wrap(~N_iter, labeller = labeller(N_iter = iter_label)) +
  geom_line(colour = "black", alpha = 0.3) +
  labs(
    x = "Number of samples used",
    title = paste0(scn_str, ": Model predictive performance across simulations") # ,
    # subtitle = "Interesting behaviour for using the first sample from each chain, otherwise very robust performance."
  ) +
  ylim(0, 1)

save_file <- paste0(gen_save_str, "_ci_model_performance_samples.png")
ggsave(save_file, plot = p2)

# if(saving) ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_ci_model_performance_samples.png") #, width = 5, height = 3.5)

p3 <- consensus_data %>%
  ggplot(aes(x = N_iter, y = Frobenius_norm, group = Simulation)) +
  facet_wrap(~N_seeds, labeller = labeller(N_seeds = seed_label)) +
  geom_line(colour = "black", alpha = 0.3) +
  labs(
    x = "Iteration used from chain",
    y = "Frobenius norm",
    title = paste0(scn_str, ": Model uncertainty across simulations")
  )

save_file <- paste0(gen_save_str, "_ci_model_uncertainty_iterations.png")
ggsave(save_file, plot = p3)

p4 <- consensus_data %>%
  ggplot(aes(x = N_seeds, y = Frobenius_norm, group = Simulation)) +
  facet_wrap(~N_iter, labeller = labeller(N_iter = iter_label)) +
  geom_line(colour = "black", alpha = 0.3) +
  labs(
    x = "Number of samples used",
    y = "Frobenius norm",
    title = paste0(scn_str, ": Model uncertainty across simulations") # ,
    # subtitle = "Interesting behaviour for using the first sample from each chain, otherwise very robust performance."
  )

save_file <- paste0(gen_save_str, "_ci_model_uncertainty_samples.png")

ggsave(save_file, plot = p4)

# consensus_data %>% glimpse()
# bayes_pooled %>% glimpse()

bayes_pooled$N_seeds <- NA


bayes_pooled <- bayes_pooled %>% dplyr::select(-Seed)

consensus_plt <- consensus_data %>%
  dplyr::filter(N_iter %in% c(1, 10, 100, 1000, 10000))

all_df <- rbind(bayes_pooled, consensus_plt, freq_data)


consensus_ind <- which(all_df$Inference == "Consensus")
models <- all_df$Inference

models[all_df$Inference == "Frequentist"] <- "Maximum likelihood (Mclust)"
models[all_df$Inference == "Bayesian"] <- "Bayesian (Pooled)"



models[consensus_ind] <- paste0(
  all_df$Inference[consensus_ind],
  " (",
  all_df$N_iter[consensus_ind],
  ", ",
  all_df$N_seeds[consensus_ind],
  ")"
)

model_order <- stringr::str_sort(unique(models), numeric = T)
n_model <- length(model_order)
model_levels <- model_order[c(2:(n_model - 1), 1, n_model)]

all_df$Model <- factor(models, # paste0(all_df$Inference, "N", all_df$N_iter, "S", all_df$N_seeds),
  levels = model_levels
)

p5 <- all_df %>%
  ggplot(aes(x = Model, y = ARI)) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip() +
  labs(
    title = paste0(scn_str, ": Predictive performance")
  ) +
  ylim(0, 1) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0))

save_file <- paste0(gen_save_str, "_all_model_performance.png")

# if(saving) ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_model_performance.png") #, width = 5, height = 3.5)

ggsave(save_file, plot = p5)

p6 <- all_df %>%
  ggplot(aes(x = Model, y = Frobenius_norm)) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip() +
  labs(
    title = paste0(scn_str, ": Uncertainty quantification") # ,
    # subtitle = "N is the length of chain, S is the number of chains combined in Consensus Inference"
  )


save_file <- paste0(gen_save_str, "_all_model_uncertainty.png")
ggsave(save_file, plot = p6)

# if(saving) ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_model_uncertainty.png") #, width = 5, height = 3.5)

p7 <- all_df %>%
  ggplot(aes(x = Model, y = K_diff)) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip()  +
  labs(
    title = paste0(scn_str, ": predicted number of clusters")
  ) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0))


save_file <- paste0(gen_save_str, "_all_predicted_K.png")
ggsave(save_file, plot = p7)

write.csv(all_df, paste0(save_dir, "/all_results.csv"))

# 
# m1 <- readMM(paste0(result_dir[1], "/Simulation1PSMN1e+06S10.txt")) %>%
#   as.matrix()
# 
# m2 <- readMM(paste0(result_dir[1], "/Simulation2PSMN1e+06S9.txt")) %>%
#   as.matrix()
# 
# t1 <- readRDS("/Users/stephen/Desktop/Work_update_2020_05_11/Input_data/base_case/cluster_IDs_1.Rds")
# t2 <- readRDS("/Users/stephen/Desktop/Work_update_2020_05_11/Input_data/base_case/cluster_IDs_2.Rds")
# 
# d1 <- read.csv("/Users/stephen/Desktop/Work_update_2020_05_11/Input_data/base_case/dataset_1.csv",
#   row.names = 1
# )
# d2 <- read.csv("/Users/stephen/Desktop/Work_update_2020_05_11/Input_data/base_case/dataset_2.csv",
#   row.names = 1
# )
# 
# mdiHelpR::compareSimilarityMatrices(m1, m2, n_row = 1)
# 
# seeds <- c(1, 10, 100)
# iter <- c(1, 10, 100)
# 
# my_cm <- list()
# n_seeds <- length(seeds)
# n_iter <- length(iter)
# for (i in 1:n_seeds) {
#   for (j in 1:n_iter) {
#     cm <- readMM(paste0(result_dir[2], "/Simulation2ConsensusMatrixN", iter[j], "S", seeds[i], ".txt")) %>%
#       as.matrix() * 1
# 
#     my_cm[[(i - 1) * n_seeds + j]] <- cm
#   }
# }
# 
# compareSimilarityMatrices(
#   matrices = my_cm,
#   title = paste0(scn_str, ": annotated consensus matrices"),
#   matrix_imposing_order = 9
# )
# # if(saving) ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_sim2_consensusmatrices.png") #, width = 5, height = 3.5)
# #
# compareSimilarityMatricesAnnotated(
#   matrices = my_cm,
#   cluster_IDs = t2,
#   title = paste0(scn_str, ": annotated consensus matrices"),
#   matrix_imposing_order = 9,
#   fill = "#D3D3D3"
# )
#
# if(saving) ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_sim2_consensusmatrices_annotated.png") #, width = 5, height = 3.5)
#
#
#
# psm1 <- readMM(paste0(result_dir[1], "/Simulation2PSMN1e+06S5.txt")) %>%
#   as.matrix() * 1
#
# # compareSimilarityMatrices(my_cm[[9]], psm1, title = "Comparison of consensus matrix and PSM")
# if(saving) ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_cm_vs_psm.png") #, width = 5, height = 3.5)
#
# cm1 <- readMM(paste0(result_dir[2], "/Simulation1ConsensusMatrixN1000S30.txt")) %>%
#   as.matrix()
#
# row.names(m1) <- colnames(m1) <- names(t1)
# row.names(m2) <- colnames(m2) <- names(t2)
# row.names(my_cm[[9]]) <- colnames(my_cm[[9]]) <- names(t2)
#
# sim_col <- simColPal()
# breaks <- defineBreaks(sim_col, lb = 0)
# annotatedHeatmap(m1, t1, col_pal = sim_col, my_breaks = breaks, main = "Sim 1")
# annotatedHeatmap(m2, t2, col_pal = sim_col, my_breaks = breaks, main = "Sim 2") #,
#                  # filename = "baseCase_sim2_psm.png")
#
# annotatedHeatmap(my_cm[[9]], t2, col_pal = sim_col, my_breaks = breaks)
# compareSimilarityMatrices(my_cm[[9]], m2)
#
#
# compareSimilarityMatricesAnnotated(my_cm[[9]], m2, cluster_IDs = t2)
#
# annotatedHeatmap(d1, t1, main = "Sim 2")
# annotatedHeatmap(d2, t2, main = "Sim 2")
#
# setMyTheme()
#
# pcaSeriesPlot(prcomp(d1)$x, t1, n_comp =6) +
#   facet_wrap(~Cluster, labeller = labeller(cluster_lables))
#
# pcaSeriesPlot(prcomp(d1)$x, t1, n_comp =5) +
#   gghighlight::gghighlight(Cluster %in% c(1,2))
#
# pcaSeriesPlot(prcomp(d2)$x, t2, n_comp =5) +
#   gghighlight::gghighlight(Cluster %in% c(1,2)) +
#   labs(
#     title = "Base case: Comparison of cluster overlap in simulation 2",
#     subtitle = "These clusters overlap quite significantly"
#   )
#
# ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_sim2_clusterOverlap.png", width = 5, height = 3.5)
#
#
# pcaSeriesPlot(prcomp(d2)$x, t2, n_comp =5) +
#   gghighlight::gghighlight(Cluster %in% c(2, 3))
#
# pcaSeriesPlot(prcomp(d2)$x, t2, n_comp =5) +
#   gghighlight::gghighlight(Cluster %in% c(3, 4))
#
#
# simple2D <- generateSimulationDataset(5, 100, 2, delta_mu = 3)
#
# annotatedHeatmap(simple2D$data, simple2D$cluster_IDs, main = "Two dimensional data")
#
# simple2D$data %>%
#   as.data.frame() %>%
#   dplyr::mutate(Cluster = factor(simple2D$cluster_IDs)) %>%
#   ggplot(aes(x = Gene_1, y = Gene_2, colour = Cluster)) +
#   geom_point() +
#   scale_color_viridis_d() +
#   labs(
#     title = "Example: 2D",
#     subtitle = "Notice that some clusters overlap",
#     x = "Gene 1",
#     y = "Gene 2"
#   )
#
# ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/simple2d_data.png", width = 5, height = 3.5)
#
#
# noStruct <- generateSimulationDataset(5, 100, 0, p_n= 2, delta_mu = 1)
#
# annotatedHeatmap(noStruct$data, noStruct$cluster_IDs, main = "No structure")
#
# noStruct$data %>%
#   as.data.frame() %>%
#   dplyr::mutate(Cluster = factor(noStruct$cluster_IDs)) %>%
#   ggplot(aes(x = Noise_1, y = Noise_2, colour = Cluster)) +
#   geom_point() +
#   scale_color_viridis_d() +
#   labs(
#     title = "Example: No structure",
#     x = "Gene 1",
#     y = "Gene 2"
#   )
#
# ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/noStruct_data.png", width = 5, height = 3.5)
#
#
# baseCase <- generateSimulationDataset(5, 100, 20, delta_mu = 1)
#
# annotatedHeatmap(baseCase$data, baseCase$cluster_IDs,
#                  main = "Base case",
#                  filename = "/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_ph.png",
#                  show_colnames = F,
#                  show_rownames = F)
#
# cluster_lables <- paste("Cluster", 1:5) %>% set_names(1:5)
#
# pcaSeriesPlot(prcomp(baseCase$data)$x, baseCase$cluster_IDs, n_comp=6) +
#   facet_wrap(~Cluster, labeller = labeller(Cluster = cluster_lables)) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 30, hjust = 1)) +
#   labs(title = "PCA series for base case data")
#
# ggsave("/Users/stephen/Desktop/Work_update_2020_05_11/Presentation/baseCase_data.png", width = 5, height = 3.5)

