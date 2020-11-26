
library(tidyverse)
library(mdiHelpR)
library(magrittr)

setMyTheme()

scenarios <- list.dirs("./Data/Simulations/", recursive = F, full.names = F)
data_files <- list.files("./Data/Simulations/", recursive = T, full.names = T, pattern = "*all_results.csv")
n_files <- length(data_files)

analysis_data <- summary_data <- vector("list", n_files) %>%
  set_names(scenarios)

models_of_interest <- c(
  "Bayesian (Pooled)",
  "Maximum likelihood (Mclust)",
  "Consensus (10, 10)",
  "Consensus (10, 50)",
  "Consensus (10, 100)",
  "Consensus (100, 10)",
  "Consensus (100, 50)",
  "Consensus (100, 100)",
  "Consensus (10000, 10)",
  "Consensus (10000, 50)",
  "Consensus (10000, 100)"
)


ari_table <- sd_table <- matrix(0, nrow = length(models_of_interest), ncol = n_files) %>% 
  set_rownames(models_of_interest) %>% 
  set_colnames(scenarios)

for (i in 1:n_files) {
  analysis_data[[i]] <- read.csv(data_files[i], row.names = 1)

  summary_data[[i]] <- curr_data <- analysis_data[[i]] %>%
    filter(Model %in% models_of_interest) %>%
    group_by(Model) %>%
    summarise(Mean_ari = mean(ARI),
              Median_ari = median(ARI),
              SD_ari = sd(ARI))
  
  ari_table[,i] <- curr_data$Median_ari[match(models_of_interest, curr_data$Model)] %>% 
    round(digits = 3)
  
  sd_table[, i] <- curr_data$SD_ari[match(curr_data$Model, models_of_interest)] %>% 
    round(digits = 3)
}

my_table <- ari_table[, c(3, 8, 9)]

colnames(my_table) <- c("Irrelevant features 100", "Simple 2D", "Small N, large P")

library(stargazer)
stargazer(my_table)

# colnames(sd_table) <- c("Irrelevant features 100", "Simple 2D", "Small N, large P")
my_table[, c(2, 3, 1)]
sd_table[, c(2, 3, 1)]


summary_data
