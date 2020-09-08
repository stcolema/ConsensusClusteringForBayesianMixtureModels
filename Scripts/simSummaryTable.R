
library(tidyverse)
library(mdiHelpR)
library(magrittr)

setMyTheme()

scenarios <- list.dirs("./Data", recursive = F, full.names = F)
data_files <- list.files("./Data", recursive = T, full.names = T)
n_files <- length(data_files)

analysis_data <- summary_data <- vector("list", n_files) %>%
  set_names(scenarios)

models_of_interest <- c(
  "Bayesian (Pooled)",
  "Maximum likelihood (Mclust)",
  "Consensus (100, 10)",
  "Consensus (100, 50)",
  "Consensus (100, 100)",
  "Consensus (10000, 10)",
  "Consensus (10, 10)"
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
  
  ari_table[,i] <- curr_data$Mean_ari[match(curr_data$Model, models_of_interest)] %>% 
    round(digits = 3)
  
  sd_table[, i] <- curr_data$SD_ari[match(curr_data$Model, models_of_interest)] %>% 
    round(digits = 3)
}

colnames(ari_table) <- c("Irrelevant features 100", "Simple 2D", "Small N, large P")
colnames(sd_table) <- c("Irrelevant features 100", "Simple 2D", "Small N, large P")
ari_table[, c(2, 3, 1)]
sd_table[, c(2, 3, 1)]


summary_data
