

library(tidyverse)
library(mdiHelpR)
library(magrittr)
library(patchwork)

setMyTheme()

scenarios <- list.dirs("./Data", recursive = F, full.names = F)
data_files <- list.files("./Data", recursive = T, full.names = T)
n_files <- length(data_files)

save_dir <- "./Images/Simulations/"

analysis_data <- summary_data <- vector("list", n_files) %>%
  set_names(scenarios)

models_of_interest <- c(
  "Bayesian (Pooled)",
  "Maximum likelihood (Mclust)",
  "Consensus (100, 50)",
  "Consensus (10000, 100)",
  "Consensus (10, 10)"
)

scenarios_pretty <- c("Irrelevant features 100", "Simple 2D", "Small N, large P")


my_table <- matrix(0, nrow = length(models_of_interest), ncol = n_files) %>% 
  set_rownames(models_of_interest) %>% 
  set_colnames(scenarios)

for (i in 1:n_files) {
  analysis_data[[i]] <- read.csv(data_files[i], row.names = 1)
  analysis_data[[i]]$Scenario <- scenarios_pretty[i]
}

my_df <- do.call("rbind", analysis_data)

my_df$Scenario <- factor(my_df$Scenario, 
                         levels = scenarios_pretty[c(2, 3, 1)]
                         )

model_labels <- my_df$Model %>% unique() %>% str_sort(numeric = T)

my_df$Model <- factor(my_df$Model, levels = model_labels[c(2:26, 1, 27)])

p1 <- my_df %>% 
  ggplot(aes(x = Model, y = ARI)) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip() +
  facet_wrap(~Scenario) +
  theme(axis.text.y=element_text(hjust=0.05)) +
  labs(title = "Predictive performance")

p1

ggsave(paste0(save_dir, "simulation_model_prediction.png"), width = 10, height = 10, plot = p1)


p2 <- my_df %>% 
  filter(Model != "Maximum likelihood (Mclust)") %>%  
  ggplot(aes(x = Model, y = Frobenius_norm)) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip() + 
  facet_wrap(~Scenario, scales = "free_x") +
  theme(axis.text.y=element_text(hjust=0.05)) +
  labs(title = "Uncertainty quantification",
       y = "Frobenius norm")


p2

ggsave(paste0(save_dir, "simulation_model_uncertainty.png"), width = 10, height = 10, plot = p2)

p3 <- p1 / p2 +
  plot_annotation(
    title = "Comparison of model performance across simulations"
  )

ggsave(paste0(save_dir, "simulation_model_performance.png"), width = 10, height = 10, plot = p3)


