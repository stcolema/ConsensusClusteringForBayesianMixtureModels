#!/usr/bin/Rscript

# Elbow plots for ensemble choice

# For similarity matrix
library(mdiHelpR)

# Pipe
library(magrittr)

# Data manipulation
library(dplyr)

# Tibbles can hold list columns
library(tibble)

# Pivot
library(tidyr)

# String manipulation
library(stringr)

# Visualisation
library(ggplot2)

set.seed(1)
setMyTheme()

# Read in CC data
alloc_data <- readRDS("./Data/Yeast/CCAllocTibble.rds")

# The datasets
datasets <- alloc_data$Dataset %>% 
  unique()

# Iterate over each dataset finding the sequential change in the consenus 
# matrices for increasing chain length
for (dataset in datasets) {
  dataset_ind <- which(alloc_data$Dataset == dataset)
  .alloc <- alloc_data[dataset_ind, ]
  s_used <- .alloc$S %>%
    unique()
  
  r_used <- .alloc$R %>%
    unique()
  
  n_r <- length(r_used)
  
  for (s in s_used) {
    .cms <- .alloc$CM[.alloc$S == s]
    
    
    for(i in 2:n_r){
      .score <- mean(sqrt((.cms[[i]] - .cms[[i-1]])**2))
      
      if(i == 2){
        scores <- .score
      } else {
        scores <- c(scores, .score)
      }
    }

    .df <- data.frame(S = s, R = r_used[-1], Score = scores, Dataset = dataset)
    if (s == s_used[1]) {
      score_df <- .df
    } else {
      score_df <- rbind(score_df, .df)
    }
  }
  if (dataset == datasets[1]) {
    plt_df <- score_df
  } else {
    plt_df <- rbind(plt_df, score_df)
  }
}

# Plot!
ensemble_elbow_plot <- plt_df %>%
  filter(S != 1) %>% #, R > 101) %>%
  ggplot(aes(x = R, y = Score, group = S, colour = as.factor(S))) +
  geom_line() +
  geom_point() +
  labs(
    title = "Ensemble choice", # "Yeast data",
    # subtitle = "Consensus clustering ensemble choice",
    x = "Chain depth",
    y = "Mean absolute difference of similarity",
    colour = "Number of chains"
  ) +
  facet_wrap(~Dataset) +
  scale_color_viridis_d() +
  theme(
    axis.text.y = element_text(hjust = 0.0, size = 10),
    axis.text.x = element_text(angle = 0, size = 10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
  ) +
  scale_x_continuous(breaks = c(0, 5000, 10000))  +
  theme(
    legend.position = c(0.85, 0.8),
    legend.direction = "vertical"
    )

ggsave("./SupplementaryMaterial/Images/Yeast/EnsembleChoicePlotAlt.png",
       plot = ensemble_elbow_plot,
       height = 6,
       width = 6
)
