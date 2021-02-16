#!/usr/bin/Rscript

# Compare the densities of the phi parameters for the Bayesian chains, 
# individually and pooled, to Consensus clustering

# Libraries
library(tibble)
library(magrittr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mdiHelpR)
setMyTheme()
set.seed(1)

# Tibbles containing samples and metadata
bayes_cont <- readRDS("./Data/Yeast/BayesContParamsTibble.rds")
cc_cont <- readRDS("./Data/Yeast/CCContParamsTibble.rds")

# Consensus clustering ensemble used
ensemble_choice <- which(cc_cont$R == max(cc_cont$R) & cc_cont$S == max(cc_cont$S))

# Extract data and shape for ggplot2
bayes_data <- bayes_cont[bayes_cont$Use_chain, ] %>% 
  apply(1, function(x){ .dat <- x$Samples; .dat$Seed <- x$Seed; .dat}) %>% 
  do.call(rbind, .) %>% 
  pivot_longer(-Seed, names_to = "Parameter")

bayes_data$Model <- paste0("Chain ", bayes_data$Seed)
bayes_data$S <- 676
bayes_data$R <- 676000

# The data from all the kept Bayesian chains
bayes_data_avg <- bayes_data
bayes_data$Model <- "Bayesian (pooled)"

# The consensus clustering samples
cc_data <- cc_cont[ensemble_choice, ]$Samples %>% 
  as.data.frame() %>% 
  pivot_longer(everything(), names_to = "Parameter")
cc_data$Seed <- NA
cc_data$S <- cc_cont[ensemble_choice, ]$S
cc_data$R <- cc_cont[ensemble_choice, ]$R
cc_data$Model <- "Consensus clustering" # paste0("CC(", cc_data$R, ", ", cc_data$S, ")")
# cc_data$Method <- "Consensus clustering"

plt_data <- rbind(bayes_data, bayes_data_avg, cc_data)

# Facet labels for continuous variables
param_labels <- c(paste0("alpha[", 1:3, "]"), paste0("phi[", c(12, 13, 23), "]"))
names(param_labels) <- unique(plt_data$Parameter)

model_labels <- c(paste0("Chain ", unique(bayes_data$Seed)), "Bayesian (pooled)", "Consensus clustering") %>% 
  set_names(c(paste0("Chain ", unique(bayes_data$Seed)), "Bayesian (pooled)", "Consensus clustering"))

# Change the model variable to a factor for plotting in an order
plt_data$Model <- factor(plt_data$Model, levels = model_labels)

# Density plots
p_param_density <- plt_data %>%
  filter(Parameter %in% c("Phi_12", "Phi_13", "Phi_23")) %>%
  ggplot(aes(x = value, fill = Model)) +
  # geom_bar(aes(y = (..count..)/sum(..count..)), alpha = 1) +#
  geom_histogram(aes(y = ..density..), binwidth = 1) +
  facet_grid(Model ~ Parameter, 
  labeller = labeller(Parameter = as_labeller(param_labels, label_parsed), 
                      Model = model_labels
  ) #,
  # scales = "free_x"
  ) +

  labs(
    # title = "Parameter density",
    x = "Value",
    y = "Density"
  ) +
  scale_fill_viridis_d()

p_param_density +
  theme(
    axis.text.y = element_text(size = 10.5),
    axis.text.x = element_text(size = 10.5),
    axis.title.y = element_text(size = 10.5),
    axis.title.x = element_text(size = 10.5),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 10.5),
    legend.text = element_text(size = 10.5)
  )  +
  theme(legend.position="bottom")

ggsave("./SupplementaryMaterial/Images/Yeast/ComparisonDensities.png",
       height = 11,
       width = 8)
