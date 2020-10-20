
library(magrittr)
library(ggplot2)
library(stringr)

mdiHelpR::setMyTheme()


data_dir <- "./Data/Simulations/"
save_dir <- "./SupplementaryMaterial/Images/Simulations/"

scenarios <- list.dirs(data_dir, recursive = F, full.names = F)

scenarios_pretty <- scenarios %>%
  str_replace_all("_", " ") %>%
  str_to_sentence()

scenarios_pretty[scenarios_pretty == "Simple 2d"] <- "2D"

# scenarios_pretty[scenarios_pretty == "Small n large p small dm"] <- bquote("Small" ~ n ~ "large" ~ p ~ "small" ~ Delta * mu)
# bquote("Eq 1:" ~ y[i] == alpha + beta * x[i] + epsilon[i] ~ "or" ) #~ .(cor2))
# bquote("Small" ~ n ~ "large" ~ p ~ "small" ~ Delta * mu)


files <- paste0(data_dir, scenarios, "/all_results.csv")

n_scn <- length(scenarios)

my_data <- list()
for (i in 1:n_scn) {
  my_data[[i]] <- read.csv(files[i], row.names = 1)
  my_data[[i]]$Scenario <- scenarios_pretty[i]
}

my_df <- do.call("rbind", my_data)

my_df$Scenario <- factor(my_df$Scenario,
  levels = scenarios_pretty[c(8, 7, 1, 11:12, 5, 6, 9, 10, 2, 4, 3)]
)


my_df$Model <- my_df$Model %>% stringr::str_replace_all("Consensus ", "CC")
my_df$Model[my_df$Model == "Maximum likelihood (Mclust)"] <- "Mclust"

model_labels <- my_df$Model %>%
  unique() %>%
  str_sort(numeric = T)
my_df$Model <- factor(my_df$Model, levels = model_labels[c(2:26, 1, 27)])

simple_2d_name <- expression(paste("2D"))
base_case_name <- expression(paste("Base case"))
no_structure_name <- expression(paste("No structure"))
irr10name <- expression(paste("Irrelevant features (", P[n], " = ", 10, ")"))
irr20name <- expression(paste("Irrelevant features (", P[n], " = ", 20, ")"))
irr100name <- expression(paste("Irrelevant features (", P[n], " = ", 100, ")"))
sd3name <- expression(paste("Large standard deviation (", sigma^2, " = ", 9, ")"))
sd9name <- expression(paste("Large standard deviation (", sigma^2, " = ", 25, ")"))
smallNlargeP <- expression(paste("Small ", N, " large ", P))
smallNlargePsmallDm <- expression(paste("Small ", N, " large ", P, " (", Delta, mu, " = ", 0.2, ")"))
varyPropName <- expression(paste("Varying proportions"))
varyPropSmallDmName <- expression(paste("Varying proportions (", Delta, mu, " = ", 0.4, ")"))

# mcclust::arandi returns NA for two partitions of a single cluster
my_df$ARI[my_df$Model == "Mclust" & my_df$Scenario == "No structure" & is.na(my_df$ARI)] <- 1.0

my_df$Scenario <- as.numeric(my_df$Scenario)

my_df$Scenario <- factor(my_df$Scenario, labels = c(
  simple_2d_name,
  no_structure_name,
  base_case_name,
  varyPropName,
  varyPropSmallDmName,
  sd3name,
  sd9name,
  smallNlargeP,
  smallNlargePsmallDm,
  irr10name,
  irr20name,
  irr100name
))


p1 <- my_df %>%
  ggplot(aes(x = Model, y = ARI)) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip() +
  facet_wrap(~Scenario, labeller = "label_parsed") +
  theme(
    axis.text.y = element_text(hjust = 0.0,  size = 13),
    axis.text.x = element_text(angle = 0, size = 13),
    axis.title.y = element_blank(),
    axis.title.x=element_text(size=15),
    plot.title = element_text(size = 22, face = "bold"),
    #       plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 13)
  ) +
  labs(title = "Predictive performance")


p1

ggsave(paste0(save_dir, "simulation_model_prediction.png"), width = 14, height = 16, plot = p1)
# ggsave(paste0(data_dir, "allmodel_performance_prediction.png"), width = 20, height = 16)


p2 <- my_df %>%
  ggplot(aes(x = Model, y = Frobenius_norm)) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip() +
  facet_wrap(~Scenario, labeller = "label_parsed") +
  theme(
    axis.text.y = element_text(hjust = 0.0,  size = 13),
    axis.text.x = element_text(angle = 0, size = 13),
    axis.title.y = element_blank(),
    axis.title.x=element_text(size=15),
    plot.title = element_text(size = 22, face = "bold"),
    #       plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 13)
  ) +
  labs(title = "Uncertainty quantification",
       y = "Frobenius norm")

ggsave(paste0(save_dir, "simulation_uncertainty.png"), width = 14, height = 16, plot = p2)

library(dplyr)

na_check <- my_df %>%
  group_by(Model, Scenario) %>%
  summarise(NAs = sum(is.na(ARI)))
# summarise(qs = function(x)
#   {if(sum(is.na(ARI) == nrow(ARI))){
#     NA
#     }
#     else {
#       quantile(ARI, c(0.25, 0.75))
#       }
#     },
#   prob = c(0.25, 0.75))

na_check[na_check$NAs == 99, ]

summ_df <- my_df %>%
  dplyr::filter(Model != "Maximum likelihood (Mclust)" | Scenario != "No structure") %>%
  group_by(Model, Scenario) %>%
  summarise(Median = median(ARI))
# summarise(qs = quantile(ARI, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75))

ir100_df <- summ_df %>% filter(Scenario == "Irrelevant features 100")

ir100_df[20:25, ]
