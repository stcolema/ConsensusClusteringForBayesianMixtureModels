
# Compare the GO term over representation for a subset of the ensembles.

library(magrittr)
library(ggplot2)
library(dplyr)

# Read in the data
goData <- read.csv("./Data/Yeast/AllGOoverRepresentationComparison.csv")


models_used <- c("Bayesian: chain 3",
                 "Bayesian: chain 5",
                 "Bayesian: chain 7",
                 "Bayesian: chain 8",
                 "Bayesian: chain 10",
                 "CC(10001,1000)"
)



plt_data <- goData %>% 
  filter(Model %in% models_used) 

plt_data$Model <- factor(plt_data$Model, models_used)
plt_data$ONTOLOGY <- factor(plt_data$ONTOLOGY, c("MF", "BP", "CC"))

p_lst <- vector("list", length = 3) %>% 
  set_names(c("MF", "BP", "CC"))

for(ont in c("MF", "BP", "CC")){
  p_lst[[ont]] <- plt_data %>% 
  filter(ONTOLOGY == ont) %>% 
  group_by(Description, Model) %>%
  mutate(Number_cluster = as.character(n())) %>%
  ggplot(aes(x = Description, y = Model)) +
  geom_point(aes(colour = log(p.adjust), size = Count)) + # , position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = Number_cluster), position = position_nudge(y = -0.3)) +
  # geom_jitter(width = 0, height = 0.4, seed = 1) +
  facet_grid(Dataset ~ ONTOLOGY, scales = "free_x") +
  scale_color_viridis_c(direction = -1) +
  labs(
    title = paste0("GO set over-representation (", ont, ")"),
    # subtitle = "Across chains",
    # x = "Descr",
    colour = "Log adjusted p-value"
  ) +
  theme(
    axis.text.y = element_text(hjust = 0.0, size = 8),
    axis.text.x = element_text(angle = 30, size = 8, vjust = 1, hjust = 1),
    axis.title.y = element_blank(),
    # axis.title.x=element_blank(),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 9)
  )
  
  ggsave(paste0("./SupplementaryMaterial/Images/Yeast/goEnrichmentComp", ont, ".png"),
         plot = p_lst[[ont]],
         height = 10,
         width = 20)
}

# p_lst[[1]]
# p_lst[[2]]
# p_lst[[3]]
#  
# dataset <- "Timecourse"
# ont <- "MF"
# plt_data %>%
#   filter(Dataset == dataset, ONTOLOGY == ont) %>%
#   ggplot(aes(x = Cluster_plt, y = Description, color = log(p.adjust), size = Count)) +
#   geom_point() +
#   facet_wrap(~ Model) +
#   scale_color_viridis_c(direction = -1) +
#   labs(
#     title = paste0("GO set over-representation"),
#     # subtitle = "Across chains",
#     x = "Cluster index",
#     colour = "Log adjusted p-value"
#   ) +
#   theme(
#     axis.text.y = element_text(hjust = 0.0, size = 8),
#     axis.text.x = element_text(angle = 0, size = 8),
#     axis.title.y = element_blank(),
#     # axis.title.x=element_blank(),
#     plot.title = element_text(size = 18, face = "bold"),
#     plot.subtitle = element_text(size = 14),
#     strip.text.x = element_text(size = 9)
#   ) +
#   facet_wrap(~Model, nrow = 2) # +
