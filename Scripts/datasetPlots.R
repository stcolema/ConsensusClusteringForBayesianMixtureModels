
library(mdiHelpR)
library(uplotter)
library(magrittr)
library(ggplot2)
library(patchwork)

set.seed(1)
setMyTheme()

K <- c(5, 5, 5)
P_s <- c(2, 500, 20)
P_n <- c(0, 0, 100)
N <- c(100, 50, 200)
dm <- c(3, 1, 1)
dataset_names <- c("Simple 2D", "Small N, large P", "Irrelevant features 100")

L <- length(N)
my_data <- list()

for (l in 1:L) {
  my_data[[l]] <- generateSimulationDataset(K[l], N[l], P_s[l], delta_mu = dm[l], p_n = P_n[l])
  my_data[[l]]$data <- scale(my_data[[l]]$data)
}

x1 <- uplotter::makeUMAPPlotData(my_data[[1]]$data, my_data[[1]]$cluster_IDs)

x1 %>% ggplot(aes(x = x, y = y, colour = as.factor(labels))) +
  geom_point() +
  labs(
    title = paste0(dataset_names[[1]], ": data"),
    x = "Feature 1",
    y = "Feature 2"
  ) +
  scale_color_viridis_d()

uplotter::visualiseUMAP(my_data[[1]]$data, as.factor(my_data[[1]]$cluster_IDs)) +
  labs(
    title = paste0(dataset_names[[1]], ": UMAP"),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_viridis_d()


uplotter::visualiseUMAP(my_data[[2]]$data, as.factor(my_data[[2]]$cluster_IDs)) +
  labs(
    title = paste0(dataset_names[[2]], ": UMAP"),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_viridis_d()

uplotter::visualiseUMAP(my_data[[3]]$data, as.factor(my_data[[3]]$cluster_IDs)) +
  labs(
    title = paste0(dataset_names[[3]], ": UMAP"),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_viridis_d()

simple2Dhm <- annotatedHeatmap(my_data[[1]]$data, my_data[[1]]$cluster_IDs,
  main = paste0(dataset_names[[1]]),
  show_rownames = F,
  show_colnames = F,
  silent = T
)

smallNlargePhm <- annotatedHeatmap(my_data[[2]]$data, my_data[[2]]$cluster_IDs,
  main = paste0(dataset_names[[2]]),
  show_rownames = F,
  show_colnames = F,
  silent = T
)

irr100hm <- annotatedHeatmap(my_data[[3]]$data, my_data[[3]]$cluster_IDs,
  main = paste0(dataset_names[[3]]),
  show_rownames = F,
  show_colnames = F,
  silent = T
)

threeSims <- wrap_plots(simple2Dhm$gtable, smallNlargePhm$gtable, irr100hm$gtable) +
  plot_annotation(
    title = "Simulated data",
    subtitle = "Random seed set to 1"
  )

ggsave(filename = "./Images/Simulations/Data.png", plot = threeSims,
       width = 13,
       height = 6)
