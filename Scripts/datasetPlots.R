
library(mdiHelpR)
library(uplotter)
library(magrittr)
library(ggplot2)
library(patchwork)

setMyTheme()



set.seed(1)

K <- c(5, 5, 5)
P_s <- c(2, 500, 20)
P_n <- c(0, 0, 100)
N <- c(100, 50, 200)
dm <- c(3, 1, 1)
dataset_names <- c("2D", "Small N, large P", "Irrelevant features")

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

row_order2D <- findOrder(my_data[[1]]$data)
col_order2D <- findOrder(t(my_data[[1]]$data))
simple2Dhm <- annotatedHeatmap(my_data[[1]]$data[row_order2D, col_order2D],
  my_data[[1]]$cluster_IDs,
  main = paste0(dataset_names[[1]]),
  show_rownames = F,
  show_colnames = F,
  cluster_rows = F,
  cluster_cols = F,
  border_color = NA,
  silent = T
)

row_order_largeP <- findOrder(my_data[[2]]$data)
col_order_largeP <- findOrder(t(my_data[[2]]$data))
smallNlargePhm <- annotatedHeatmap(my_data[[2]]$data[row_order_largeP, col_order_largeP],
  my_data[[2]]$cluster_IDs,
  main = paste0(dataset_names[[2]]),
  show_rownames = F,
  show_colnames = F,
  cluster_rows = F,
  cluster_cols = F,
  silent = T
)

row_order_irr100 <- findOrder(my_data[[3]]$data)
col_order_irr100 <- findOrder(t(my_data[[3]]$data))
irr100hm <- annotatedHeatmap(my_data[[3]]$data[row_order_irr100, ],
  my_data[[3]]$cluster_IDs,
  main = paste0(dataset_names[[3]]),
  show_rownames = F,
  show_colnames = F,
  cluster_rows = F,
  cluster_cols = F,
  silent = T
)

threeSims <- wrap_plots(simple2Dhm$gtable, smallNlargePhm$gtable, irr100hm$gtable,
                        nrow = 3
                        ) +
  plot_annotation(
    title = "Simulated data",
    subtitle = "Random seed set to 1",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  ) 

threeSims

ggsave(filename = "./Images/Simulations/Data.png", plot = threeSims,
       width = 6,
       height = 6)
