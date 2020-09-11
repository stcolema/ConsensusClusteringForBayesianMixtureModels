
library(tibble)
library(ggplot2)
library(mdiHelpR)
library(magrittr)
library(pheatmap)
library(patchwork)

mdiHelpR::setMyTheme()
sim_col <- simColPal()
breaks <- defineBreaks(sim_col, lb = 0, ub = 1)

# 
cc_dir <- "C:/Users/stephen/Documents/PhD/Year_1/Consensus_inference/Consensus_inference_gen/Analysis/MDI_yeast_dataset/CMDLineMDI/ConsensusClustering/"

R <- c(100, 500)
S <- c(50, 100)
all_comb <- expand.grid(R, S)

cc_tibs <- paste0(cc_dir, "R", all_comb$Var1, "S", all_comb$Var2, "/compare_tibble.rds")

my_tibs <- lapply(cc_tibs[-1], readRDS)

cmR500S50 <- my_tibs[[1]]$similarity_matrix[[1]]
cmR100S100 <- my_tibs[[2]]$similarity_matrix[[1]]
cmR500S100 <- my_tibs[[3]]$similarity_matrix[[1]]



cm_order <- findOrder(cmR500S100)
ph_lst <- list(cmR500S50, cmR100S100, cmR500S100) %>% 
lapply(function(x) {
  pheatmap(x[cm_order, cm_order],
           silent = T,
           color = sim_col,
           breaks = breaks,
           legend = F,
           cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           show_colnames = F
  )
}) %>%
  lapply(extract2, "gtable")
my_list <- list(plot_spacer(), ph_lst[1:3])
layout <- "
##AA
BBCC
"
row.names(all_comb) <- c("R", "S")

my_lst <- list(gridExtra::tableGrob(all_comb[2:4, ]), unlist(ph_lst))
p1 <- patchwork::wrap_plots(ph_lst) + 
  plot_layout(design = layout) +
  plot_annotation(
    title = "Consensus matrices for time course data",
    theme= theme(axis.text.y=element_text(size = 10.5),
          axis.text.x=element_text(size = 10.5),
          axis.title.y=element_text(size = 10.5),
          axis.title.x=element_text(size = 10.5),
          plot.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 14),
          strip.text.x = element_text(size = 10.5),
          legend.text = element_text(size = 10.5)
  )
  )

p1

ggsave("./Images/Yeast/CMs.png", plot = p1,
       width = 6,
       height = 6)
