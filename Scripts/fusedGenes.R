
# Look at fused genes across datasets

# Where the data lives
library(tibble)

# Data manipulation
library(tidyr)

# BiocManager::install("clusterProfiler")
library(clusterProfiler)

# BiocManager::install("biomaRt")
library(biomaRt)

# Annotaiton database for yeast
# library(AnnotationHub)
library(org.Sc.sgd.db)

# Visualisation
library(ggplot2)
library("scales")
library(ggthemes)
library(patchwork)

# String manipulation
library(stringr)

# Semantic similarity of GO terms
library(GOSemSim)

# Heatmap colouring
library(mdiHelpR)
library(pheatmap)

# Pipe and associated functions
library(magrittr)

# LaTex tables
library(stargazer)

# Find genes that fused across the three datasets
findFusedYeast <- function(samples, fusion_threshold = 0.5) {
  S <- nrow(samples[[1]])

  genes_fused <- (colSums(samples[[1]] == samples[[2]] &
    samples[[1]] == samples[[3]] &
    samples[[2]] == samples[[3]]) > (S * fusion_threshold)) %>%
    which() %>%
    names() %>%
    stringr::str_remove_all("Dataset1_")

  genes_fused
}

findFusedPairwise <- function(samples, fusion_threshold = 0.5) {
  S <- nrow(samples[[1]])

  genes_fused <- (colSums(samples[[1]] == samples[[2]]) > (S * fusion_threshold)) %>%
    which() %>%
    names() %>%
    stringr::str_remove_all("Dataset1_")

  genes_fused
}

fusedTable <- function(genes, cl,
                       descriptions = org.Sc.sgdDESCRIPTION,
                       gene_names = org.Sc.sgdGENENAME) {
  cl_fused <- cl[match(genes, names(cl))]

  ## Bimap interface:
  # descriptions <- org.Sc.sgdDESCRIPTION

  # Convert to a list
  fused_descriptions <- as.list(descriptions[genes])

  fused_df <- fused_descriptions %>%
    unlist() %>%
    data.frame(Description = .) %>%
    cbind(data.frame(Gene = genes, Cluster = cl_fused))

  # Get the gene names that are mapped to an ORF identifier
  gene_names <- org.Sc.sgdGENENAME
  fused_df$Name <- as.list(gene_names[genes]) %>%
    unlist()

  fused_df <- fused_df[order(fused_df$Cluster), ]

  fused_df$Cluster <- as.numeric(as.factor(fused_df$Cluster))

  # Reorder columns
  as.tibble(fused_df[, c(2, 4, 3, 1)])
}

# === Setup ====================================================================

setMyTheme()

# Pheatmap visuals
col_pal <- simColPal()
breaks <- defineBreaks(col_pal, lb = 0)

# random seed
set.seed(1)

# Original data modelled
data_dir <- "./Data/Yeast/Original_data/"
datasets <- c("Timecourse", "ChIP-chip", "PPI")
data_files <- list.files(data_dir, full.names = T)[c(2, 1, 3)]
orig_data <- data_files %>%
  lapply(read.csv, row.names = 1) %>%
  set_names(datasets)

# The number of items in each dataset
N <- nrow(orig_data[[1]])

# The number of datasets
L <- length(datasets)

# Genes present
gene_names <- row.names(orig_data[[1]])

# Read in the CC and Bayes tibbles.
alloc_data <- readRDS("./Data/Yeast/CCAllocTibble.rds")
continuous_data <- readRDS("./Data/Yeast/CCContParamsTibble.rds")

fusion_threshold <- 0.5

S_final <- max(alloc_data$S)
R_final <- max(alloc_data$R)

cc_results <- alloc_data[which(alloc_data$S == S_final & alloc_data$R == R_final), ]

datasets <- cc_results$Dataset %>%
  unique()

L <- length(datasets)

fused_genes_pairwise <- list()
count <- 0

for (l in 1:(L - 1)) {
  for (m in (l + 1):L) {
    count <- count + 1
    fused_genes_pairwise[[count]] <- .genes <- findFusedPairwise(cc_results$Samples[c(l, m)])
    if (count == 1) {
      fused_count <- length(.genes)
    } else {
      fused_count <- c(fused_count, length(.genes))
    }
  }
}

# === All datasets =============================================================

genes_fused <- findFusedYeast(cc_results$Samples)

fused_df <- fusedTable(genes_fused,
  cc_results$Cl[[1]],
  descriptions = org.Sc.sgdDESCRIPTION,
  gene_names = org.Sc.sgdGENENAME
) %>%
  dplyr::arrange(Cluster, Name)

stargazer(fused_df, summary = FALSE, rownames = FALSE) %>%
  dplyr::arrange(Cluster, Name)

p1 <- orig_data$Timecourse[match(fused_df$Gene, row.names(orig_data$Timecourse)), ] %>%
  rownames_to_column("Gene") %>%
  mutate(Cluster = as.factor(fused_df$Cluster)) %>%
  # filter(Cluster %in% plt_clusters) %>%
  pivot_longer(-c(Gene, Cluster), names_to = "Timepoint", values_to = "Expression") %>%
  mutate(Timepoint = as.numeric(str_remove_all(Timepoint, "X"))) %>%
  ggplot(aes(x = Timepoint, y = Expression, group = Gene, colour = Cluster)) +
  geom_line() +
  facet_wrap(~Cluster, ncol = 1)

p2 <- orig_data$`ChIP-chip`[match(fused_df$Gene, row.names(orig_data$Timecourse)), colSums(orig_data$`ChIP-chip`) > 5] %>%
  rownames_to_column("Gene") %>%
  mutate(Cluster = as.factor(fused_df$Cluster)) %>%
  pivot_longer(-c(Gene, Cluster), names_to = "TF", values_to = "Interaction") %>%
  # filter(Cluster %in% plt_clusters) %>%
  ggplot(aes(x = TF, y = Gene, fill = Interaction)) +
  geom_tile() +
  # coord_flip()+
  scale_fill_gradient(low = "white", high = "#146EB4") +
  facet_wrap(~Cluster, ncol = 1, scales = "free_y") +
  theme(
    # axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank() # ,
    # axis.title.y = element_text(size = 10.5),
    # axis.title.x = element_text(size = 10.5),
    # plot.title = element_text(size = 18, face = "bold"),
    # plot.subtitle = element_text(size = 14),
    # strip.text.x = element_text(size = 10.5),
    # legend.text = element_text(size = 10.5)
    # strip.background = element_blank(),
    # strip.text.x = element_blank()
  ) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))

p3 <- orig_data$PPI[match(fused_df$Gene, row.names(orig_data$Timecourse)), colSums(orig_data$PPI) > 10] %>%
  rownames_to_column("Gene") %>%
  mutate(Cluster = as.factor(fused_df$Cluster)) %>%
  pivot_longer(-c(Gene, Cluster), names_to = "Protein", values_to = "Interaction") %>%
  # filter(Cluster %in% c(9, 12)) %>%
  ggplot(aes(x = Protein, y = Gene, fill = Interaction)) +
  geom_tile() +
  # coord_flip()+
  scale_fill_gradient(low = "white", high = "#146EB4") +
  facet_wrap(~Cluster, ncol = 1, scales = "free_y") +
  theme(
    # axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank() # ,
    # axis.title.y = element_text(size = 10.5),
    # axis.title.x = element_text(size = 10.5),
    # plot.title = element_text(size = 18, face = "bold"),
    # plot.subtitle = element_text(size = 14),
    # strip.text.x = element_text(size = 10.5),
    # legend.text = element_text(size = 10.5)
    # strip.background = element_blank(),
    # strip.text.x = element_blank()
  ) # +
# theme(axis.text.y = element_blank())

p_patch <- p1 + p2 + p3 +
  plot_layout(guides = "collect")
p_patch

# === Timecourse + ChIP-chip ===================================================

fused_timecourse_chipchip <- findFusedPairwise(cc_results$Samples[1:2])

timecourse_chipchip_table <- fusedTable(fused_timecourse_chipchip,
  cc_results$Cl[[1]],
  descriptions = org.Sc.sgdDESCRIPTION,
  gene_names = org.Sc.sgdGENENAME
) %>%
  dplyr::arrange(Cluster, Name)

# LaTex table for the fused clusters
stargazer(timecourse_chipchip_table, summary = FALSE, rownames = FALSE)

# Trnasform the data to ggplot form
long_chipchip_data <- orig_data$`ChIP-chip`[match(timecourse_chipchip_table$Gene, row.names(orig_data$Timecourse)), ] %>%
  rownames_to_column("Gene") %>%
  mutate(Cluster = as.factor(timecourse_chipchip_table$Cluster)) %>%
  pivot_longer(-c(Gene, Cluster), names_to = "TF", values_to = "Interaction")

long_timecourse_data <- orig_data$Timecourse[match(timecourse_chipchip_table$Gene, row.names(orig_data$Timecourse)), ] %>%
  rownames_to_column("Gene") %>%
  mutate(Cluster = as.factor(timecourse_chipchip_table$Cluster)) %>%
  pivot_longer(-c(Gene, Cluster), names_to = "Timepoint", values_to = "Expression")

# Look at the global pattern
long_timecourse_data %>%
  mutate(Timepoint = as.numeric(str_remove_all(Timepoint, "X"))) %>%
  ggplot(aes(x = Timepoint, y = Expression, group = Gene)) +
  geom_line()

long_chipchip_data %>%
  ggplot(aes(x = TF, y = Gene, fill = Interaction)) +
  geom_tile() +
  coord_flip() +
  scale_fill_gradient(low = "white", high = "#146EB4")

# Cluster with any interactions in the ChIP-chip data
# plt_clusters <- c(1, 2, 5, 9, 11, 12, 16, 17, 18, 20, 25, 26)

# Separate out the clusters into 4 plots to view the combinations of the
# Timecourse and ChIP-chip in a legible way
p1 <- long_timecourse_data %>%
  filter(Cluster %in% c(1:7)) %>%
  mutate(Timepoint = as.numeric(str_remove_all(Timepoint, "X"))) %>%
  ggplot(aes(x = Timepoint, y = Expression, group = Gene, colour = Cluster)) +
  geom_line() +
  facet_wrap(~Cluster, ncol = 1)

p2 <- long_timecourse_data %>%
  filter(Cluster %in% c(8:14)) %>%
  mutate(Timepoint = as.numeric(str_remove_all(Timepoint, "X"))) %>%
  ggplot(aes(x = Timepoint, y = Expression, group = Gene, colour = Cluster)) +
  geom_line() +
  facet_wrap(~Cluster, ncol = 1)

p3 <- long_timecourse_data %>%
  filter(Cluster %in% c(15:21)) %>%
  mutate(Timepoint = as.numeric(str_remove_all(Timepoint, "X"))) %>%
  ggplot(aes(x = Timepoint, y = Expression, group = Gene, colour = Cluster)) +
  geom_line() +
  facet_wrap(~Cluster, ncol = 1)

p4 <- long_timecourse_data %>%
  filter(Cluster %in% c(22:26)) %>%
  mutate(Timepoint = as.numeric(str_remove_all(Timepoint, "X"))) %>%
  ggplot(aes(x = Timepoint, y = Expression, group = Gene, colour = Cluster)) +
  geom_line() +
  facet_wrap(~Cluster, ncol = 1)


# Now the same for the ChIP-chip data
# This could have been done in two lists, but not worth redoing
p5 <- long_chipchip_data %>%
  filter(Cluster %in% c(1:7)) %>%
  ggplot(aes(x = TF, y = Gene, fill = Interaction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  facet_wrap(~Cluster, ncol = 1, scales = "free_y") +
  theme(
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))

p6 <- long_chipchip_data %>%
  filter(Cluster %in% c(8:14)) %>%
  ggplot(aes(x = TF, y = Gene, fill = Interaction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  facet_wrap(~Cluster, ncol = 1, scales = "free_y") +
  theme(
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))


p7 <- long_chipchip_data %>%
  filter(Cluster %in% c(15:21)) %>%
  ggplot(aes(x = TF, y = Gene, fill = Interaction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  facet_wrap(~Cluster, ncol = 1, scales = "free_y") +
  theme(
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))

p8 <- long_chipchip_data %>%
  filter(Cluster %in% c(22:26)) %>%
  ggplot(aes(x = TF, y = Gene, fill = Interaction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  facet_wrap(~Cluster, ncol = 1, scales = "free_y") +
  theme(
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))

# Look at these fused clusters across both datasets
p1 + p5 + plot_layout(guides = "collect")
p2 + p6 + plot_layout(guides = "collect")
p3 + p7 + plot_layout(guides = "collect")
p4 + p8 + plot_layout(guides = "collect")

# Keep only those with some interactions in the TF data and not singletons
cl_to_plt <- c(1:2, 5, 9, 11:12, 16:17, 20, 26)

# The timecourse data for these clusters
p9 <- long_timecourse_data %>%
  filter(Cluster %in% cl_to_plt) %>%
  mutate(Timepoint = as.numeric(str_remove_all(Timepoint, "X"))) %>%
  ggplot(aes(x = Timepoint, y = Expression, group = Gene, colour = Cluster)) +
  geom_line() +
  facet_wrap(~Cluster, ncol = 1) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    legend.text = element_text(size = 13)
  )

# We also want to highlight the TFs that define the strucutre in the ChIp-chip
# dataset
tfs_important <- c(
  "ACE2", "DIG1", "FKH1", "FKH2", "MBP1", "MCM1", "NDD1",
  "STE12", "SWI4", "SWI5", "SWI6", "TEC1", "YOX1"
)

# Find descriptions of these TFs
descriptions <- org.Sc.sgdDESCRIPTION

# Convert to a list
orfs <- as.list(org.Sc.sgdCOMMON2ORF)

# Remove probes that do not map in COMMON2ORF
tf_orfs <- orfs[match(tfs_important, names(orfs))] %>% 
  unlist()

tf_descriptions <- as.list(descriptions[tf_orfs])

tf__df <- tf_descriptions %>%
  unlist() %>%
  data.frame(Description = .) %>%
  cbind(data.frame(Gene = tfs_important, ORF = tf_orfs)) %>%
  dplyr::arrange(Gene)

stargazer(tf__df[, c(3, 2, 1)], summary = FALSE, rownames = FALSE)

# Now we do an awful hack to separate out these on the x-axis
tfs_pres <- long_chipchip_data$TF %>%
  unique()

# Rearranged to have space between the TFs of interest
tfs_pres_rearranged <- tfs_pres[c(1:3, 17, 4:9, 19, 10:16, 23, 18, 20:22, 24:26, 29, 27:28, 30:33, 34, 35:40, 66, 41:42, 44:50, 83, 51:58, 84, 59:65, 82, 67:73, 85, 75:81, 74, 86, 88:93, 87, 94:102, 107, 103:106, 108:110, 43, 111:117)]

# We need a continuous variable for this trick to work
long_chipchip_data$TF_num <- as.numeric(factor(long_chipchip_data$TF, levels = tfs_pres_rearranged))
TF_breaks_used <- long_chipchip_data$TF_num[match(tfs_important, long_chipchip_data$TF)]

# Now plot the ChIP-chip data with the important proteins highlighted
p10 <- long_chipchip_data %>%
  filter(Cluster %in% cl_to_plt) %>%
  ggplot(aes(x = TF_num, y = Gene, fill = Interaction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#146EB4") +
  facet_wrap(~Cluster, ncol = 1, scales = "free_y") +
  theme(
    axis.text.x = element_text(size = 11.5),
    axis.ticks = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(labels = tfs_important, breaks = TF_breaks_used) +
  theme(axis.text.y = element_blank()) +
  labs(x = "Transcription Factor") +
  geom_vline(
    xintercept = long_chipchip_data$TF_num[match(tfs_important, long_chipchip_data$TF)],
    colour = "red",
    lty = 2, 
    alpha = 0.4
  )

p9 + p10 +
  plot_layout(guides = "collect", widths = c(2, 5)) +
  plot_annotation(
    title = "Integrated genes",  #"Consensus clustering",
    subtitle = "Timecourse and ChIP-chip datasets",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 16)
    )
  )

ggsave("./SupplementaryMaterial/Images/Yeast/timecourseChIPchipFused.png",
  height = 10,
  width = 14
)


long_timecourse_data %>%
  filter(Cluster %in% c(1, 2, 5, 9, 11, 12, 16, 17, 20, 26)) %>%
  mutate(Timepoint = as.numeric(str_remove_all(Timepoint, "X"))) %>%
  ggplot(aes(x = Timepoint, y = Expression, group = Gene)) + #, colour = Cluster)) +
  geom_line() +
  facet_wrap(~Cluster, ncol = 1) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 21),
    plot.subtitle = element_text(size = 18)
  ) +
  geom_vline(
    xintercept = c(9, 18, 21, 32, 35),
    colour = "grey"
  ) +
  geom_hline(yintercept = c(-1, 0, 1), lty = 2, colour = "dark grey") +
  labs(
    subtitle = "Integrated genes",
    title = "Timecourse data"
  )

ggsave("./SupplementaryMaterial/Images/Yeast/timecourseOverlappingClusters.png",
       height = 12, 
       width = 7)


cl_5_and_12 <- timecourse_chipchip_table$Gene[timecourse_chipchip_table$Cluster %in% c(5, 12)]

inds <- match(cl_5_and_12, row.names(alloc_data$CM[[1]]))

alloc_data$CM[[118]][inds, ] %>% 
  extract(, colSums(.) > 0.5) %>% 
  pheatmap()

alloc_data$CM[[119]][inds, ] %>% 
  extract(, colSums(.) > 0.5) %>% 
  pheatmap()


cl18 <- alloc_data$CM[[118]][c(1:3, which(row.names(alloc_data$CM[[1]]) == "YLR209C")), ] %>% 
  magrittr::extract(, colSums(.) > 0.05)
  
pheatmap(cl18[, 15:45], cluster_cols = T)

# === GO term attempts =========================================================

# GO term over-representation in the fused clusters
cluster_comp <- cc_results$Cl[[1]][match(genes_fused, names(cc_results$Cl[[1]]))] %>%
  doClusterComparison(
    orig_data$Timecourse,
    geneTable,
    universe,
    ont = "ALL",
    drop_na = drop_na,
    min_cluster_size = 2
  )

cluster_comp <- cc_results$Cl[[1]][match(fused_timecourse_chipchip, names(cc_results$Cl[[1]]))] %>%
  doClusterComparison(
    orig_data$Timecourse,
    geneTable,
    universe,
    ont = "ALL",
    drop_na = drop_na,
    min_cluster_size = 2
  )

bayes_alloc_data <- readRDS("./Data/Yeast/BayesAllocTibble.rds")

n_chains <- bayes_alloc_data$Seed %>%
  unique() %>%
  length()

findFusedYeast(bayes_alloc_data[1:3, ]$Samples)

bayes_fused <- list()
n_bayes <- 0
for (i in 1:n_chains) {
  curr_ind <- 1 + (i - 1) * L
  if (bayes_alloc_data$Use_chain[curr_ind]) {
    curr_fused <- findFusedYeast(bayes_alloc_data$Samples[curr_ind:(curr_ind + (L - 1))])
    n_bayes <- n_bayes + 1
    bayes_fused[[n_bayes]] <- curr_fused
  }
}
