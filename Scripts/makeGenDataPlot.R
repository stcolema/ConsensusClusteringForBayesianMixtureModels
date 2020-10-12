
# Visualise the generated data for the simulations

library(mdiHelpR)
library(ggplot2)
library(magrittr)
library(patchwork)

set.seed(1)

scn_table <- data.frame(
  Scenario = c("Simple 2D", "No structure", "Base Case", rep("Large N, small P", 3), rep("Large standard deviation", 3), rep("Irrelevant features", 5), rep("Small N, large P", 2), "Varying proportions"),
  N = c(100, 100, 2e2, 1e4, 1e4, 1e4, 2e2, 2e2, 2e2, 2e2, 2e2, 2e2, 2e2, 2e2, 50, 50, 200),
  P_s = c(2, 0, 20, 4, 4, 4, 20, 20, 20, 20, 20, 20, 20, 20, 500, 500, 20),
  P_n = c(0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 10, 20, 100, 200, 0, 0, 0),
  K = c(5, 1, 5, 5, 50, 50, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
  Delta_mu = c(3, 0, 1, 1, 1, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.2, 1),
  sigma2 = c(1, 1, 1, 1, 1, 1, 3, 5, 10, 1, 1, 1, 1, 1, 1, 1, 1),
  Pi = c(rep("vec(1/K)", 16), "(0.5, 0.25, 0.125, 0.0675, 0.0675)")
)

colnames(scn_table) <- c(
  "Scenario",
  "$N$",
  "$P_s$",
  "$P_n$",
  "$K$",
  "$\\Delta_{\\mu}$",
  "$\\sigma^2$",
  "$\\pi$"
)

sims_used <- scn_table[c(1:3, 7, 8, 11:13, 17, 15:16), ] %>%
  set_rownames(1:nrow(.))

scn <- "Simple 2D"
row_ind <- which(sims_used$Scenario == scn)
N <- sims_used[["$N$"]][row_ind]
P_s <- sims_used[["$P_s$"]][row_ind]
P_n <- sims_used[["$P_n$"]][row_ind]
K <- sims_used[["$K$"]][row_ind]
dm <- sims_used[["$\\Delta_{\\mu}$"]][row_ind]
s2 <- sims_used[["$\\sigma^2$"]][row_ind]
pi <- sims_used[["$\\pi$"]][row_ind]

if (pi == "vec(1/K)") {
  pi <- rep(1 / K, K)
} else {
  pi <- "(0.5, 0.25, 0.125, 0.0675, 0.0675)" %>%
    stringr::str_remove_all("[() ]") %>%
    strsplit(",") %>%
    unlist() %>%
    as.numeric()
}

# Generate dataset
simple2D_gen <- generateSimulationDataset(K, N, P_s,
  delta_mu = dm,
  p_n = P_n,
  cluster_sd = s2,
  pi = pi
)

col_order <- findOrder(t(noStructure_gen$data))

simple_2d_ph <- annotatedHeatmap(scale(simple2D_gen$data[, col_order]), simple2D_gen$cluster_IDs,
  main = "2D", # paste0(scn, " (seed 1)"),
  show_rownames = F,
  show_colnames = F,
  cluster_cols = F,
  silent = T
)

scn <- "No structure"
row_ind <- which(sims_used$Scenario == scn)
N <- sims_used[["$N$"]][row_ind]
P_s <- sims_used[["$P_s$"]][row_ind]
P_n <- sims_used[["$P_n$"]][row_ind]
K <- sims_used[["$K$"]][row_ind]
dm <- sims_used[["$\\Delta_{\\mu}$"]][row_ind]
s2 <- sims_used[["$\\sigma^2$"]][row_ind]
pi <- sims_used[["$\\pi$"]][row_ind]

if (pi == "vec(1/K)") {
  pi <- rep(1 / K, K)
} else {
  pi <- "(0.5, 0.25, 0.125, 0.0675, 0.0675)" %>%
    stringr::str_remove_all("[() ]") %>%
    strsplit(",") %>%
    unlist() %>%
    as.numeric()
}

# Generate dataset
noStructure_gen <- generateSimulationDataset(K, N, P_s,
  delta_mu = dm,
  p_n = P_n,
  cluster_sd = s2,
  pi = pi
)


col_order <- findOrder(t(noStructure_gen$data))

no_structure_ph <- annotatedHeatmap(scale(noStructure_gen$data[, col_order]), noStructure_gen$cluster_IDs,
  main = scn, # paste0(scn, " (seed 1)"),
  show_rownames = F,
  show_colnames = F,
  cluster_cols = F,
  silent = T
)

scn <- "Base Case"
row_ind <- which(sims_used$Scenario == scn)[1]
N <- sims_used[["$N$"]][row_ind]
P_s <- sims_used[["$P_s$"]][row_ind]
P_n <- sims_used[["$P_n$"]][row_ind]
K <- sims_used[["$K$"]][row_ind]
dm <- sims_used[["$\\Delta_{\\mu}$"]][row_ind]
s2 <- sims_used[["$\\sigma^2$"]][row_ind]
pi <- sims_used[["$\\pi$"]][row_ind]

if (pi == "vec(1/K)") {
  pi <- rep(1 / K, K)
} else {
  pi <- "(0.5, 0.25, 0.125, 0.0675, 0.0675)" %>%
    stringr::str_remove_all("[() ]") %>%
    strsplit(",") %>%
    unlist() %>%
    as.numeric()
}


# Generate dataset
baseCase_gen <- generateSimulationDataset(K, N, P_s,
  delta_mu = dm,
  p_n = P_n,
  cluster_sd = s2,
  pi = pi
)

col_order <- findOrder(t(baseCase_gen$data))

base_case_ph <- annotatedHeatmap(scale(baseCase_gen$data[, col_order]),
  baseCase_gen$cluster_IDs,
  main = scn, # paste0(scn, " (seed 1)"),
  show_rownames = F,
  show_colnames = F,
  cluster_cols = F,
  silent = T
)

scn <- "Large standard deviation"
row_ind <- which(sims_used$Scenario == scn)[1]
N <- sims_used[["$N$"]][row_ind]
P_s <- sims_used[["$P_s$"]][row_ind]
P_n <- sims_used[["$P_n$"]][row_ind]
K <- sims_used[["$K$"]][row_ind]
dm <- sims_used[["$\\Delta_{\\mu}$"]][row_ind]
s2 <- sims_used[["$\\sigma^2$"]][row_ind]
pi <- sims_used[["$\\pi$"]][row_ind]

if (pi == "vec(1/K)") {
  pi <- rep(1 / K, K)
} else {
  pi <- "(0.5, 0.25, 0.125, 0.0675, 0.0675)" %>%
    stringr::str_remove_all("[() ]") %>%
    strsplit(",") %>%
    unlist() %>%
    as.numeric()
}


# Generate dataset
largeStadDev_gen <- generateSimulationDataset(K, N, P_s,
  delta_mu = dm,
  p_n = P_n,
  cluster_sd = s2,
  pi = pi
)

col_order <- findOrder(t(largeStadDev_gen$data))

large_std_dev_ph <- annotatedHeatmap(scale(largeStadDev_gen$data[, col_order]),
  largeStadDev_gen$cluster_IDs,
  main = scn, # paste0(scn, " (seed 1)"),
  show_rownames = F,
  show_colnames = F,
  cluster_cols = F,
  silent = T
)

scn <- "Small N, large P"
row_ind <- which(sims_used$Scenario == scn)[2]
N <- sims_used[["$N$"]][row_ind]
P_s <- sims_used[["$P_s$"]][row_ind]
P_n <- sims_used[["$P_n$"]][row_ind]
K <- sims_used[["$K$"]][row_ind]
dm <- sims_used[["$\\Delta_{\\mu}$"]][row_ind]
s2 <- sims_used[["$\\sigma^2$"]][row_ind]
pi <- sims_used[["$\\pi$"]][row_ind]

if (pi == "vec(1/K)") {
  pi <- rep(1 / K, K)
} else {
  pi <- "(0.5, 0.25, 0.125, 0.0675, 0.0675)" %>%
    stringr::str_remove_all("[() ]") %>%
    strsplit(",") %>%
    unlist() %>%
    as.numeric()
}


# Generate dataset
smallNlargeP_gen <- generateSimulationDataset(K, N, P_s,
  delta_mu = dm,
  p_n = P_n,
  cluster_sd = s2,
  pi = pi
)

col_order <- findOrder(t(smallNlargeP_gen$data))


small_n_large_p_ph <- annotatedHeatmap(scale(smallNlargeP_gen$data[, col_order]),
  smallNlargeP_gen$cluster_IDs,
  main = scn, # paste0(scn, " (seed 1)"),
  show_rownames = F,
  show_colnames = F,
  cluster_cols = F,
  silent = T
)


scn <- "Irrelevant features"
row_ind <- which(sims_used$Scenario == scn)[3]
N <- sims_used[["$N$"]][row_ind]
P_s <- sims_used[["$P_s$"]][row_ind]
P_n <- sims_used[["$P_n$"]][row_ind]
K <- sims_used[["$K$"]][row_ind]
dm <- sims_used[["$\\Delta_{\\mu}$"]][row_ind]
s2 <- sims_used[["$\\sigma^2$"]][row_ind]
pi <- sims_used[["$\\pi$"]][row_ind]

if (pi == "vec(1/K)") {
  pi <- rep(1 / K, K)
} else {
  pi <- "(0.5, 0.25, 0.125, 0.0675, 0.0675)" %>%
    stringr::str_remove_all("[() ]") %>%
    strsplit(",") %>%
    unlist() %>%
    as.numeric()
}


# Generate dataset
irrelevantFeatures_gen <- generateSimulationDataset(K, N, P_s,
  delta_mu = dm,
  p_n = P_n,
  cluster_sd = s2,
  pi = pi
)

col_order <- findOrder(t(irrelevantFeatures_gen$data))


irr_features_ph <- annotatedHeatmap(scale(irrelevantFeatures_gen$data[, col_order]),
  irrelevantFeatures_gen$cluster_IDs,
  main = scn, # main = paste0("Irrelevant features (seed 1, ", P_n, " irrelevant features)"),
  show_rownames = F,
  show_colnames = F,
  cluster_cols = F,
  silent = T
)

scn <- "Varying proportions"
row_ind <- which(sims_used$Scenario == scn)[1]
N <- sims_used[["$N$"]][row_ind]
P_s <- sims_used[["$P_s$"]][row_ind]
P_n <- sims_used[["$P_n$"]][row_ind]
K <- sims_used[["$K$"]][row_ind]
dm <- sims_used[["$\\Delta_{\\mu}$"]][row_ind]
s2 <- sims_used[["$\\sigma^2$"]][row_ind]
pi <- sims_used[["$\\pi$"]][row_ind]

if (pi == "vec(1/K)") {
  pi <- rep(1 / K, K)
} else {
  pi <- "(0.5, 0.25, 0.125, 0.0675, 0.0675)" %>%
    stringr::str_remove_all("[() ]") %>%
    strsplit(",") %>%
    unlist() %>%
    as.numeric()
}


# Generate dataset
varyingProportions_gen <- generateSimulationDataset(K, N, P_s,
  delta_mu = dm,
  p_n = P_n,
  cluster_sd = s2,
  pi = pi
)


col_order <- findOrder(t(varyingProportions_gen$data))

varying_prop_ph <- annotatedHeatmap(scale(varyingProportions_gen$data[, col_order]),
  varyingProportions_gen$cluster_IDs,
  main = scn, # paste0(scn, " (seed 1)"),
  show_rownames = F,
  show_colnames = F,
  cluster_cols = F,
  silent = T
)


# simple_2d_ph
# no_structure_ph
# base_case_ph
# large_std_dev_ph
# small_n_large_p_ph
# irr_features_ph
# varying_prop_ph

gen_data_plt <- (wrap_elements(panel = simple_2d_ph$gtable, clip = FALSE) +
  wrap_elements(panel = no_structure_ph$gtable, clip = FALSE) +
  wrap_elements(panel = base_case_ph$gtable, clip = FALSE) +
  wrap_elements(panel = large_std_dev_ph$gtable, clip = FALSE) +
  wrap_elements(panel = small_n_large_p_ph$gtable, clip = FALSE) +
  wrap_elements(panel = irr_features_ph$gtable, clip = FALSE) +
  wrap_elements(panel = varying_prop_ph$gtable, clip = FALSE)) + plot_annotation(
  title = "Generated data",
  subtitle = "Seed set to 1"
)

ggsave("./Large_documents/First_year_report/Images/Simulations/generated_datasets.png", gen_data_plt,
  width = 12, height = 10
)
