
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


# Find genes that fused across the three datasets
findFusedYeast <- function(samples, fusion_threshold = 0.5){
  S <- nrow(samples[[1]])
  
  genes_fused <- (colSums(samples[[1]] == samples[[2]] &
                            samples[[1]] == samples[[3]] & 
                            samples[[2]] == samples[[3]]
  ) > (S * fusion_threshold)) %>% 
    which() %>% 
    names() %>% 
    stringr::str_remove_all("Dataset1_")
  
  genes_fused
}

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

(colSums(cc_results$Samples[[1]] == cc_results$Samples[[2]]) > (S_final * fusion_threshold)) %>%
  sum()

(colSums(cc_results$Samples[[1]] == cc_results$Samples[[3]]) > (S_final * fusion_threshold)) %>%
  sum()

(colSums(cc_results$Samples[[2]] == cc_results$Samples[[3]]) > (S_final * fusion_threshold)) %>%
  sum()

genes_fused <- findFusedYeast(cc_results$Samples)


## Bimap interface:
x <- org.Sc.sgdDESCRIPTION

# Convert to a list
fused_descriptions <- as.list(x[genes_fused])

y <- fused_descriptions %>% 
  unlist() %>% 
  data.frame() %>% 
  cbind(data.frame(Gene = genes_fused, Cluster = cl_fused)) 

y <- y[order(y$Cluster),]

y$Cluster <- as.numeric(as.factor(y$Cluster))


library(stargazer)

stargazer(y[, c(2, 3, 1)], summary=FALSE, rownames=FALSE)


cc_results$Cl[[1]][match(genes_fused, names(cc_results$Cl[[1]]))]
cc_results$Cl[[2]][match(genes_fused, names(cc_results$Cl[[2]]))]
cc_results$Cl[[3]][match(genes_fused, names(cc_results$Cl[[3]]))]

unique(cc_results$Cl[[2]][match(genes_fused, names(cc_results$Cl[[2]]))]) %>% length()

      
cluster_comp <- cc_results$Cl[[1]][match(genes_fused, names(cc_results$Cl[[1]]))] %>%
  doClusterComparison(
         orig_data[[dataset]],
         geneTable,
         universe,
         ont = ontology,
         drop_na = drop_na
  )

bayes_alloc_data <- readRDS("./Data/Yeast/BayesAllocTibble.rds")

n_chains <- bayes_alloc_data$Seed %>% 
  unique() %>% 
  length()

findFusedYeast(bayes_alloc_data[1:3, ]$Samples)

bayes_fused <- list()
n_bayes <- 0
for(i in 1:n_chains){
  curr_ind <- 1 + (i - 1) * L
  if(bayes_alloc_data$Use_chain[curr_ind]){
    curr_fused <- findFusedYeast(bayes_alloc_data$Samples[curr_ind:(curr_ind + (L - 1))])
    n_bayes <- n_bayes + 1
    bayes_fused[[n_bayes]] <- curr_fused
  }
  
}

