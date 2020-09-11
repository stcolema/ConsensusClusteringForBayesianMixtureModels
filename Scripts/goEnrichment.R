#!/bin/Rscript

# GO term over-representation

# === R setup ==================================================================

#
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

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




clusterList <- function(cl, geneTable, items_to_drop = NULL, exclude_singletons = T) {
  cl_lst <- list()
  labels_present <- unique(cl)
  K <- length(labels_present)
  
  list_ind <- 0
  
  for (i in 1:K) {
    k <- labels_present[i]
    cl_k <- cl == k
    
    if ((!exclude_singletons) || sum(cl_k) > 1) {
      list_ind <- list_ind + 1
      
      # Fidn the genes
      genes_in_cluster <- names(cl)[cl_k]
      
      # Drop the gene that cannot be mapped to an ENTREZID
      genes_in_cluster <- genes_in_cluster[genes_in_cluster != items_to_drop] # drop_na]
      
      # Find relevant ENTREZ IDs
      cluster_IDs <- geneTable$GeneID[match(genes_in_cluster, geneTable$Locus)]
      
      cl_lst[[list_ind]] <- cluster_IDs
    }
  }
  names(cl_lst) <- paste0("K", 1:list_ind)
  cl_lst
}

# ggplot theme
mdiHelpR::setMyTheme()

# Pheatmap visuals
col_pal <- simColPal()
breaks <- defineBreaks(col_pal, lb = 0)

# random seed
set.seed(1)

# Find directory
platform <- .Platform

linux <- T
if (platform$OS.type[1] == "windows") {
  linux <- F
}

if (linux) {
  phd_dir <- "/mnt/c/Users/stephen/Documents/PhD/"
} else {
  phd_dir <- "C:/Users/stephen/Documents/PhD/"
}

# Directory to work within
home_dir <- paste0(phd_dir, "Year_1/Consensus_inference/Consensus_inference_gen/Analysis/MDI_yeast_dataset/CMDLineMDI/")

cc_dir <- "ConsensusClustering/"
cc_ver <- "R500S100/"

bayes_dir <- "BayesianYeastN676000T1000Seed1K275/"
bayes_dir_2 <- "BayesianYeastN676000T1000Seed2K275/"
bayes_dir_3 <- "BayesianYeastN676000T1000Seed3K275/"
bayes_dir_4 <- "BayesianYeastN676000T1000Seed4K275/"

tib_name <- "compare_tibble.rds"

cc_tf <- paste0(home_dir, cc_dir, cc_ver, tib_name)

bayes_tf <- paste0(home_dir, bayes_dir, tib_name)
bayes_tf_2 <- paste0(home_dir, bayes_dir_2, tib_name)
bayes_tf_3 <- paste0(home_dir, bayes_dir_3, tib_name)
bayes_tf_4 <- paste0(home_dir, bayes_dir_4, tib_name)

# Save results
save_dir <- "./Images/Yeast/"

# Read in tibble containing PSM, predicted clustering, etc.
cc_tib <- readRDS(cc_tf)

bayes_tib <- readRDS(bayes_tf)
bayes_tib_2 <- readRDS(bayes_tf_2)
bayes_tib_3 <- readRDS(bayes_tf_3)
bayes_tib_4 <- readRDS(bayes_tf_4)

# === Yeast data ===============================================================

# Download and read in the Yeast genome
temp <- tempfile()
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz", temp)
x1 <- unz(temp, "GCF_000146045.2_R64_genomic.gff")
Gff2GeneTable(temp)
unlink(temp)

# The Gff2GeneTable function automatically saves to the working directory an
# object called "geneTable.rda". Load this.
load(paste0(getwd(), "/geneTable.rda"))

# Gff2GeneTable(paste0(my_d, "GCF_000146045.2_R64_genomic.gff"))
# load(paste0(my_d, "geneTable.rda"))


# Genes present in the dataset
genes_in_dataset <- match(row.names(timecourse_data), geneTable$Locus)

# Missing from database
missing <- which(is.na(genes_in_dataset))

if (length(missing) > 0) {
  cat("\nGene(s) present in analysed data not present in database. Missing:\n")
  cat(row.names(timecourse_data)[missing], "\n\n")
}

# The values used in filter the biomaRt goMap
eg <- geneTable$GeneID

# Find the dataset to use
ensembl <- useMart("ensembl")
available_datasets <- listDatasets(ensembl)

# Find the dataset name
dataset_name <- available_datasets %>%
  dplyr::filter(grepl("Saccharomyces|Yeast", description)) %>%
  dplyr::select(dataset) %>%
  unlist()

# Generate the dataset
yeast_ensemble <- useMart("ensembl", dataset = dataset_name)

attributes <- listAttributes(yeast_ensemble)

# Some attributes of interest (probably overboard)
attributes_to_use <- attributes %>%
  dplyr::filter(description %in% c(
    "GO term accession",
    "GO term name",
    "GO term definition",
    "GO term evidence code",
    "KEGG Pathway and Enzyme ID"
  )) %>%
  dplyr::select(name) %>%
  unlist()

# === GO annotation ============================================================

# Create GO map
gomap <- getBM(
  attributes = c("entrezgene_id",
                 "ensembl_gene_id",
                 "external_gene_name",
                 "go_id",
                 "name_1006",
                 "definition_1006",
                 "go_linkage_type",
                 "description"
  ),
  filters = "ensembl_gene_id",
  values = row.names(timecourse_data),
  # filters = "entrezgene_id",
  # values = eg[match(row.names(timecourse_data), geneTable$Locus)],
  mart = yeast_ensemble
)

# buld GO map saves a bunch of .rda files to the working directory
bgoMap <- buildGOmap(gomap)


# head(geneTable)
genes_present <- geneTable$GeneID[match(row.names(timecourse_data), geneTable$Locus)]

# === Consensus clustering =====================================================

enrichPlot <- function(cl, geneTable, universe,
                       drop = NULL, 
                       exclude_singletons = T,
                       ...){
  
  clustered_genes <- clusterList(cl, geneTable,
                                 items_to_drop = drop,
                                 exclude_singletons = exclude_singletons
  )
  
  ck <- compareCluster(
    geneCluster = clustered_genes,
    fun = "enrichGO",
    OrgDb = "org.Sc.sgd.db",
    universe = universe,
    ...
  )
  
  p <- dotplot(ck)
  
  p <- p$data %>%
    dplyr::mutate(Cluster_str = factor(str_replace_all(Cluster, "\\n", " "),
                                       levels = str_replace_all(levels(p$data$Cluster), "\\n", " "),
                                       ordered = T
    )) %>%
    ggplot(aes(x = Cluster_str, y = Description)) +
    geom_point(aes(colour = p.adjust, size = Count)) #  +
    # theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1.)) +
    # labs(
    #   x = "Predicted cluster (consensus clustering)",
    #   y = "GO description",
    #   title = "Timecourse",
    #   subtitle = "Enrichment across clusters"
    # )
  p
}

ontology <- "MF"

cc_cl <- mcclust::maxpear(cc_tib$similarity_matrix[[1]], max.k = 250)$cl

drop_na <- row.names(timecourse_data)[missing]

p1 <- enrichPlot(cc_cl, geneTable, genes_present[-missing],
                       drop = drop_na, 
                       exclude_singletons = T,
                 keyType = "ENTREZID",
                 ont = ontology,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.20)

p1 +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1.)) +
  labs(
    x = "Predicted cluster (Bayesian)",
    y = "GO description",
    title = "Timecourse",
    subtitle = "Enrichment across clusters"
  )

ggsave(paste0(save_dir, "/goEnrichmentCC.png"), p1, width = 6, height = 6)


# === Bayesian inference =======================================================

bayes_pred_timecourse <- bayes_tib$pred_allocation[[1]]

p2 <- enrichPlot(bayes_pred_timecourse, geneTable, genes_present[-missing],
                 drop = drop_na, 
                 exclude_singletons = T,
                 keyType = "ENTREZID",
                 ont = ontology,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.20)

p2 +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1.)) +
  labs(
    x = "Predicted cluster (Bayesian)",
    y = "GO description",
    title = "Timecourse",
    subtitle = "Enrichment across clusters"
  )

ggsave(paste0(save_dir, "/goEnrichmentBayes1.png"), p2, width = 6, height = 6)


# === Comparison ===============================================================

cc_data <- p1$data
bayes_data <- p2$data

cc_data$Inference <- "CC"
bayes_data$Inference <- "Bayesian"

plt_data <- rbind(cc_data, bayes_data)
plt_data$Cluster_plt <- stringr::str_extract(plt_data$Cluster_str, "\\d+") %>% as.numeric() 


p3 <- plt_data %>% ggplot(aes(x = Cluster_plt, y = Description)) +
  geom_point(aes(colour = p.adjust, size = Count)) +
  labs(
    title = "GO set over-representation",
    x = "Cluster index",
    colour = "Adjusted p-value"
  ) +
  theme(axis.text.y=element_text(hjust=0.0, size = 8),
        axis.text.x=element_text(angle=0, size = 8),
        axis.title.y=element_blank(),
        # axis.title.x=element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 10.5)
  ) +
  facet_wrap(~Inference) +
  xlim(1, 55)
p3

ggsave(paste0(save_dir, "/GOoverRepresentationBayesCC.png"),
       p3,
       width = 8,
       height = 6
)

# === Compare chains GO ========================================================

bayes_pred_timecourse_2 <- bayes_tib_2$pred_allocation[[1]]
bayes_pred_timecourse_3 <- bayes_tib_3$pred_allocation[[1]]
bayes_pred_timecourse_4 <- bayes_tib_4$pred_allocation[[1]]

p_lst <- lapply(list(bayes_pred_timecourse_2,
                                         bayes_pred_timecourse_3,
                                         bayes_pred_timecourse_4),
                                    enrichPlot,
                                    geneTable, genes_present[-missing],
                                     drop = drop_na, 
                                     exclude_singletons = T,
                                     keyType = "ENTREZID",
                                     ont = ontology,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.20)



dat_lst <- list()
dat_lst[[1]] <- bayes_data
dat_lst[[1]]$Chain <- 1
for(i in 2:4){
  
  dat_lst[[i]] <- p_i$data
  dat_lst[[i]]$Inference <- "Bayesian"
  dat_lst[[i]]$Chain <- i
}

bayes_data<- do.call(rbind, dat_lst)

bayes_data$Model <- paste0(bayes_data$Inference, ": chain ", bayes_data$Chain)

cc_data$Chain <- NA
cc_data$Model <- "Consensus Clustering"

new_plt_data <- rbind(cc_data, bayes_data)
new_plt_data$Cluster_plt <- stringr::str_extract(new_plt_data$Cluster_str, "\\d+") %>% as.numeric()


p4 <- new_plt_data %>% ggplot(aes(x = Cluster_plt, y = Description)) +
  geom_point(aes(colour = p.adjust, size = Count)) +
  labs(
    title = "GO set over-representation",
    # subtitle = "Across chains",
    x = "Cluster index",
    colour = "Adjusted p-value"
  ) +
  theme(axis.text.y=element_text(hjust=0.0, size = 8),
        axis.text.x=element_text(angle=0, size = 8),
        axis.title.y=element_blank(),
        # axis.title.x=element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 9)
  ) +
  facet_wrap(~Model, nrow = 1) # +
# xlim(1, 55)
p4

ggsave(paste0(save_dir, "/GOoverRepresentation.png"),
       p4,
       width = 12,
       height = 8
)
