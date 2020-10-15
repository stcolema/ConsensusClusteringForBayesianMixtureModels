#!/bin/Rscript

# GO term over-representation in the predicted clusters for each dataset for 
# the long chains and a number of ensembles.

# === Packages =================================================================

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

# Pipe and associated functions
library(magrittr)

# === Functions ================================================================

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



goOverRep <- function(cl, X, geneTable, universe = row.names(X), drop_na = NULL, ont = "ALL", save_dir = NULL) {
  labels_present <- unique(cl)
  n_labels <- length(labels_present)

  ego <- list()

  for (i in 1:n_labels) {
    k <- labels_present[i]
    N_k <- sum(cl == k)

    # Fidn the genes
    genes_in_cluster <- row.names(X)[cl == k]

    # Drop the gene that cannot be mapped to an ENTREZID
    if (!is.null(drop_na)) {
      genes_in_cluster <- genes_in_cluster[genes_in_cluster != drop_na]
    }

    # Find relevant ENTREZ IDs
    cluster_IDs <- geneTable$GeneID[match(genes_in_cluster, geneTable$Locus)]

    # Enrichment using chosen ontology
    ego[[i]] <- enrichGO(
      gene = cluster_IDs,
      OrgDb = org.Sc.sgd.db,
      universe = universe, # [-missing],
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.15
    )

    # png(paste0(home_dir, cc_dir, cc_ver, "goEnrichmentCluster", k, ".png"))
    p <- dotplot(ego[[i]],
      showCategory = 10,
      title = paste0("Timecourse: Enriched genes for cluster ", k)
    )
    # p <- dotplot(ego[[38]], showCategory = 10, title = paste0("Timecourse: Enriched genes for cluster ", 38))




    if (sum(ego[[i]]@result$p.adjust < 0.05) > 0) {
      p2 <- p$data %>%
        mutate(Description_sent = stringr::str_to_sentence(Description)) %>%
        ggplot(aes(x = GeneRatio, y = Description)) +
        geom_point(aes(colour = p.adjust, size = Count)) +
        labs(
          x = "Gene ratio in cluster",
          y = "Description",
          title = "Timecourse data",
          subtitle = paste0("Enriched genes for cluster ", k),
          caption = paste0("Cluster membership is ", N_k, " items.")
        )


      p_cc_upset <- enrichplot::upsetplot(ego[[i]], n = 10) +
        labs(
          title = "Timecourse data",
          subtitle = paste0("Enriched genes for cluster ", k),
          caption = paste0("Cluster membership is ", N_k, " items.")
        )

      if (!is.null(save_dir)) {
        ggsave(paste0(save_dir, "/goEnrichmentCluster", k, ".png"), p2,
          width = 8, height = 4
        )
        ggsave(paste0(save_dir, "/UpsetPlotCluster", k, ".png"), p_cc_upset,
          width = 8, height = 4
        )
      }
    }

    cc_new <- ego[[i]]@result
    cc_new$K <- NULL
    cc_new$N <- NULL
    if (nrow(cc_new) > 0) {
      cc_new$K <- k
      cc_new$N <- sum(cl == k)
    }

    if (i == 1) {
      cc_go <- cc_new
    } else {
      cc_go <- rbind(cc_go, cc_new)
    }
  }
  cc_go
}



doClusterComparison <- function(cl, X, geneTable, universe, ont = "MF", drop_na = NULL) {
  clustered_genes <- list()
  labels_present <- unique(cl)
  n_labels <- length(labels_present)

  for (i in 1:n_labels) {
    k <- labels_present[i]

    # Fidn the genes
    genes_in_cluster <- row.names(X)[cl == k]

    # Drop the gene that cannot be mapped to an ENTREZID
    genes_in_cluster <- genes_in_cluster[genes_in_cluster != drop_na]

    # Find relevant ENTREZ IDs
    cluster_IDs <- geneTable$GeneID[match(genes_in_cluster, geneTable$Locus)]

    clustered_genes[[i]] <- cluster_IDs
  }

  names(clustered_genes) <- paste0("K", labels_present)

  clustered_genes <- clusterList(cl, geneTable,
    items_to_drop = drop_na,
    exclude_singletons = T
  )

  cl_comp <- compareCluster(
    geneCluster = clustered_genes,
    fun = "enrichGO",
    OrgDb = "org.Sc.sgd.db",
    universe = universe,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )

  p <- dotplot(cl_comp)


  p2 <- p$data %>%
    dplyr::mutate(Cluster_str = factor(str_replace_all(Cluster, "\\n", " "),
      levels = str_replace_all(levels(p$data$Cluster), "\\n", " "),
      ordered = T
    )) %>%
    ggplot(aes(x = Cluster_str, y = Description)) +
    geom_point(aes(colour = log(p.adjust), size = Count)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1.))
  p2
}

# === Setup ====================================================================

# ggplot theme
mdiHelpR::setMyTheme()

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

# Cluster analysis tibbles
tib_dir <- "./Data/Yeast/"
tib_files <- list.files(tib_dir, pattern = "AllocTibble.rds$")
tib_details <- stringr::str_match(tib_files, "(\\S*)Alloc")
tibbles <- list.files(tib_dir, pattern = "AllocTibble.rds$", full.names = T) %>%
  lapply(readRDS)

n_tibbles <- length(tibbles)

names(tibbles) <- tib_details[, 2]

# GO ontologies to investigate
ontologies <- c("MF", "BP", "CC")

# Save results
go_dir <- "Go_enrichment"

# === Yeast data ===============================================================

# Download and read in the Yeast genome if it is not already saved
if(! file.exists(paste0(getwd(), "/geneTable.rda"))){
temp <- tempfile()
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz", temp)
x1 <- unz(temp, "GCF_000146045.2_R64_genomic.gff")
Gff2GeneTable(temp)
unlink(temp)
}

# The Gff2GeneTable function automatically saves to the working directory an
# object called "geneTable.rda". Load this.
load(paste0(getwd(), "/geneTable.rda"))

# Gff2GeneTable(paste0(my_d, "GCF_000146045.2_R64_genomic.gff"))
# load(paste0(my_d, "geneTable.rda"))


# Genes present in the dataset
genes_in_dataset <- match(gene_names, geneTable$Locus)

# Missing from database
missing <- which(is.na(genes_in_dataset))

if (length(missing) > 0) {
  cat("\nGene(s) present in analysed data not present in database. Missing:\n")
  cat(gene_names[missing], "\n\n")
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

# === GO map ===================================================================

# Create GO map
gomap <- getBM(
  attributes = c(
    attributes_to_use,
    "entrezgene_id",
    "ensembl_gene_id",
    "external_gene_name",
    # "name_1006",
    # "definition_1006",
    "description"
  ),
  filters = "ensembl_gene_id",
  values = gene_names,
  mart = yeast_ensemble
)

# buld GO map saves a bunch of .rda files to the working directory
bgoMap <- buildGOmap(gomap)


# head(geneTable)
genes_present <- geneTable$GeneID[match(gene_names, geneTable$Locus)]

# === Random partition =========================================================

R_max <- max(tibbles$CC$R)
S_max <- max(tibbles$CC$S)

cc_tib <- tibbles$CC
K_cc <- cc_tib$Cl[which(cc_tib$S == S_max & cc_tib$R == R_max)] %>% 
  lapply(function(x){length(unique(x))}) %>% 
  unlist() %>% 
  set_names(datasets)

rand_cl <- K_cc %>% lapply(function(k){
  sample(1:k, size = N, replace = T)
})

drop_na <- gene_names[missing]
universe <- geneTable$GeneID[match(gene_names, geneTable$Locus)]

rand_go <- list()
for(l in 1:L){
rand_go[[l]] <- tryCatch(
    {
      doClusterComparison(rand_cl[[l]], orig_data[[l]], geneTable, universe, 
                          ont = "ALL",
                          drop_na = drop_na)
    },
    error=function(cond) {
      message(paste("No enriched GO terms found."))
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("Clusterign cause a warning?"))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
}

rand_time_cl <- doClusterComparison(rand_cl[[1]], orig_data[[1]], geneTable, universe, ont = "ALL", drop_na = drop_na)
rand_time_cl
# Error in compareCluster(geneCluster = clustered_genes, fun = "enrichGO",  : 
# No enrichment found in any of gene cluster, please check your input... 

# === Dataset choice ===========================================================

my_data <- NULL
R_used <- c(501, 1001,5001,10001)
S_used <- c(100, 500, 1000)

for (dataset in datasets) {
  save_dir <- paste0("./Images/Yeast/", dataset, "/")
  dir.create(save_dir)

  drop_na <- gene_names[missing]
  universe <- geneTable$GeneID[match(gene_names, geneTable$Locus)]

  ontology <- "ALL"
  # for (ont in ontologies) {
    for (i in 1:n_tibbles) {

      # Select the tibble
      .tib <- tibbles[[i]]

      # Indices for analyses matching the dataset
      curr_ind <- which(.tib$Dataset == dataset)

      # Meta data
      method <- "Consensus clustering"
      if (names(tibbles)[i] == "Bayes") {
        curr_ind <- which(.tib$Dataset == dataset & .tib$Use_chain)
        method <- "Bayesian"
      } else {
        curr_ind <- which(.tib$Dataset == dataset & .tib$R %in% R_used & .tib$S %in% S_used)
      }

      n_ind <- length(curr_ind)

      # === Consensus clustering =====================================================

      cluster_comp <- .tib[curr_ind, ]$Cl %>%
        lapply(doClusterComparison,
          orig_data[[dataset]],
          geneTable,
          universe,
          ont = ontology,
          drop_na = drop_na
        )

      # cc_p <- lapply(cc_tib, function(x) {
      #   cl <- mcclust::maxpear(x$similarity_matrix[[curr_ind]], max.k = 250)$cl
      #   p <- doClusterComparison(cl, orig_data[[dataset_name]], geneTable, universe,
      #                            ont = ont,
      #                            drop_na = drop_na
      #   )
      #   p
      # })

      for (j in 1:n_ind) {
        if (names(tibbles)[i] == "Bayes") {
          chain_num <- .tib$Seed[curr_ind][j]
          r <- .tib$R[curr_ind][j]
          subtitle <- paste0("Enrichment across clusters, chain ", chain_num)
          save_name <- paste0(save_dir, dataset_name, "goEnrichmentAllClustersBayes", chain_num, ".png")
          s <- NA
        } else {
          r <- .tib$R[curr_ind][j]
          s <- .tib$S[curr_ind][j]
          subtitle <- paste0("Enrichment across clusters, Consensus(", r, ",", s, ")")
          save_name <- paste0(save_dir, dataset_name, "goEnrichmentAllClustersCCR", r, "S", s, ".png")
          chain_num <- NA
        }

        # p_curr <- cluster_comp[[j]] +
        #   labs(
        #     x = "Predicted cluster",
        #     y = "GO description",
        #     title = dataset,
        #     subtitle = subtitle
        #   ) +
        #   scale_color_viridis_c(direction = -1)

        # ggsave(save_name, p_curr, width = 14, height = 10)


        curr_data <- cluster_comp[[j]]$data
        curr_data$Inference <- method
        curr_data$R <- r
        curr_data$S <- s
        curr_data$Chain <- chain_num
        curr_data$Dataset <- dataset
        # curr_data$Ontology <- ont

        if (method == "Bayesian") {
          curr_data$Model <- paste0("Bayesian: chain ", chain_num)
        } else {
          curr_data$Model <- paste0("CC(", r, ",", s, ")")
        }

        curr_data$Cluster_plt <- stringr::str_extract(curr_data$Cluster_str, "\\d+") %>%
          as.numeric()
        # my_data$Model <- factor(my_data$Model, levels = unique(my_data$Model))
        
        if (is.null(my_data)) {
          my_data <- curr_data
        } else {
          my_data <- rbind(my_data, curr_data)
        }
      }
    }

    # ego <- goOverRep(cc_cl, timecourse_data, geneTable, universe = universe,
    #                  drop_na = drop_na, ont = ont, save_dir = NULL)

    # Save annotation data.frame
    # write.csv(cc_go, paste0(save_dir[2], "/goEnrichmentCC.csv"))


    # Compare all clusters in a single plot
    # p1 <- doClusterComparison(cc_cl, timecourse_data, geneTable, universe,
    #   ont = ont,
    #   drop_na = drop_na
    # )
    # p1 +
    #   labs(
    #     x = "Predicted cluster (consensus clustering)",
    #     y = "GO description",
    #     title = dataset_name,
    #     subtitle = "Enrichment across clusters"
    #   ) +
    #   scale_color_viridis_c(direction = -1)
    #
    # # ggsave(paste0(save_dir, "/goEnrichmentAllClustersCC.png"), p1, width = 14, height = 10)
    #

    # === Bayesian inference =======================================================

    # bayes_p <- lapply(bayes_tibs, function(x) {
    #   cl <- x$pred_allocation[[curr_ind]]
    #   p <- doClusterComparison(cl, orig_data[[dataset_name]], geneTable, universe,
    #     ont = ont,
    #     drop_na = drop_na
    #   )
    #   p
    # })
    #
    # for (i in 1:10) {
    #   curr_data <- bayes_p[[i]]$data
    #   curr_data$Inference <- "Bayesian"
    #   curr_data$Chain <- i
    #   curr_data$Model <- paste0("Bayesian: chain ", i)
    #
    #   if (i == 1) {
    #     bayes_data <- curr_data
    #   } else {
    #     bayes_data <- rbind(bayes_data, curr_data)
    #   }
    # }
    #
    # # === Comparison ===============================================================
    #
    # # cc_data <- p1$data
    # # cc_data$Inference <- "Consensus clustering"
    # cc_data$Chain <- NA
    # # cc_data$Model <- "Consensus Clustering"
    #
    # bayes_data$R <- 676000
    # bayes_data$S <- 676
    #
    # # Bind the data
    # new_plt_data <- rbind(cc_data, bayes_data)

  for(ont in ontologies){
  
    # Plot
    p_comp <- my_data %>%
      filter(Dataset == dataset, ONTOLOGY == ont) %>%
      ggplot(aes(x = Cluster_plt, y = Description, color = log(p.adjust), size = Count)) +
      geom_point() +
      facet_wrap(~Model) +
      scale_color_viridis_c(direction = -1) +
      labs(
        title = paste0(dataset, ": GO set over-representation (", ont, ")"),
        # subtitle = "Across chains",
        x = "Cluster index",
        colour = "Log adjusted p-value"
      ) +
      theme(
        axis.text.y = element_text(hjust = 0.0, size = 8),
        axis.text.x = element_text(angle = 0, size = 8),
        axis.title.y = element_blank(),
        # axis.title.x=element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 9)
      ) +
      facet_wrap(~Model, nrow = 2) # +

    # xlim(1, 55)
    # p_comp

    p_comp_2 <- p_comp$data %>%
      group_by(Description, Model) %>%
      mutate(Number_cluster = as.character(n())) %>%
      ggplot(aes(x = Description, y = Model)) +
      geom_point(aes(colour = log(p.adjust), size = Count)) + # , position = position_stack(reverse = TRUE)) +
      geom_text(aes(label = Number_cluster), position = position_nudge(y = -0.3)) +
      # geom_jitter(width = 0, height = 0.4, seed = 1) +
      scale_color_viridis_c(direction = -1) +
      labs(
        title = paste0(dataset, ": GO set over-representation (", ont, ")"),
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


    ggsave(paste0(save_dir, dataset, "GOoverRepresentationComparison", ont, ".png"),
      p_comp,
      width = 14,
      height = 14
    )

    ggsave(paste0(save_dir, dataset, "GOoverRepresentationComparison", ont, "Alt.png"),
      p_comp_2,
      width = 14,
      height = 8
    )
    }
    # write.csv(p_comp$data, paste0(save_dir, dataset, "GOoverRepresentationComparison", ont, ".csv"))

    # p1 <- new_plt_data %>%
    #   filter(Chain %in% c(1:5, NA)) %>%
    #   ggplot(aes(x = Cluster_plt, y = Description, color = log(p.adjust), size = Count)) +
    #   geom_point() +
    #   facet_wrap(~Model) +
    #   scale_color_viridis_c(direction = -1) +
    #   labs(
    #     title = "GO set over-representation",
    #     # subtitle = "Across chains",
    #     x = "Cluster index",
    #     colour = "Log adjusted p-value"
    #   ) +
    #   theme(axis.text.y=element_text(hjust=0.0, size = 8),
    #         axis.text.x=element_text(angle=0, size = 8),
    #         axis.title.y=element_blank(),
    #         # axis.title.x=element_blank(),
    #         plot.title = element_text(size = 18, face = "bold"),
    #         plot.subtitle = element_text(size = 14),
    #         strip.text.x = element_text(size = 9)
    #   ) +
    #   facet_wrap(~Model, nrow = 1) # +
    # p1
    #
    # p2 <- new_plt_data %>%
    #   filter(Chain %in% c(6:10, NA)) %>%
    #   ggplot(aes(x = Cluster_plt, y = Description, color = log(p.adjust), size = Count)) +
    #   geom_point() +
    #   facet_wrap(~Model) +
    #   scale_color_viridis_c(direction = -1) +
    #   labs(
    #     title = "GO set over-representation",
    #     # subtitle = "Across chains",
    #     x = "Cluster index",
    #     colour = "Log adjusted p-value"
    #   ) +
    #   theme(axis.text.y=element_text(hjust=0.0, size = 8),
    #         axis.text.x=element_text(angle=0, size = 8),
    #         axis.title.y=element_blank(),
    #         # axis.title.x=element_blank(),
    #         plot.title = element_text(size = 18, face = "bold"),
    #         plot.subtitle = element_text(size = 14),
    #         strip.text.x = element_text(size = 9)
    #   ) +
    #   facet_wrap(~Model, nrow = 1) # +
    # p2
    #
    # p1 / p2  +
    #   plot_layout(guides = 'collect')
}

my_data$Model <- factor(my_data$Model, levels = unique(my_data$Model))

write.csv(my_data, "./Data/Yeast/AllGOoverRepresentationComparison.csv")
