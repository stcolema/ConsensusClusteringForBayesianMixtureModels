#!/bin/Rscript

# GO term over-representation

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

# Read in analysed data
data_dir <- paste0(phd_dir, "Year_1/Consensus_inference/Consensus_inference_gen/Data/YeastData/")

harbison_data <- read.csv(paste0(
  data_dir,
  "harbison_ppi_marina_multinomial.csv"
),
row.names = 1
)

timecourse_data <- read.csv(paste0(
  data_dir,
  "marina_ppi_harbison.csv"
),
row.names = 1
)

ppi_data <- read.csv(paste0(
  data_dir,
  "ppi_marina_harbisonreduced2.csv"
),
row.names = 1
)

N <- nrow(harbison_data)

# Directory to work within
home_dir <- paste0(phd_dir, "Year_1/Consensus_inference/Consensus_inference_gen/Analysis/MDI_yeast_dataset/CMDLineMDI/")

all_dirs <- list.dirs(home_dir, recursive = F)

bayesian_dirs <- all_dirs[grep(all_dirs, pattern = "*BayesianYeastN*")]

cc_dir <- "ConsensusClustering/"
cc_ver <- "R500S100/"

tib_name <- "/compare_tibble.rds"

# bayes_tf <- paste0(home_dir, bayes_dir, tib_name)
bayes_tf <- paste0(bayesian_dirs, tib_name)
cc_tf <- paste0(home_dir, cc_dir, cc_ver, tib_name)

# bayes_tf_2 <- paste0(home_dir, bayes_dir_2, tib_name)
# bayes_tf_3 <- paste0(home_dir, bayes_dir_3, tib_name)
# bayes_tf_4 <- paste0(home_dir, bayes_dir_4, tib_name)

# Save results
go_dir <- "Go_enrichment"

# save_dir <- paste0(c(paste0(home_dir, bayes_dir), paste0(home_dir, cc_dir, cc_ver)), go_dir)
#
# for (d in save_dir) {
#   dir.create(d)
# }


# Read in tibble containing PSM, predicted clustering, etc.
cc_tib <- readRDS(cc_tf)
bayes_tibs <- lapply(bayes_tf, readRDS)

datasets <- cc_tib$dataset
dataset_names <- c("Timecourse", "ChIP-ChIP", "PPI")
ontologies <- c("MF", "BP", "CC", "ALL")

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
  values = row.names(timecourse_data),
  # attributes = c(attributes_to_use, "entrezgene_id"),
  # filters = "entrezgene_id",
  # values = eg[match(row.names(timecourse_data), geneTable$Locus)],
  mart = yeast_ensemble
)

# buld GO map saves a bunch of .rda files to the working directory
bgoMap <- buildGOmap(gomap)


# head(geneTable)
genes_present <- geneTable$GeneID[match(row.names(timecourse_data), geneTable$Locus)]


# === Dataset choice ===========================================================

for (dataset_used in datasets) {
  curr_ind <- match(dataset_used, datasets)
  dataset_name <- dataset_names[curr_ind]

  save_dir <- paste0("./Images/Yeast/", dataset_name, "/")
  dir.create(save_dir)

  # ont <- "ALL"

  for (ont in ontologies) {

    # === Consensus clustering =====================================================

    cc_cl <- mcclust::maxpear(cc_tib$similarity_matrix[[curr_ind]], max.k = 250)$cl

    labels_present <- unique(cc_cl)
    n_labels <- length(labels_present)

    drop_na <- row.names(timecourse_data)[missing]
    universe <- geneTable$GeneID[match(row.names(timecourse_data), geneTable$Locus)]
    # ego <- goOverRep(cc_cl, timecourse_data, geneTable, universe = universe,
    #                  drop_na = drop_na, ont = ont, save_dir = NULL)

    # Save annotation data.frame
    # write.csv(cc_go, paste0(save_dir[2], "/goEnrichmentCC.csv"))


    # Compare all clusters in a single plot
    p1 <- doClusterComparison(cc_cl, timecourse_data, geneTable, universe,
      ont = ont,
      drop_na = drop_na
    )
    p1 +
      labs(
        x = "Predicted cluster (consensus clustering)",
        y = "GO description",
        title = dataset_name,
        subtitle = "Enrichment across clusters"
      ) +
      scale_color_viridis_c(direction = -1)

    # ggsave(paste0(save_dir, "/goEnrichmentAllClustersCC.png"), p1, width = 14, height = 10)


    # === Bayesian inference =======================================================

    bayes_p <- lapply(bayes_tibs, function(x) {
      cl <- x$pred_allocation[[curr_ind]]
      p <- doClusterComparison(cl, timecourse_data, geneTable, universe,
        ont = ont,
        drop_na = drop_na
      )
      p
    })

    for (i in 1:10) {
      curr_data <- bayes_p[[i]]$data
      curr_data$Inference <- "Bayesian"
      curr_data$Chain <- i
      curr_data$Model <- paste0("Bayesian: chain ", i)

      if (i == 1) {
        bayes_data <- curr_data
      } else {
        bayes_data <- rbind(bayes_data, curr_data)
      }
    }

    # === Comparison ===============================================================

    cc_data <- p1$data
    cc_data$Inference <- "Consensus clustering"
    cc_data$Chain <- NA
    cc_data$Model <- "Consensus Clustering"

    # Bind the data
    new_plt_data <- rbind(cc_data, bayes_data)
    new_plt_data$Cluster_plt <- stringr::str_extract(new_plt_data$Cluster_str, "\\d+") %>% as.numeric()
    new_plt_data$Model <- factor(new_plt_data$Model, levels = unique(new_plt_data$Model))
    new_plt_data$Ontology <- ont

    # Plot
    p_comp <- new_plt_data %>%
      ggplot(aes(x = Cluster_plt, y = Description, color = log(p.adjust), size = Count)) +
      geom_point() +
      facet_wrap(~Model) +
      scale_color_viridis_c(direction = -1) +
      labs(
        title = paste0(dataset_name, ": GO set over-representation (", ont, ")"),
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
    p_comp

    p_comp_2 <- p_comp$data %>%
      group_by(Description, Model) %>% 
      mutate(Number_cluster = as.character(n())) %>% 
      ggplot(aes(x = Description, y = Model)) +
      geom_point(aes(colour = log(p.adjust), size = Count)) + #, position = position_stack(reverse = TRUE)) +
      geom_text(aes(label = Number_cluster), position = position_nudge(y = -0.3)) +
      # geom_jitter(width = 0, height = 0.4, seed = 1) +
      scale_color_viridis_c(direction = -1) +
      labs(
        title = paste0(dataset_name, ": GO set over-representation (", ont, ")"),
        # subtitle = "Across chains",
        # x = "Descr",
        colour = "Log adjusted p-value"
      ) +
      theme(
        axis.text.y = element_text(hjust = 0.0, size = 8),
        axis.text.x = element_text(angle = 30, size = 8, vjust = 1, hjust=1),
        axis.title.y = element_blank(),
        # axis.title.x=element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text.x = element_text(size = 9)
      )
    
    
    ggsave(paste0(save_dir, "GOoverRepresentationComparison", ont, ".png"),
      p_comp,
      width = 14,
      height = 14
    )
    
    ggsave(paste0(save_dir, "GOoverRepresentationComparison", ont, "Alt.png"),
           p_comp_2,
           width = 14,
           height = 8
    )

    write.csv(p_comp$data, paste0(save_dir, "GOoverRepresentationComparison", ont, ".csv"))

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
}
