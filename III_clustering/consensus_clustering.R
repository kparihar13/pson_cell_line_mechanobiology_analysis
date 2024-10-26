# Aim: To perform consensus clustering on the ECDF-based distance matrix of the features

# load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(fpc)
  library(ConsensusClusterPlus)
})

# set the working directory same as the location of the script -----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# global parameters -------
features <- c("area", "circularity", "aspect_ratio",
              "cell_stiffness", "motility")
# define the metric to be used for clustering
# "wass1" (1-Wasserstein) or "koln_simr" (Kolmogrov-Smirnov)
metric_of_interest <- "wass1"

# load distance matrix -----------
metric <- setNames(vector("list", length(features)), features)
for (f in features) {
  metric[[f]] <- read.csv(paste(f, "/", metric_of_interest, ".csv", sep=""))
}

# Helper function for clustering -------
# ConsensusClusteringPlus output object, distance object, 
# max number of clusters K -> ggplots
# produces line plots for cluster statistics (PAC, CHI) vs number of clusters
cluster.statistics <- function(ccout, metric.dist, K){
  # named list to store cluster statistics 
  cluster.metrics <- setNames(vector("list", 3), 
                              c("k","PAC","CHI"))
  
  # total number of cell line-substrate pairs
  N <- nrow(as.matrix(metric.dist))
  
  for (k in 2:K) {
    cluster.metrics[["k"]] <- c(cluster.metrics[["k"]], k)
    
    # PAC calculation based on definition in
    # https://www.nature.com/articles/s41598-020-58766-1#Sec11
    CDF01 <- sum(ccout[[k]]$consensusMatrix[which(lower.tri(ccout[[k]]$consensusMatrix))] <= 0.1)/(N*(N - 1)/2)
    CDF09 <- sum(ccout[[k]]$consensusMatrix[which(lower.tri(ccout[[k]]$consensusMatrix))] <= 0.9)/(N*(N - 1)/2)
    cluster.metrics[["PAC"]] <- c(cluster.metrics[["PAC"]], 
                                  round(CDF09 - CDF01, 3))

    # CHI
    cluster.metrics[["CHI"]] <- c(cluster.metrics[["CHI"]], 
                                  round(cluster.stats(metric.dist, 
                                                      ccout[[k]]$consensusClass)$ch, 2))
    
  }

  cluster.metric <- cluster.metrics %>% 
    as_tibble()
  
  plt.pac <- cluster.metric %>% 
    ggplot(aes(x = k, y = PAC)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of Clusters",
         y = "Proportion of ambigous clustering (PAC)") +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))

  plt.chi <- cluster.metric %>%
    ggplot(aes(x = k, y = CHI)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of Clusters",
         y = "Calinski-Harabasz Index (CHI)") +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  
  return(list(plt.pac, plt.chi))
}

# Consensus based clustering ----------
# Partitioning Around Medoids (PAM) clustering: for subsampled set of distance matrix
# pairwise consensus values, the proportion that two items occupied
# the same cluster out of the number of times they occurred in the same subsample,
# are calculated and stored in a symmetrical consensus matrix for each k.
# Heirarichical clustering is then performed on the consensus matrix

# to store the results from consensus clustering
ccout <- setNames(vector("list", length(features)), features)
icl <- setNames(vector("list", length(features)), features)

# iterate through each feature
for (f in features) {
  
  print(paste0("Performing consensus clustering on ", f))
  
  # folder where plots produced internally by ConsensusClusterPlus will be saved
  folder <- paste0(f, "/", metric_of_interest, "_consensus_clustering")
  
  # ConsensusClusterPlus: determining cluster number and class membership by stability evidence
  ccout[[f]] <- ConsensusClusterPlus(d = as.dist(metric[[f]]), 
                                     # max number to clusters to evaluate
                                     maxK = 9, 
                                     # number of times to subsample
                                     reps = 1000,
                                     # proportion of items to sample
                                     # i.e. 80% of the cell line-substrate pairs are sampled
                                     # for creating each subsample set
                                     pItem = 0.8,
                                     # clustering algorithm to use on the 
                                     # subsampled distance matriux
                                     clusterAlg = "pam",
                                     # seed for reproducibility
                                     seed = 10,
                                     # folder to save plots 
                                     title = folder,
                                     # file type for saving plots
                                     plot = "png")
  
  # calculate cluster-consensus and item-consensus
  icl[[f]] <- calcICL(ccout[[f]])
  
  # plot cluster statistics, PAC and CHI
  plts <- cluster.statistics(ccout[[f]], as.dist(metric[[f]]), 9)
  
  # save the plots
  if (metric_of_interest == "wass1")
    ggsave(paste0("../Figures/Supplementary_Figures6_and_7/", f, "_pac_chi.png"), 
           plot = plot_grid(plotlist = plts, ncol = 2, nrow = 1), 
           dpi = 300, width = 2000, height = 1000, units = "px")
  else
    ggsave(paste0("../Figures/Supplementary_Figure8/", f, "_pac_chi.png"), 
           plot = plot_grid(plotlist = plts, ncol = 2, nrow = 1), 
           dpi = 300, width = 2000, height = 1000, units = "px")
}

# save ccout and icl in <feature_name>.RData
save(ccout, icl, file = paste0(metric_of_interest, "_consensus_clustering.RData"))

