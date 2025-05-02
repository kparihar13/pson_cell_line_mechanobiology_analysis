# Aim: to identify boundary cases in the consensus clustering results
# For a particular feature, boundary cases correspond to cell line-substrate pairs
# that have itemConsensus values < 0.8 for the cluster they are assigned to
# by heirarchical clustering of the consensus matrix

# load packages ---------
suppressPackageStartupMessages({
  library(tidyverse)
})

# set working directory to the location of this script -----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# global parameters ------
features <- c(
  "area", "aspect_ratio", "circularity",
  "cell_stiffness", "motility"
)
# define the metric to be used for clustering
# "wass1" (1-Wasserstein) or "kolm_smir" (Kolmogrov-Smirnov)
metric_of_interest <- "kolm_smir"

# load the results from consensus clustering -----
load(paste0(metric_of_interest, "_consensus_clustering.RData"))

# Note down the boundary cases for each feature -------
weak.cases <- setNames(vector("list", length(features)), features)
# optimal number of clusters determined based on PAC and CHI
if (metric_of_interest == "wass1") {
  k.clusters <- list(
    "area" = 5, "circularity" = 3, "aspect_ratio" = 5,
    "cell_stiffness" = 4, "motility" = 2
  )
} else {
  k.clusters <- list(
    "area" = 4, "circularity" = 6, "aspect_ratio" = 4,
    "cell_stiffness" = 3, "motility" = 2
  )
}

for (f in features) {
  # get the itemConsensus for each cell line-substrate pair in each of the clusters
  df.wider <- icl[[f]]$itemConsensus %>%
    filter(k == k.clusters[[f]]) %>%
    select(-k) %>%
    pivot_wider(names_from = cluster, values_from = itemConsensus)

  # get the cluster id for each cell line-substrate pair
  # ConsensusClusterPlus alloted it based on heirarchical clustering of the consensus matrix
  clusterid.df <- as.data.frame(ccout[[f]][[k.clusters[[f]]]]$consensusClass) %>%
    rownames_to_column(var = "item")
  colnames(clusterid.df) <- c("item", "clusterid")

  # add assigned cluster id to the wider dataframe containing itemConsensus values
  df.wider <- df.wider %>%
    inner_join(., clusterid.df, by = "item")

  # cases assigned to cluster 1 and have itemConsensus < 0.8
  weak.cases[[f]] <- df.wider %>%
    filter(clusterid == 1 & `1` < 0.8)
  # cases assigned to cluster i and have itemConsensus < 0.8
  for (i in 2:k.clusters[[f]]) {
    temp.df <- df.wider %>%
      filter(clusterid == i)
    if (sum(temp.df[, i + 1] < 0.8) > 0 &
      # also check if the cluster has more than one cell line pair
      nrow(temp.df) > 1) {
      weak.cases[[f]] <- rbind(
        weak.cases[[f]],
        temp.df[which(temp.df[, i + 1] < 0.8), ]
      )
    }
  }
}

save(ccout, icl, weak.cases,
  file = paste0(metric_of_interest, "_consensus_clustering.RData")
)
