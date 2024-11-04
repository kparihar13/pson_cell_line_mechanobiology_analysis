# Aim: To plot the heatmap showing the cluster labels for each cell-sub pair 
# and the KDEs for each cluster identified by the consensus clustering

# Load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(pracma)
})

# set working directory same as the current script location -------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# global parameters --------
features <- c("area", "circularity", "aspect_ratio",
              "cell_stiffness", "motility")

# set the cutoff for min. number of observations needed for 
# a cell-substrate pair to be included in the analysis
min_cutoff <- 25

# define tissue respective colors
cell.line.label.colors = read.csv("../Figures/cell_line_label_colors.csv", 
                                  stringsAsFactors = FALSE)  %>%
  # change cell_id to factor with levels as the current order of cell lines
  # done to ensure the order of cell lines in the plots
  mutate(cl_id = factor(cl_id, levels = unique(cl_id)))

# color for clusters based on the number of clusters
cluster_colors <- list(c("#0000ff","#ff0000"), 
                       c("#0000ff","#000000","#ff0000"), 
                       c("#0000ff","#006400","#000000","#ff0000"),
                       c("#0000ff","#006400","orange4","#000000","#ff0000"),
                       c("#0000ff","magenta","#006400","orange4","#000000","#ff0000"))

# metric of interest
# "wass1" (1-Wasserstein) or "kolm_smir" (Kolmogrov-Smirnov)
metric_of_interest <- "wass1"

# load the data --------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)  
for (f in features){
  data[[f]] <- read_tsv(paste('../data/',f,'.tsv',sep = ''), 
                        show_col_types = FALSE) %>% 
    as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>% 
    mutate_at(c('cl_id', 'sub_id'), as.factor)
}

# remove cases with less than minimum number of data points -------
# find the cell-subs pairs with data points < min_cutoff
data_overview <- setNames(vector("list", length(features)), features)
for (f in features) {
  data_overview[[f]] <- data[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>% 
    group_by(cl_sub) %>% 
    summarise(totalcount = n(), .groups = "drop") %>% 
    as_tibble() %>% 
    filter(totalcount < min_cutoff)
}
# remove the cell-subs pairs with data points < min_cutoff
for (f in features) {
  data[[f]] <- data[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>%
    filter(!(cl_sub %in% data_overview[[f]]$cl_sub)) 
}


# load consensus cluster .RData object and merge labels with data --------

if (metric_of_interest == "wass1") {
  load("wass1_consensus_clustering.RData")
} else {
  load("kolm_smir_consensus_clustering.RData")
}

# optimal number of clusters determined based on PAC and CHI
if (metric_of_interest == "wass1") {
  k.clusters <- list("area" = 5, "circularity" = 3, "aspect_ratio" = 5,
                     "cell_stiffness" = 4, "motility" = 2)
} else {
  k.clusters <- list("area" = 4, "circularity" = 6, "aspect_ratio" = 4,
                     "cell_stiffness" = 3, "motility" = 2)
}

# to store cluster mapping for each cell-sub pair
cluster.map <- setNames(vector("list", length(features)), features)
# to store the boundary cases
boundary.cases <- setNames(vector("list", length(features)), features)
for (f in features){
  # get the cluster ids for each cell line-substrate pair
  cluster.map[[f]] <- ccout[[f]][[k.clusters[[f]]]]$consensusClass %>% 
    as_tibble()
  colnames(cluster.map[[f]]) <- c("label")
  # add the cell line and substrate names
  cluster.map[[f]] <- cluster.map[[f]] %>% 
    add_column(., cl_sub = names(ccout[[f]][[k.clusters[[f]]]]$consensusClass)) %>%
    separate_wider_delim(., cols = cl_sub, delim = "_", names = c("cl_id", "sub_id"))
  
  # replace all "." with "-" in cell line names
  cluster.map[[f]]$cl_id <- gsub("\\.", "-", cluster.map[[f]]$cl_id) 
  # replace all "." with " " in substrate names
  cluster.map[[f]]$sub_id <- gsub("\\.", " ", cluster.map[[f]]$sub_id)
  # for some reason there is an "X" before 22Rv1, remove it 
  cluster.map[[f]]$cl_id <- gsub("X22Rv1", "22Rv1", cluster.map[[f]]$cl_id)
  
  # Reorder columns: cell and subs ids first, followed by cluster label
  cluster.map[[f]] <- cluster.map[[f]][,c(2:3,1)]
  
  # add cluster labels to the data
  data[[f]] <- data[[f]] %>% 
    inner_join(., cluster.map[[f]], by = c("cl_id", "sub_id"))


  # note the two possible clusters for each boundary case in weak.cases dataframe
  # boundary cases: cell line-substrate pairs with itemConsensus < 0.8 in the 
  # cluster assigned based heirarchical clustering of the consensus matrix
  label1 <- c()
  label2 <- c()
  # iterate through the boundary cases
  for (i in 1:nrow(weak.cases[[f]])){
    # label1: cluster with max itemConsensus value (though it is < 0.8)
    label1 <- c(label1, order(as.numeric(weak.cases[[f]][i,2:(k.clusters[[f]]+1)]), 
                              decreasing = TRUE)[1])
    # label2: cluster with 2nd highest itemConsensus value
    label2 <- c(label2, order(as.numeric(weak.cases[[f]][i,2:(k.clusters[[f]]+1)]), 
                              decreasing = TRUE)[2])
  }
  boundary.cases[[f]] <- weak.cases[[f]] %>% 
    cbind(label1) %>% 
    cbind(label2) %>%
    separate(item, into = c("cl_id","sub_id"), sep="_", remove = TRUE)
  
  # replace all "." with "-" in cell line names
  boundary.cases[[f]]$cl_id <- gsub("\\.", "-", boundary.cases[[f]]$cl_id) 
  # replace all "." with " " in substrate names
  boundary.cases[[f]]$sub_id <- gsub("\\.", " ", boundary.cases[[f]]$sub_id)
  # for some reason there is an "X" before 22Rv1, remove it 
  boundary.cases[[f]]$cl_id <- gsub("X22Rv1", "22Rv1", boundary.cases[[f]]$cl_id)
  
  boundary.cases[[f]] <- boundary.cases[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep=""))
  
  # remove the clusterid (which was allotted to it based on heirarichical clustering 
  # of the consensus matrix) column, which corresponds to the highest itemConsensus value
  # for that cell-substrate pair
  # Again note that this is a boundary case with highest itemConsensus < 0.8
  boundary.cases[[f]] <- boundary.cases[[f]][,-(k.clusters[[f]]+3)]
}

# Determine cluster specific median values -----------
median.clusters <- setNames(vector("list", length(features)), features)
for (f in features){
  median.clusters[[f]] = data[[f]] %>% 
    # remove the boundary cases
    filter(!(cl_sub %in% boundary.cases[[f]]$cl_sub)) %>%
    # group by the cluster label
    group_by(label) %>% 
    # compute the median for each cluster
    summarise(median_value = median(feature_value)) %>% 
    ungroup() %>% 
    # arrange the clusters in descending order of median value
    arrange(median_value) %>% 
    # allot color to the clusters
    mutate(colors = cluster_colors[[k.clusters[[f]] - 1]]) %>% 
    # re-arrange back based on the cluster label (1,2,3,4,5)
    # Note that the cluster labels (1, 2, 3, 4, 5) from consensus clustering 
    # do not necessarily correspond to the order of the clusters based on median values 
    arrange(label)
}

# Compute substrate-specific cell line median values ----------
medians <- setNames(vector("list", length(features)), features) 
median.decimals <- list("area" = 0, "circularity" = 2, "aspect_ratio" = 2, 
                        "cell_stiffness" = 0, "motility" = 1)

for (f in features) {
  temp <- list()
  # get substrate specifc median values for each cell line
  for (subs in unique(data[[f]]$sub_id)){
    temp[[subs]] <- data[[f]] %>% 
      filter(sub_id == subs) %>%
      group_by(cl_id) %>% 
      summarise(median = round(median(feature_value),
                               median.decimals[[f]]), 
                .groups = "drop_last") %>% 
      select(cl_id,median) %>%
      ungroup() 
  }
  
  # combine all the substrates into a single dataframe
  medians[[f]] <- bind_rows(temp, .id = "sub_id") %>% 
    pivot_wider(names_from = sub_id, values_from = median) %>% 
    # Need to order the cell lines according the order in 
    # in cell.line.label.colors for consistency in plotting
    mutate(cl_id = factor(cl_id, levels = cell.line.label.colors$cl_id)) %>%
    arrange(cl_id) %>%
    column_to_rownames("cl_id")
}

# Heatmaps showing cluster label for cell-substrate pairs -----------
for (f in features){
  # keep only the cell lines for which median values are available
  # in case any of the features does not have all the 30 cell lines
  cell.line.label.colors.feature <- cell.line.label.colors %>% 
    filter(cl_id %in% unique(cluster.map[[f]]$cl_id))
  
  cluster.map.pivot <- cluster.map[[f]] %>% 
    pivot_wider(names_from = sub_id, values_from = label)
  
  cluster.map.pivot <- cluster.map.pivot %>% 
    # order in desired order of cell lines
    mutate(cl_id = factor(cl_id, levels = cell.line.label.colors.feature$cl_id)) %>%
    arrange(cl_id) %>%
    column_to_rownames('cl_id')
  
  # get the colors for the clusters
  colors <- structure(median.clusters[[f]]$colors, 
                      names = as.character(median.clusters[[f]]$label))
  
  # Get the indices (i,j) for boundary.cases in cluster.map.pivot
  # column for storing row index
  boundary.cases[[f]]$row_index <- rep(0,nrow(boundary.cases[[f]]))
  # column for storing column index
  boundary.cases[[f]]$col_index <- rep(0,nrow(boundary.cases[[f]]))
  # looping through the boundary.cases matrix
  for (i in 1:nrow(boundary.cases[[f]])) {
    boundary.cases[[f]][i,"row_index"] <- which(rownames(cluster.map.pivot) == 
                                                  boundary.cases[[f]][i,"cl_id"])
    boundary.cases[[f]][i,"col_index"] <- which(colnames(cluster.map.pivot) == 
                                                  boundary.cases[[f]][i,"sub_id"])
  }
  
  
  ht = Heatmap(as.matrix(cluster.map.pivot),
               col = colors,
               na_col = "white",
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               row_dend_reorder = FALSE,
               column_dend_reorder = FALSE,
               show_row_dend = FALSE,
               show_column_dend = FALSE,
               row_names_side = 'left',
               column_names_side = 'top',
               column_names_rot = 45,
               row_labels = cell.line.label.colors.feature$label, 
               column_labels = colnames(cluster.map.pivot),
               row_names_gp = gpar(fontsize = 10, 
                                   col = cell.line.label.colors.feature$tissue_col),
               column_names_gp = gpar(fontsize = 10),
               # !!!
               # uncomment if PC-3 is removed from the motility data
               # row_split = ifelse((f != "motility"), 
               #                    rep(1:8, c(4,5,3,3,3,4,6,2)),
               #                    rep(1:8, c(4,4,3,3,3,4,6,2))), 
               # splitting rows in terms of tissue types for aesthetics
               row_split = rep(1:8, c(4,5,3,3,3,4,6,2)),
               row_title = NULL, 
               row_gap = unit(1.2, "mm"),
               width = unit(9, "cm"), 
               height = unit(18, "cm"),
               rect_gp = gpar(col = "white", lwd = 2),
               # adapted from ComplexHeatmap documentation
               cell_fun = function(j, i, x, y, w, h, col) { 
                 # if (i,j) are boundary.cases then use 2 grid.rect():
                 # 1st half of the width and color corresponding to label1
                 # 2nd half of the width and color corresponding to label2 
                 for (k in 1:nrow(boundary.cases[[f]])) {
                   if (i == boundary.cases[[f]][k,"row_index"] & 
                       j == boundary.cases[[f]][k,"col_index"]) {
                     # Set the color of half the cell based on label1 
                     grid.rect(x = x-w/4, y = y, width = w/2, height = h, 
                               gp = gpar(col = "white", lwd = 2,
                                         fill = colors[[boundary.cases[[f]][k,"label1"]]]))
                     # Set the color of remaining half the cell based on label2 
                     grid.rect(x = x+w/4, y = y, width = w/2, height = h, 
                               gp = gpar(col = "white", lwd = 2,
                                         fill = colors[[boundary.cases[[f]][k,"label2"]]]))
                   }
                 }
                 
                 # add median values to each cell
                 grid.text(as.matrix(medians[[f]])[i, j], x, y,
                           gp = gpar(fontsize = 10, col = 'white')) 
               },
               show_heatmap_legend = FALSE)  
  
  # save the heatmaps
  if (metric_of_interest == "wass1") {
    if (f != "circularity"){
      png(paste0("../Figures/Figure5/", f,"_cluster_heatmap.png"), 
          res = 300, width = 1500, height = 2500, units = "px")
    } else {
      png(paste0("../Figures/Supplementary_Figure8/", f,"_cluster_heatmap.png"), 
          res = 300, width = 1500, height = 2500, units = "px")
    }
  } else {
    if (f %in% c("area", "cell_stiffness")){
      png(paste0("../Figures/Supplementary_Figure9/", f,"_cluster_heatmap.png"), 
          res = 300, width = 1500, height = 2500, units = "px")
    } else {
      png(paste0(f, "/", metric_of_interest,"_consensus_clustering/cluster_heatmap.png"), 
          res = 300, width = 1500, height = 2500, units = "px")
    }
  }

  draw(ht)
  dev.off()
}

# helper function to compute KDEs ------------
# Adapted from
# https://thirdorderscientist.org/homoclinic-orbit/2013/10/24/kernel-density-estimation-for-random-variables-with-bounded-support-mdash-the-transformation-trick
# https://medium.com/mlearning-ai/density-estimation-for-bounded-variables-7d68f633e772
# https://github.com/Aurelien-Pelissier/Medium/blob/main/Density%20estimation%20for%20bounded%20variables/Density_estimation.py

# list, numeric, numeric, numeric -> dataframe
# produces kernel density estimate for the data in list between range (from, to)
# and eps (numeric) is used to stabilize the KDE at the boundary value 0
KDE_func <- function(data, from, to, eps){

  # eps: area, cell_stiffness -> 0; motility -> 5; aspect_ratio -> -0.9
  # for stability at boundary value 0
  # minimum value is pretty far from 0 for area and cell stiffness
  data = data + eps
  # transform the data using log(x) from (0,Inf) -> (-Inf, Inf)
  x.transformed = log(data) 
  # compute the kernel density estimate in log scale
  density.transformed = density(x.transformed, 
                                from = log(from), 
                                to = log(to), 
                                n = 2^16)
  
  # Back transform to original scale
  x.vals = exp(density.transformed$x)
  density.vals = density.transformed$y/x.vals
  
  # return dataframe containing the KDE
  return(data.frame(x = x.vals - eps, y = density.vals)) 
}

# Vector, numeric, numeric -> dataframe
# produces kernel density estimate for the data in list between range (from, to)
# logit(x) transformation is used for circularity as it is bounded between 0 and 1
KDE_circularity <- function(data, from, to){
  
  # for boundary data of value 1, to ensure it isn't removed during logit transformation
  data = data - 0.05
  # for stability at boundary value 0; 
  # only 6 cases out of 20563 had value < 0.06 which will be removed due to following condition
  data = data[which(data > 0.01)]
  
  # transform the data using logit(x) from (0,1] -> (-Inf, Inf)
  # logit(x) = log(x/(1-x))
  x.transformed = logit(data) 
  # compute the kernel density estimate in logit scale
  density.transformed = density(x.transformed, 
                                from = logit(from), 
                                to = logit(to), 
                                n = 2^16)
  
  # Back transform to original scale
  x.vals = exp(density.transformed$x)/(1 + exp(density.transformed$x))
  density.vals = density.transformed$y/(x.vals*(1-x.vals))
  
  # return dataframe containing the KDE
  return(data.frame(x = x.vals + 0.05, y = density.vals)) 
}

# Plot KDEs for each cluster ----------
from_values <- list("area" = 10, "aspect_ratio" = 0.1, 
                    "cell_stiffness" = 100, "motility" = 5)
# negative for aspect ratio to approx shift from [1, Inf) -> (0, Inf)
eps <- list("area" = 0, "aspect_ratio" = -0.9, 
            "cell_stiffness" = 0, "motility" = 5)

for (f in features){
  KDE <- list()
  for (i in 1:k.clusters[[f]]){
    
    temp <- data[[f]] %>% 
      # remove the boundary.cases
      filter(!(cl_sub %in% boundary.cases[[f]]$cl_sub)) %>%
      # keep only the cluster of interest
      filter(label == i)
    
    if (f == "circularity") 
      { KDE_temp <- KDE_circularity(temp$feature_value, from = 1e-3, to = 0.95) }
    else { 
      KDE_temp <- KDE_func(temp$feature_value, 
                           from = from_values[[f]], 
                           to = max(data[[f]]$feature_value) + eps[[f]], 
                           eps[[f]]) }
    
    KDE[[i]] = data.frame(x = KDE_temp$x, y = KDE_temp$y) 
  }
  
  # save the KDEs
  if (metric_of_interest == "wass1") {
    if (f != "circularity"){
      png(paste0("../Figures/Figure5/", f,"_cluster_KDEs.png"), 
          res = 300, width = 1800, height = 1200, units = "px")
    } else {
      png(paste0("../Figures/Supplementary_Figure8/", f,"_cluster_KDEs.png"), 
          res = 300, width = 1800, height = 1200, units = "px")
    }
  } else {
    png(paste0(f, "/", metric_of_interest,"_consensus_clustering/cluster_KDEs.png"), 
        res = 300, width = 1800, height = 1200, units = "px")
  }

  KDE_labels <- ggplot() + theme_classic()
  
  # Plotting based in order of colors 
  for (c in cluster_colors[[k.clusters[[f]]-1]]){
    KDE_labels <- KDE_labels + 
      geom_line(data = KDE[[which(median.clusters[[f]]$colors == c)]], 
                aes(x = x, y = y), 
                col = c) +
      geom_area(data = KDE[[which(median.clusters[[f]]$colors == c)]], 
                aes(x = x, y = y), 
                fill = c, 
                alpha = 0.5) 
  }
  
  # Add axis labels
  if (f == "area") 
  { KDE_labels <- KDE_labels + 
      xlim(0,4000) + 
    labs(x = expression("Area (" * mu*"m)"^2), 
         y = 'density') }
  else if (f == "circularity") 
  {KDE_labels <- KDE_labels + 
      labs(x = "Circularity", 
           y = 'density')} 
  else if (f == "aspect_ratio") 
  {KDE_labels <- KDE_labels + 
      xlim(1,4) + 
    labs(x = "Aspect Ratio", 
         y = 'density')} 
  else if (f == "cell_stiffness") 
  {KDE_labels <- KDE_labels + 
      xlim(0,18000) + 
    labs(x = "Cell Stiffness (Pa)", 
         y = 'density') } 
  else if(f == "motility") 
  {KDE_labels <- KDE_labels + 
      xlim(0,120) + 
    labs(x = expression("Speed (" * mu*"m/hr)"), 
         y = 'density')}
  
  print(KDE_labels)
  dev.off()
}
