# Packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(pracma)
})

# load the data --------
features <- c("area", "circularity", "aspect_ratio","stiffness", "motility")

data <- setNames(vector("list", length(features)), features)  # define an empty named list for storing data
for (f in features){
  data[[f]] <- read_tsv(paste('./data/',f,'_data.csv',sep = ''), show_col_types = FALSE) %>% as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>% mutate_at(c('cl_id', 'sub_id'), as.factor)
}

# remove cases with less than minimum number of data points -------
min_cutoff <- 25  # set the cutoff for min. # of data points required

# find the cell-subs pairs with data points < min_cutoff
data_overview <- setNames(vector("list", length(features)), features)
for (f in features) {
  data_overview[[f]] <- data[[f]] %>% mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>% group_by(cl_sub) %>% 
    summarise(totalcount = n(), .groups = "drop") %>% as_tibble() %>% filter(totalcount < min_cutoff)
}

for (f in features) {
  data[[f]] <- data[[f]] %>% mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>%
    filter(!(cl_sub %in% data_overview[[f]]$cl_sub)) #%>% select(-cl_sub)
}

# # create data frame to note which cell-subs pairs have (or do not have) the feature values --------
# cell_subs <- data[[1]] %>% group_by(cl_id,sub_id) %>% select(cl_id,sub_id) %>% unique() %>% add_column(!!features[1] := 1)
# for (f in features[2:5]){
#   cell_subs <- cell_subs %>% left_join(.,(data[[f]] %>% group_by(cl_id,sub_id) %>% 
#                                             select(cl_id,sub_id) %>% 
#                                             unique() %>% add_column(!!f := 1)), 
#                                        by = c("cl_id", "sub_id"))
# }

# define tissue and respective colors ------------
cell.lines.labels.color = read.csv("cell_lines_labels_color_updated.csv", stringsAsFactors = FALSE)
subs_short = list("500Pa Coll"= "500pacol","500Pa FN"= "500pafn",
                  "30kPa Coll"= "30kpacol", "30kPa FN"= "30kpafn",
                  "HA Coll" = "hacol", "HA FN"="hafn", "Glass" = "glass")
cluster_colors <- list(c("#0000ff","#ff0000"), c("#0000ff","#000000","#ff0000"), 
                       c("#0000ff","#006400","#000000","#ff0000"),
                       c("#0000ff","#006400","orange4","#000000","#ff0000"),
                       c("#0000ff","magenta","#006400","orange4","#000000","#ff0000"))

#folder <- "Figure2_wass1"
folder <- "Figure2_ks"

# load consensus cluster RData object and merge labels with data -------------
if (folder == "Figure2_wass1") {
  load("consensus_cluster_wass1_121723.RData")
} else if (folder == "Figure2_ks") {
  load("consensus_cluster_ks_090424.RData")
}

# number of determined based on PAC and CHI
if (folder == "Figure2_wass1") {
  k.clusters <- list("area" = 5, "circularity" = 3, "aspect_ratio" = 5,
                     "stiffness" = 4, "motility" = 2)
} else if (folder == "Figure2_ks") {
  k.clusters <- list("area" = 4, "circularity" = 6, "aspect_ratio" = 4,
                     "stiffness" = 3, "motility" = 2)
}

# define cluster mapping for each cell-sub pair
cluster.map <- setNames(vector("list", length(features)), features)
for (f in features){
  cluster.map[[f]] <- ccout[[f]][[k.clusters[[f]]]]$consensusClass %>% as_tibble()
  colnames(cluster.map[[f]]) <- c("label")
  cluster.map[[f]] <- cluster.map[[f]] %>% add_column(., cl_sub = names(ccout[[f]][[k.clusters[[f]]]]$consensusClass)) %>%
    separate_wider_delim(., cols = cl_sub, delim = "_", names = c("cl_id", "sub_id"))
  
  # replace all "." with "-"
  cluster.map[[f]]$cl_id <- gsub("\\.", "-", cluster.map[[f]]$cl_id) 
  # replace all "." with " "
  cluster.map[[f]]$sub_id <- gsub("\\.", " ", cluster.map[[f]]$sub_id)
  # for some reason there is an "X" before 22Rv1; remove it 
  cluster.map[[f]]$cl_id <- gsub("X22Rv1", "22Rv1", cluster.map[[f]]$cl_id)
  
  # putting cell and subs id columns first, followed by cluster label
  cluster.map[[f]] <- cluster.map[[f]][,c(2:3,1)]
  
  # merge cluster labels with the expt. data
  data[[f]] <- data[[f]] %>% inner_join(.,cluster.map[[f]], by = c("cl_id", "sub_id"))
}

# note the two possible clusters for each case in weak.cases dataframe
for (f in features){
  label1 <- c()
  label2 <- c()
  for (i in 1:nrow(weak.cases[[f]])){
    label1 <- c(label1, order(as.numeric(weak.cases[[f]][i,2:(k.clusters[[f]]+1)]), 
                              decreasing = TRUE)[1])
    label2 <- c(label2, order(as.numeric(weak.cases[[f]][i,2:(k.clusters[[f]]+1)]), 
                              decreasing = TRUE)[2])
  }
  weak.cases[[f]] <- weak.cases[[f]] %>% 
    cbind(label1) %>% cbind(label2) %>%
    # select(all_of(c("item","label1","label2"))) %>%
    separate(item, into = c("cl_id","sub_id"), sep="_", remove = TRUE)
  
  # replace all "." with "-"
  weak.cases[[f]]$cl_id <- gsub("\\.", "-", weak.cases[[f]]$cl_id) 
  # replace all "." with " "
  weak.cases[[f]]$sub_id <- gsub("\\.", " ", weak.cases[[f]]$sub_id)
  # for some reason there is an "X" before 22Rv1; remove it 
  weak.cases[[f]]$cl_id <- gsub("X22Rv1", "22Rv1", weak.cases[[f]]$cl_id)
  
  weak.cases[[f]] <- weak.cases[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep=""))
  
  # weak.cases[[f]] <- weak.cases[[f]][,c(5,1:4)]
  #weak.cases[[f]] <- weak.cases[[f]][,c(1:(k.clusters[[f]]+2),(k.clusters[[f]]+4),(k.clusters[[f]]+5))]
  weak.cases[[f]] <- weak.cases[[f]][,-(k.clusters[[f]]+3)]
}

# Determine cluster specific median values -----------
medians <- setNames(vector("list", length(features)), features)
for (f in features){
  medians[[f]] = data[[f]] %>% 
    ##group by the cluster label
    group_by(label) %>% 
    ##compute the median for each cluster
    summarise(median_value = median(feature_value)) %>% 
    ungroup() %>% 
    ##arrange the clusters in descending order of median value
    arrange(median_value) %>% 
    ##allot color to the clusters
    mutate(colors = cluster_colors[[k.clusters[[f]] - 1]]) %>% 
    ##re-arrange back based on the cluster label (1,2,3,4,5)
    arrange(label)
}

# Heatmaps showing cluster label for cell-sub pairs -----------

cell.line.order <- c("SK-MEL-2", "A375", "WM266-4", "MeWo",
                     "RWPE-1", "22Rv1", "LnCaP", "DU145", "PC-3",
                     "hTERT-HPNE", "Panc-1", "Capan-1",
                     "SKOV-3", "Caov-3", "OVCAR-3",
                     "NL20", "NCI-H2126", "NCI-H2087",
                     "HCT116", "HT29", "SW480", "SW620",
                     "hTERT-HME1", "MCF10A-JSB", "T-47D", "MCF7", "MDA-MB-231", "HCC1937",
                     "U-87", "T98G")

for (f in features){
  cluster.map.pivot <- cluster.map[[f]] %>% 
    pivot_wider(names_from = sub_id, values_from = label)
  
  cluster.map.pivot <- cluster.map.pivot %>% 
    # order in desired order of cell lines
    mutate(cl_id = factor(cl_id, levels = cell.line.order)) %>%
    arrange(cl_id) %>%
    column_to_rownames('cl_id')
  
  cell.lines.labels.color.feature <- left_join(data.frame(cl_id = rownames(cluster.map.pivot)),
                                               cell.lines.labels.color, by = 'cl_id')
  
  colors <- structure(medians[[f]]$colors, names = as.character(medians[[f]]$label))
  
  ##Get the indices (i,j) for weak.cases in cluster.map.pivot
  ##column for storing row index
  weak.cases[[f]]$row_index <- rep(0,nrow(weak.cases[[f]]))
  ##column for storing column index
  weak.cases[[f]]$col_index <- rep(0,nrow(weak.cases[[f]]))
  ##looping through the weak.cases matrix
  for (i in 1:nrow(weak.cases[[f]])) {
    weak.cases[[f]][i,"row_index"] <- which(rownames(cluster.map.pivot) == 
                                              weak.cases[[f]][i,"cl_id"])
    weak.cases[[f]][i,"col_index"] <- which(colnames(cluster.map.pivot) == 
                                              weak.cases[[f]][i,"sub_id"])
  }
  
  if (f != "motility"){
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
                 #show_column_names = FALSE,
                 column_names_rot = 45,
                 row_labels = cell.lines.labels.color.feature$label, 
                 column_labels = colnames(cluster.map.pivot),
                 row_names_gp = gpar(fontsize = 10, col = cell.lines.labels.color.feature$tissue_col),
                 column_names_gp = gpar(fontsize = 10),
                 row_split = rep(1:8, c(4,5,3,3,3,4,6,2)), 
                 row_title = NULL, row_gap = unit(1.2, "mm"),
                 width = unit(5, "cm"), height = unit(12, "cm"),
                 rect_gp = gpar(col = "white", lwd = 2),
                 ##Add cell_fun with 2 grid.rect() for weak.cases: 
                 ##1. for half the width and color corresponding to label1
                 ##2. for half the width and color corresponding to label2
                 cell_fun = function(j, i, x, y, w, h, col) { 
                   ##check if i and j are in the weak.cases matrix
                   for (k in 1:nrow(weak.cases[[f]])) {
                     if (i == weak.cases[[f]][k,"row_index"] & j == weak.cases[[f]][k,"col_index"]) {
                       ##Set the color of half the cell based on label1 
                       grid.rect(x = x-w/4, y = y, width = w/2, height = h, 
                                 gp = gpar(col = "white", lwd = 2,
                                           fill = colors[[weak.cases[[f]][k,"label1"]]]))
                       ##Set the color of remaining half the cell based on label2 
                       grid.rect(x = x+w/4, y = y, width = w/2, height = h, 
                                 gp = gpar(col = "white", lwd = 2,
                                           fill = colors[[weak.cases[[f]][k,"label2"]]]))
                     }
                   }
                 },
                 show_heatmap_legend = FALSE) } else {
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
                 #show_column_names = FALSE,
                 column_names_rot = 45,
                 row_labels = cell.lines.labels.color.feature$label, 
                 column_labels = colnames(cluster.map.pivot),
                 row_names_gp = gpar(fontsize = 10, col = cell.lines.labels.color.feature$tissue_col),
                 column_names_gp = gpar(fontsize = 10),
                 row_split = rep(1:8, c(4,4,3,3,3,4,6,2)), 
                 row_title = NULL, row_gap = unit(1.2, "mm"),
                 width = unit(5, "cm"), height = unit(12, "cm"),
                 rect_gp = gpar(col = "white", lwd = 2),
                 ##Add cell_fun with 2 grid.rect() for weak.cases: 
                 ##1. for half the width and color corresponding to label1
                 ##2. for half the width and color corresponding to label2
                 cell_fun = function(j, i, x, y, w, h, col) { 
                   ##check if i and j are in the weak.cases matrix
                   for (k in 1:nrow(weak.cases[[f]])) {
                     if (i == weak.cases[[f]][k,"row_index"] & j == weak.cases[[f]][k,"col_index"]) {
                       ##Set the color of half the cell based on label1 
                       grid.rect(x = x-w/4, y = y, width = w/2, height = h, 
                                 gp = gpar(col = "white", lwd = 2,
                                           fill = colors[[weak.cases[[f]][k,"label1"]]]))
                       ##Set the color of remaining half the cell based on label2 
                       grid.rect(x = x+w/4, y = y, width = w/2, height = h, 
                                 gp = gpar(col = "white", lwd = 2,
                                           fill = colors[[weak.cases[[f]][k,"label2"]]]))
                     }
                   }
                 },
                 show_heatmap_legend = FALSE) 
  }
  
  png(paste("./",f,"/",folder,"/heatmap_", 
            k.clusters[[f]],"clusters_weak_cases_MCF10A_order_changed.png",sep=''), 
      res = 300, width = 1500, height = 2000)
  draw(ht)
  dev.off()
}
rm(ht, temp)

# functions to compute KDEs ------------
##Based on https://thirdorderscientist.org/homoclinic-orbit/2013/10/24/kernel-density-estimation-for-random-variables-with-bounded-support-mdash-the-transformation-trick
##https://medium.com/mlearning-ai/density-estimation-for-bounded-variables-7d68f633e772
##https://github.com/Aurelien-Pelissier/Medium/blob/main/Density%20estimation%20for%20bounded%20variables/Density_estimation.py

##eps: area, stiffness -> 0; motility -> 5; aspect_ratio -> -0.9
KDE_func <- function(data, from, to, eps){
  ##transform the data using log to (0,Inf) -> (-Inf, Inf)
  
  ###No need to worry about boundary cases as minimum value is pretty far from 0 unlike aspect ratio
  data = data + eps
  x.transformed = log(data) 
  density.transformed = density(x.transformed, from = log(from), to = log(to), n = 2^16)
  
  ##Back transform
  x.vals = exp(density.transformed$x)
  density.vals = density.transformed$y/x.vals
  
  ##return dataframe
  return(data.frame(x = x.vals - eps, y = density.vals)) 
}

KDE_circularity <- function(data, from, to){
  ##transform the data using logit to (0,1] -> (-Inf, Inf)
  ##logit(x) = log(x/(1-x))
  
  ###for boundary data of value 1, to ensure it isn't removed during logit transformation
  data = data - 0.05
  ###for stability at boundary value 0; only 6 cases out of 20563 had value < 0.06 which will be removed due to following condition
  data = data[which(data > 0.01)]
  
  x.transformed = logit(data) 
  density.transformed = density(x.transformed, from = logit(from), to = logit(to), n = 2^16)
  
  ##Back transform
  x.vals = exp(density.transformed$x)/(1 + exp(density.transformed$x))
  density.vals = density.transformed$y/(x.vals*(1-x.vals))
  
  ##return dataframe
  return(data.frame(x = x.vals + 0.05, y = density.vals)) 
}

# Plot KDEs for each cluster -------------------
from_values <- list("area" = 10, "aspect_ratio" = 0.1, "stiffness" = 100, "motility" = 5)
eps <- list("area" = 0, "aspect_ratio" = -0.9, "stiffness" = 0, "motility" = 5)
for (f in features){
  KDE <- list()
  for (i in 1:k.clusters[[f]]){
    temp <- data[[f]] %>% filter(label == i)
    if (f == "circularity") { KDE_temp <- KDE_circularity(temp$feature_value, from = 1e-3, to = 0.95) }
    else { 
      KDE_temp <- KDE_func(temp$feature_value, from = from_values[[f]], to = max(data[[f]]$feature_value) + eps[[f]], eps[[f]]) }
    
    KDE[[i]] = data.frame(x = KDE_temp$x, y = KDE_temp$y) 
  }
  
  png(paste("./",f,"/",folder,"/kde_",k.clusters[[f]],"clusters.png",sep = ''), res = 300, width = 1800, height = 1200)
  KDE_labels <- ggplot() + theme_classic()
  
  ###Plotting based in order of colors
  for (c in cluster_colors[[k.clusters[[f]]-1]]){
    KDE_labels <- KDE_labels + geom_line(data = KDE[[which(medians[[f]]$colors == c)]], aes(x = x, y = y), col = c) +
      geom_area(data = KDE[[which(medians[[f]]$colors == c)]], aes(x = x, y = y), fill = c, alpha = 0.5)
  }
  
  
  if (f == "area") { KDE_labels <- KDE_labels + xlim(0,4000) + labs(x = expression("Area (" * mu*"m)"^2), y = 'density') }
  else if (f == "circularity") {KDE_labels <- KDE_labels + labs(x = "Circulrity", y = 'density')} 
  else if (f== "aspect_ratio") {KDE_labels <- KDE_labels + xlim(1,4) + labs(x = "Aspect Ratio", y = 'density')} 
  else if (f == "stiffness") {KDE_labels <- KDE_labels + xlim(0,18000) + labs(x = "Cell Stiffness (Pa)", y = 'density') } 
  else if(f == "motility") {KDE_labels <- KDE_labels + xlim(0,120) + labs(x = expression("Speed (" * mu*"m/hr)"), y = 'density')}
  
  print(KDE_labels)
  dev.off()
}
rm(KDE_temp)