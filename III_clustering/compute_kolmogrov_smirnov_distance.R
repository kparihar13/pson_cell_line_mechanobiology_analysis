# Aim: Compute the pariwise Kolmogrov-Smirnov (KS) distance between the empirical CDFs of
# feature each cell-substrate pair.

# Load packages --------------
suppressPackageStartupMessages({
  library(tidyverse)
  # geostats is a package for computing the 1-Wasserstein distance
  library(geostats)
})

# set the working directory same as the script location ------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# global parameters ------
features <- c("area", "circularity", "aspect_ratio",
              "cell_stiffness", "motility")

# set the cutoff for min. number of observations needed for 
# a cell-substrate pair to be included in the analysis
min_cutoff <- 25

# load the data --------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)  
for (f in features){
  data[[f]] <- read_tsv(paste('../../data/',f,'.tsv',sep = ''), 
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

for (f in features) {
  data[[f]] <- data[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>%
    filter(!(cl_sub %in% data_overview[[f]]$cl_sub)) 
}

# create list of cell-substrate pair data vectors for KS calculations ---------
data.list <- setNames(vector("list", length(features)), features)
for (f in features) {
  # for each feature, create a list of data vectors for each cell-sub pair
  data.list[[f]] <- setNames(vector("list", length(unique(data[[f]]$cl_sub))), 
                             unique(data[[f]]$cl_sub))
  # interate through all cell-sub pairs
  for (cs in unique(data[[f]]$cl_sub)){
    data.list[[f]][[cs]] <- filter(data[[f]], cl_sub == cs)$feature_value
  }
}

# compute pairwise KS distance -------
for (f in 1:length(features)) {
  ks <- as.matrix(ksdist(data.list[[f]]))
  # save the KS distance matrix to a csv file
  write.csv(as.data.frame(ks), paste(features[[f]],".csv",sep = ''), 
            row.names = FALSE)
}


