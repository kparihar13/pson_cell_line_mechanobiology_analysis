# Aim: To perform correlation analysis correlation analysis using median values 

# load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
  library(Hmisc)
})

# set the working directory as the location of the script -------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# global parameters -------
features <- c("area", "circularity", "aspect_ratio",
              "cell_stiffness", "motility")

# load the data --------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)  
for (f in features){
  data[[f]] <- read_tsv(paste('../data/',f,'.tsv',sep = ''), 
                        show_col_types = FALSE) %>% as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>% 
    mutate_at(c('cl_id', 'sub_id'), as.factor)
}

# remove cases with less than minimum number of data points -------
# set the cutoff for min. number of data points required for each cell-subs pair
min_cutoff <- 25  

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
    filter(!(cl_sub %in% data_overview[[f]]$cl_sub)) %>% 
    select(-cl_sub)
}

# median value of features for each cell line-substrate pair -------
medians <- setNames(vector("list", length(features)), features) 
for (f in features){
  medians[[f]] <- data[[f]] %>% 
    group_by(cl_id, sub_id) %>% 
    summarise(!!f := median(feature_value), .groups = "drop") %>% 
    as_tibble()
}

# collect the medians for each feature into one dataframe
# first define a dataframe with cl_id and sub_id columns
# such that each element in unique(data[[1]]$cl_id) x unique(data[[1]]$sub_id) 
# is present. This is done to ensure that all the cell line-substrate pairs 
# are present in the final dataframe even if few of the cell line-substrate pairs
# may have been removed above due to less number of data points (<25)
all.features <- expand.grid(unique(data[[1]]$cl_id), unique(data[[1]]$sub_id))
colnames(all.features) <- c("cl_id", "sub_id")

for (f in features) {
  all.features <- all.features %>% 
    left_join(.,medians[[f]], by = c("cl_id", "sub_id"))
}

colnames(all.features)[3:7] <- c("Area","Circularity","Aspect Ratio",
                                   "Cell Stiffness","Cell Speed")

# correlation analysis using median values -------
corr.spearman <- rcorr(as.matrix(all.features %>% 
                                 select(-all_of(c("cl_id","sub_id")))), 
                     type = "spearman")

# plot the correlation matrix -------
png(paste("../Figures/supplementary_figure1.png",sep=""),
    res = 300, width = 1200, height = 1200)

plt = corrplot(corr.spearman$r,
  p.mat = corr.spearman$P,
  method = "color",
  type = "lower",
  insig = "blank",
  addCoef.col = "white",
  addgrid.col = "black",
  diag = FALSE,
  tl.col = "black",
  cl.align = "l",
  cl.cex = 1,
  cl.length = 5,
  col = colorRampPalette(c("blue", "white", "red"))(100)
)

dev.off()
