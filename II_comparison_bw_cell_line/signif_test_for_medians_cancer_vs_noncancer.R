# Aim: For a particular tissue type, perform a permutation test 
# to compare the median values of the features between the cancer and non-cancer cell lines

# load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(rcompanion)
  library(boot)
  library(foreach)
  library(doParallel)
})
registerDoParallel(8)

# set the working directory same as the location of the script -----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Global parameters ------
features <- c("area", "circularity", "aspect_ratio", 
              "cell_stiffness", "motility")
subs <- c("HA Coll", "HA FN", "500Pa Coll", "500Pa FN", 
          "30kPa Coll", "30kPa FN", "Glass")

# define the tissue of interest
tissues_of_interest <- c("Prostate", "Pancreas", "Lung", "Breast")

# get the cell lines for the tissue of interest
cell_lines <-  read.csv("cell_line_label_colors.csv", 
                        stringsAsFactors = FALSE) %>%
  # keep only relevant columns
  select(cl_id, tissue) %>%
  # filter to keep only the cell lines of interest
  filter(tissue %in% tissues_of_interest)

# minimum cutoff for number of observations to be included in the analysis
min_cutoff <- 25

# load data -----------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)  
data_overview <- setNames(vector("list", length(features)), features)

for (f in features){
  data[[f]] <- read_tsv(paste('../data/',f,'.tsv',sep = ''), 
                        show_col_types = FALSE) %>% 
    as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>%
    # add column for cancer and non-cancer status
    mutate(status = ifelse(cl_id %in% c("RWPE-1","hTERT-HPNE","NL20",
                                        "hTERT-HME1","MCF10A-JSB"), 
                           "Non-cancer", "Cancer")) %>%
    # filter to only the cell lines of interest
    filter(cl_id %in% cell_lines$cl_id) %>%
    # fix the order of substrates for plotting purposes
    mutate(sub_id = fct_relevel(sub_id, subs)) %>%
    # fix the order for status
    mutate(status = fct_relevel(status, c("Non-cancer", "Cancer")))
  
  # find cases with number of observations < min_cutoff
  data_overview[[f]] <- data[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>% 
    group_by(cl_sub) %>% 
    summarise(totalcount = n(), .groups = "drop") %>% 
    as_tibble() %>% 
    filter(totalcount < min_cutoff)
  
  # remove cases with <min_cutoff values
  data[[f]] <- data[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>%
    filter(!(cl_sub %in% data_overview[[f]]$cl_sub)) %>% 
    select(-cl_sub)
}

# perform permutation test ------------
for (t in tissues_of_interest) {
  
  # cell lines for the tissue of interest
  cell_lines_tissue <- cell_lines %>% 
    filter(tissue == t) %>% 
    pull(cl_id)
  
  # stat.test with respect to each of the non-cancer cell lines
  stat.test <- setNames(vector("list", length(features)), features)
  # to store the adjusted pvalues after BH correction
  stat.test.bh <- setNames(vector("list", length(features)), features)
  
  for (f in features) {

    nc_cells <- data[[f]] %>% 
      # keep only the cell lines of interest
      filter(cl_id %in% cell_lines_tissue) %>%
      # get the non-cancer cell lines
      filter(status == "Non-cancer") %>% 
      pull(cl_id) %>% 
      unique()
    
    c_cells <- data[[f]] %>% 
      # keep only the cell lines of interest
      filter(cl_id %in% cell_lines_tissue) %>%
      # get the cancer cell lines
      filter(status == "Cancer") %>% 
      pull(cl_id) %>% 
      unique()
    
    # make sure that there are non-cancer cell lines
    if (length(nc_cells) > 0) {
      stat.test[[f]] <- setNames(vector("list", length(nc_cells)), nc_cells)
      stat.test.bh[[f]] <- setNames(vector("list", length(nc_cells)), nc_cells)
      
      for (nc in nc_cells) {
        stat.test[[f]][[nc]] <- matrix(NA, length(c_cells), length(subs))
        stat.test[[f]][[nc]] <- as.data.frame(stat.test[[f]][[nc]])
        colnames(stat.test[[f]][[nc]]) <- subs
        rownames(stat.test[[f]][[nc]]) <- c_cells
        
        # perform permutation tests
        for (cancer in c_cells) {
          # temp.data.pre <- 
          
          # Note that output is returned in order because .inorder = TRUE (default in foreach)
          # can just simply substitute values in stat.test[[f]]
          pval <- foreach(s=subs, .combine=c) %dopar% {
            temp.data <- data[[f]] %>% 
              # filter to keep only the nc and c cell lines
              filter(cl_id %in% c(nc,cancer)) %>%
              # filter to keep only the substrate of interest
              filter(sub_id == s)
            # check there is one cancer and one non-cancer cell line
            # keep only the cell id and status column followed by 
            # removing the duplicates
            if (nrow(temp.data %>% 
                     select(cl_id, status) %>% 
                     distinct()) == 2) {
              # significance test
              # Adapted from https://rcompanion.org/handbook/F_15.html
              t <- percentileTest(feature_value ~ status,
                                  data = temp.data,
                                  test = "median",
                                  r    = 5000)
              round(t$Result$p.value, 3)
            }
            # if condition not met than just return NA
            else { NA }
          }
          # store the pvalues  
          stat.test[[f]][[nc]][cancer, ] <- pval 
        }
        
        # adjust for multiple testing using BH correction
        p <- as.vector(as.matrix(stat.test[[f]][[nc]]))
        stat.test.bh[[f]][[nc]] <- matrix(p.adjust(p, "BH"), 
                                          nrow(stat.test[[f]][[nc]]), 
                                          ncol(stat.test[[f]][[nc]]))
        rownames(stat.test.bh[[f]][[nc]]) <- rownames(stat.test[[f]][[nc]])
        colnames(stat.test.bh[[f]][[nc]]) <- colnames(stat.test[[f]][[nc]])
      }
    }
  }
  
  save(stat.test.bh, file = paste0(t, "_BH_adj_pvals.RData"))
}
