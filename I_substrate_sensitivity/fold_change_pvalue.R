# Aim: To calculate the fold change and p-values for each of the physical features

# load packages ---------
suppressPackageStartupMessages({
  library(tidyverse)
  library(foreach)
  library(doParallel)
  registerDoParallel(cores = 8)
})

# set the working directory to the location of the script ------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# global parameters ---------
features <- c("area", "circularity", "aspect_ratio", 
              "cell_stiffness", "motility")

# define fold change (to, from)
folds <- list(
  "30k-500Pa Coll" = c("30kPa Coll", "500Pa Coll"),
  "30kPa-500Pa FN" = c("30kPa FN", "500Pa FN"),
  "HA-500Pa Coll" = c("HA Coll", "500Pa Coll"),
  "HA-500Pa FN" = c("HA FN", "500Pa FN"),
  "30kPa-HA Coll" = c("30kPa Coll", "HA Coll"),
  "30kPa-HA FN" = c("30kPa FN", "HA FN"),
  "HA Coll-FN" = c("HA Coll", "HA FN"),
  "500Pa Coll-FN" = c("500Pa Coll", "500Pa FN"),
  "30kPa Coll-FN" = c("30kPa Coll", "30kPa FN"),
  "Glass-30kPa Coll" = c("Glass", "30kPa Coll"),
  "Glass-30kPa FN" = c("Glass", "30kPa FN")
)

# load the data --------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)
for (f in features) {
  data[[f]] <- read_tsv(paste("../data/", f, ".tsv", sep = ""), 
                        show_col_types = FALSE) %>%
    as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
}

# define required custom functions -----------
ratio <- function(x) {
  return(x[1] / x[2])
}

# run two-side permutation test ------------
# minimum number of observations in each group
min.cutoff <- 25

# iterate through features
for (f in features) {
  print(f)
  cell_lines <- unique(data[[f]]$cl_id)
  subs <- unique(data[[f]]$sub_id)
  # number of repetitions in the permutation test
  reps <- 1000
  fold_value <- data.frame(
    col1 = double(), col2 = double(), col3 = double(),
    col4 = double(), col5 = double(), col6 = double(),
    col7 = double(), col8 = double(), col9 = double(),
    col10 = double(), col11 = double(),
    stringsAsFactors = FALSE
  )
  fold_pvalue <- data.frame(
    col1 = double(), col2 = double(), col3 = double(),
    col4 = double(), col5 = double(), col6 = double(),
    col7 = double(), col8 = double(), col9 = double(),
    col10 = double(), col11 = double(),
    stringsAsFactors = FALSE
  )
  
  # iterate through cell lines
  for (i in 1:length(cell_lines)) {
    print(paste("Cell line:", cell_lines[[i]]))

    pvalue_temp <- list()
    ratio_temp <- list()
    # iterate through folds
    for (fc in folds) {
      print(paste("fold:", fc))
      # filter to keep only the data for feature, cell line, and substrates of interest
      temp <- data[[f]] %>% 
        filter((cl_id == cell_lines[[i]]) & (sub_id %in% fc))
  
      # calculate the median values for each of the substrates
      obs_ratio <- temp %>%
        group_by(sub_id) %>%
        summarise(median_value = median(feature_value))
      
      # check if the number of observations in each group is greater than the minimum cutoff
      if ((nrow(filter(temp, sub_id == fc[1])) >= min.cutoff) &&
        (nrow(filter(temp, sub_id == fc[2])) >= min.cutoff)) {
        # calculate the observed ratio of medians
        obs_ratio <- obs_ratio %>%
          # to make sure the order of substrates are in correct order
          # 1st row -> to and 2nd row -> from
          arrange(factor(sub_id, levels = fc)) %>%
          summarise(ratio(median_value)) %>%
          pull()
        
        # Null hypothesis: same median, i.e. log2(ratio) = 0
        medianratio <- foreach(i = 1:reps, .combine = c) %dopar% {
          oneratio <- temp %>%
            # permutate the substrate labels, keeping feature values intact
            mutate(permsubs = sample(sub_id)) %>%
            group_by(permsubs) %>%
            summarise(median_value = median(feature_value)) %>%
            arrange(factor(permsubs, levels = fc)) %>%
            summarise(log2(ratio(median_value))) %>%
            pull()
        }

        # all the values in permutation that are outside the extremum defined by log2(obs_ratio)
        # two-sided test, thats why abs() is used
        # Note: +1 in numerator and denominator to ensure finite sample type-I error control
        # +1 basically corresponds to including the observed value in the permutation distribution
        pvalue_temp <- c(pvalue_temp, 
                         (1 + sum(abs(medianratio) >= abs(log2(obs_ratio)))) / (1 + reps))

        ratio_temp <- c(ratio_temp, obs_ratio)
      } else {
        # if the number of observations in either group is less than the minimum cutoff
        pvalue_temp <- c(pvalue_temp, NA)
        ratio_temp <- c(ratio_temp, NA)
      }
    }
    fold_value[i, ] <- ratio_temp
    fold_pvalue[i, ] <- pvalue_temp
  }

  rownames(fold_value) <- cell_lines
  colnames(fold_value) <- names(folds)
  rownames(fold_pvalue) <- cell_lines
  colnames(fold_pvalue) <- names(folds)
  
  # flatten fold_pvalue for Benjamini-Hochberg correction
  pvalue_bh_corrected <- fold_pvalue %>%
    rownames_to_column(var = "cell_line") %>%
    # convert to long format
    gather(key = "fold", value = "pvalue", -cell_line) %>%
    mutate(cell_line_fold = paste(cell_line, fold, sep = "_")) %>%
    select(cell_line_fold, pvalue) %>%
    column_to_rownames(var = "cell_line_fold") %>%
    as.matrix()
  # convert to named vector
  cell_line_fold <- rownames(pvalue_bh_corrected)
  pvalue_bh_corrected <- as.vector(pvalue_bh_corrected)
  names(pvalue_bh_corrected) <- cell_line_fold
  
  # perform Benjamini-Hochberg correction
  pvalue_bh_corrected <- p.adjust(pvalue_bh_corrected, method = "BH")
  
  # convert corrected values back into data frame
  pvalue_bh_corrected <- as.matrix(pvalue_bh_corrected)
  colnames(pvalue_bh_corrected) <- "pvalue"
  pvalue_bh_corrected <- as.data.frame(pvalue_bh_corrected) %>%
    rownames_to_column(var = "cell_line_fold") %>%
    mutate(cell_line = gsub("_.*", "", cell_line_fold),
           fold = gsub(".*_", "", cell_line_fold)) %>%
    select(cell_line, fold, pvalue) %>%
    # convert to wide format
    spread(key = "fold", value = "pvalue") %>%
    # order cell lines same as the names in the fold_pvalue
    mutate(cell_line = factor(cell_line, levels = rownames(fold_pvalue))) %>%
    arrange(cell_line) %>%
    column_to_rownames(var = "cell_line")
  
  # order the columns of pvalue_bh_corrected
  pvalue_bh_corrected <- pvalue_bh_corrected[, colnames(fold_pvalue)]
  
  # save the fold values and p-values
  save(fold_value, pvalue_bh_corrected, file = paste(f,"_folds.RData", sep = ""))
}
