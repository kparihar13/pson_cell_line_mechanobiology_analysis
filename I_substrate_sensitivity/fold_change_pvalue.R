# Aim: To calculate the fold change and p-values for each of the physical features

# load packages ---------
suppressPackageStartupMessages({
  library(tidyverse)
  library(foreach)
  library(doParallel)
  registerDoParallel(cores = 8)
})

# load the data --------
features <- c("area", "circularity", "aspect_ratio", "cell_stiffness", "motility")

# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)
for (f in features) {
  data[[f]] <- read_tsv(paste("../data/", f, ".csv", sep = ""), 
                        show_col_types = FALSE) %>% as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
}

# define fold change (to, from) ----------
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
          ## to make sure the order of substrates are in correct order
          ## 1st row -> to and 2nd row -> from
          arrange(factor(sub_id, levels = fc)) %>%
          summarise(ratio(median_value)) %>%
          pull()
        
        # Null hypothesis: same median, i.e. log2(ratio) = 0
        medianratio <- foreach(i = 1:reps, .combine = c) %dopar% {
          oneratio <- temp %>%
            ## permutate the substrate labels, keeping feature values intact
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
                         (sum(abs(medianratio) >= abs(log2(obs_ratio))) + 1) / (1 + reps))

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
  
  # save the fold values and p-values
  save(fold_value, fold_pvalue, file = paste(f,"_folds.RData", sep = ""))
}
