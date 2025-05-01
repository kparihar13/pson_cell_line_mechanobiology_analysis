# Aim: For a particular tissue type, perform a permutation test
# to compare the median values of the features between the cancer and non-cancer cell lines

# load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(rcompanion)
  library(foreach)
  library(doParallel)
})
registerDoParallel(8)

# set the working directory same as the location of the script -----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Global parameters ------
features <- c(
  "area", "circularity", "aspect_ratio",
  "cell_stiffness", "motility"
)
subs <- c(
  "HA Coll", "HA FN", "500Pa Coll", "500Pa FN",
  "30kPa Coll", "30kPa FN", "Glass"
)

# define the tissue of interest
tissues_of_interest <- c("Prostate", "Pancreas", "Lung", "Breast")

# get the cell lines for the tissue of interest
cell_lines <- read.csv("../Figures/cell_line_label_colors.csv",
  stringsAsFactors = FALSE
) %>%
  # keep only relevant columns
  select(cl_id, tissue) %>%
  # filter to keep only the cell lines of interest
  filter(tissue %in% tissues_of_interest)

# minimum cutoff for number of observations to be included in the analysis
min_cutoff <- 25

# number of permutations
reps <- 1000

# load data -----------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)
data_overview <- setNames(vector("list", length(features)), features)

for (f in features) {
  data[[f]] <- read_tsv(paste("../data/", f, ".tsv", sep = ""),
    show_col_types = FALSE
  ) %>%
    as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>%
    # add column for cancer and non-cancer status
    mutate(status = ifelse(cl_id %in% c(
      "RWPE-1", "hTERT-HPNE", "NL20",
      "hTERT-HME1", "MCF10A-JSB"
    ),
    "Non-cancer", "Cancer"
    )) %>%
    # filter to only the cell lines of interest
    filter(cl_id %in% cell_lines$cl_id) %>%
    # fix the order of substrates for plotting purposes
    mutate(sub_id = fct_relevel(sub_id, subs)) %>%
    # fix the order for status
    mutate(status = fct_relevel(status, c("Non-cancer", "Cancer")))

  # find cases with number of observations < min_cutoff
  data_overview[[f]] <- data[[f]] %>%
    mutate(cl_sub = paste(cl_id, "_", sub_id, sep = "")) %>%
    group_by(cl_sub) %>%
    summarise(totalcount = n(), .groups = "drop") %>%
    as_tibble() %>%
    filter(totalcount < min_cutoff)

  # remove cases with <min_cutoff values
  data[[f]] <- data[[f]] %>%
    mutate(cl_sub = paste(cl_id, "_", sub_id, sep = "")) %>%
    filter(!(cl_sub %in% data_overview[[f]]$cl_sub)) %>%
    select(-cl_sub)
}

# perform permutation test ------------

# function to perform permutation test
perm_test <- function(data, reps) {
  # get the values for the two groups
  group1 <- data %>%
    filter(status == "Non-cancer") %>%
    pull(feature_value)
  group2 <- data %>%
    filter(status == "Cancer") %>%
    pull(feature_value)

  # calculate the observed difference in medians
  obs_diff <- median(group1) - median(group2)

  # store the permuted differences
  perm_diff <- vector("numeric", reps)

  # perform the permutation test
  for (i in 1:reps) {
    # shuffle the data
    shuffled_data <- sample(c(group1, group2), replace = FALSE)
    # calculate the permuted difference
    perm_diff[i] <- median(shuffled_data[1:length(group1)]) -
      median(shuffled_data[(length(group1) + 1):length(shuffled_data)])
  }

  # calculate the pvalue
  # two-sided test, thats why abs() is used
  # Note: +1 in numerator and denominator to ensure finite sample type-I error control
  # +1 basically corresponds to including the observed value in the permutation distribution
  pval <- (1 + sum(abs(perm_diff) >= abs(obs_diff))) / (1 + reps)

  return(pval)
}

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
      # define an empty dataframe for storing the pvalues
      stat.test[[f]] <- data.frame()
      stat.test.bh[[f]] <- data.frame()

      for (nc in nc_cells) {
        temp.test <- matrix(NA, length(c_cells), length(subs))
        temp.test <- as.data.frame(temp.test)
        colnames(temp.test) <- subs
        rownames(temp.test) <- c_cells

        # perform permutation tests
        for (cancer in c_cells) {
          # Note that output is returned in order because .inorder = TRUE (default in foreach)
          # can just simply substitute values in temp.test
          pval <- foreach(s = subs, .combine = c) %dopar% {
            temp.data <- data[[f]] %>%
              # filter to keep only the nc and c cell lines
              filter(cl_id %in% c(nc, cancer)) %>%
              # filter to keep only the substrate of interest
              filter(sub_id == s)
            # check there is one cancer and one non-cancer cell line
            # keep only the cell id and status column followed by
            # removing the duplicates
            if (nrow(temp.data %>%
              select(cl_id, status) %>%
              distinct()) == 2) {
              # perform permutation test
              perm_test(temp.data, reps)
            }
            # if condition not met than just return NA
            else {
              NA
            }
          }
          # store the pvalues
          temp.test[cancer, ] <- pval
        }

        temp.test <- temp.test %>%
          rownames_to_column(var = "cancer") %>%
          # add the non-cancer cell line name
          mutate(non_cancer = ifelse(str_equal(nc, "MCF10A-JSB"), "MCF10A", nc)) %>%
          pivot_longer(
            cols = -c(non_cancer, cancer),
            names_to = "sub_id",
            values_to = "pval"
          ) %>%
          mutate(test_name = paste0(non_cancer, "_", cancer, "_", sub_id)) %>%
          select(test_name, pval)

        stat.test[[f]] <- rbind(stat.test[[f]], temp.test)
      }
    }

    stat.test[[f]] <- stat.test[[f]] %>%
      column_to_rownames(var = "test_name")

    # adjust for multiple testing using BH correction
    # prepare the data for BH correction
    stat.test.bh[[f]] <- as.vector(as.matrix(stat.test[[f]]))
    names(stat.test.bh[[f]]) <- rownames(stat.test[[f]])
    # perform BH correction
    stat.test.bh[[f]] <- p.adjust(stat.test.bh[[f]], method = "BH")
    # convert the vector back to a dataframe
    stat.test.bh[[f]] <- as.matrix(stat.test.bh[[f]])
    colnames(stat.test.bh[[f]]) <- "pval"
    stat.test.bh[[f]] <- as.data.frame(stat.test.bh[[f]]) %>%
      rownames_to_column(var = "test_name") %>%
      separate(test_name, c("non_cancer", "cancer", "sub_id"), sep = "_") %>%
      # pivot data to wide format
      pivot_wider(names_from = sub_id, values_from = pval)

    # order the columns
    stat.test.bh[[f]] <- stat.test.bh[[f]][, c("non_cancer", "cancer", subs)]
  }

  save(stat.test.bh, file = paste0(t, "_BH_adj_pvals.RData"))
}
