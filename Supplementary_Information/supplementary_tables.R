# Aim: Count the number of measurements for each physical feature for each of the
# cell line-substrate pairs

# load the packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(kableExtra)
})

# set the working directory as the location of the script -------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# global parameters -------
features <- c(
  "area", "circularity", "aspect_ratio",
  "cell_stiffness", "motility"
)

cell.line.names <- read.csv("../Figures/cell_line_label_colors.csv",
  stringsAsFactors = FALSE
) %>%
  select(cl_id, label)

subs <- c(
  "500Pa Coll", "500Pa FN",
  "30kPa Coll", "30kPa FN",
  "HA Coll", "HA FN",
  "Glass"
)

# load the data --------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)
for (f in features) {
  data[[f]] <- read_tsv(paste("../data/", f, ".tsv", sep = ""),
    show_col_types = FALSE
  ) %>% as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>%
    mutate_at(c("cl_id", "sub_id"), as.factor)
}

# Count the number of measurements -----
counts <- setNames(vector("list", length(features)), features)

for (f in features) {
  counts[[f]] <- data[[f]] %>%
    group_by(cl_id, sub_id) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    # pivot to wider
    pivot_wider(names_from = sub_id, values_from = n) %>%
    # join with cell.line.names to get the cell line labels
    left_join(cell.line.names, by = "cl_id") %>%
    # make cl_id factor with levels in the order of cell.line.names
    mutate(cl_id = factor(cl_id, levels = cell.line.names$cl_id)) %>%
    arrange(cl_id) %>%
    select(-cl_id) %>%
    # make cl_id to rownames
    column_to_rownames(var = "label")
}

# Compare the counts for area, aspect ratio, and circularity ----
# expect them to be the same for all cell line-substrate pairs

if (sum(counts[["area"]] != counts[["aspect_ratio"]]) != 0 |
  sum(counts[["area"]] != counts[["circularity"]]) != 0) {
  message("The counts for area, aspect ratio, and circularity are NOT same for all cell line-substrate pairs")
} else {
  message("The counts for area, aspect ratio, and circularity are same for all cell line-substrate pairs")
}

# save the counts as a latex table -------
for (f in c("area", "cell_stiffness", "motility")) {
  if (f == "area") {
    temp <- "Number of spread area, aspect ratio, and circularity measurements for each cell line-substrate combination."
  } else if (f == "cell_stiffness") {
    temp <- "Number of cell stiffness measurements for each cell line-substrate combination."
  } else if (f == "motility") {
    temp <- "Number of cell speed measurements for each cell line-substrate combination."
  }

  # change NA to 0
  counts[[f]][is.na(counts[[f]])] <- 0

  # order the columns
  counts[[f]] <- counts[[f]][, c(
    "500Pa Coll", "500Pa FN",
    "30kPa Coll", "30kPa FN",
    "HA Coll", "HA FN",
    "Glass"
  )]

  # Change 500Pa and 30kPa to 500 Pa and 30 kPa
  colnames(counts[[f]]) <- gsub("500Pa", "500 Pa", colnames(counts[[f]]))
  colnames(counts[[f]]) <- gsub("30kPa", "30 kPa", colnames(counts[[f]]))

  counts[[f]] |>
    kable(
      format = "latex",
      align = "c",
      # caption differs based on feature, so can't define here
      caption = temp,
      label = paste0("count-", f),
      booktabs = TRUE,
      linesep = "",
      escape = FALSE
    ) |>
    kable_styling(latex_options = "hold_position") |>
    save_kable(paste0("supplementary_table_", f, ".tex"))
}
