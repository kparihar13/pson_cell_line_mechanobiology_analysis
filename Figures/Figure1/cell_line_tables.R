# Aim: to plot cell line tables using gt

# load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(flextable)
  library(officer)
})

# set working directory same as the current script location -------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load cell line data -------
cell_lines <- read.csv("../cell_line_label_colors.csv",
  stringsAsFactors = FALSE
) |>
  select(tissue, label)

colnames(cell_lines) <- c("Tissue", "Cell Line")

# plot cell line tables using flextable ---------
tab1 <- cell_lines |>
  filter(Tissue %in% c("Skin", "Prostate") | (Tissue == "Brain" & `Cell Line` == "U-87")) |>
  flextable() |>
  # merge cells in the first column corresponding to same tissue
  merge_v(j = "Tissue") |>
  theme_box() |>
  # center align for all columns
  align(j = c("Tissue", "Cell Line"), align = c("center", "center"), part = "all") |>
  # remove bold from header
  bold(j = c("Tissue", "Cell Line"), bold = FALSE, part = "header") |>
  # set font size
  fontsize(j = c("Tissue", "Cell Line"), size = 6, part = "all") |>
  # set background color
  bg(j = c("Tissue", "Cell Line"), bg = "white", part = "all") |>
  # set border color
  border(j = c("Tissue", "Cell Line"), border = fp_border(color = "black", width = 0.5), part = "all") |>
  # set width of columns
  width(j = "Tissue", width = 10, unit = "mm")

save_as_image(tab1, "cell_lines_table1.png", res = 300)

tab2 <- cell_lines |>
  filter(Tissue %in% c("Pancreas", "Ovary", "Lung") | (Tissue == "Brain" & `Cell Line` == "T98G")) |>
  flextable() |>
  # merge cells in the first column corresponding to same tissue
  merge_v(j = "Tissue") |>
  theme_box() |>
  # center align for all columns
  align(j = c("Tissue", "Cell Line"), align = c("center", "center"), part = "all") |>
  # remove bold from header
  bold(j = c("Tissue", "Cell Line"), bold = FALSE, part = "header") |>
  # set font size
  fontsize(j = c("Tissue", "Cell Line"), size = 6, part = "all") |>
  # set background color
  bg(j = c("Tissue", "Cell Line"), bg = "white", part = "all") |>
  # set border color
  border(j = c("Tissue", "Cell Line"), border = fp_border(color = "black", width = 0.5), part = "all") |>
  # set width of columns
  width(j = "Tissue", width = 10, unit = "mm")

save_as_image(tab2, "cell_lines_table2.png", res = 300)

tab3 <- cell_lines |>
  filter(Tissue %in% c("Colon", "Breast")) |>
  flextable() |>
  # merge cells in the first column corresponding to same tissue
  merge_v(j = "Tissue") |>
  theme_box() |>
  # center align for all columns
  align(j = c("Tissue", "Cell Line"), align = c("center", "center"), part = "all") |>
  # remove bold from header
  bold(j = c("Tissue", "Cell Line"), bold = FALSE, part = "header") |>
  # set font size
  fontsize(j = c("Tissue", "Cell Line"), size = 6, part = "all") |>
  # set background color
  bg(j = c("Tissue", "Cell Line"), bg = "white", part = "all") |>
  # set border color
  border(j = c("Tissue", "Cell Line"), 
         border = fp_border(color = "black", width = 0.5), 
         part = "all") |>
  # set width of columns
  width(j = "Tissue", width = 10, unit = "mm")

save_as_image(tab3, "cell_lines_table3.png", res = 300)
