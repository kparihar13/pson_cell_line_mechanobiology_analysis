# Aim: create horizontal and vertically aligned Tissue legend for figures

# load packages 
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
})

# set the working directory same as the script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load the tissue color
tissue.colors = read.csv('cell_line_label_colors.csv', 
                         stringsAsFactors = FALSE) %>%
  select(tissue, tissue_col) %>%
  distinct()

# create the horizontal legend
png('tissue_legend_horizontal.png', res = 300, width = 1200, height = 400)
lgd = Legend(labels = tissue.colors$tissue, 
             title = "Tissue", 
             title_gp = gpar(fontsize = 12), 
             labels_gp = gpar(fontsize = 12), 
             legend_gp = gpar(fill = unique(cell_lines_labels_colors$tissue_col)), 
             ncol = 4, 
             title_position = "lefttop")
draw(lgd)
dev.off()

# create the vertical legend
png('tissue_legend_vertical.png', res = 300, width = 400, height = 1200)
lgd = Legend(labels = tissue.colors$tissue, 
             title = "Tissue", 
             title_gp = gpar(fontsize = 12), 
             labels_gp = gpar(fontsize = 12), 
             legend_gp = gpar(fill = unique(cell_lines_labels_colors$tissue_col)), 
             ncol = 1, 
             title_position = "topcenter")
draw(lgd)
dev.off()
