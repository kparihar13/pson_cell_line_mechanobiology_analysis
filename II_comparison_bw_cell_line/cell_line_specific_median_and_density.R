# Aim: plot cell line specific median values for each feature along with density plot

# load required libraries -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
})

# set the working directory as the source file location ------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Global parameters ------
features <- c("area", "circularity", "aspect_ratio","cell_stiffness", "motility")
# tissue and respective colors for text color in heatmaps
cell.line.label.colors = read.csv("../Figures/cell_line_label_colors.csv", 
                                  stringsAsFactors = FALSE) %>%
  # change cell_id to factor with levels as the current order of cell lines
  # done to ensure the order of cell lines in the plots
  mutate(cl_id = factor(cl_id, levels = unique(cl_id)))

# load the data --------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)
for (f in features) {
  data[[f]] <- read_tsv(paste("../data/", f, ".tsv", sep = ""), 
                        show_col_types = FALSE) %>% as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>% 
    mutate_at(c("cl_id", "sub_id"), as.factor)
}

# Determine cell line specific median values -----------
medians <- setNames(vector("list", length(features)), features)
# number of digits after decimal for each feature 
decimal <- list("area"=0,"circularity"=2,"aspect_ratio"=2,
                "cell_stiffness"=0,"motility"=1)
for (f in features){
  medians[[f]] = data[[f]] %>% 
    group_by(cl_id) %>% 
    summarise(median = round(median(feature_value),decimal[[f]])) %>% 
    ungroup() %>%
    column_to_rownames(var = "cl_id")
}

# draw heatmap ------------
# range within which the median values are shown in the heatmap 
# with "black" text color and outside the range with "white" text color
# done for visibility of text in heatmap
xlimit.text <- list("area" = c(700,2500), 
                    "circularity" = c(0.5,0.85),
                    "aspect_ratio" = c(1.27,2), 
                    "cell_stiffness" = c(2600,9000),
                    "motility" = c(15,45))
# x-limit for density plots to be used in annotation
xlimit.ha <- list("area" = c(0,3000), 
                  "circularity" = c(0,1),
                  "aspect_ratio" = c(0.95,3), 
                  "cell_stiffness" = c(0,15000),
                  "motility" = c(0,125))
# scale for joyplot
j.scale <- list("area" = 3, 
                "circularity" = 5,
                "aspect_ratio" = 5, 
                "cell_stiffness" = 3,
                "motility" = 3)
# x-axis label for density plots
xlab <- list("area" = expression("Area (" * mu*"m)"^2), 
             "circularity" = "Circularity",
             "aspect_ratio" = "Aspect Ratio", 
             "cell_stiffness" = "Cell Stiffness (kPa)",
             "motility" = expression("Speed (" * mu*"m/hr)"))
# tick marks for the legends
at <- list("area" = c(0,1000,2000,3000), 
           "circularity" = c(0.2,0.4,0.6,0.8,1),
           "aspect_ratio" = c(1,1.5,2,2.5), 
           "cell_stiffness" = c(1000,4000,7000,10000),
           "motility" = c(0,20,40,60))

for (f in features) {
  # keep only the cell lines for which median values are available
  # in case any of the features does not have all the 30 cell lines
  cell.line.label.colors.temp <- cell.line.label.colors %>% 
    filter(cl_id %in% rownames(medians[[f]]))
  
  # get the median values in the required order
  temp.matrix <- medians[[f]] %>%
    rownames_to_column(var = 'cl_id') %>%
    # order in desired order of cell lines
    mutate(cl_id = factor(cl_id, levels = cell.line.label.colors.temp$cl_id)) %>%
    arrange(cl_id) %>%
    select(-cl_id) %>%
    as.matrix()
  
  # remove column name so that it does come in figure
  colnames(temp.matrix) <- " "
  
  # create a list of feature values for each cell line for density plots
  data.list <- setNames(vector("list", length(cell.line.label.colors.temp$cl_id)), 
                        cell.line.label.colors.temp$cl_id)
  for (cell in cell.line.label.colors.temp$cl_id) {
    data.list[[cell]] <- filter(data[[f]], cl_id == cell)$feature_value
  }
  
  # needed to cell_stiffness separately as it has much larger scale needing manual
  # setting of x-axis tick marks for the density plots
  if (f != "cell_stiffness") {
    ha <- rowAnnotation(dens = anno_density(data.list, 
                                            type = "line", 
                                            xlim = xlimit.ha[[f]],
                                            joyplot_scale = j.scale[[f]],
                                            border = FALSE, 
                                            gp = gpar(fill = "#CCCCCC50"),
                                            axis_param = list(gp = gpar(fontsize = 12), 
                                                              labels_rot = 0),
                                            width = unit(4, "cm")), 
                        show_annotation_name = FALSE)
    ht <- Heatmap(temp.matrix,
                  col = hcl.colors(12, "Blue-Red"),
                  cluster_rows = FALSE,  
                  show_row_dend = FALSE,
                  show_row_names = TRUE,
                  row_labels = cell.line.label.colors.temp$label,
                  row_names_gp = gpar(fontsize = 12, 
                                      col = cell.line.label.colors.temp$tissue_col),
                  row_names_side = 'left',
                  row_names_centered = TRUE,
                  width = unit(2, "cm"),
                  height = unit(20, "cm"),
                  right_annotation = ha,
                  row_split = c(1:nrow(temp.matrix)),
                  row_title = NULL, row_gap = unit(1.7, "mm"),
                  cell_fun = function(j, i, x, y, w, h, col) { 
                    # add text (median value) to each grid with black or white color
                    if(temp.matrix[i, j] > xlimit.text[[f]][2] || 
                       temp.matrix[i, j] < xlimit.text[[f]][1]){
                      grid.text(temp.matrix[i, j], x, y,
                                gp = gpar(fontsize = 12, col = 'white'))
                    }
                    else {
                      grid.text(temp.matrix[i, j], x, y,
                                gp = gpar(fontsize = 12)) }
                  },
                  show_heatmap_legend = FALSE,
                  heatmap_legend_param = list(title = "", 
                                              #title_position = "bottomcenter", 
                                              #title_gp = gpar(fontsize = 12),
                                              labels_gp = gpar(fontsize = 11), 
                                              border = "black", 
                                              direction = "horizontal",
                                              grid_width = unit(3.2, "cm"), 
                                              grid_height = unit(0.4, "cm"))
    )
  } else {
    ha <- rowAnnotation(dens = anno_density(data.list, type = "line", xlim = xlimit.ha[[f]],
                                            joyplot_scale = j.scale[[f]], border = FALSE, 
                                            gp = gpar(fill = "#CCCCCC50"),
                                            axis_param = list(gp = gpar(fontsize = 12), 
                                                              labels_rot = 0,
                                                              at = c(0,5000,10000,15000), 
                                                              labels = c(0,5,10,15)),
                                            width = unit(4, "cm")), 
                        show_annotation_name = FALSE)
    ht <- Heatmap(temp.matrix,
                  col = hcl.colors(12, "Blue-Red"),
                  cluster_rows = FALSE,  
                  show_row_dend = FALSE,
                  show_row_names = TRUE,
                  row_labels = cell.line.label.colors.temp$label,
                  row_names_gp = gpar(fontsize = 12, 
                                      col = cell.line.label.colors.temp$tissue_col),
                  row_names_side = 'left',
                  row_names_centered = TRUE,
                  width = unit(2, "cm"),
                  height = unit(20, "cm"),
                  right_annotation = ha,
                  row_split = c(1:nrow(temp.matrix)),
                  row_title = NULL, row_gap = unit(1.7, "mm"),
                  cell_fun = function(j, i, x, y, w, h, col) { 
                    # add text (median value) to each grid with black or white color
                    if(temp.matrix[i, j] > xlimit.text[[f]][2] || 
                       temp.matrix[i, j] < xlimit.text[[f]][1]){
                      grid.text(temp.matrix[i, j], x, y,
                                gp = gpar(fontsize = 12, col = 'white'))
                    }
                    else {
                      grid.text(temp.matrix[i, j], x, y,
                                gp = gpar(fontsize = 12)) }
                  },
                  show_heatmap_legend = FALSE,
                  heatmap_legend_param = list(title = "",
                                              # title_position = "topcenter", title_gp = gpar(fontsize = 12),
                                              labels_gp = gpar(fontsize = 11),
                                              border = "black", 
                                              direction = "horizontal",
                                              at = c(1000,4000,7000,10000), 
                                              labels = c(1000,4000,7000,10000),
                                              grid_width = unit(3.2, "cm"), 
                                              grid_height = unit(0.4, "cm"))
                  )
  }
  
  # for plots with heatmap legend
  # png(paste("legends/",f,".png",sep=""),res = 300, width = 2000, height = 3000)
  
  # save the heatmap
  if (f %in% c("area", "cell_stiffness", "motility")) {
    png(paste("../Figures/Figure3/",f,".png",sep=""),
        res = 300, width = 2000, height = 3500)
  } else {
    png(paste("../Figures/Supplementary_Figure4/",f,".png",sep=""),
        res = 300, width = 2000, height = 3500)
  }
  
  draw(ht, heatmap_legend_side = "bottom")
  
  # add x-axis label for the density plots and "Median Value" label for heatmap
  decorate_annotation("dens", { 
    grid.text(xlab[[f]], 
              x = unit(0.5,"npc"),
              y = unit(-1.79*nrow(temp.matrix),"line"),
              # !!!
              # y = ifelse(f == "motility", 
              #            unit(-1.85*nrow(temp.matrix),"line"),
              #            unit(-1.79*nrow(temp.matrix),"line")),
              gp = gpar(fontsize = 12))

    # !!!
    # if (f != "motility") {grid.text("Median\nvalue", 
    #                                 x = unit(-0.27,"npc"), 
    #                                 y = unit(-1.76*nrow(temp.matrix),"line"),
    #                                 gp = gpar(fontsize = 12))}
    # else {grid.text("Median\nvalue", x
    #                 = unit(-0.27,"npc"), 
    #                 y = unit(-1.82*nrow(temp.matrix),"line"),
    #                 gp = gpar(fontsize = 12))}
    grid.text("Median\nvalue", 
              x = unit(-0.27,"npc"),
              y = unit(-1.76*nrow(temp.matrix),"line"),
              gp = gpar(fontsize = 12))
  })
  dev.off()
  
  # plot the legends for the heatmap
  col_fun = colorRamp2(seq(min(temp.matrix), max(temp.matrix), length.out=12), 
                       hcl.colors(12, "Blue-Red"))
  if (f %in% c("area", "cell_stiffness", "motility")) {
    png(paste("../Figures/Figure3/",f,"_legend.png",sep=""),
        res = 300, width = 1000, height = 300)
  } else {
    png(paste("../Figures/Supplementary_Figure4/",f,"_legend.png",sep=""),
        res = 300, width = 1000, height = 300)
  }
  
  draw(Legend(col_fun = col_fun, 
              title = "",
              labels_gp = gpar(fontsize = 11),
              border = "black", 
              direction = "horizontal",
              at = at[[f]], 
              labels = at[[f]],
              grid_width = unit(3.2, "cm"), 
              grid_height = unit(0.4, "cm")))
  
  dev.off()
}

