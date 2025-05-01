# Aim: To calculate an estimate for decorrelation time for each cell line and substrate
# and plot a heatmap of the (decorrelation time/total time) ratio.

# Note: we assume that the decorrelation time is the time at which the cos(theta) goes below 0.1

# Load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(ggsignif)
})

# set working directory same as the current script location -------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# parameters -------
cell.line.label.colors <- read.csv("../Figures/cell_line_label_colors.csv",
  stringsAsFactors = FALSE
) |>
  # remove Caov-3 since it is not used in the motility analysis
  filter(cl_id != "Caov-3")

# substrates
subs <- c(
  "500pacol" = "500 Pa Coll",
  "500pafn" = "500 Pa FN",
  "30kpacol" = "30 kPa Coll",
  "30kpafn" = "30 kPa FN",
  "hacol" = "HA Coll",
  "hafn" = "HA FN",
  "glass" = "Glass"
)

decorr.cutoff <- 0.2

# Decorrelation time -------
decorr.time <- matrix(NA, nrow = nrow(cell.line.label.colors), ncol = length(subs))
# row and column names
rownames(decorr.time) <- cell.line.label.colors$cl_id
colnames(decorr.time) <- names(subs)

# track cases where autocorrelation stays > 0.2 for whole of cell tracking time
boundary.cases <- matrix(0, nrow = nrow(cell.line.label.colors), ncol = length(subs))
# row and column names
rownames(boundary.cases) <- cell.line.label.colors$cl_id
colnames(boundary.cases) <- names(subs)

for (c in cell.line.label.colors$cl_id) {
  # load data
  da.data <- read.table(paste0("directional_autocorrelation/", c, ".tsv"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  for (s in names(subs)) {
    print(paste("Processing", c, "on", s))

    # check whether s is present in da.data
    # substrates with < 25 tracks are not included in the analysis
    if (s %in% da.data$substrate) {
      # filter da.data for s
      temp <- da.data |>
        filter(substrate == s)

      # find the first row with mean <= decorr_cutoff
      decorr.row <- which(temp$mean <= decorr.cutoff)[1]
      if (!is.na(decorr.row)) {
        # interpolate decorr_row time with previous row to get an estimate for decorrelation time
        # slope of the line is (y2 - y1) / (x2 - x1)
        m <- (temp$mean[decorr.row] - temp$mean[decorr.row - 1]) /
          (temp$dt[decorr.row] - temp$dt[decorr.row - 1])
        # intercept of the line is y - mx
        b <- temp$mean[decorr.row] - m * temp$dt[decorr.row]
        # at y = 0.1, x = (y - b) / m; y: cos(theta) and x: time
        decorr.time[c, s] <- round((decorr.cutoff - b) / m, 1)
      } else {
        # if there is no such row, that implies that the mean cos(theta) never goes below 0.1
        # set decorrelation time as the time in the last row
        decorr.time[c, s] <- round(tail(temp$dt, 1), 1)
        boundary.cases[c, s] <- 1
      }
    } else {
      print(paste0("Not enough data for ", c, " on ", s))
    }
  }
}

# change column names for heatmap
colnames(decorr.time) <- subs[colnames(decorr.time)]

# plot heatmap for decorr.time using complexHeatmap
ht <- Heatmap(decorr.time,
  col = colorRamp2(
    breaks = c(5, 40, 100),
    hcl_palette = "PuBu",
    reverse = TRUE
  ),
  na_col = "white",
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_labels = cell.line.label.colors$label,
  row_names_gp = gpar(
    fontsize = 12,
    col = cell.line.label.colors$tissue_col
  ),
  row_split = rep(1:8, c(4, 5, 3, 2, 3, 4, 6, 2)),
  row_gap = unit(1.2, "mm"),
  row_title = NULL,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 2),
  width = unit(10, "cm"),
  height = unit(20, "cm"),
  heatmap_legend_param = list(
    title = "",
    labels_gp = gpar(fontsize = 12),
    at = c(0, 30, 60, 90, 120),
    labels = c(0, 30, 60, 90, 120),
    legend_height = unit(4, "cm"),
    border = "black"
  ),
  cell_fun = function(j, i, x, y, w, h, col) {
    # add text decorrelation time to each grid with black or white color
    if (!is.na(decorr.time[i, j])) {
      if (decorr.time[i, j] > 40) {
        if (boundary.cases[i, j] == 0) {
          grid.text(decorr.time[i, j], x, y,
            gp = gpar(fontsize = 12, col = "white")
          )
        } else {
          grid.text(paste0(decorr.time[i, j], "*"), x, y,
            gp = gpar(fontsize = 12, col = "white")
          )
        }
      } else {
        if (boundary.cases[i, j] == 0) {
          grid.text(decorr.time[i, j], x, y,
            gp = gpar(fontsize = 12)
          )
        } else {
          grid.text(paste0(decorr.time[i, j], "*"), x, y,
            gp = gpar(fontsize = 12)
          )
        }
      }
    }
  }
)

png(paste0("../Figures/Supplementary_Figure16/persistence_heatmap.png"),
  res = 300, width = 1800, height = 2700, units = "px"
)
draw(ht)
dev.off()
