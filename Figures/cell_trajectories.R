# Aim: To plot the cell trajectories for each cell line and substrate
# that are used to calculate speed and directional autocorrelation.

# Load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(celltrackR)
})

# set working directory same as the current script location -------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# parameters -------
# cell line names
cell_lines <- read.csv("cell_line_label_colors.csv",
  stringsAsFactors = FALSE
) |>
  # only need the cl_id column
  select(cl_id) |>
  # remove Caov-3 since it is not used in the motility analysis
  filter(cl_id != "Caov-3") |>
  pull(cl_id)

# substrates
subs <- c(
  "30kpacol", "30kpafn", "500pacol", "500pafn",
  "hacol", "hafn", "glass"
)

# set the cutoff for min. number of observations needed for
# a cell-substrate pair to be included in the analysis
min_cutoff <- 25

# Trajectory (Rose) plots -------
for (c in cell_lines) {
  for (s in subs) {
    # check if the directory exists
    if (dir.exists(paste0("../raw_tracking_data_motility/", c, "/processed_data/", s))) {
      data <- read.tracks.csv(paste0("../raw_tracking_data_motility/", c, "/processed_data/", s, "/cell_tracks.tsv"),
        id.column = "ID",
        time.column = "time",
        pos.columns = c("x", "y"),
        header = TRUE,
        sep = "\t"
      )

      if (length(data) >= min_cutoff) {
        png(paste0("cell_trajectories/", c, "_", s, ".png"),
          width = 10, height = 10, units = "in", res = 300
        )
        plot(normalizeTracks(data), xlab = "x", ylab = "y", pch.start = NULL, cex = 0.5)
        dev.off()
      } else {
        print(paste("Not enough data for", c, "on", s))
      }
    }
  }
}
