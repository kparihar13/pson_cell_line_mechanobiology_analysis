# Aim: To calculate and plot directional autocorrelation for each cell line and substrate

# Load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(celltrackR)
})

# set working directory same as the current script location -------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# parameters -------
# cell line names
cell_lines <- read.csv("../Figures/cell_line_label_colors.csv",
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

# Calculate directional autocorrelation --------

# Returns TRUE if cell was not stationary in first and last step of a track
# Used to filter out pair of steps (s1, s2) where in either one of them the
# cell was stationary
check <- function(x) {
  all(sapply(list(head(x, 2), tail(x, 2)), trackLength) > 0.0)
}

for (c in cell_lines) {
  da.data <- data.frame()

  for (s in subs) {
    print(paste("Processing", c, "on", s))

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
        # calculate the directional autocorrelation using celltrackR function
        # adapted from https://ingewortel.github.io/celltrackR/vignettes-out/ana-methods.html
        temp <- aggregate(data, overallNormDot,
          FUN = "mean.se",
          filter.subtracks = check
        )
        temp$substrate <- s
        temp$dt <- temp$i * timeStep(data)
        # remove case with NA in upper or lower
        # these correspond to situations when there is only one value for that dt
        temp <- temp[!is.na(temp$upper) & !is.na(temp$lower), ]

        # add to da dataframe
        da.data <- rbind(da.data, temp)
      } else {
        print(paste0("Not enough data for ", c, " on ", s))
      }
    }
  }

  # save the data as tsv
  write.table(da.data,
    file = paste0("directional_autocorrelation/", c, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  # plot directional autocorrelation
  plt <- ggplot(da.data, aes(x = dt, y = mean, color = substrate, fill = substrate)) +
    geom_hline(yintercept = 0) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line() +
    labs(
      title = c,
      x = expression(paste(Delta, " t (min)")),
      y = expression(paste("cos(", theta, " )"))
    ) +
    scale_color_manual(values = c(
      "30kpacol" = "#E41A1C",
      "30kpafn" = "#377EB8",
      "500pacol" = "#4DAF4A",
      "500pafn" = "#FF7F00",
      "hacol" = "#FFFF33",
      "hafn" = "#A65628",
      "glass" = "#999999"
    ), label = c(
      "30kpacol" = "30 kPa Coll",
      "30kpafn" = "30 kPa FN",
      "500pacol" = "500 Pa Coll",
      "500pafn" = "500 Pa FN",
      "hacol" = "HA Coll",
      "hafn" = "HA FN",
      "glass" = "Glass"
    )) +
    scale_fill_manual(values = c(
      "30kpacol" = "#E41A1C",
      "30kpafn" = "#377EB8",
      "500pacol" = "#4DAF4A",
      "500pafn" = "#FF7F00",
      "hacol" = "#FFFF33",
      "hafn" = "#A65628",
      "glass" = "#999999"
    ), label = c(
      "30kpacol" = "30 kPa Coll",
      "30kpafn" = "30 kPa FN",
      "500pacol" = "500 Pa Coll",
      "500pafn" = "500 Pa FN",
      "hacol" = "HA Coll",
      "hafn" = "HA FN",
      "glass" = "Glass"
    )) +
    theme_bw() +
    # no legend
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.x = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      plot.title = element_text(size = 12, hjust = 0.5, color = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 12, color = "black"),
    ) +
    # legend in one row
    guides(
      fill = guide_legend(nrow = 1, byrow = TRUE),
      color = guide_legend(nrow = 1, byrow = TRUE)
    )

  # save plot
  ggsave(
    filename = paste0("../Figures/Supplementary_Figures17_21/", c, ".png"),
    plot = plt,
    width = 6,
    height = 4,
    dpi = 300
  )
}
