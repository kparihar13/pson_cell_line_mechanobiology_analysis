# Aim: density plots for hTERT-HPNE cell area and Capan-1 cell stiffness

# load required packages -------
suppressPackageStartupMessages({
  library(tidyverse)
})

# set the working directory same as the location of this script -------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Global parameters -------
features <- c("area", "cell_stiffness")

# load the data --------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)
for (f in features) {
  data[[f]] <- read_tsv(paste("../data/", f, ".tsv", sep = ""),
    show_col_types = FALSE
  ) %>%
    as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>%
    mutate_at(c("cl_id", "sub_id"), as.factor)
}

# function to calculate binwidth for density plots -------
# Dataframe -> Numeric
# Calculate binwidth for density plots using Freedman-Diaconis rule
bw.function <- function(data) {
  bw <- 2 * IQR(data) / length(data)^(1 / 3)
  return(bw)
}

# plot KDE for hTERT-HPNE area -----------
df <- data[["area"]] %>%
  filter(cl_id == "hTERT-HPNE")

plt_htert <- df %>%
  ggplot(aes(feature_value)) +
  geom_histogram(aes(y = after_stat(density)),
    colour = 1,
    fill = "white",
    binwidth = bw.function(df$feature_value)
  ) +
  geom_density() +
  labs(
    x = expression("Area (" * mu * " m)"^2),
    y = "Density"
  ) +
  theme_bw() +
  theme(
    # eliminates gridlines, and chart border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # draw x and y axis lines
    axis.line = element_line(color = "black"),
    # set axis labels and text size
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  )

# add text "hTERT-HPNE across all substrates" to the plot right upper corner
plt_htert <- plt_htert +
  annotate("text",
    x = Inf,
    y = Inf,
    hjust = 1,
    vjust = 1,
    label = "hTERT-HPNE    \nacross all substrates",
    size = 5,
    colour = "black"
  )

ggsave("../Figures/Supplementary_Figure7/hTERTHPNE_area_KDE.png",
  plt_htert,
  dpi = 300, width = 1500, height = 900, units = "px"
)

# Plot KDE for Capan-1 cell stiffness ------------
df <- data[["cell_stiffness"]] %>%
  filter(cl_id == "Capan-1")

plt_capan1 <- df %>% ggplot(aes(feature_value)) +
  geom_histogram(aes(y = after_stat(density)),
    colour = 1,
    fill = "white",
    binwidth = bw.function(df$feature_value)
  ) +
  geom_density() +
  labs(
    x = "Cell Stiffness (Pa)",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    # eliminates gridlines, and chart border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # draw x and y axis lines
    axis.line = element_line(color = "black"),
    # set axis labels and text size
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  )

# add text "Capan-1 across all substrates" to the plot right upper corner
plt_capan1 <- plt_capan1 +
  annotate("text",
    x = Inf,
    y = Inf,
    hjust = 1,
    vjust = 1,
    label = "Capan-1         \nacross all substrates",
    size = 5,
    colour = "black"
  )

ggsave("../Figures/Supplementary_Figure7/Capan1_cell_stiffness_KDE.png",
  plt_capan1,
  dpi = 300, width = 1500, height = 900, units = "px"
)
