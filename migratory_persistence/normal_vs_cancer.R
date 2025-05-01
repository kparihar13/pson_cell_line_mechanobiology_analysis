# Aim: To generate barplots for each cell line in a tissue of interest

# load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
})

# set the working directory same as the location of the script -----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Global parameters ------
subs <- c(
  "500 Pa Coll", "500 Pa FN",
  "30 kPa Coll", "30 kPa FN",
  "HA Coll", "HA FN",
  "Glass"
)

# colors for boxplots
cell.line.colors <- c(
  "#9B8ED8", "#CC6677", "#DDCC77",
  "#88CCEE", "#117733", "#EFEDE8"
)

# define the tissue of interest
tissues_of_interest <- c("Breast", "Lung", "Pancreas", "Prostate")

# get the cell lines for the tissue of interest
cell_lines <- read.csv("../Figures/cell_line_label_colors.csv",
  stringsAsFactors = FALSE
) %>%
  # filter to keep only the cell lines of interest
  filter(tissue %in% tissues_of_interest)

# load data ---------
data <- read_tsv(paste("directional_persistence.tsv", sep = ""),
  show_col_types = FALSE
) %>%
  as_tibble() %>%
  # filter to only the cell lines of interest
  filter(cl_id %in% cell_lines$cl_id) %>%
  # fix the order of substrates for plotting purposes
  mutate(sub_id = fct_relevel(sub_id, subs)) %>%
  # add column for status and label (for plotting purposes)
  inner_join(., cell_lines %>% select(cl_id, label, status), by = "cl_id") %>%
  # fix the order for status and labels
  mutate(
    status = fct_relevel(status, c("Non-cancer", "Cancer")),
    label = fct_relevel(label, cell_lines$label)
  )

# Barplots  ----------
for (t in tissues_of_interest) {
  # cell lines for the tissue of interest
  cell_lines_tissue <- cell_lines %>%
    filter(tissue == t) %>%
    pull(cl_id)

  # get the non-cancer cell lines
  nc_cells <- cell_lines %>%
    # keep only the cell lines of interest
    filter(tissue == t) %>%
    # get the non-cancer cell lines
    filter(status == "Non-cancer") %>%
    select(all_of(c("cl_id", "label")))

  c_cells <- cell_lines %>%
    # keep only the cell lines of interest
    filter(tissue == t) %>%
    # get the cancer cell lines
    filter(status == "Cancer") %>%
    pull(cl_id)

  # make sure that there are non-cancer cell lines
  if (nrow(nc_cells) > 0) {
    # order of cell lines in the plot
    cell_order <- c(nc_cells$label, c_cells)
    # color for each of the cell line
    manual.colors <- setNames(vector("list", length(cell_order)), cell_order)
    manual.colors <- cell.line.colors[1:length(manual.colors)]

    # keep only subs with more than one cell line having
    subs.to.keep <- data %>%
      select(all_of(c("sub_id", "cl_id"))) %>%
      distinct() %>%
      # only keep cell lines of interest
      filter(cl_id %in% cell_lines_tissue) %>%
      group_by(sub_id) %>%
      summarise(count = n()) %>%
      filter(count != 1)

    # create the barplot
    plt <- data %>%
      # filter to keep only the substrates of interest and
      # cell lines of interest
      filter(
        (sub_id %in% subs.to.keep$sub_id),
        (cl_id %in% cell_lines_tissue)
      ) %>%
      ggplot(.) +
      aes(
        x = sub_id,
        y = persistence,
        fill = factor(label,
          levels = cell_order
        )
      ) +
      geom_bar(
        stat = "identity",
        position = "dodge",
        width = 0.7,
        colour = "black"
      ) +
      labs(
        x = "",
        y = "Migratory Persistence (min)",
        fill = paste0(t, ": ")
      ) +
      # Specify colours
      scale_fill_manual(values = manual.colors) +
      theme_bw() +
      # set the font size
      theme(
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title.y = element_text(size = 13, colour = "black"),
        legend.text = element_text(size = 13, colour = "black"),
        legend.title = element_text(size = 13, colour = "black")
      ) +
      theme(legend.position = "top") +
      guides(fill = guide_legend(nrow = 1))

    # save the plot
    png(paste0("../Figures/Supplementary_Figure15/", t, ".png"),
      res = 300, width = 2700, height = 1500
    )
    print(plt)
    dev.off()
  }
}
