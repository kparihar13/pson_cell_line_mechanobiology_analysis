# Aim: Plots for Figure 2 (change in physical feature with change in substrate)
# Except Circularity, which is saved as Supplementary Figure 2

# load packages ---------
suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(gdata)
})

# define some useful lists ----------

features <- c("area", "circularity", "aspect_ratio", "cell_stiffness", "motility")

feature_names <- list(
  "area" = "Area", "circularity" = "Circularity",
  "aspect_ratio" = "Aspect Ratio",
  "cell_stiffness" = "Cell Stiffness", "motility" = "Motility"
)

folds <- list(
  "30k-500Pa Coll" = "colfold",
  "30kPa-500Pa FN" = "fnfold",
  "HA-500Pa Coll" = "ha500colfold",
  "HA-500Pa FN" = "ha500fnfold",
  "30kPa-HA Coll" = "ha30kcolfold",
  "30kPa-HA FN" = "ha30kfnfold",
  "HA Coll-FN" = "hyalfold",
  "500Pa Coll-FN" = "fn500colfold",
  "30kPa Coll-FN" = "fn30kcolfold",
  "Glass-30kPa Coll" = "glasscolfold",
  "Glass-30kPa FN" = "glassfnfold"
)

# folds of interest
fold_names <- list(
  "colfold" = "30k-500Pa Coll",
  "fnfold" = "30k-500Pa FN",
  "ha500colfold" = "HA-500Pa Coll",
  "ha500fnfold" = "HA-500Pa FN",
  "glasscolfold" = "Glass-30kPa Coll",
  "glassfnfold" = "Glass-30kPa FN"
)

# text color cell line names
cell_line_label_colors_main <- read.csv("cell_lines_labels_color.csv",
  stringsAsFactors = FALSE
)

# Plot fold change ---------------

# Dataframe data, String column_name -> Dataframe with color
# specifies the fold change based color for bar chart
get_color <- function(data, column_name) {
  data <- data %>%
    # blue if > 0, else red
    mutate(
      anno_bar = !!sym(column_name),
      color = ifelse(!!sym(column_name) > 0, "steelblue", "brown")
    ) %>%
    select(anno_bar, color)

  # if fold change is to high or to low, then just take it as 2 or -2 for plotting purposes
  data[which(data$anno_bar > 2), 1] <- 2
  data[which(data$anno_bar < -2), 1] <- -2

  return(data)
}

for (f in features) {
  # load the fold change values and respective pvalues
  load(paste(f, "_folds.RData", sep = ""))

  fold_value <- fold_value %>%
    rownames_to_column(var = "cl_id") %>%
    # order in desired order of cell lines
    mutate(cl_id = factor(cl_id, levels = cell_line_label_colors_main$cl_id)) %>%
    arrange(cl_id) %>%
    column_to_rownames("cl_id")

  fold_pvalue <- fold_pvalue %>%
    rownames_to_column(var = "cl_id") %>%
    # order in desired order of cell lines
    mutate(cl_id = factor(cl_id, levels = cell_line_label_colors_main$cl_id)) %>%
    arrange(cl_id) %>%
    column_to_rownames("cl_id")

  # delete row named PC-3 from fold_value and fold_pvalue
  # all fold change values are NA for this due to lack of data (< min.cutoff of 25 cells)
  # in at least one of the substrates being compared in each fold
  fold_value <- fold_value[-which(rownames(fold_value) == "PC-3"), ]
  fold_pvalue <- fold_pvalue[-which(rownames(fold_pvalue) == "PC-3"), ]

  # log2 transform the ratio values
  fold_value <- log2(fold_value)

  # change the colnames for fold_value and fold_pvalue for plotting convenience
  # used in annotation creation
  for (fc in names(folds)) {
    colnames(fold_value)[which(colnames(fold_value) == fc)] <- folds[[fc]]
    colnames(fold_pvalue)[which(colnames(fold_pvalue) == fc)] <- folds[[fc]]
  }

  # get the colors for cell lines
  # done in case some cell line is not present in all the features
  cell_line_label_colors <- left_join(data.frame(cl_id = rownames(fold_value)),
    cell_line_label_colors_main,
    by = "cl_id"
  )

  ht_opt(RESET = TRUE)
  # list to store the heatmaps for different folds
  ht_list <- list()

  # 1st fold: 500Pa Coll -> 30kPa Coll (colfold)
  # specify the fold change based color for bar chart
  bar_colfold <- get_color(fold_value, "colfold")
  # define the heatmap annotation
  ha_colfold <- rowAnnotation(
    anno_emp = anno_empty(width = unit(0.1, "cm"), border = FALSE),
    colfold = anno_barplot(unlist(bar_colfold$anno_bar),
      baseline = 0, ylim = c(-2, 2),
      axis_param = list(
        side = "bottom",
        labels_rot = 0,
        at = seq(-2, 2, 0.5),
        labels = c(-2, "", -1, "", 0, "", 1, "", 2),
        gp = gpar(fontsize = 9)
      ),
      border = FALSE, bar_width = 0.7,
      gp = gpar(fill = unlist(bar_colfold$color))
    ),
    show_annotation_name = FALSE, width = unit(2.6, "cm")
  )

  # define the heatmap: empty matrix, only interested in the annotation
  ht_list[["colfold"]] <- Heatmap(matrix(nc = 0, nr = nrow(cell_line_label_colors)),
    rect_gp = gpar(type = "none"),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_labels = cell_line_label_colors$label,
    row_names_gp = gpar(
      fontsize = 10,
      col = cell_line_label_colors$tissue_col
    ),
    row_names_side = "left",
    left_annotation = ha_colfold,
    height = unit(13, "cm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  # 2nd fold: 500Pa FN -> 30kPa FN (fnfold)
  # specify the fold change based color for bar chart
  bar_fnfold <- get_color(fold_value, "fnfold")
  # define the heatmap annotation
  ha_fnfold <- rowAnnotation(
    fnfold = anno_barplot(unlist(bar_fnfold$anno_bar),
      baseline = 0, ylim = c(-2, 2),
      axis_param = list(
        side = "bottom", labels_rot = 0,
        at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
        labels = c(-2, "", -1, "", 0, "", 1, "", 2),
        gp = gpar(fontsize = 9)
      ), border = FALSE,
      bar_width = 0.7, gp = gpar(fill = unlist(bar_fnfold$color))
    ),
    show_annotation_name = FALSE, width = unit(2.5, "cm")
  )
  # define the heatmap: empty matrix, only interested in the annotation
  ht_list[["fnfold"]] <- Heatmap(matrix(nc = 0, nr = nrow(cell_line_label_colors)),
    rect_gp = gpar(type = "none"),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_labels = cell_line_label_colors$label,
    row_names_gp = gpar(
      fontsize = 10,
      col = cell_line_label_colors$tissue_col
    ),
    row_names_side = "left",
    left_annotation = ha_fnfold,
    height = unit(13, "cm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  # 3rd fold: 500Pa Coll -> HA Coll (ha500colfold)
  # specify the fold change based color for bar chart
  bar_ha500colfold <- get_color(fold_value, "ha500colfold")
  # define the heatmap annotation
  ha_ha500colfold <- rowAnnotation(
    ha500colfold = anno_barplot(unlist(bar_ha500colfold$anno_bar),
      baseline = 0, ylim = c(-2, 2),
      axis_param = list(
        side = "bottom", labels_rot = 0,
        at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
        labels = c(-2, "", -1, "", 0, "", 1, "", 2),
        gp = gpar(fontsize = 9)
      ), border = FALSE,
      bar_width = 0.7, gp = gpar(fill = unlist(bar_ha500colfold$color))
    ),
    show_annotation_name = FALSE, width = unit(2.5, "cm")
  )
  # define the heatmap: empty matrix, only interested in the annotation
  ht_list[["ha500colfold"]] <- Heatmap(matrix(nc = 0, nr = nrow(cell_line_label_colors)),
    rect_gp = gpar(type = "none"),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_labels = cell_line_label_colors$label,
    row_names_gp = gpar(
      fontsize = 10,
      col = cell_line_label_colors$tissue_col
    ),
    row_names_side = "left",
    left_annotation = ha_ha500colfold,
    height = unit(13, "cm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  # 4th fold: 500Pa FN -> HA FN (ha500fnfold)
  # specify the fold change based color for bar chart
  bar_ha500fnfold <- get_color(fold_value, "ha500fnfold")
  # define the heatmap annotation
  ha_ha500fnfold <- rowAnnotation(
    ha500fnfold = anno_barplot(unlist(bar_ha500fnfold$anno_bar),
      baseline = 0, ylim = c(-2, 2),
      axis_param = list(
        side = "bottom", labels_rot = 0,
        at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
        labels = c(-2, "", -1, "", 0, "", 1, "", 2),
        gp = gpar(fontsize = 9)
      ), border = FALSE,
      bar_width = 0.7, gp = gpar(fill = unlist(bar_ha500fnfold$color))
    ),
    show_annotation_name = FALSE, width = unit(2.5, "cm")
  )
  # define the heatmap: empty matrix, only interested in the annotation
  ht_list[["ha500fnfold"]] <- Heatmap(matrix(nc = 0, nr = nrow(cell_line_label_colors)),
    rect_gp = gpar(type = "none"),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_labels = cell_line_label_colors$label,
    row_names_gp = gpar(fontsize = 10, col = cell_line_label_colors$tissue_col),
    row_names_side = "left",
    left_annotation = ha_ha500fnfold,
    height = unit(13, "cm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  # 5th fold: 30kPA Coll -> Glass (glasscolfold)
  # specify the fold change based color for bar chart
  bar_glasscolfold <- get_color(fold_value, "glasscolfold")
  ## define the heatmap annotation
  ha_glasscolfold <- rowAnnotation(
    glasscolfold = anno_barplot(unlist(bar_glasscolfold$anno_bar),
      baseline = 0, ylim = c(-2, 2),
      axis_param = list(
        side = "bottom",
        labels_rot = 0,
        at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
        labels = c(-2, "", -1, "", 0, "", 1, "", 2),
        gp = gpar(fontsize = 9)
      ),
      border = FALSE, bar_width = 0.7,
      gp = gpar(fill = unlist(bar_glasscolfold$color))
    ),
    show_annotation_name = FALSE, width = unit(2.5, "cm")
  )
  # define the heatmap: empty matrix, only interested in the annotation
  ht_list[["glasscolfold"]] <- Heatmap(matrix(nc = 0, nr = nrow(cell_line_label_colors)),
    rect_gp = gpar(type = "none"),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_labels = cell_line_label_colors$label,
    row_names_gp = gpar(
      fontsize = 10,
      col = cell_line_label_colors$tissue_col
    ),
    row_names_side = "left",
    left_annotation = ha_glasscolfold,
    height = unit(13, "cm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  # 6th fold: 30kPa FN -> Glass (glassfnfold)
  # specify the fold change based color for bar chart
  bar_glassfnfold <- get_color(fold_value, "glassfnfold")
  # define the heatmap annotation
  ha_glassfnfold <- rowAnnotation(
    glassfnfold = anno_barplot(unlist(bar_glassfnfold$anno_bar),
      baseline = 0, ylim = c(-2, 2),
      axis_param = list(
        side = "bottom", labels_rot = 0,
        at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
        labels = c(-2, "", -1, "", 0, "", 1, "", 2),
        gp = gpar(fontsize = 9)
      ), border = FALSE,
      bar_width = 0.7, gp = gpar(fill = unlist(bar_glassfnfold$color))
    ),
    show_annotation_name = FALSE, width = unit(2.5, "cm")
  )
  # define the heatmap: empty matrix, only interested in the annotation
  ht_list[["glassfnfold"]] <- Heatmap(matrix(nc = 0, nr = nrow(cell_line_label_colors)),
    rect_gp = gpar(type = "none"),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_labels = cell_line_label_colors$label,
    row_names_gp = gpar(
      fontsize = 10,
      col = cell_line_label_colors$tissue_col
    ),
    row_names_side = "left",
    left_annotation = ha_glassfnfold,
    height = unit(13, "cm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  ha_gap <- rowAnnotation(empty = anno_empty(border = FALSE, width = unit(0.2, "mm")))

  if (f != "circularity") {
    png(paste("../Figures/Figure2/", f, ".png", sep = ""),
      res = 300, width = 3000, height = 2000
    )
  } else {
    png(paste("../Figures/Supplementary_Figure2/", f, "_A.png", sep = ""),
      res = 300, width = 3000, height = 2000
    )
  }

  draw(ht_list[["colfold"]] + ha_gap + ha_gap + ht_list[["fnfold"]] +
    ha_gap + ha_gap + ht_list[["ha500colfold"]] +
    ha_gap + ha_gap + ht_list[["ha500fnfold"]] +
    ha_gap + ha_gap + ht_list[["glasscolfold"]] +
    ha_gap + ha_gap + ha_gap + ht_list[["glassfnfold"]], auto_adjust = FALSE)

  for (k in 1:length(names(fold_names))) {
    decorate_annotation(names(fold_names)[k], {
      pvalue_temp <- fold_pvalue %>%
        select(names(fold_names)[k]) %>%
        pull()
      fold_temp <- fold_value %>%
        select(names(fold_names)[k]) %>%
        pull()

      # place *, **, *** based on p-values
      if (sum(is.na(pvalue_temp)) > 0) {
        pvalue_temp[which(is.na(pvalue_temp))] <- 1000
      }
      eps <- ifelse(f == "motility" && k == 1, 0.95, 1)
      for (i in 1:length(pvalue_temp)) {
        y_pvalue <- unit((length(pvalue_temp) - i + 0.4) * 130 / length(pvalue_temp), "mm")
        if (pvalue_temp[[i]] <= 0.001) {
          if (abs(fold_temp[[i]]) > 2) {
            x.temp <- sign(fold_temp[[i]]) * 2
            grid.text(paste(sprintf(fold_temp[[i]], fmt = "%#.1f"), "***", sep = ""),
              x = unit(0.5 + 0.08 * sign(fold_temp[[i]]), "npc") +
                unit(eps * x.temp * 1.25 / 1.75, "cm"),
              y = y_pvalue + unit(0.25, "mm"),
              gp = gpar(fontsize = 8)
            )
          } else if (abs(fold_temp[[i]]) >= 1) {
            grid.text("***",
              x = unit(0.5 + 0.08 * sign(fold_temp[[i]]), "npc") +
                unit(eps * fold_temp[[i]] * 1.25 / 2.1, "cm"),
              y = y_pvalue, gp = gpar(fontsize = 8)
            )
          } else {
            grid.text("***",
              x = unit(0.5 + 0.08 * sign(fold_temp[[i]]), "npc") +
                unit(eps * fold_temp[[i]] * 1.25 / 2, "cm"),
              y = y_pvalue,
              gp = gpar(fontsize = 8)
            )
          }
        } else if (pvalue_temp[[i]] <= 0.01) {
          if (abs(fold_temp[[i]]) > 2) {
            x.temp <- sign(fold_temp[[i]]) * 2
            grid.text(paste(sprintf(fold_temp[[i]], fmt = "%#.1f"), "**", sep = ""),
              x = unit(0.5 + 0.06 * sign(fold_temp[[i]]), "npc") +
                unit(eps * x.temp * 1.25 / 1.8, "cm"),
              y = y_pvalue + unit(0.25, "mm"),
              gp = gpar(fontsize = 8)
            )
          } else if (abs(fold_temp[[i]]) >= 1) {
            grid.text("**",
              x = unit(0.5 + 0.06 * sign(fold_temp[[i]]), "npc") +
                unit(eps * fold_temp[[i]] * 1.25 / 2.15, "cm"),
              y = y_pvalue, gp = gpar(fontsize = 8)
            )
          } else {
            grid.text("**",
              x = unit(0.5 + 0.06 * sign(fold_temp[[i]]), "npc") +
                unit(eps * fold_temp[[i]] * 1.25 / 2.05, "cm"),
              y = y_pvalue, gp = gpar(fontsize = 8)
            )
          }
        } else if (pvalue_temp[[i]] <= 0.05) {
          grid.text("*",
            x = unit(0.5 + 0.05 * sign(fold_temp[[i]]), "npc") +
              unit(eps * fold_temp[[i]] * 1.25 / 2.2, "cm"),
            y = y_pvalue, gp = gpar(fontsize = 8)
          )
        }
      }

      # add the fold names
      grid.text(fold_names[[names(fold_names)[k]]],
        x = unit(0.5, "npc"),
        y = unit(1, "npc") + unit(0.85, "line"), gp = gpar(fontsize = 10)
      )

      if (fold_names[[names(fold_names)[k]]] == "30k-500Pa Coll") {
        # add the feature name on top
        grid.text(feature_names[[f]],
          x = unit(3.75, "npc"),
          y = unit(1.08, "npc"), gp = gpar(fontsize = 12)
        )
        # x-axis label
        grid.text(expression("log"[2] * "(ratio of medians)"),
          x = unit(3.75, "npc"),
          y = unit(-0.07, "npc"), gp = gpar(fontsize = 12)
        )
      }
    })
  }

  dev.off()
  keep("folds", "fold_names", "features", "get_color",
    "feature_names", "cell_line_label_colors_main",
    sure = TRUE
  )
}
