# Aim: Plots for Supplementary Figure 4 (change in physical feature with change in substrate)

# load packages ---------
suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ComplexHeatmap)
})

# set directory as the source file location ----------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# global parameters -----------
features <- c("circularity", "aspect_ratio")

feature_names <- list(
  "circularity" = "Circularity",
  "aspect_ratio" = "Aspect Ratio"
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

# text color cell line names
cell_line_label_colors_main <- read.csv("../Figures/cell_line_label_colors.csv",
  stringsAsFactors = FALSE
) %>%
  # change cell_id to factor with levels as the current order of cell lines
  # done to ensure the order of cell lines in the plots
  mutate(cl_id = factor(cl_id, levels = unique(cl_id)))

# plot fold change ---------------

# Dataframe, String -> Dataframe
# specifies the fold change based color for bar chart
get_color <- function(data, column_name) {
  data <- data %>%
    # blue if > 0, else red
    mutate(
      anno_bar = !!sym(column_name),
      color = ifelse(!!sym(column_name) > 0, "#005AB5", "#DC3220")
    ) %>%
    select(anno_bar, color)

  # if fold change is to high or to low, then just take it as 2 or -2 for plotting purposes
  data[which(data$anno_bar > 2), 1] <- 2
  data[which(data$anno_bar < -2), 1] <- -2

  return(data)
}

# Dataframe obtained from get_color() -> ComplexHeatmap barplot annotation object
# produces a barplot annotation object for the heatmap
get_bar_plot <- function(df) {
  row_bar <- anno_barplot(unlist(df$anno_bar),
    baseline = 0,
    ylim = c(-2, 2),
    axis_param = list(
      side = "bottom",
      labels_rot = 0,
      at = seq(-2, 2, 0.5),
      labels = c(-2, "", -1, "", 0, "", 1, "", 2),
      gp = gpar(fontsize = 9)
    ),
    border = FALSE,
    bar_width = 0.7,
    gp = gpar(fill = unlist(df$color))
  )

  return(row_bar)
}

# Row annotation oject, boolean, dataframe -> Heatmap object
# produces a heatmap object with the row annotation, boolean specifies whether
# to include rownames or not, and dataframe contains cell line names and tissue colors
# used as text color for rownames
create_heatmap <- function(row_anno, include_rownames, cell_line_label_colors) {
  # empty matrix as only interested in the annotation
  ht <- Heatmap(matrix(nc = 0, nr = nrow(cell_line_label_colors)),
    rect_gp = gpar(type = "none"),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = include_rownames,
    show_column_names = FALSE,
    row_labels = cell_line_label_colors$label,
    row_names_gp = gpar(
      fontsize = 10,
      col = cell_line_label_colors$tissue_col
    ),
    row_names_side = "left",
    left_annotation = row_anno,
    height = unit(13, "cm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  return(ht)
}

# folds of interest
fold_names <- list(
  "colfold" = "30k-500Pa Coll",
  "fnfold" = "30k-500Pa FN",
  "ha500colfold" = "HA-500Pa Coll",
  "ha500fnfold" = "HA-500Pa FN",
  "hyalfold" = "HA Coll-FN"
)

for (f in features) {
  # load the fold change values and respective pvalues
  load(paste(f, "_folds.RData", sep = ""))

  fold_value <- fold_value %>%
    rownames_to_column(var = "cl_id") %>%
    # order in desired order of cell lines
    mutate(cl_id = factor(cl_id, levels = cell_line_label_colors_main$cl_id)) %>%
    arrange(cl_id) %>%
    column_to_rownames("cl_id")

  pvalue_bh_corrected <- pvalue_bh_corrected %>%
    rownames_to_column(var = "cl_id") %>%
    # order in desired order of cell lines
    mutate(cl_id = factor(cl_id, levels = cell_line_label_colors_main$cl_id)) %>%
    arrange(cl_id) %>%
    column_to_rownames("cl_id")

  # log2 transform the ratio values
  fold_value <- log2(fold_value)

  # change the colnames for fold_value and pvalue_bh_corrected for plotting convenience
  # used in annotation creation
  for (fc in names(folds)) {
    colnames(fold_value)[which(colnames(fold_value) == fc)] <- folds[[fc]]
    colnames(pvalue_bh_corrected)[which(colnames(pvalue_bh_corrected) == fc)] <- folds[[fc]]
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
    colfold = get_bar_plot(bar_colfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["colfold"]] <- create_heatmap(
    ha_colfold, TRUE,
    cell_line_label_colors
  )

  # 2nd fold: 500Pa FN -> 30kPa FN (fnfold)
  # specify the fold change based color for bar chart
  bar_fnfold <- get_color(fold_value, "fnfold")
  # define the heatmap annotation
  ha_fnfold <- rowAnnotation(
    fnfold = get_bar_plot(bar_fnfold),
    show_annotation_name = FALSE, width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["fnfold"]] <- create_heatmap(
    ha_fnfold, FALSE,
    cell_line_label_colors
  )

  # 3rd fold: 500Pa Coll -> HA Coll (ha500colfold)
  # specify the fold change based color for bar chart
  bar_ha500colfold <- get_color(fold_value, "ha500colfold")
  # define the heatmap annotation
  ha_ha500colfold <- rowAnnotation(
    ha500colfold = get_bar_plot(bar_ha500colfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["ha500colfold"]] <- create_heatmap(
    ha_ha500colfold, FALSE,
    cell_line_label_colors
  )

  # 4th fold: 500Pa FN -> HA FN (ha500fnfold)
  # specify the fold change based color for bar chart
  bar_ha500fnfold <- get_color(fold_value, "ha500fnfold")
  # define the heatmap annotation
  ha_ha500fnfold <- rowAnnotation(
    ha500fnfold = get_bar_plot(bar_ha500fnfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["ha500fnfold"]] <- create_heatmap(
    ha_ha500fnfold, FALSE,
    cell_line_label_colors
  )

  # 5th fold: HA FN -> HA Coll (ha500fnfold)
  # specify the fold change based color for bar chart
  bar_hyalfold <- get_color(fold_value, "hyalfold")
  # define the heatmap annotation
  ha_hyalfold <- rowAnnotation(
    hyalfold = get_bar_plot(bar_hyalfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["hyalfold"]] <- create_heatmap(
    ha_hyalfold, FALSE,
    cell_line_label_colors
  )

  # for empty gap between plots
  ha_gap <- rowAnnotation(empty = anno_empty(
    border = FALSE,
    width = unit(0.2, "mm")
  ))

  # file to save the plot
  if (f == "aspect_ratio") {
    png(paste("../Figures/Supplementary_Figure4/", f, "_1.png", sep = ""),
      res = 300, width = 2200, height = 2000
    )
  } else {
    png(paste("../Figures/Supplementary_Figure5/", f, "_1.png", sep = ""),
        res = 300, width = 2200, height = 2000
    )
  }
  # combine the plots
  draw(
    ht_list[["colfold"]] +
      ha_gap + ht_list[["fnfold"]] +
      ha_gap + ht_list[["ha500colfold"]] +
      ha_gap + ht_list[["ha500fnfold"]] +
      ha_gap + ha_gap + ht_list[["hyalfold"]],
    auto_adjust = FALSE
  )

  # add the significance level markers, fold names, feature name and x-axis label
  for (k in 1:length(names(fold_names))) {
    decorate_annotation(names(fold_names)[k], {
      pvalue_temp <- pvalue_bh_corrected %>%
        select(names(fold_names)[k]) %>%
        pull()
      fold_temp <- fold_value %>%
        select(names(fold_names)[k]) %>%
        pull()

      # place *, **, *** based on p-values
      if (sum(is.na(pvalue_temp)) > 0) {
        pvalue_temp[which(is.na(pvalue_temp))] <- 1000
      }
      eps <- 1
      for (i in 1:length(pvalue_temp)) {
        y_pvalue <- unit((length(pvalue_temp) - i + 0.4) * 130 / length(pvalue_temp), "mm")
        if (pvalue_temp[[i]] <= 0.01) {
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
        } else if (pvalue_temp[[i]] <= 0.05) {
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
        } else if (pvalue_temp[[i]] <= 0.1) {
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
        y = unit(1, "npc") + unit(0.85, "line"),
        gp = gpar(fontsize = 10)
      )

      if (fold_names[[names(fold_names)[k]]] == "30k-500Pa Coll") {
        # add the feature name on top
        grid.text(feature_names[[f]],
          x = unit(2.97, "npc"),
          y = unit(1.08, "npc"),
          gp = gpar(fontsize = 12)
        )
        # x-axis label
        grid.text(expression("log"[2] * "(ratio of medians)"),
          x = unit(2.97, "npc"),
          y = unit(-0.07, "npc"),
          gp = gpar(fontsize = 12)
        )
      }
    })
  }

  dev.off()
}


# folds of interest
fold_names <- list(
  "fn500colfold" = "500Pa Coll-FN",
  "fn30kcolfold" = "30kPa Coll-FN",
  "ha30kcolfold" = "30kPa-HA Coll",
  "ha30kfnfold" = "30kPa-HA FN",
  "glasscolfold" = "Glass-30kPa Coll",
  "glassfnfold" = "Glass-30kPa FN"
)

for (f in features) {
  # load the fold change values and respective pvalues
  load(paste(f, "_folds.RData", sep = ""))

  fold_value <- fold_value %>%
    rownames_to_column(var = "cl_id") %>%
    # order in desired order of cell lines
    mutate(cl_id = factor(cl_id, levels = cell_line_label_colors_main$cl_id)) %>%
    arrange(cl_id) %>%
    column_to_rownames("cl_id")

  pvalue_bh_corrected <- pvalue_bh_corrected %>%
    rownames_to_column(var = "cl_id") %>%
    # order in desired order of cell lines
    mutate(cl_id = factor(cl_id, levels = cell_line_label_colors_main$cl_id)) %>%
    arrange(cl_id) %>%
    column_to_rownames("cl_id")

  # log2 transform the ratio values
  fold_value <- log2(fold_value)

  # change the colnames for fold_value and pvalue_bh_corrected for plotting convenience
  for (fc in names(folds)) {
    colnames(fold_value)[which(colnames(fold_value) == fc)] <- folds[[fc]]
    colnames(pvalue_bh_corrected)[which(colnames(pvalue_bh_corrected) == fc)] <- folds[[fc]]
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

  # 1st fold: 500Pa FN -> 500Pa Coll (fn500colfold)
  # specify the fold change based color for bar chart
  bar_fn500colfold <- get_color(fold_value, "fn500colfold")
  # define the heatmap annotation
  ha_fn500colfold <- rowAnnotation(
    anno_emp = anno_empty(width = unit(0.1, "cm"), border = FALSE),
    fn500colfold = get_bar_plot(bar_fn500colfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["fn500colfold"]] <- create_heatmap(
    ha_fn500colfold, TRUE,
    cell_line_label_colors
  )

  # 2nd fold: 30kPa FN -> 30kPa Coll (fn30kcolfold)
  # specify the fold change based color for bar chart
  bar_fn30kcolfold <- get_color(fold_value, "fn30kcolfold")
  # define the heatmap annotation
  ha_fn30kcolfold <- rowAnnotation(
    fn30kcolfold = get_bar_plot(bar_fn30kcolfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["fn30kcolfold"]] <- create_heatmap(
    ha_fn30kcolfold, FALSE,
    cell_line_label_colors
  )

  # 3rd fold: HA Coll -> 30kPa Coll (ha30kcolfold)
  # specify the fold change based color for bar chart
  bar_ha30kcolfold <- get_color(fold_value, "ha30kcolfold")
  # define the heatmap annotation
  ha_ha30kcolfold <- rowAnnotation(
    ha30kcolfold = get_bar_plot(bar_ha30kcolfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["ha30kcolfold"]] <- create_heatmap(
    ha_ha30kcolfold, FALSE,
    cell_line_label_colors
  )

  # 4th fold: HA FN -> 30kPa FN (ha30kfnfold)
  # specify the fold change based color for bar chart
  bar_ha30kfnfold <- get_color(fold_value, "ha30kfnfold")
  # define the heatmap annotation
  ha_ha30kfnfold <- rowAnnotation(
    ha30kfnfold = get_bar_plot(bar_ha30kfnfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["ha30kfnfold"]] <- create_heatmap(
    ha_ha30kfnfold, FALSE,
    cell_line_label_colors
  )

  # 5th fold: 30kPA Coll -> Glass (glasscolfold)
  # specify the fold change based color for bar chart
  bar_glasscolfold <- get_color(fold_value, "glasscolfold")
  ## define the heatmap annotation
  ha_glasscolfold <- rowAnnotation(
    glasscolfold = get_bar_plot(bar_glasscolfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["glasscolfold"]] <- create_heatmap(
    ha_glasscolfold, FALSE,
    cell_line_label_colors
  )

  # 6th fold: 30kPa FN -> Glass (glassfnfold)
  # specify the fold change based color for bar chart
  bar_glassfnfold <- get_color(fold_value, "glassfnfold")
  # define the heatmap annotation
  ha_glassfnfold <- rowAnnotation(
    glassfnfold = get_bar_plot(bar_glassfnfold),
    show_annotation_name = FALSE,
    width = unit(2.5, "cm")
  )
  # define the heatmap
  ht_list[["glassfnfold"]] <- create_heatmap(
    ha_glassfnfold, FALSE,
    cell_line_label_colors
  )

  # for empty gap between heatmaps
  ha_gap <- rowAnnotation(empty = anno_empty(
    border = FALSE,
    width = unit(0.2, "mm")
  ))

  # file to save the plot
  if (f == "aspect_ratio") {
    png(paste("../Figures/Supplementary_Figure4/", f, "_2.png", sep = ""),
      res = 300, width = 2500, height = 1800
    )
  } else {
    png(paste("../Figures/Supplementary_Figure5/", f, "_2.png", sep = ""),
        res = 300, width = 2500, height = 1800
    )
  }
  # combine the plots
  draw(
    ht_list[["fn500colfold"]] + ha_gap +
      ht_list[["fn30kcolfold"]] + ha_gap +
      ht_list[["ha30kcolfold"]] + ha_gap +
      ht_list[["ha30kfnfold"]] + ha_gap +
      ht_list[["glasscolfold"]] + ha_gap + ha_gap +
      ht_list[["glassfnfold"]],
    auto_adjust = FALSE
  )


  # add the significance level markers, fold names, feature name and x-axis label
  for (k in 1:length(names(fold_names))) {
    decorate_annotation(names(fold_names)[k], {
      pvalue_temp <- pvalue_bh_corrected %>%
        select(names(fold_names)[k]) %>%
        pull()
      fold_temp <- fold_value %>%
        select(names(fold_names)[k]) %>%
        pull()

      # place *, **, *** based on p-values
      if (sum(is.na(pvalue_temp)) > 0) {
        pvalue_temp[which(is.na(pvalue_temp))] <- 1000
      }
      eps <- 1
      for (i in 1:length(pvalue_temp)) {
        y_pvalue <- unit((length(pvalue_temp) - i + 0.4) * 130 / length(pvalue_temp), "mm")
        if (pvalue_temp[[i]] <= 0.01) {
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
        } else if (pvalue_temp[[i]] <= 0.05) {
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
        } else if (pvalue_temp[[i]] <= 0.1) {
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
        y = unit(1, "npc") + unit(0.85, "line"),
        gp = gpar(fontsize = 10)
      )

      if (fold_names[[names(fold_names)[k]]] == "500Pa Coll-FN") {
        # add the x-axis label
        grid.text(expression("log"[2] * "(ratio of medians)"),
          x = unit(3.6, "npc"),
          y = unit(-0.07, "npc"),
          gp = gpar(fontsize = 12)
        )
      }
    })
  }

  dev.off()
}
