# Aim: To perform PCA using median values for the features

# load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# set the working directory as the location of the script -------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# global features -------
features <- c("area", "circularity", "aspect_ratio",
              "cell_stiffness", "motility")

cell.line.tissue <- read.csv("../Figures/cell_line_label_colors.csv",
                             stringsAsFactors = FALSE) 

tissue.color <- cell.line.tissue %>%
  select(tissue, tissue_col) %>%
  distinct() %>%
  pull(tissue_col)
names(tissue.color) <- cell.line.tissue %>%
  select(tissue, tissue_col) %>%
  distinct() %>%
  pull(tissue)

substrate.color <- c("steelblue","darkorange","forestgreen",
                     "brown1","mediumpurple","sienna","orchid")
names(substrate.color) <- c("500Pa Coll", "500Pa FN", "30kPa Coll", 
                            "30kPa FN", "HA FN", "HA Coll", "Glass")

# load the data --------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)  
for (f in features){
  data[[f]] <- read_tsv(paste('../data/',f,'.tsv',sep = ''), 
                        show_col_types = FALSE) %>% as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>% 
    mutate_at(c('cl_id', 'sub_id'), as.factor)
}

# remove cases with less than minimum number of data points -------
# set the cutoff for min. number of data points required for each cell-subs pair
min_cutoff <- 25  

# find the cell-subs pairs with data points < min_cutoff
data_overview <- setNames(vector("list", length(features)), features)
for (f in features) {
  data_overview[[f]] <- data[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>% 
    group_by(cl_sub) %>% 
    summarise(totalcount = n(), .groups = "drop") %>% 
    as_tibble() %>% 
    filter(totalcount < min_cutoff)
}

for (f in features) {
  data[[f]] <- data[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>%
    filter(!(cl_sub %in% data_overview[[f]]$cl_sub)) %>% 
    select(-cl_sub)
}

# median value of features for each cell line-substrate pair -------
medians <- setNames(vector("list", length(features)), features) 
for (f in features){
  medians[[f]] <- data[[f]] %>% 
    group_by(cl_id, sub_id) %>% 
    summarise(!!f := median(feature_value), .groups = "drop") %>% 
    as_tibble()
}

# collect the medians for each feature into one dataframe
# first define a dataframe with cl_id and sub_id columns
# such that each element in unique(data[[1]]$cl_id) x unique(data[[1]]$sub_id) 
# is present. This is done to ensure that all the cell line-substrate pairs 
# are present in the final dataframe even if few of the cell line-substrate pairs
# may have been removed above due to less number of data points (<25)
medians.merged <- expand.grid(unique(data[[1]]$cl_id), unique(data[[1]]$sub_id))
colnames(medians.merged) <- c("cl_id", "sub_id")

for (f in features) {
  medians.merged <- medians.merged %>% 
    left_join(.,medians[[f]], by = c("cl_id", "sub_id"))
}

colnames(medians.merged)[3:7] <- c("Area","Circularity","Aspect Ratio",
                                   "Cell Stiffness","Speed")

# Principal component analysis (PCA) -------------
# prepare the data
pca.data <- medians.merged %>% 
  mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>%
  column_to_rownames(var = "cl_sub") %>%
  select(-cl_id, -sub_id) %>% 
  # remove rows with missing values
  # corresponding to case with number of observations < min_cutoff
  drop_na() %>%
  # convert to matrix
  as.matrix() 

# perform PCA
pca.res <- prcomp(pca.data, center = TRUE, scale. = TRUE)
# Prints variance summary for all principal components.
summary(pca.res) 

# sdev^2 captures these eigenvalues from the PCA result
pc.var <- pca.res$sdev^2 
# the eigenvalues can be used to calculate the percentage variance explained by each PC
pc.per <- round(pc.var/sum(pc.var)*100, 1) 

# plot showing variance explained by each PC
plt.variance <- data.frame(PC = c("PC1","PC2","PC3","PC4","PC5"),
           Variance_explained = pc.per) %>% 
  ggplot(.) +
  aes(x = PC, y = Variance_explained) +
  geom_bar(stat = "identity") + 
  labs(x = "Principal Component", 
       y = "Variance Explained (%)") +
  scale_y_continuous(
    # don't expand y scale at the lower end
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        # eliminate vertical grid lines
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # eliminate panel border
        panel.border = element_blank(),
        # draw x and y axis lines
        axis.line = element_line(color = 'black'))

# Plot loadings for each feature in each PC
plt.loadings <- as.data.frame(pca.res$rotation) %>% 
  rownames_to_column(var = "feature") %>% 
  pivot_longer(cols = PC1:PC4,
               names_to = "PC",
               values_to = "loadings") %>% 
  mutate(loadings = abs(loadings)) %>% 
  ggplot(.) +
  aes(x = feature, 
      y = loadings, 
      group = PC) +
  geom_bar(sta = "identity") +  
  facet_wrap(~PC, ncol = 2) +
  labs(x = "Feature", 
       y = "|Loadings|") +
  scale_y_continuous(
    # don't expand y scale at the lower end
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        # tilt the y axis labels by 30 degrees
        axis.text.x = element_text(angle = 30, hjust = 1),
        # eliminate vertical grid lines
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

# combine variance and loadings plots
ggsave("../Figures/Supplementary_Figure6/pca_variance_loadings.png", 
       plot_grid(plt.variance, plt.loadings, 
                 nrow = 1, ncol = 2,
                 labels = "AUTO", 
                 rel_widths = c(1, 2)), 
       dpi = 300, units = "px", 
       width = 4000, height = 2000)

# 2D PCA plots: one PC vs another PC 
# prepare the data for plotting
pca.res.df <- pca.res$x[,1:4] %>% 
  as_tibble(rownames = NA) %>% 
  add_column(cl_sub = rownames(.)) %>% 
  mutate(cl_id = mapply(function(x) strsplit(x,"_")[[1]][1], cl_sub), 
         sub_id = mapply(function(x) strsplit(x,"_")[[1]][2], cl_sub)) %>%
  # add tissue type
  left_join(., cell.line.tissue, by = "cl_id") %>%
  # keep only the required columns
  select(tissue, cl_id, sub_id, status, PC1, PC2, PC3, PC4) %>% 
  mutate_at(c('tissue','cl_id', 'sub_id'), as.factor)

# dataframe, list, list, string -> 2D scatter ggplot
# plots and saves 2D scatter plots for data points in the df on each pair of PCA components
# points are colored based on the category (eg. tissue) using the manual color list
# variance explained for each PC is displayed on the axes
plotpca <- function(df, var.explained, manual.color, category) {
  
  # PC pairs
  pc.pairs <- combn(colnames(df)[grepl("PC", colnames(df))], 
                    2, simplify = FALSE)
  
  # to store plot for each pair of PCs
  plts <- list()
  
  for (i in 1:length(pc.pairs)) {
    xaxis <- pc.pairs[[i]][1]
    yaxis <- pc.pairs[[i]][2]
    
    plts[[i]] <- ggplot(df) +
      aes(x = !!sym(xaxis), 
          y = !!sym(yaxis), 
          color = !!sym(category), 
          shape = status) +
      geom_point(size=3) +
      scale_color_manual(values = manual.color) +
      xlab(paste0(xaxis," (",var.explained[as.numeric(gsub("PC","",xaxis))],"%",")")) + 
      ylab(paste0(yaxis," (",var.explained[as.numeric(gsub("PC","",yaxis))],"%",")")) +
      coord_fixed() +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            # eliminate rid lines
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none")
  }
  
  combined <- plot_grid(plotlist = plts, 
                        nrow = 3, ncol = 2, 
                        align = "hv")
  
  legend <- get_plot_component(
    plts[[1]] +
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 10),
        legend.title = element_blank()
      ),
    "guide-box-top",
    return_all = TRUE
  )
  
  return(plot_grid(legend, combined, 
                   nrow = 2, rel_heights = c(0.1, 1)))
}

tissue.pca.plots <- plotpca(pca.res.df, pc.per, tissue.color, "tissue")     

ggsave("../Figures/Supplementary_Figure6/pca_tissue.png", 
       tissue.pca.plots, 
       dpi = 300, units = "px", width = 2500, height = 3500)

subs.pca.plots <- plotpca(pca.res.df, pc.per, substrate.color, "sub_id")     

ggsave("../Figures/Supplementary_Figure6/pca_substrate.png", 
       subs.pca.plots, 
       dpi = 300, units = "px", width = 2500, height = 3500)





