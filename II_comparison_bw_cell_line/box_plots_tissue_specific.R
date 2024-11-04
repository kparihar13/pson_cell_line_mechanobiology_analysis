# Aim: To generate boxplots for each cell line in a tissue of interest
# marked with significance levels based on BH adjusted p-values
# calculated in the signif_test_for_cancer_vs_noncancer.R script

# load packages -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ggsignif)
})

# set the working directory same as the location of the script -----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Global parameters ------
features <- c("area", "circularity", "aspect_ratio","cell_stiffness", "motility")
subs <- c("HA Coll", "HA FN", "500Pa Coll", "500Pa FN", 
          "30kPa Coll", "30kPa FN", "Glass")

# colors for boxplots
cell.line.colors <- c("firebrick2", "orange2", "forestgreen", "mediumorchid", 
                      "turquoise3", "palevioletred1")

# define the tissue of interest
tissues_of_interest <- c("Prostate", "Pancreas", "Lung", "Breast")

# get the cell lines for the tissue of interest
cell_lines <-  read.csv("../Figures/cell_line_label_colors.csv", 
                        stringsAsFactors = FALSE) %>%
  # filter to keep only the cell lines of interest
  filter(tissue %in% tissues_of_interest)

# minimum cutoff for number of observations to be included in the analysis
min_cutoff <- 25

# load data ---------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)  
data_overview <- setNames(vector("list", length(features)), features)

for (f in features){
  data[[f]] <- read_tsv(paste('../data/',f,'.tsv',sep = ''), 
                        show_col_types = FALSE) %>% 
    as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>%
    # filter to only the cell lines of interest
    filter(cl_id %in% cell_lines$cl_id) %>%
    # fix the order of substrates for plotting purposes
    mutate(sub_id = fct_relevel(sub_id, subs)) %>%
    # add column for status and label (for plotting purposes)
    inner_join(., cell_lines %>% select(cl_id, label, status), by = "cl_id") %>%
    # fix the order for status and labels
    mutate(status = fct_relevel(status, c("Non-cancer", "Cancer")),
           label = fct_relevel(label, cell_lines$label))
  
  # find cases with number of observations < min_cutoff
  data_overview[[f]] <- data[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>% 
    group_by(cl_sub) %>% 
    summarise(totalcount = n(), .groups = "drop") %>% 
    as_tibble() %>% 
    filter(totalcount < min_cutoff)
  
  # remove cases with < min_cutoff values
  data[[f]] <- data[[f]] %>% 
    mutate(cl_sub = paste(cl_id,"_",sub_id,sep='')) %>%
    filter(!(cl_sub %in% data_overview[[f]]$cl_sub)) %>% 
    select(-cl_sub)
}

# Boxplots with p-values  ----------

# Numeric -> String
# converts numeric pvalue to string 
# - '***' if pvalue < 0.001
# - '**' if pvalue < 0.01 
# - '*' if pvalue < 0.05
# - 'ns' otherwise
pval.signif <- function(x) {
  if (x < 0.001) { return("***") }
  else if (x < 0.01) { return("**") }
  else if (x < 0.05) { return("*") }
  else { return("ns") }
}

# vectorize the function
pval.signif.vec <- Vectorize(pval.signif)

# y-axis labels
y.labels <- list("area" = expression("Area (" * mu*"m)"^2), 
                 "circularity" = "Circularity",
                 "aspect_ratio" = "Aspect Ratio", 
                 "cell_stiffness" = "Cell Stiffness (Pa)",
                 "motility" = expression("Speed (" * mu*"m/hr)"))

# for pvalue bars above the boxplots
# dependent on tissue of interest due to differing ranges of y-axis 
offset <- list("Prostate" = list("area" = 130, 
                                 "circularity" = 0.04, 
                                 "aspect_ratio" = 0.3, 
                                 "cell_stiffness" = 500,
                                 "motility" = 5.5),
               "Pancreas" = list("area" = 300, 
                                 "circularity" = 0.04, 
                                 "aspect_ratio" = 0.15,
                                 "cell_stiffness" = 900,
                                 "motility" = 5.5),
               "Lung" = list("area" = 85, 
                             "circularity" = 0.04, 
                             "aspect_ratio" = 0.065, 
                             "cell_stiffness" = 250,
                             "motility" = 4),
               "Breast" = list("area" = 105, 
                               "circularity" = 0.04, 
                               "aspect_ratio" = 0.2, 
                               "cell_stiffness" = 850,
                               "motility" = 10))

for (t in tissues_of_interest) {
  
  # cell lines for the tissue of interest
  cell_lines_tissue <- cell_lines %>% 
    filter(tissue == t) %>% 
    pull(cl_id)
  
  # load the pre-computed BH corrected p-values
  load(paste0(t,"_BH_adj_pvals.RData"))

  # create annotation dataframe for adding p-values to the boxplots
  # adapted from https://const-ae.github.io/ggsignif/
  #pvalue.df <- setNames(vector("list", length(features)), features)
  
  for (f in features) {
    
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
      manual.colors <- setNames(vector("list",length(cell_order)), cell_order)
      manual.colors <- cell.line.colors[1:length(manual.colors)]
      
      # keep only subs with more than one cell line having 
      # number of observations > min_cutoff
      subs.to.keep <- data[[f]] %>%
        select(all_of(c("sub_id", "cl_id"))) %>%
        distinct() %>%
        # only keep cell lines of interest
        filter(cl_id %in% cell_lines_tissue) %>%
        group_by(sub_id) %>%
        summarise(count = n()) %>%
        filter(count != 1)
      
      # create the boxplot
      plt <- data[[f]] %>% 
        # filter to keep only the substrates of interest and 
        # cell lines of interest
        filter((sub_id %in% subs.to.keep$sub_id),
               (cl_id %in% cell_lines_tissue)) %>%
        ggplot(.) +
        aes(x=sub_id, 
            y=feature_value, 
            fill=factor(label, 
                        levels = cell_order)) + 
        geom_boxplot(position="dodge", 
                     outliers = F, 
                     staplewidth = 0.5) +
        labs(x = "Substrate", 
             y = y.labels[[f]], 
             fill = paste0(t,": ")) +
        # Specify colours
        scale_fill_manual(values=manual.colors) +
        theme_bw() +
        theme(legend.position = "top") + 
        guides(fill = guide_legend(nrow = 1))
      
      # extract plot details as dataframe to add
      # p-value annotation on the plots
      gg <- ggplot_build(plt)
      
      # create dataframe for pvalue annotation
      anno.df <- data[[f]] %>% 
        # keep only cell lines of interest
        filter(cl_id %in% cell_lines_tissue) %>%
        # only need cl_id, label, sub_id, status
        select(all_of(c("cl_id", "sub_id", "label", "status"))) %>% 
        # remove duplicate rows
        distinct() %>%
        # filter to keep only substrates of interest
        filter(sub_id %in% subs.to.keep$sub_id) %>%
        # make cols factor for ordering
        mutate(sub_id = factor(sub_id, levels = subs),
               label = factor(label, levels = cell_order)) %>%
        arrange(sub_id, label)
      
      # get the x position (center) of all the boxplots
      anno.df$x <- gg$data[[1]]$x
      # get the topmost point in the boxplots
      anno.df$ymax <- gg$data[[1]]$ymax
      # determine the topmost point within a substrate type
      y.temp <- anno.df %>% 
        group_by(sub_id) %>% 
        summarise(y.pos = max(ymax))
      # add back to the anno dataframe
      anno.df <- anno.df %>% 
        left_join(., y.temp, by = "sub_id")
      
      # create an empty dataframe with same columns as anno.df
      # this is needed for case of more than one non-cancer cell line
      anno.df.ext <- anno.df[FALSE,]
      # iterate substrate wise to fill the dataframe
      for (s in subs) {
        # now iterate non-cancer cell line wise
        for (nc in nc_cells$label) {
          temp <- anno.df %>%
            # keep only s, one of the non-cancer cell lines, and all cancer cell lines
            filter(sub_id == s & label %in% c(nc, c_cells))
          
          # Determine the start and end of brackets
          # start of the bracket is the non-cancer cell line
          start.temp <- temp %>% 
            filter(status == "Non-cancer") %>%
            mutate(start = x, status = "Cancer") %>%
            select(all_of(c("sub_id","status","start")))
          
          temp <- temp %>% 
            # add the start position in the anno dataframe
            left_join(., start.temp, by = c("sub_id","status")) %>%
            # end is just the center of boxplots 
            mutate(end = x) %>%
            # add a column called noncan for later joining with pvalues
            add_column(noncan = rep(nc, nrow(.)))
          
          # add to anno.df.ext
          anno.df.ext <- rbind(anno.df.ext, temp)
        }
      }
      
      # create dataframe containing pvalues, to be merged with anno.df.ext
      pvalue.df <- NULL
      for (nc in 1:nrow(nc_cells)) {
        # extarct the adj pval for a non-cancer cell line
        temp.anno <- stat.test.bh[[f]][[nc_cells$cl_id[nc]]] %>%
          t() %>% 
          as.data.frame() %>%
          rownames_to_column(var = "sub_id") %>%
          # pivot to longer format dataframe
          pivot_longer(cols = c_cells, 
                       names_to = "cl_id", 
                       values_to = "adj.pval") %>%
          # remove cases with NA which correspond to number of observations < min_cutoff
          drop_na() 
        
        # add column to specify which nc cell
        temp.anno$noncan <- rep(nc_cells$label[nc], nrow(temp.anno))
        
        # add temp_anno to apvalue.df
        # if pvalue.df is empty then assign temp.anno
        if(is.null(pvalue.df)){ pvalue.df <- temp.anno }
        # else concatenate the temp.anno to pvalue.df
        else { pvalue.df <- rbind(pvalue.df, temp.anno) }
      }
      
      # change numeric pvalues to ***, **, *, ns
      pvalue.df <- pvalue.df %>%
        mutate(adj.pval.anno = pval.signif.vec(adj.pval))
      
      anno.df.ext <- anno.df.ext %>% 
        # join to annotation_df to get the pvalues
        left_join(., pvalue.df, by = c("cl_id","sub_id","noncan")) %>%
        # drop rows with NA, basically corresponds to non-cancer cell line rows
        # as pvalue is for non-cancer vs cancer cell lines 
        # so cl_id column in pvalue.df only has cancer cell lines
        drop_na() %>%
        # remove non-relevant columns
        # "label" can be removed as the dataframe only has cancer cell lines
        # because the pvalue.df has only cancer cell lines
        select(-all_of(c("label","x", "ymax", "status", "adj.pval"))) %>%
        # order sub_id, noncan, cl_id factor for plotting purposes
        mutate(sub_id = factor(sub_id, levels = subs),
               noncan = factor(noncan, levels = nc_cells$label),
               cl_id = factor(cl_id, levels = c_cells)) %>%
        # Show signif for one nc and then next: 1st by substrate, 2nd by nc_cells, 3rd by c_cells
        # Show signif alternatingly for nc: 1st by substrate, 2nd by c_cells, 3rd by nc_cells
        arrange(sub_id, cl_id, noncan)
        
      # need a counter column for defining the offset of brackets
      anno.df.ext$subs_count <- rep(0, nrow(anno.df.ext))
      for (s in subs) {
        # just save (0,1,2,...till number of cancer lines in substrate s)
        anno.df.ext$subs_count[which(anno.df.ext$sub_id == s)] <- 1:table(anno.df.ext$sub_id)[s]
      }
      anno.df.ext <- anno.df.ext %>%
        # define the y position of the bracket using the counter 
        mutate(y.pos = y.pos + offset[[t]][[f]]*subs_count)
      
      # add pvalue brackets to the plot
      plt <- plt +
        geom_signif(annotations = anno.df.ext$adj.pval.anno,
                    y_position = anno.df.ext$y.pos,
                    xmin = anno.df.ext$start,
                    xmax = anno.df.ext$end,
                    tip_length = 0, 
                    # to bring text closer to the bracket
                    vjust = 0.4,
                    # text size
                    textsize = 2.75)
      
      # save the plot
      if (f %in% c("circularity", "aspect_ratio")){
        png(paste0("../Figures/Supplementary_Figure5/",t, "_", f, ".png"), 
            res = 300, width = 2300, height = 1500)
      } else {
        png(paste0("../Figures/Figure4/",t, "_", f, ".png"), 
            res = 300, width = 2300, height = 1500)
      }
      print(plt)
      dev.off()
    }
  }
}