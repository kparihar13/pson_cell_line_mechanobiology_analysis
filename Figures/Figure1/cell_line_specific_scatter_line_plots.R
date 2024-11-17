# Had to reorient "Glass" label in Gimp after hiding tick marks
# before replotting, shift "Glass" to the left so the "a" is aligned 
# with third tick mark

# Aim: plot scatter plots with lines for cell line specific median values

# load packages ------
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(plotrix)
  library(cowplot)
})

# set the working directory same as the location of the script ------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# some global variables ---------
features <- c("area", "circularity", "aspect_ratio",
              "cell_stiffness", "motility")
cell_line <- 'T98G'
subs <- data.frame('sub_id' = c('30kPa Coll', '30kPa FN', 
                                '500Pa Coll', '500Pa FN',
                                'HA Coll', 'HA FN', 'Glass'),
                   'subs_kPa' = c(30, 30, 0.5, 0.5,0.3, 0.3, 500),
                   'subs_type_I' = as.factor(c('Collagen', 'Fibronectin', 
                                               'Collagen', 'Fibronectin',
                                               'Hyaluronic acid', 'Hyaluronic acid', 
                                               'Glass')),
                   'subs_type_II' = as.factor(c('Collagen', 'Fibronectin', 
                                                'Collagen', 'Fibronectin',
                                                'Collagen', 'Fibronectin', 
                                                'Glass'))) 

# load the data --------
# define an empty named list for storing data
data <- setNames(vector("list", length(features)), features)  
for (f in features){
  data[[f]] <- read_tsv(paste0('../../data/',f,'.tsv'), 
                        show_col_types = FALSE) %>% 
    as_tibble()
  colnames(data[[f]]) <- c("cl_id", "sub_id", "feature_value")
  data[[f]] <- data[[f]] %>% 
    mutate_at(c('cl_id', 'sub_id'), as.factor)
}

# remove cases with less than minimum number of data points -------
min_cutoff <- 25  # set the cutoff for min. # of data points required

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

# Keep only the cell line specific data and calculate substrate wise median ------
for (f in features) {
  data[[f]] <- data[[f]] %>% 
    filter(cl_id == cell_line) %>%
    group_by(sub_id) %>%
    summarize(median = median(feature_value)) %>%
    inner_join(., subs, by = 'sub_id')
}

# Plot with axis breaks using base R ------

# from https://stackoverflow.com/questions/47890742/logarithmic-scale-plot-in-r
# function to create tick marks on x-axis
log10Tck <- function(side, type){
  # switch() to select one of a list of alternatives
  lim <- switch(side, 
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceiling(lim[2])
  return(switch(type, 
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}

plot_feature <- function(data, f, margins, offset) {
  x.axis.limits <- c(0.15, max(data[[f]]$subs_kPa)*1.25)
  y.axis.limits <- c(min(data[[f]]$median)*0.997, max(data[[f]]$median)*1.003)
  
  if (f == "area") { 
    ylabel <- bquote(.(yLab) ~ "Area (" * mu*"m"^2*")") } else if (f == "circularity") {
      ylabel <- paste(yLab,"Circularity")} else if (f== "aspect_ratio") {
        ylabel <- paste(yLab,"Aspect Ratio")} else if (f == "cell_stiffness") {
          ylabel <- paste(yLab,"Cell Stiffness (Pa)")} else if (f == "motility") {
            ylabel <- bquote(.(yLab) ~ "Speed (" * mu*"m/hr)")}
  
  par(mar=margins)
  
  # plot 500Pa and 30kPa Coll with line
  temp <- data[[f]] %>% filter(grepl("Pa Coll",sub_id))
  plot(temp$subs_kPa, temp$median, 
       type = "b", col = "blue", pch = 19, cex = 2.3, lwd = 2.5, cex.lab=1.5, 
       log = 'x', axes = F, 
       xlim = x.axis.limits, ylim = y.axis.limits,
       xlab= 'Substrate Stiffness (kPa)', ylab=ylabel)
  # add grid lines
  grid(lty = 'solid')
  
  # plot 500Pa and 30kPa FN with line
  temp <- data[[f]] %>% filter(grepl("Pa FN",sub_id))
  points(temp$subs_kPa, temp$median, 
         type = "b", col = "red", pch = 19, cex = 2.3, lwd = 2.5,
         xlim = x.axis.limits, ylim = y.axis.limits)
  
  # plot HA Coll
  temp <- data[[f]] %>% filter(grepl("HA Coll",sub_id))
  points(temp$subs_kPa, temp$median, 
         col = "blue", pch = 17, cex = 2.3, 
         xlim = x.axis.limits, ylim = y.axis.limits)
  
  # plot HA FN
  temp <- data[[f]] %>% filter(grepl("HA FN",sub_id))
  points(temp$subs_kPa, temp$median, 
         col = "red", pch = 17, cex = 2.3, 
         xlim = x.axis.limits, ylim = y.axis.limits)
  
  # plot Glass
  temp <- data[[f]] %>% filter(grepl("Glass",sub_id))
  points(temp$subs_kPa, temp$median, 
         col = "black", pch = 15, cex = 2.3, 
         xlim = x.axis.limits, ylim = y.axis.limits)
  
  # add x-axis
  axis(1, at=log10Tck('x','major'), tcl= 0.75,
       labels = round(log10Tck('x','major')), cex.axis=1.2) 
  axis(1, at=log10Tck('x','minor'), tcl= 0.4, labels = NA) 
  # add y-axis
  axis(2, tcl= 0.75, cex.axis=1.2) 
  # add box around the plot
  box()
  # add axis break for Glass on x-axis
  axis.break(1, 200, style="slash", bgcol="white", brw=0.04) 
  text(y=par()$usr[3]-offset, x=1200, "Glass", 
       xpd=TRUE, adj=1, cex=1.2)
}

# plotting median
yLab <- "Median"

png(paste0(cell_line,".png"),
    res = 300, width = 2500, height = 2500)
par(mfrow=c(2,2), mgp=c(2.2,0.5,0))
plot_feature(data, "area", c(2,4,4,2), 45)
plot_feature(data,"aspect_ratio", c(2,4,4,2), 0.1)
plot_feature(data,"cell_stiffness", c(4,4,2,2), 135)
plot_feature(data,"motility", c(4,4,2,2), 2)
dev.off()

# Get the legend -----------
# remove factor for subs_type_II
temp <- data[[f]] %>% 
  select(-subs_type_I) %>%
  mutate(subs_type_II = as.character(subs_type_II))

temp$subs_type_II[which(temp$sub_id == "HA Coll")] <- "HA Collagen"
temp$subs_type_II[which(temp$sub_id == "HA FN")] <- "HA Fibronectin"

lgd <- temp %>%
  # change subs_type_II to factor
  mutate(subs_type_II = as.factor(subs_type_II)) %>%
  ggplot(.) +
  aes(x = subs_kPa, y = median, shape = subs_type_II, color = subs_type_II) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("Collagen" = 19, "Fibronectin" = 19,
                                "HA Collagen" = 17, "HA Fibronectin" = 17, 'Glass' = 15)) +
  scale_color_manual(values = c("Collagen" = 'blue', "Fibronectin" = 'red',
                                "HA Collagen" = 'blue', "HA Fibronectin" = 'red',
                                'Glass' = 'black')) +
  theme_bw()

# extract legend from ggplot using cowplot 
legend <- get_plot_component(
  lgd +
    theme(
      legend.text = element_text(size = 12),
      legend.title = element_blank()
    ),
  "guide-box-right"
)

# save the legend
ggsave("substrate_legend.png",
       plot_grid(legend, ncol = 1),
       device = "png", 
       width = 600, 
       height = 480,
       units = "px"
)
