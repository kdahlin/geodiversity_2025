#Script to plot metrics against increasing scale
library(ggplot2)
library(openxlsx2)
library(stringr) #Allows python-style string formatting
library(dplyr)

#####################################################################
#Set working directory
setwd("C:/Users/hallerdi/Documents/Dillon_Classes/geodiversity_2025")
#####################################################################

# read in the long data
in.data <- read.csv("./results/all_metrics_long.csv")

# split out scene names
files <- in.data$file
parts <- strsplit(files, "[_\\.]")
scenes <- unlist(lapply(parts, function(x) {paste(x[1],x[6], sep = "_")}))
in.data$scene <- scenes

#Get unique scene and variable names
unique_scenes <- unique(scenes)
unique_variables <- unique(in.data$func)

#Need to create normalized values for each scene
#n_value = (value - 1mvalue)/1mvalue
norm.data <- in.data %>%
  group_by(scene, func) %>%
  mutate(
    norm_value = (value - value[res == 1]) / value[res == 1]
  ) %>%
  ungroup()

#Different form of normalization. We normalize between all scenes at the same scale for the same function
scale.norm.data <- in.data %>%
  group_by(res, func) %>%
  mutate(
    max = max(value),
    min = min(value),
    norm_value = (value - min) / (max - min)
  ) %>%
  ungroup() 


#Iterate over variables, saving out a raw value and a normalized value for each metric
for (variable in unique_variables){
  #Subset data to only that currently of interest
  subset <- norm.data[norm.data$func==variable,]
  ggplot(data = subset)+
    geom_line(mapping = aes(x = res, y = norm_value*100, color = site, group = scene))+
    labs(x = "Pixel size (meters)", y = "% change relative to 1-meter pixel size", title = variable)
  ggsave(str_glue("./sandbox/results/scales/{variable}_normalized.png"), width = 5, height = 3.5, units = "in", limitsize = FALSE)
  ggplot(data = subset)+
    geom_line(mapping = aes(x = res, y = value, color = site, group = scene))+
    labs(x = "Pixel size (meters)", y = "Value", title = variable)
  ggsave(str_glue("./sandbox/results/scales/{variable}_rawvalue.png"), width = 5, height = 3.5, units = "in", limitsize = FALSE)
  subset <- scale.norm.data[scale.norm.data$func == variable,]
  ggplot(data = subset)+
    geom_line(mapping = aes(x = res, y = norm_value, color = site, group = scene))+
    labs(x = "Pixel size (meters)", y = "Relative value", title = variable)
  ggsave(str_glue("./sandbox/results/scales/{variable}_scale_normalized.png"), width = 5, height = 3.5, units = "in", limitsize = FALSE)
}

####################
#Summarized per-site
####################

#Summarize by mean, min and max for each site
grouped.data <- in.data %>% group_by(func, site, res) %>% 
  summarise(mean = mean(value, na.rm = TRUE)) %>%
  mutate(res = as.numeric(res))
norm.grouped.data <- grouped.data %>%
  group_by(res, func) %>%
  mutate(
    max = max(mean),
    min = min(mean),
    norm_value = (mean - min) / (max - min)
  ) %>%
  ungroup() 


#Massive plot of all metrics
ggplot(grouped.data, aes(x = res, y = mean, color = site, group=site)) +
  geom_point() +
  geom_line() +
  #geom_errorbar(aes(x = res, ymin = min, ymax = max), width = .2)+
  facet_wrap(~ func, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(grouped.data$res)))+
  theme_bw() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "OAES" = "aquamarine4", "CLBJ" = "aquamarine4", 
                                "ORNL" = "royalblue4", "MLBS" = "royalblue4", 
                                "RMNP"="purple3", "TEAK"="purple3", "WREF"="purple3", 
                                "WOOD"="darkorange3", "OSBS"="darkorange3", "UNDE"="darkorange3"))+
  labs(x = "", y = "", title = "Mean change over resolutions") +
  theme(title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 20))
ggsave("./sandbox/results/scales/MeanChangeResolutions.png", width = 70 , height = 40 , units = "cm", bg = "white")
ggplot(norm.grouped.data, aes(x = res, y = norm_value, color = site, group=site)) +
  geom_point() +
  geom_line() +
  #geom_errorbar(aes(x = res, ymin = min, ymax = max), width = .2)+
  facet_wrap(~ func, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(norm.grouped.data$res)))+
  theme_bw() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "OAES" = "aquamarine4", "CLBJ" = "aquamarine4", 
                                "ORNL" = "royalblue4", "MLBS" = "royalblue4", 
                                "RMNP"="purple3", "TEAK"="purple3", "WREF"="purple3", 
                                "WOOD"="darkorange3", "OSBS"="darkorange3", "UNDE"="darkorange3"))+
  labs(x = "", y = "", title = "Relative mean change over resolutions") +
  theme(title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 20))
ggsave("./sandbox/results/scales/NormMeanChangeResolutions.png", width = 70 , height = 40 , units = "cm", bg = "white")
