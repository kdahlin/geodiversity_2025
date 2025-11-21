#Script to plot metrics against increasing scale
library(ggplot2)
library(openxlsx2)
library(stringr) #Allows python-style string formatting
library(dplyr)

#Set working directory
setwd("C:/Users/hallerdi/Documents/Dillon_Classes/geodiversity_2025")

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

#Iterate over variables
for (variable in unique_variables){
  #Subset data to only that currently of interest
  subset <- norm.data[norm.data$func==variable,]
  ggplot(data = subset)+
    geom_line(mapping = aes(x = res, y = norm_value*100, color = site, group = scene))+
    labs(x = "Pixel size (meters)", y = "% change relative to 1-meter pixel size", title = variable)
  ggsave(str_glue("./sandbox/results/scales/{variable}.png"), width = 5, height = 3.5, units = "in", limitsize = FALSE)
}
