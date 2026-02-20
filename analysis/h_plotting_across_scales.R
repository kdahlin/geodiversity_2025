#Script to plot metrics against increasing scale
library(ggplot2)
library(openxlsx2)
library(stringr) #Allows python-style string formatting
library(dplyr)

# read in the long data
in.data <- read.csv("./results/all_metrics_long.csv")

# rename functions so they are pretty

old.names <- unique(in.data$func)

names <- c("Sa", "Sq", "s10z", "Sdq", "Sdq6", "Sdr", "SBI", "SCI", "Ssk",
  "Sku", "Sds", "Sfd", "Srw1", "Srw2", "Srw3", "Std1", "Std2", 
  "Svi", "Stxr1", "Stxr2", "Ssc", "Sv", "Sph", "Sk", "Smean", 
  "Spk", "Svk", "Scl1", "Scl2", "Sdc","TPIv1", "TRIv1", "VRM", 
  "SAR", "rEnt", "Slope", "nor", "east", "TPIv2", "TRIv2", 
  "TRIriley", "TRIrmsd", "rough", "stdv", "MAD", "nmodes", 
  "range", "CuPro", "CuPl", "CuTot")

for (i in 1:length(names)) {
  in.data$func[in.data$func == old.names[i]] <- names[i]
}

remove.metrics <- c("Sa", "Sq", "s10z", "Sdq6", "Sph", "Svi", "Scl1",
                    "TRIv1", "TRIriley", "TRIrmsd", "rough",
                    "Std1", "TPIv1", "TPIv1", "TRIv2", "srw1", "Scl2")

check <- !(in.data$func %in% remove.metrics)
Metrics.sub <- filter(in.data, check)


grouped.data <- Metrics.sub %>% group_by(func, site, res) %>% 
  summarise(mean = mean(value, na.rm = TRUE)) %>%
  mutate(res = as.numeric(res))

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

i <- 48
subset <- in.data[in.data$func==unique_variables[i],]
subset %>%
  group_by(res, site) %>%
  mutate(mean.value = mean(value),
         plus.sd = mean(mean.value) + sd(value),
         minus.sd = mean(mean.value) - sd(value)) %>%
  ggplot(., aes(group = site, fill = site)) +
    geom_ribbon(mapping = aes(x = res, ymin = minus.sd, 
                            ymax = plus.sd), alpha = 0.2) +
  geom_line(aes(x = res, y = mean.value, color = site), 
            linewidth = 2) +
  labs(x = "Pixel size (meters)", y = "Actual Value", 
       title = unique_variables[i]) +
  theme_bw()+
  facet_grid(site~.)
  

ggplot(data = subset) +
  geom_line(mapping = aes(x = res, y = value, color = site, group = scene)) +
  labs(x = "Pixel size (meters)", y = "Actual Value", 
       title = unique_variables[i])


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
