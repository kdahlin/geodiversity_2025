######Plotting######################
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(ggrepel)
library(ggbeeswarm)

Metrics <- read_csv("./results/all_metrics_wide.csv") 
View(Metrics)


remove.metrics <- c("sa_1", "sq_1", "s10z_1", "sdq6_1", "sph_1", "svi_1", "scl_1",
                    "TRI_1", "TRIriley_1", "TRIrmsd_1", "roughness_1",
                    "std_1", "TPI_1", "tri_1")

check <- !(Metrics$func %in% remove.metrics)
Metrics.sub <- filter(Metrics, check)

# do some renaming
i <- which(Metrics.sub$func == "srw_1")
Metrics.sub$func[i] <- "srw1"

i <- which(Metrics.sub$func == "srw_2")
Metrics.sub$func[i] <- "srw2"

i <- which(Metrics.sub$func == "srw_3")
Metrics.sub$func[i] <- "srw3"

i <- which(Metrics.sub$func == "curvature_profile_1")
Metrics.sub$func[i] <- "Profile Curvature"

i <- which(Metrics.sub$func == "curvature_planform_1")
Metrics.sub$func[i] <- "Planform Curvature"

i <- which(Metrics.sub$func == "curvature_total_1")
Metrics.sub$func[i] <- "Total Curvature"

i <- which(Metrics.sub$func == "stxr_1")
Metrics.sub$func[i] <- "stxr1"

i <- which(Metrics.sub$func == "stxr_2")
Metrics.sub$func[i] <- "stxr2"

i <- which(Metrics.sub$func == "scl_2")
Metrics.sub$func[i] <- "scl2"

#Pivot data to format for plotting
Metrics_Long<- Metrics.sub %>%
  mutate(func = gsub("_.*", "", func)) %>%  # remove anything after the first "_" from func column
  pivot_longer(
    cols = -func,               # except first row
    names_to = "Raster",        # new column for raster
    values_to = "value"         # new column for the metric value
  ) %>%
  rename(metrics = func) %>%       # rename func to metrics
  select(metrics, value, Raster) %>%
  mutate(raster_id = substr(Raster, 1, 4),   # add new column with locations, ID first 4 letters in raster
         resolution = sub(".*_(\\d+)m\\.tif$", "\\1", Raster),     # extract number before "m.tif" 
         Tile = sub(".*_(\\d+)_\\d+m\\.tif$", "\\1", Raster)) %>%  # second-to-last number   
  drop_na()  #remove NA values

#filter resolutions
Metrics_LongRes <- Metrics_Long %>% filter(resolution==1)

#Identify Outliers
findoutlier <- function(x) {    #Define function to identify outlier
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
  }
Metrics_LongRes <- Metrics_LongRes %>%    #Include new columns that shows if observation is an outlier
  group_by(metrics, raster_id) %>%
  mutate(outlier = ifelse(findoutlier(value), value, NA))

#Beeswarm plot for all
ggplot(Metrics_LongRes, aes(x = raster_id, y = value, color = raster_id)) +
  geom_beeswarm() +
  facet_wrap(~ metrics, scales = "free_y") +  # one panel per metric
  scale_color_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3"))+
  labs(x = "Raster ID",
       y = "Value",
       title = "Metrics Across 4 different surfaces") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey80"))

ggsave("BeeswarmAll.png", width = 70 , height = 40 , units = "cm", bg = "white")

#Plot subset for individual types
#Scale metrics so that can be plotted on the same 0-1 scale
Metrics_LongResScaled <- Metrics_LongRes %>%
  group_by(metrics) %>%
  mutate(value_scaled = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

#Rerun individual types
#filter data per metrics
type1 <- Metrics_LongResScaled %>%
  filter(metrics %in% c("s10z", "sdq6", "smean", "sph", "sv", "svk")) %>% mutate(type = "Type 1")
type2 <- Metrics_LongResScaled %>%
  filter(metrics %in% c("mad", "range","roughness", 
                        "sdc", "sdq", "sk", "slope", "spk", "stdv",
                        "tri", "TRI", "TRIriley", "TRIrmsd", "vrm", "sa", "sq")) %>% mutate(type = "Type 2")
type3 <- Metrics_LongResScaled %>%
  filter(metrics %in% c("scl", "sku", "srw", "raster.entropy", "stxr")) %>% mutate(type = "Type 3")
type4 <- Metrics_LongResScaled %>%
  filter(metrics %in% c("sfd", "sci", "svi","tpi", "TPI")) %>% mutate(type = "Type 4")
type5 <- Metrics_LongResScaled %>%
  filter(metrics %in% c("sbi", "ssk")) %>% mutate(type = "Type 5")
type6 <- Metrics_LongResScaled %>%
  filter(metrics %in% c("northness", "sdr","curvature", "sds")) %>% mutate(type = "Type 6")
type7 <- Metrics_LongResScaled %>%
  filter(!metrics %in% c("s10z", "sdq6", "smean", "sph", "sv", "svk",
                         "AdjSD", "RIE", "roughness", "sdc", "sdq", "sk", "slope", "spk", 
                         "tri", "TRI", "TRIriley", "TRIrmsd", "vrm", "sa", "sq",
                         "scl", "sku", "srw", "raster.entropy", "stxr", 
                         "sfd", "sci", "svi","tpi", "TPI",
                         "sbi", "ssk",
                         "northness", "sdr", "curvature", "sds"))
All <- bind_rows(type1, type2, type3, type4, type5, type6)

###Plot individual types
#Boxplot
ggplot(type6, aes(x = raster_id, y = value_scaled, color = raster_id)) +
  geom_boxplot() +
  theme_minimal() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3"))+
  labs(x = "", y = "", title = "Type 6") +
  theme(legend.position="none", 
        axis.text = element_text(size=6))
ggsave("Boxplot_type.png", width = 6 , height = 6 , units = "cm")

#Boxplots (with individual metrics)
ggplot(type1, aes(x = raster_id, y = value_scaled, color = raster_id, fill=metrics)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_brewer(palette="Set3") +
  scale_color_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3"))+
  labs(x = "", y = "", title = "Metrics Across 4 different surfaces") +
  theme(axis.text = element_text(size=6))

#Beeswarm
ggplot(type1, aes(x = raster_id, y = value_scaled, color = raster_id, shape = metrics)) +
  geom_beeswarm() +
  theme_minimal() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3"))+
  labs(
    x = "Raster ID",
    y = "Value",
    title = "Metrics Across 4 different surfaces")

###Plot All
ggplot(All, aes(x = raster_id, y = value_scaled, color = raster_id)) +
  geom_beeswarm() +
  facet_wrap(~ type, scales = "free_y", nrow=2) +  # metrics by type
  theme_bw() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3")) +
  labs(x = "", y = "", title = "") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "honeydew2"))

ggsave("Types.png", width = 30 , height = 20 , units = "cm")



###Plot beeswarm for all variables
Metrics_LongRes_all <- Metrics_LongRes %>%
  mutate(raster_id = factor(raster_id, levels = c("CPER", "OAES", "CLBJ",
                                                  "WOOD", "OSBS","UNDE",
                                                  "ORNL", "MLBS",  "RMNP", 
                                                  "TEAK", "WREF")),
         metrics = factor(metrics,
                          levels = c("range", "smean", "svk", "spk", "mad",
                                     "slope", "northness", "eastnesss", "sdq",
                                     "sdr", "sar",
                                     "sbi", "sci", "ssk", "sku", "sds", "svi",
                                     "sdc", "nmodes",
                                     "ssc", "Profile Curvature", "Planform Curvature",
                                     "Total Curvature",
                                     "sk", "stdv", "vrm", "tpi",
                                     "sfd", "srw1", "srw2", "srw3", "raster.entropy", 
                                     "std", "stxr1", "stxr2", "scl2")
         ))

ggplot(Metrics_LongRes_all, aes(x = raster_id, y = value, color = raster_id, shape = Tile)) +
  geom_beeswarm(size=3) +
  facet_wrap(~ metrics, scales = "free_y", nrow=6) +  # metrics by type
  theme_bw() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "OAES" = "aquamarine4", "CLBJ" = "aquamarine4", 
                                "ORNL" = "royalblue4", "MLBS" = "royalblue4", 
                                "RMNP"="purple3", "TEAK"="purple3", "WREF"="purple3", 
                                "WOOD"="darkorange3", "OSBS"="darkorange3", "UNDE"="darkorange3"), guide="none") +
  scale_shape_manual(name = "Tile",
                     values = c("1" = 1, "2" = 2, "3" = 0, "4" = 5, "5" = 16, "6" = 17, "7" = 18, "8" = 15, "9" = 4))+
  labs(x = "", y = "", title = "")


+
  theme(legend.position = "right",
        legend.text = element_text(size=25),
        legend.title = element_text(size=30),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 25),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill = "honeydew2"))




###Plot beeswarm for top 12 variables from PCA
Metrics_LongRes_top <- Metrics_LongRes %>%
  filter(metrics %in% c("range", "curvature.total", "sbi", "stxr", "curvature.profile", "svk",
                        "sku", "eastness", "northness", "nmodes", "raster.entropy", "srw")) %>%
  mutate(raster_id = factor(raster_id, levels = c("CPER", "OAES", "CLBJ",
                                                  "WOOD", "OSBS","UNDE",
                                                  "ORNL", "MLBS",  
                                                  "RMNP", "TEAK", "WREF")),
         metrics = factor(metrics,
                          levels = c("range", "curvature.total", "sbi", "stxr", "curvature.profile", "svk",
                                     "sku", "eastness", "northness", "nmodes", "raster.entropy", "srw"),
                          labels = c("Range", "Total curvature", "Surface bearing index", "Texture aspect ratio", 
                                     "Profile curvature", "Reduced valley depth", "Surface kurtosis",
                                     "Aspect (easting)", "Aspect (northing)",
                                     "No. of modes", "Raster entropy", "Dominant radial wavelength")))
ggplot(Metrics_LongRes_top, aes(x = raster_id, y = value, color = raster_id, shape = Tile)) +
  geom_beeswarm(size=3) +
  facet_wrap(~ metrics, scales = "free_y", nrow=3) +  # metrics by type
  theme_bw() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "OAES" = "aquamarine4", "CLBJ" = "aquamarine4", 
                                "ORNL" = "royalblue4", "MLBS" = "royalblue4", 
                                "RMNP"="purple3", "TEAK"="purple3", "WREF"="purple3", 
                                "WOOD"="darkorange3", "OSBS"="darkorange3", "UNDE"="darkorange3"), guide="none") +
  scale_shape_manual(name = "Tile",
                     values = c("1" = 1, "2" = 2, "3" = 0, "4" = 5, "5" = 16, "6" = 17, "7" = 18, "8" = 15, "9" = 4))+
  labs(x = "", y = "", title = "") +
  theme(legend.position = "right",
        legend.text = element_text(size=25),
        legend.title = element_text(size=30),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 25),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill = "honeydew2"))
ggsave("Top12.png", width = 50 , height = 30 , units = "cm")

###Data Exploration###
Metrics_LongRes_test <- Metrics_LongRes %>% filter(metrics %in% c("range")) %>% filter(raster_id != "RMNP")
ggplot(Metrics_LongRes_test, aes(x = raster_id, y = value, color = raster_id)) + geom_beeswarm(size=3)























################### Resolutions ####################
#Grouped and calculate mean for 9 tiles
Metrics_Long_met <- Metrics_Long %>% group_by(metrics, raster_id, resolution) %>% 
  summarise(mean = mean(value, na.rm = TRUE)) %>%
  mutate(resolution = as.numeric(resolution))
Metrics_Long_sa <- Metrics_Long_met %>% filter(metrics=="range")

#Mean change over resolutions
ggplot(Metrics_Long_sa, aes(x = resolution, y = mean, color = raster_id, group=raster_id)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = sort(unique(Metrics_Long_sa$resolution)))+
  scale_color_manual(values = c("CPER" = "aquamarine4", "OAES" = "aquamarine4", "CLBJ" = "aquamarine4", 
                                "ORNL" = "royalblue4", "MLBS" = "royalblue4", 
                                "RMNP"="purple3", "TEAK"="purple3", "WREF"="purple3", 
                                "WOOD"="darkorange3", "OSBS"="darkorange3", "UNDE"="darkorange3"))

ggplot(Metrics_Long_met, aes(x = resolution, y = mean, color = raster_id, group=raster_id)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ metrics, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(Metrics_Long_sa$resolution)))+
  theme_bw() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "OAES" = "aquamarine4", "CLBJ" = "aquamarine4", 
                                "ORNL" = "royalblue4", "MLBS" = "royalblue4", 
                                "RMNP"="purple3", "TEAK"="purple3", "WREF"="purple3", 
                                "WOOD"="darkorange3", "OSBS"="darkorange3", "UNDE"="darkorange3"))+
  labs(x = "", y = "", title = "Mean change over resolutions") +
  theme(title = element_text(size=60),
        legend.text = element_text(size=35),
        legend.title = element_text(size=40),
        axis.text = element_text(size = 35),
        strip.text = element_text(size = 40))
ggsave("MeanChangeResolutions.png", width = 70 , height = 40 , units = "cm", bg = "white")

#Percent change relative to 1 m resolution
Metrics_Long_Percent <- Metrics_Long %>%
  group_by(metrics, raster_id, resolution) %>%
  summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(resolution = as.numeric(resolution)) %>%
  group_by(metrics, raster_id) %>%
  mutate(
    mean_1m = mean[resolution == 1],
    pct_change = ((mean - mean_1m) / mean_1m) * 100
  ) %>%
  ungroup()

ggplot(Metrics_Long_Percent, aes(x = resolution, y = pct_change, color = raster_id, group=raster_id)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ metrics, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(Metrics_Long_sa$resolution)))+
  theme_bw() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "OAES" = "aquamarine4", "CLBJ" = "aquamarine4", 
                                "ORNL" = "royalblue4", "MLBS" = "royalblue4", 
                                "RMNP"="purple3", "TEAK"="purple3", "WREF"="purple3", 
                                "WOOD"="darkorange3", "OSBS"="darkorange3", "UNDE"="darkorange3"))+
  labs(x = "", y = "", title = "% mean change relative to 1m resolution") +
  theme(title = element_text(size=60),
        legend.text = element_text(size=35),
        legend.title = element_text(size=40),
        axis.text = element_text(size = 35),
        strip.text = element_text(size = 40))

ggsave("MeanChangeResolutions_Percent.png", width = 70 , height = 40 , units = "cm", bg = "white")
