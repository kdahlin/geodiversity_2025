######Plotting######################
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(ggrepel)
library(ggbeeswarm)
library(egg)

Metrics <- read_csv("./results/all_metrics_wide.csv") 
View(Metrics)


remove.metrics <- c("sa_1", "sq_1", "s10z_1", "sdq6_1", "sph_1", "svi_1", "scl_1",
                    "TRI_1", "TRIriley_1", "TRIrmsd_1", "roughness_1",
                    "std_1", "TPI_1", "tri_1", "srw_1", "scl_2")

check <- !(Metrics$func %in% remove.metrics)
Metrics.sub <- filter(Metrics, check)

# do some renaming

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

###Plot beeswarm for all variables
Metrics_LongRes_all <- Metrics_LongRes %>%
  mutate(raster_id = factor(raster_id, levels = c("CPER", "OAES", "CLBJ",
                                                  "WOOD", "OSBS","UNDE",
                                                  "ORNL", "MLBS",  "RMNP", 
                                                  "TEAK", "WREF")),
         metrics = factor(metrics,
                          levels = c("range", "smean", "svk", "spk", "mad",
                                     "slope", "northness", "eastness", "sdq",
                                     "sdr", "sar",
                                     "sbi", "sci", "ssk", "sku", "sds", "sv",
                                     "sdc", "nmodes",
                                     "ssc", "Profile Curvature", "Planform Curvature",
                                     "Total Curvature",
                                     "sk", "stdv", "vrm", "tpi",
                                     "sfd", "srw1", "srw2", "srw3", "raster.entropy", 
                                     "std", "stxr1", "stxr2", "scl2")
         ))

x <- ggplot(Metrics_LongRes_all, aes(x = raster_id, y = value, color = raster_id, shape = Tile)) +
  geom_beeswarm(size=1) +
  facet_wrap(~ metrics, scales = "free_y", ncol=4) +  # metrics by type
  theme_bw() +
  scale_color_manual(values = c("CPER" = "#a6cee3", "OAES" = "#a6cee3", "CLBJ" = "#a6cee3", 
                                "ORNL" = "#b2df8a", "MLBS" = "#b2df8a", 
                                "RMNP"="#33a02c", "TEAK"="#33a02c", "WREF"="#33a02c", 
                                "WOOD"="#1f78b4", "OSBS"="#1f78b4", "UNDE"="#1f78b4"), guide="none") +
  scale_shape_manual(name = "Scene",
    values = c("1" = 1, "2" = 2, "3" = 0, "4" = 5, "5" = 16, 
                                "6" = 17, "7" = 18, "8" = 15, "9" = 4))+
  labs(x = "", y = "", title = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 5)
  ) +
  theme(legend.position = "none")

tag_facet(x, tag_pool = c(letters, LETTERS), size = 3)

ggsave("./results/all_metrics_beewarm.png", width = 6.5 , height = 8 , units = "in", 
       dpi = 300)
