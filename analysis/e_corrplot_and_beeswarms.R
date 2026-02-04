######Plotting######################
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(ggrepel)
library(ggbeeswarm)
library(egg)

Metrics <- read.csv("./results/all_metrics_wide.csv") 
View(Metrics)

# rename to match manuscript & store
Metrics$func <- c("Sa", "Sq", "s10z", "Sdq", "Sdq6", "Sdr", "SBI", "SCI", "Ssk",
                  "Sku", "Sds", "Sfd", "Srw1", "Srw2", "Srw3", "Std1", "Std2", 
                  "Svi", "Stxr1", "Stxr2", "Ssc", "Sv", "Sph", "Sk", "Smean", 
                  "Spk", "Svk", "Scl1", "Scl2", "Sdc","TPIv1", "TRIv1", "VRM", 
                  "SAR", "rEnt", "Slope", "nor", "east", "TPIv2", "TRIv2", 
                  "TRIriley", "TRIrmsd", "rough", "stdv", "MAD", "nmodes", 
                  "range", "CuPro", "CuPl", "CuTot")

write.csv(Metrics[,-1], "./results/all_metrics_wide.csv", row.names = Metrics$func)

################################################################################
### data exploration at 1 m resolution 
in.data <- read.csv("./results/all_metrics_wide.csv", row.names = 1) 

# transpose data so each row is a subsite and each column is a metric
in.data.t <- as.data.frame(t(in.data))

# get row names (now site x scene x resolution)
samples <- rownames(in.data.t)

# split out resolutions
sample.res <- strsplit(samples, "[_\\.]")
res <- unlist(lapply(sample.res, function(x) {x[9]}))

# rename rows so they aren't so clunky
rnames <- row.names(in.data.t)
split.rnames <- strsplit(rnames, "_")

rnames.new <- paste0(split.rnames[[1]][1], "_", split.rnames[[1]][6],
                     "_", res[1])

for (i in 2:length(rnames)) {
  rnames.new <- rbind(rnames.new, 
                      paste0(split.rnames[[i]][1], "_", split.rnames[[i]][6],
                             "_", res[i]))
}

# attach the better names
row.names(in.data.t) <- rnames.new[,1]

# get just 1m resolution
in.data.1m <- subset(in.data.t, res == "1m")

#remove the res from the row names
rownames(in.data.1m) <- gsub("_1m", "", rownames(in.data.1m))

# look at column variances (to remove columns with no variance)
variances <- sapply(in.data.1m, var)

# get the column names with no variance
which.var <- which(variances == 0)

# remove columns with no variance
in.data.1m <- in.data.1m[,-(which.var)]

# look at correlations
cors <- cor(in.data.1m)

#write.csv(cors, "./results/correlation_matrix_1m.csv", row.names = TRUE)

png(filename = "./results/correlation_matrix.png", 
    width = 6, height = 6, units = "in", res = 300, bg = "white")
# visualize with corrplot
corrplot(cors, order = "FPC", method = "circle", tl.cex = 0.6)
dev.off()

########## Remove highly correlated ones ###########################

remove.metrics <- c("Sa", "Sq", "s10z", "Sdq6", "Sph", "Svi", "Scl1",
                    "TRIv1", "TRIriley", "TRIrmsd", "rough",
                    "Std1", "TPIv1", "TPIv1", "TRIv2", "srw_1", "Scl2")

check <- !(Metrics$func %in% remove.metrics)
Metrics.sub <- filter(Metrics, check)

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
                                     "sfd", "srw2", "srw3", "raster.entropy", 
                                     "std", "stxr1", "stxr2")
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
