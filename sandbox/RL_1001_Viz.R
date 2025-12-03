#install packages
library(terra)
library(geodiv)
library(mapdata)
library(ggplot2)
library(cowplot)
library(tidyr)
library(sf)
library(rasterVis)

setwd("C://Users/rache/Documents/geodiversity_2025/sandbox/")

ORNL <- "C://Users/rache/Documents/geodiversity_2025/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_2025-11-14_1_1m.tif"
CPER <- "C://Users/rache/Documents/geodiversity_2025/processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-11-14_9_1m.tif"
RMNP <- "C://Users/rache/Documents/geodiversity_2025/processed_tifs/RMNP/RMNP_2020_DEM_mosaic_2025-11-14_2_1m.tif"
WOOD <- "C://Users/rache/Documents/geodiversity_2025/processed_tifs/WOOD_2020_DEM_mosaic_20251005.tif"

r1<-rast(CPER)
plot(r1)
summary(r1)
min <- global(r1, fun = "min", na.rm = TRUE)
min
max <- global(r1, fun = "max", na.rm = TRUE)
max

ORNL_slope <- terrain(r1, v = "slope", neighbors = 8, unit = "degrees")
plot(ORNL_slope)
mean_slope <- global(ORNL_slope, fun = "mean", na.rm = TRUE)
mean_slope

ORNL_aspect <- terrain(r1, v="aspect", v = "aspect", neighbors = 8, unit = "radians")
plot(ORNL_aspect)
mean_aspect <- global(ORNL_aspect, fun = "mean", na.rm = TRUE)
mean_aspect

r2<-rast(RMNP)
plot(r2)

###Calulate global surface gradient metrics
#Average roughness
ORNL_roughness <- sa(r1)
ORNL_roughness

#root mean square roughness
ORNL_rms_roughness <- sq(r1)
ORNL_rms_roughness

#ten-point height
ORNL_10ptht <- s10z(r1)
ORNL_10ptht

#root mean square slope of surface, 2-point method
ORNL_rms_slope_2pt <- sdq(r1)
ORNL_rms_slope_2pt

#root mean square slope of surface, 7-point method
ORNL_rms_slope_7pt <- sdq6(r1)
ORNL_rms_slope_7pt

#surface area ratio
ORNL_area <- sdr(r1)
ORNL_area

#Surface bearing index
ORNL_surfacebearing <- sbi(r1)
ORNL_surfacebearing

#core fluid retention index
ORNL_fluid <- sci(r1)
ORNL_fluid

#skewness
ORNL_skewness <- ssk(r1)
ORNL_skewness

#kurtosis
ORNL_kurtosis <- sku(r1)
ORNL_kurtosis

#summit density
ORNL_summit <- sds(r1)
ORNL_summit

#3d fractal dimension
ORNL_fractal <- sfd(r1)
ORNL_fractal

#dominant radial wavelength, radial wavelength index, mean half wavelength
ORNL_radial <- srw(r1)
ORNL_radial

#angle of dominating texture, texture direction index
ORNL_texturedirection <- std(r1)
ORNL_texturedirection

#valley fluid retention index
ORNL_valleyfluid <- svi(r1)
ORNL_valleyfluid

#texture aspect ratio
ORNL_textureaspect <- stxr(r1)
ORNL_textureaspect

#mean summit curvature
ORNL_curvature <- ssc(r1)
ORNL_curvature

#maximum valley depth
ORNL_valleydepth <- sv(r1)
ORNL_valleydepth

#maximum peak height
ORNL_peak<- sph(r1)
ORNL_peak

#core roughness depth
ORNL_roughnessdepth<- sk(r1)
ORNL_roughnessdepth

#mean peak height
ORNL_peakheight<- smean(r1)
ORNL_peakheight

#reduced valley depth
ORNL_reducedvalley<- svk(r1)
ORNL_reducedvalley

#reduced peak height
ORNL_reducedpeak<- spk(r1)
ORNL_reducedpeak

#correlation length
ORNL_corrlength<- scl(r1)
ORNL_corrlength

#bearing area curve height interval
ORNL_bearing<- sdc(r1, low=0, high=0.05)
ORNL_bearing

###################Test####################

# create 11x11 grid over raster extent
grid_sf <- st_make_grid(
  st_as_sfc(st_bbox(r1)),   # convert raster extent to sf geometry
  n = c(11, 11),             # 11x11 grid
  what = "polygons"
) |> st_as_sf()            # return as sf object
st_crs(grid_sf) <- st_crs(r1) #coord system

# assign IDs
grid_sf$id <- seq_len(nrow(grid_sf))

# calculate metrics (ORNL) for for each grid
for (i in seq_len(nrow(grid_sf))) {
  sub_r <- crop(r1, vect(grid_sf[i, ]))
  mat <- as.matrix(sub_r, wide = TRUE)
  grid_sf$ORNL_curvature[i] <- geodiv::ssc(mat)
}

# Plot for each grid
ggplot() +
  geom_sf(data = grid_sf, aes(fill = ORNL_curvature), color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  labs(
    title = "ORNL Curvature (ssc) per 11×11 Grid Cell",
    fill = "Curvature"
  ) +
  theme_minimal()

grid_sf$ORNL_curvature

###################MultiscaleDTM####################
#https://cran.r-project.org/web/packages/MultiscaleDTM/MultiscaleDTM.pdf
#Calculates multi-scale geomorphometric terrain attributes from regularly gridded digital terrain models using a variable focal windows size
library(MultiscaleDTM)

#fits a quadratic surface and can be used to calculate slope, aspect, curvatures, and provide a map of discrete landform classes
#default 3x3 window
qfit1<-Qfit(r1, w=3, metrics =c('profc'), slope_tolerance = 2, force_center=TRUE, include_scale=TRUE, na.rm=TRUE)
profc <- qfit1$profc_3x3
summary(profc)
plot(profc)

grid_vect <- vect(grid_sf) # Convert grid_sf to terra vector

# Extract profc for each polygon directly
# This computes the mean of all raster cells inside each polygon
vals_dtm <- terra::extract(qfit1$profc_3x3, grid_vect, fun = mean, na.rm = TRUE)

# Attach to grid_sf
grid_sf$profc.dtm <- vals_dtm[,2] 

ggplot() +
  geom_sf(data = grid_sf, aes(fill = profc.dtm), color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", limits = c(-0.005, 0.005)) +
  labs(
    title = "Profile Curvature using MultiScaleDTM Package",
    fill = "Profile curvature"
  ) +
  theme_minimal()

#Calculating roughness via adjusted standard deviation
DTM_adjsd1<-AdjSD(r1, include_scale=TRUE)
mean_adjsd1 <- global(DTM_adjsd1, fun = "mean", na.rm = TRUE)
mean_adjsd1
plot(DTM_adjsd1)

#Roughness Index-Elevation
DTM_RoughnessIndexEle<-RIE(r1, include_scale=TRUE)
mean_rie<- global(DTM_RoughnessIndexEle, fun = "mean", na.rm = TRUE)
mean_rie
plot(DTM_RoughnessIndexEle)

######################SpatialEco###########################
library(spatialEco)

SpaEco_curv<-curvature(r1, type=c("profile", "planform", "total"))  #profile curvature from SpaEco package
summary(values(SpaEco_curv))

#SpaEco_curv_clip <- clamp(SpaEco_curv, lower=-0.01, upper=0.01) #set lower and upper limit/ remove outliers
#plot(SpaEco_curv_clip, zlim=c(-0.01, 0.01))
# Extract profc for each polygon directly
vals_spaeco <- terra::extract(SpaEco_curv, grid_vect, fun = mean, na.rm = TRUE)

# Attach to grid_sf
grid_sf$profc.spaeco <- vals_spaeco[,2] 

#plot
ggplot() +
  geom_sf(data = grid_sf, aes(fill = profc.spaeco), color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", limits = c(-0.005, 0.005)) +
  labs(
    title = "Profile Curvature using SpatialEco Package",
    fill = "Profile curvature"
  ) +
  theme_minimal()

#topographic position
SpaEco_tpi <- tpi(r1, scale = 3, win = "rectangle", normalize = FALSE, zero.correct = FALSE)
plot(SpaEco_tpi)
mean_tpi <- global(SpaEco_tpi, fun = "mean", na.rm = TRUE)
mean_tpi

#terrain ruggedness
SpaEco_tri <- tri(r1, s = 3, exact = TRUE)
mean_tri <- global(SpaEco_tri, fun = "mean", na.rm = TRUE)
mean_tri
plot(SpaEco_tri)

#Vector ruggedness measure
SpaEco_vrm<-VRM(r1, include_scale=TRUE)
mean_vrm<- global(SpaEco_vrm, fun = "mean", na.rm = TRUE)
mean_vrm
plot(SpaEco_vrm)

#surface area ratio
SpaEco_sar <- sar(r1, s = NULL, scale = TRUE)
mean_sar <- global(SpaEco_sar, fun = "mean", na.rm = TRUE)
mean_sar
plot(SpaEco_sar)

#Raster entropy
rEnt <- raster.entropy(r1, d=3, categorical = FALSE, global = TRUE)
median_rEnt <- global(rEnt, fun=median, na.rm=TRUE)
median_rEnt
plot(rEnt, limits=c(2.1962, 2.1973))
rEnt

#Calculate focal statistics for raster
library(moments)
moments_skew <- raster.moments(r1, type = "skewness", s = 3)
moments_skew
plot(moments_skew)
mean_skew <- global(moments_skew, fun = "mean", na.rm = TRUE)
mean_skew

moments_kurtosis <- raster.moments(r1, type = "kurtosis", s = 3)
moments_kurtosis
plot(moments_kurtosis)
mean_kurtosis <- global(moments_skew, fun = "mean", na.rm = TRUE)
mean_kurtosis

###################Others####################

#texture direction metrics
texturedirection <- std(r1, create_plot = FALSE, option=1)
texturedirection

###Texture image creation
window <- matrix(1, nrow=7, ncol=7)
system.time(
  output_raster<-focal_metrics(r1, window, metrics=list('sa'), progress=TRUE))
print(output_raster)
rasterVis::levelplot(output_raster[[1]], margin=F, par.settings=eviTheme, 
          ylab=NULL, xlab=NULL, main='Sa')

######Plotting######################
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(ggrepel)
library(ggbeeswarm)
Metrics <- read_csv("C:/Users/rache/Documents/geodiversity_2025/results/all_metrics_wide.csv") 
View(Metrics)

#Pivot data to format for plotting
Metrics_Long<- Metrics %>%
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

#Boxplot for all
ggplot(Metrics_LongRes, aes(x = raster_id, y = value, color = raster_id)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 1) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.8) +
  facet_wrap(~ metrics, scales = "free_y") +  # one panel per metric
  scale_color_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3"))+
  labs(x = "Raster ID",
       y = "Value",
       title = "Metrics Across 4 different surfaces") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey80"))
ggsave("BoxplotAll.png", width = 70 , height = 40 , units = "cm", bg = "white")

#Identify Outliers
findoutlier <- function(x) {    #Define function to identify outlier
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
  }
Metrics_LongRes <- Metrics_LongRes %>%    #Include new columns that shows if observation is an outlier
  group_by(metrics, raster_id) %>%
  mutate(outlier = ifelse(findoutlier(value), value, NA))

#filter data per metrics
type1 <- Metrics_LongRes %>%
  filter(metrics %in% c("s10z", "sdq6", "smean", "sph", "sv", "svk"))
type2 <- Metrics_LongRes %>%
  filter(metrics %in% c("mad", "range","roughness", 
                        "sdc", "sdq", "sk", "slope", "spk", "stdv",
                        "tri", "TRI", "TRIriley", "TRIrmsd", "vrm", "sa", "sq"))
type3 <- Metrics_LongRes %>%
  filter(metrics %in% c("scl", "sku", "srw", "raster.entropy", "stxr"))
type4 <- Metrics_LongRes %>%
  filter(metrics %in% c("sfd", "sci", "svi","tpi", "TPI"))
type5 <- Metrics_LongRes %>%
  filter(metrics %in% c("sbi", "ssk"))
type6 <- Metrics_LongRes %>%
  filter(metrics %in% c("northness", "sdr","curvature", "sds"))
type7 <- Metrics_LongRes %>%
  filter(!metrics %in% c("s10z", "sdq6", "smean", "sph", "sv", "svk",
                         "mad", "range", "roughness", "sdc", "sdq", "sk", "slope", "spk", "stdv",
                         "tri", "TRI", "TRIriley", "TRIrmsd", "vrm", "sa", "sq",
                         "scl", "sku", "srw", "raster.entropy", "stxr", 
                         "sfd", "sci", "svi","tpi", "TPI",
                         "sbi", "ssk",
                         "northness", "sdr", "curvature", "sds"))

#Boxplots by types
ggplot(type6, aes(x = raster_id, y = value, color = raster_id)) +
  geom_boxplot() +
  geom_text_repel(data = subset(type6, !is.na(outlier)),
            aes(label = index),
            hjust = -0.8) +
  facet_wrap(~ metrics, scales = "free_y") +  # one panel per metric
  theme_minimal() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3"))+
  labs(x = "", y = "", title = "") +
  theme(legend.position="none",
        axis.text = element_text(size=6))

ggsave("Boxplot6.png", width = 13.5 , height = 10 , units = "cm")
ggsave("Boxplot5.png", width = 13.5 , height = 6 , units = "cm")

#Violin plot for all
ggplot(Metrics_Long, aes(x = raster_id, y = value, fill = raster_id)) +
  geom_violin(trim=FALSE, alpha=0.5) +
  geom_boxplot(width=0.1) +
  facet_wrap(~ metrics, scales = "free_y") +  # one panel per metric
  scale_fill_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3"))+
  labs(x = "Raster ID",
       y = "Value",
       title = "Metrics Across 4 different surfaces") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey80"))

ggsave("ViolinAll.png", width = 70 , height = 40 , units = "cm", bg = "white")

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

###Plot beeswarm for top 6 variables from PCA
Metrics_LongRes_top <- Metrics_LongRes %>%
  filter(metrics %in% c("range", "slope", "sku", "stdv", "curvature.profile", "sfd" ,"stxr")) %>%
  mutate(raster_id = factor(raster_id, levels = c("CPER", "OAES", "CLBJ",
                                                  "WOOD", "OSBS","UNDE",
                                                  "ORNL", "MLBS", 
                                                  "RMNP", "TEAK", "WREF")),
         metrics = factor(metrics,
                     levels = c("range", "slope", "sku", "stdv", "curvature.profile", "sfd" ,"stxr"),
                     labels = c("Range", "Slope", "Surface kurtosis", "Standard deviation", 
                                "Profile curvature", "Fractal dimension", "Texture aspect ratio")))
ggplot(Metrics_LongRes_top, aes(x = raster_id, y = value, color = raster_id, shape = Tile)) +
  geom_beeswarm(size=3) +
  facet_wrap(~ metrics, scales = "free_y", nrow=2) +  # metrics by type
  theme_bw() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "OAES" = "aquamarine4", "CLBJ" = "aquamarine4", 
                                "ORNL" = "royalblue4", "MLBS" = "royalblue4", 
                                "RMNP"="purple3", "TEAK"="purple3", "WREF"="purple3", 
                                "WOOD"="darkorange3", "OSBS"="darkorange3", "UNDE"="darkorange3"), guide="none") +
  scale_shape_manual(name = "Tile",
    values = c("1" = 1, "2" = 2, "3" = 0, "4" = 5, "5" = 16, "6" = 17, "7" = 18, "8" = 15, "9" = 4))+
  labs(x = "", y = "", title = "") +
  theme(legend.position = c(0.8, 0.22),
        legend.text = element_text(size=25),
        legend.title = element_text(size=30),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 25),
        strip.text = element_text(size = 45),
        strip.background = element_rect(fill = "honeydew2"))
ggsave("Top6.png", width = 45 , height = 15 , units = "cm")

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
