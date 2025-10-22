#install packages
library(terra)
library(geodiv)
library(raster)
library(mapdata)
library(ggplot2)
library(cowplot)
library(sf)
library(rasterVis)

setwd("C://Users/rache/Documents/geodiversity_2025/sandbox/")

ORNL <- "C://Users/rache/Documents/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif"
CPER <- "C://Users/rache/Documents/geodiversity_2025/processed_tifs/CPER_2020_DEM_mosaic_20251005.tif"
RMNP <- "C://Users/rache/Documents/geodiversity_2025/processed_tifs/RMNP_2020_DEM_mosaic_20251005.tif"
WOOD <- "C://Users/rache/Documents/geodiversity_2025/processed_tifs/WOOD_2020_DEM_mosaic_20251005.tif"

r1<-rast(ORNL)
plot(r1)

ORNL_slope <- terrain(r1, v="slope")
ORNL_aspect <- terrain(r1, v="aspect")
plot(ORNL_slope)
plot(ORNL_aspect)

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
adjsd1<-AdjSD(r1, include_scale=TRUE)
plot(adjsd1)

######################SpatialEco###########################
library(spatialEco)

SpaEco_curv<-curvature(r1, type="profile")  #profile curvature from SpaEco package
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
plot(SpaEco_tri)

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
