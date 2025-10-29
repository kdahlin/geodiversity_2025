library(terra)
library(geodiv)
library(aqp)
#get path to mosaic data

mosaic.dem <- "C:\\Users\\rutgers1\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\ORNL_2018_DEM_mosaic_20250925.tif"
r1 <- rast(mosaic.dem)
plot(r1)

metrics <- terrain(r1)
plot(metrics)

aspect <- terrain(r1, v = "aspect")
plot(aspect)
roughness <- terrain(r1, v = "roughness")
plot(roughness)
units(roughness) <- "m"
plot(roughness)
units(roughness) <- "cm"
plot(roughness)

install.packages("geodiv")
library(terra)
mosaic.dem <- "C:\\Users\\rutgers1\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\WOOD_2020_DEM_mosaic_20251005.tif"
r1 <- rast(mosaic.dem)
plot(r1)
library(geodiv)
# import raster image
data("normforest")
normforest <- terra::unwrap(normforest)

#calculate the fractal dimension
sfd <- sfd(normforest)
sfd(normforest)

# import raster image
data(normforest)
normforest <- terra::unwrap(normforest)

# calculate Scl20, the minimum distance to an autocorrelation value of 0.2 in the AACF
Scl20 <- scl(normforest)[1]

# import raster image
data(normforest)
normforest <- terra::unwrap(normforest)

# determine the 0-5% height interval of the
# bearing area curve
val <- sdc(normforest, 0.0, 0.5)
sdc(normforest, 0, 0.5)

mosaic.dem <- "C:\\Users\\rutgers1\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\WOOD_2020_DEM_mosaic_20251005.tif"
r1 <- rast(mosaic.dem)
data(r1)
r1 <- terra::unwrap(r1)

# import raster image
data(r1)
r1 <- terra::unwrap(r1)
# calculate the fractal dimension
Sfd <- sfd(r1)
sfd(r1)

# import raster image
data(r1)
r1 <- terra::unwrap(r1)
# determine the core roughness depth
Sk <- sk(normforest)
sk(r1)

WOOD <- "C:\\Users\\rutgers1\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\WOOD_2020_DEM_mosaic_20251005.tif"
r1 <- rast(WOOD)
plot(r1)

WOOD_slope <- terrain(r1, v = "slope")
WOOD_aspect <- terrain(r1, v = "aspect")
plot(WOOD_slope)
plot(WOOD_aspect)

#core roughness depth
WOOD_roughnessdepth <- sk(r1)
WOOD_roughnessdepth

#bearing area curve height interval
WOOD_bearing <- sdc(r1, low=0, high=0.05)
WOOD_bearing

#maximum peak height
WOOD_peak<- sph(r1)
WOOD_peak

#texture aspect ratio
WOOD_texture <- stxr(r1, threshold = 0.2)
WOOD_texture

#average roughness
WOOD_roughness <- sa(r1)
WOOD_roughness

WOOD <- "C:\\Users\\rutgers1\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\WOOD_2020_DEM_mosaic_20251005.tif"
r1 <- rast(WOOD)
plot(r1)

# import raster image
r1 <- terra::unwrap(r1)
# find the surface roughness
roughness <- sa(r1)

WOOD_roughness <- sa(r1)
WOOD_roughness
#root mean square roughness
WOOD_rms_roughness <- sq(r1)
WOOD_rms_roughness

#dominant texture direction
WOOD_texturedir <- std(r1, create_plot = FALSE, option = c(1,2))
WOOD_texturedir

WOOD_texturedirtest <- std(r1)

install.packages("devtools")
devtools::install_github("bioXgeo/geodiv")
library(raster)
library(rasterVis)
library(mapdata)
library(maptools)
library(rgeos)
library(ggplot2)
library(tidyverse)
library(parallel)
library(sf)
library(rasterVis)
library(ggmap)
library(corrplot)
library(gridExtra)
library(cowplot)
library(factoextra)
library(cluster)

#NOT geodiv package- trying to measure curvature using spatialEco package
# step 1 install package
install.packages("spatialEco")
library(spatialEco)
library(terra) #already did it above but just to be sure

# step 2 load DEM
elev<-rast("C:\\Users\\rutgers1\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\WOOD_2020_DEM_mosaic_20251005.tif")
# step 3 calculate curvature types

profile_curv<-curvature(elev,type="profile")  #profile curvature

planform_curv<-curvature(elev,type="planform") # plan curvature

total_curv<-curvature(elev,type="total")  #total/sum curvature

# step 4 plot
plot(profile_curv,main="CPER Profile curvature")
plot(planform_curv,main="CPER Plan curvature")
plot(total_curv,main="CPER total curvature")
print(total_curv)
plot(total_curv, main="WOOD Total Curvature")

#get path to mosaic data
wood.dem <- "C:\\Users\\rutgers1\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\WOOD_2020_DEM_mosaic_20251005.tif"
r4 <- rast(wood.dem)
plot(r4)

library(spatialEco)
library(terra) #already did it above but just to be sure
#step 2 load DEM
elev<-rast("C:\\Users\\rutgers1\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\WOOD_2020_DEM_mosaic_20251005.tif")
# step 3 calculate curvature types
profile_curv<-curvature(elev,type="profile")  #profile curvature

planform_curv<-curvature(elev,type="planform") # plan curvature

total_curv<-curvature(elev,type="total")  #total/sum curvature

# step 4 plot

plot(profile_curv,main="WOOD Profile curvature")
plot(planform_curv,main="WOOD Plan curvature")
plot(total_curv,main="WOOD total curvature")

print(total_curv)

library(geodiv)

#average roughness of the surface
roughness <- sa(r4)
roughness

rtmean <- sq(r4)
rtmean

tenpoint <- s10z(r4)
tenpoint

rms2pt <- sdq(r4)
rms2pt

rms7pt <- sdq6(r4)
rms7pt

#surface area ratio
surfarearatio <- sdr(r4)
surfarearatio

#surface bearing index
surfbi <- sbi(r4)
surfbi

#height bearing area curve
heightbac <- sdc(r4, low=0, high=0.05)
heightbac

#correlation length
corrlength <- scl(r4)
corrlength

#reduced valley depth
redvaldepth <- svk(r4)
redvaldepth

redpk <- spk(r4)
redpk

meanpkht <- smean(r4)
meanpkht

coreroughnessdepth <- sk(r4)
coreroughnessdepth

maxpeakh <- sph(r4)
maxpeakh

valleydep

metrics<-terrain(r1)
plot(metrics)
help(terrain)

aspect<-terrain(r1, v="aspect")
plot(aspect)

install.packages("MultiscaleDTM")
install.packages("rasterVis")
install.packages("mapdata")

library(MultiscaleDTM)
library(terra)
library(geodiv)
library(raster)
library(mapdata)
library(ggplot2)
library(cowplot)
library(sf)
library(rasterVis)
library(lattice)
#fits a quadratic surface and can be used to calculate slope, aspect, curvatures, and provide a map of discrete landform classes

#default 3x3 window
qfit1<-Qfit(r4)
plot(qfit1)



#11x11 window
qfit2<-Qfit(r4, w=11, metrics =c('qslope', 'qaspect', 'qeastness', 'qnorthness',
                                 'profc', 'planc', 'twistc', 'meanc', 'features'),
            slope_tolerance = 2, force_center=TRUE, include_scale=TRUE, na.rm=TRUE)
plot(qfit2)
meanc <- qfit2$meanc_11x11
meanc
plot(meanc)

features <- qfit2$features_11x11
features
plot(features)

# Convert grid_sf to terra vector
grid_vect <- vect(grid_sf)

# create 11x11 grid over raster extent
grid_sf <- st_make_grid(
  st_as_sfc(st_bbox(r4)),   # convert raster extent to sf geometry
  n = c(11, 11),             # 3x3 grid
  what = "polygons"
) |> st_as_sf()            # return as sf object
st_crs(grid_sf) <- st_crs(r4) #coord system

# assign IDs
grid_sf$id <- seq_len(nrow(grid_sf))

# calculate metrics (ORNL) for for each grid
for (i in seq_len(nrow(grid_sf))) {
  sub_r <- crop(r4, vect(grid_sf[i, ]))
  mat <- as.matrix(sub_r, wide = TRUE)
  grid_sf$WOOD_curvature[i] <- geodiv::ssc(mat)
}

# Plot for each grid
ggplot() +
  geom_sf(data = grid_sf, aes(fill = WOOD_curvature), color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  labs(
    title = "WOOD Curvature (ssc) per 11×11 Grid Cell",
    fill = "Curvature"
  ) +
  theme_minimal()

# assign IDs
grid_sf$id <- seq_len(nrow(grid_sf))

# Convert grid_sf to terra vector
grid_vect <- vect(grid_sf)

# Extract meanc for each polygon directly
# This computes the mean of all raster cells inside each polygon
vals <- terra::extract(qfit2$meanc_11x11, grid_vect, fun = mean, na.rm = TRUE)

# Attach to grid_sf
grid_sf$meanc <- vals[,2] 

ggplot() +
  geom_sf(data = grid_sf, aes(fill = meanc), color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  labs(
    title = "Local Mean Curvature (aggregated by 11×11 grid)",
    fill = "Mean curvature"
  ) +
  theme_minimal()
