library(terra)

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
WOOD_texturedirtest
