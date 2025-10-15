install.packages("cowplot")
library(terra)
library(sf)
library(geodiv)
library(ggplot2)
library(cowplot)

mosaic.dem <- "C:\\Users\\courtney\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\RMNP_2020_DEM_mosaic_20251005.tif"
r1 <- rast(mosaic.dem)
plot(r1)

metrics <- terrain(r1)
plot(metrics)

aspect <- terrain(r1, v = "aspect")
plot(aspect)
