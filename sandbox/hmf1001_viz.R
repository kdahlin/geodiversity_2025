library(terra)
library(raster)
library(tidyverse)
library(geodiv)

## get path to mosaic:: 
#note to self, use below code for creating
mosaic_dem <-"processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif"
r1 <- rast(mosaic_dem)
plot(r1)
metrics <- terrain(r1)
plot(metrics)

aspect <- terrain (r1, v = "aspect")
aspect

#Oak Ridge (stripey) 
# Prairie Potholes (lumpy) 
# Central Plains (flat) 
# Rocky Mountains (high relief) 
# what are 




