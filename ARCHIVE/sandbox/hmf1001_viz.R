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


library(terra)
library(sf)
library(dplyr)
library(readr)
library(stringr)


uno <-("processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-10-22_1.tif") 
dos <- "processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-10-22_2.tif"
tres <- "processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-10-22_3.tif"
quatro <- "processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-10-22_4.tif"
cinco <- "processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-10-22_5.tif"
seis<- "processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-10-22_6.tif"
siete <- "processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-10-22_7.tif"
ocho <-"processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-10-22_8.tif"
nueve <-"processed_tifs/CPER/CPER_2020_DEM_mosaic_2025-10-22_9.tif"
rast1 <- rast(uno)
rast2 <- rast(dos) 
rast3 <- rast(tres)
rast4 <- rast(quatro)
rast5 <- rast(cinco)
rast6 <- rast(seis)
rast7 <- rast(siete)
rast8 <- rast(ocho)
rast9 <- rast(nueve)


plot(rast1)
plot(rast2)
plot(rast3)
plot(rast4)
plot(rast5)
plot(rast6)
plot(rast7)
plot(rast8)
plot(rast9)

#merging all of the rasters for CPER
x <- list(rast1, rast2, rast3, rast4, rast5, rast6, rast7, rast8, rast9)
m <- do.call(merge, x)
plot(m)
plot(rast5)
m

utm_point <- data.frame(
  x = 521000,   # Easting (m)
  y = 4518000   # Northing (m)
)

utm_sf <- st_as_sf(utm_point, coords = c("x", "y"), crs = 32613)



                    
                    
                    
                    
         