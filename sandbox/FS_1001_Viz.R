library(terra)
#Get the path to MOSAIC dataset 
mosaic_dem <- "F:/Semester Third/GEO 871/GitHub/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif"
r1 <- rast(mosaic_dem)
plot(r1)

metrics <- terrain(r1)
plot(metrics)

aspect <- terrain(r1)
plot(metrics)



