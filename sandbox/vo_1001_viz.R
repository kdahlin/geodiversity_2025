library(terra)

#Get path to mosaic data
mosaic_dem <- "G:\\Y3S1\\GEO 871 Seminar in Physical Geography\\Git_R\\geodiversity_2025-main\\processed_tifs\\ORNL_2018_DEM_mosaic_20250925.tif"

r1 <- rast(mosaic_dem)
plot(r1)

metrics <- terrain(r1)
plot(metrics)

#Aspect calculation and visualization
aspect <- terrain(r1, v = "aspect")
plot(aspect)
