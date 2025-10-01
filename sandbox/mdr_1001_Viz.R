library(terra)

#get path to mosaic data

mosaic.dem <- "C:\\Users\\rutgers1\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\ORNL_2018_DEM_mosaic_20250925.tif"
r1 <- rast(mosaic.dem)
plot(r1)

metrics <- terrain(r1)
plot(metrics)

aspect <- terrain(r1, v = "aspect")
plot(aspect)
