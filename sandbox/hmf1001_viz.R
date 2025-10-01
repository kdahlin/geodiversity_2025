library(terra)
library()
## get path to mosaic:: 

mosaic_dem <- "c:\\Users\\harryfuess\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\ORNL_2018_DEM_mosaic20250925.tif"
r1 <- rast(mosaic_dem)
plot(r1)
metrics <- terrain(r1)
plot(metrics)

aspect <- terrain (r1, v = "aspect")
