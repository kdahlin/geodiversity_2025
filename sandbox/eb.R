library(terra)

#get path to mosaic data
mosaic_dem<- "C:/Users/Emma/871/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif"

r1<-rast(mosaic_dem)
plot(r1)

metrics<-terrain(r1)
plot(metrics)
help(terrain)

aspect<-terrain(r1, v="aspect")
plot(aspect)
