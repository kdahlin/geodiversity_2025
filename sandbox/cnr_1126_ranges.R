library(terra)

#Get path to mosaic data
mosaic_dem <- "C:\\Users\\courtney\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\WREF\\WREF_2023_DEM_mosaic_2025-11-23_9_1m.tif"

r1 <- rast(mosaic_dem)
plot(r1)

#dem_extent <- ext(r1)
#print(dem_extent)
#center_x <- (dem_extent$xmin + dem_extent$xmax) / 2
#center_y <- (dem_extent$ymin + dem_extent$ymax) / 2

#center_coords_utm <- c(center_x, center_y)
#print(center_coords_utm)

min_elevation <- min(r1[])
max_elevation <- max(r1[])
print(c(min_elevation, max_elevation))

