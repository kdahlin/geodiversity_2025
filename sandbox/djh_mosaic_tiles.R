# -----
# Step 7: Mosaic AOP tiles!
# -----
library(sf)
library(terra)

#Only part that should change between runs
sitename = "CPER"
year = "2020" #Year of data acquisition, not current year


date = toString(Sys.Date())

NEON.path <- paste0("./NEON_data/",sitename,"/")

dem.files <- list.files(NEON.path, pattern = "DTM.tif$", recursive = TRUE)

bli = 33 #Index of the tile that forms the bottom left of the mosaic (bottom-left index)

tile1 <- rast(paste0(NEON.path, dem.files[bli]))
tile2 <- rast(paste0(NEON.path, dem.files[bli+1]))
tile3 <- rast(paste0(NEON.path, dem.files[bli+7]))
tile4 <- rast(paste0(NEON.path, dem.files[bli+8]))

mosaic.dem <- mosaic(tile1, tile2, tile3, tile4)

plot(mosaic.dem) # so pretty!

#Make output directory
dir.create(paste0("./processed_tifs/", sitename))
out.mosaic.name <- paste0("./processed_tifs/", sitename, "/",
                          sitename, "_", year, "_DEM_mosaic_",
                          date, "_9.tif")

writeRaster(mosaic.dem, out.mosaic.name,
            overwrite = FALSE)

# -----
# all done!
# -----