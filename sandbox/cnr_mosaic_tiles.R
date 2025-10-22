# -----
# Step 7: Mosaic AOP tiles!
# -----

NEON.path <- paste0("./NEON_data/ORNL/DP3.30024.001/neon-aop-products/2018/",
                    "FullSite/D07/2018_ORNL_4/L3/DiscreteLidar/DTMGtif/")

dem.files <- list.files(NEON.path)

tile1 <- rast(paste0(NEON.path, dem.files[1]))
tile2 <- rast(paste0(NEON.path, dem.files[2]))
tile3 <- rast(paste0(NEON.path, dem.files[3]))
tile4 <- rast(paste0(NEON.path, dem.files[4]))

mosaic.dem <- mosaic(tile1, tile2, tile3, tile4)

plot(mosaic.dem) # so pretty!

out.mosaic.name <- paste0("./processed_tifs/",
                          site, "_", year, "_DEM_mosaic_",
                          date, ".tif")

writeRaster(mosaic.dem, out.mosaic.name,
            overwrite = FALSE)

# -----
# all done!
# -----