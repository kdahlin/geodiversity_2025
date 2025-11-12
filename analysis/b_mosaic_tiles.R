# -----
# Step 7: Mosaic AOP tiles!
# -----
library(sf)
library(terra)

#####Only part that should change between runs############
sitename = "WOOD"
year = "2019" #Year of data acquisition, not current year
#setwd("C:/...")#Uncomment this line if your working directory is not already set to the repo!
##########################################################

date = toString(Sys.Date())

NEON.path <- paste0("./NEON_data/",sitename,"/")#Find path to raw NEON data

dem.files <- list.files(NEON.path, pattern = "DTM.tif$", recursive = TRUE)

#Index of the tile that forms the bottom left of the mosaic (bottom-left index)
#Why these numbers? The raw data is identified by northing/easting. That means
#that when the files are listed alphabetically, they are listed first from South
#to North and then from West to East. Since the raw downloads were 7x7, that means
#that incrementing the index by 1 will pick the next tile up and incrementing by
#7 will pick the next tile over. Since we're doing 2x2 mosaics, we need to start
#from 1, then go up by 2 to 3, then go up by two to five to get the three mosaics
#in the first column, then go up by 14 to proceed two tiles over to start the next
#column of mosaics.
blis <- c(1,3,5,15,17,19,29,31,33)

for (i in 1:length(blis)){
  bli <- blis[i]
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
                            date, "_", bli, "_1m.tif")
  
  writeRaster(mosaic.dem, out.mosaic.name,
              overwrite = FALSE)
}

# -----
# all done!
# -----