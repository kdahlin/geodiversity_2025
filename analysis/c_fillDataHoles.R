#This script fills holes in the mosaic raster by averaging the values of adjacent
#pixels (Queen's case)
library(terra)
library(tools)

#setwd(C:/...)#Uncomment this line if needed to set working directory to repo

#Define NEON site location codes (subdirectory names)
locations <- c("ORNL",
               "RMNP",
               "CPER",
               "WOOD")

#iterate over each NEON site  as a subdir of processed_tifs to get tif files
tifs <- c()
for (loc in locations) {
  dir_path <- file.path("processed_tifs", loc)
  tif_files <- list.files(path = dir_path, pattern = "\\d.tif$", full.names = TRUE)
  tifs <- c(tifs, tif_files)
}
#Iterate over tifs
for (tif in tifs){
  r1 <- rast(tif)
  #Apply a focal average of a 3x3 window to any NA values
  f1 <- focal(r1, w=3, fun = "mean", na.policy = "only", na.rm=T)
  #Wouldn't normally recommend this, but I am overwriting the original mosaics to
  #save space
  writeRaster(f1, tif, overwrite = T)
}
