#This script generates histograms of each mosaic's heights and slopes
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
  tif_files <- list.files(path = dir_path, pattern = "1m.tif$", full.names = TRUE)
  tifs <- c(tifs, tif_files)
}

#Iterate over tifs
for (tif in tifs) {
  r1 <- rast(tif)
  p1 <- plot(r1, pax = list(labels = FALSE, tick = FALSE))
  rh1 <- hist(r1, maxcell = 4000000)
  s1 <- terrain(r1)
}
