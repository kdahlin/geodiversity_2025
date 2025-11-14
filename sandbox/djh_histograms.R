#This script generates histograms of each mosaic's heights and slopes
library(terra)
library(tools)
library(ggplot2)

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
  tifname <- strsplit(tif, "/")[[1]][3]
  tifname <- strsplit(tifname, "\\.")[[1]][1]
  
  #Plots for elevation
  r1 <- rast(tif)
  png(paste0("sandbox/results/", tifname, "_elev_map.png"))
  rp1 <- plot(r1, pax = list(labels = FALSE, tick = FALSE))
  dev.off()
  png(paste0("sandbox/results/", tifname, "_elev_hist.png"))
  rh1 <- hist(r1, maxcell = 4000000)
  dev.off()
  
  #Plots for slope
  s1 <- terrain(r1)
  png(paste0("sandbox/results/", tifname, "_slope_map.png"))
  sp1 <- plot(s1, pax = list(labels = FALSE, tick = FALSE))
  dev.off()
  png(paste0("sandbox/results/", tifname, "_slope_hist.png"))
  sh1 <- hist(s1, maxcell = 4000000)
  dev.off()
}
