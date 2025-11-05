
#Batch Processing Script to Load Resampled Rasters at Multiple Resolutions

#Author : LAB
#Date   : 11/05/2025

# Load required libraries

library(terra)
library(tidyverse)

#################################################################
# ONLY CHANGE THESE SETTINGS IF NEEDED

# Set working directory
setwd("/Users/leobaldiga/Documents/GitHub/geodiversity_2025/") 

#Define NEON site location codes (processed_tifs/subdirectory names)
locations <- c("ORNL", "RMNP", "CPER", "WOOD")

# Define target resolutions in meters
target_resolutions <- c(5, 10, 15, 30, 50, 100)

################################################################## 

#load all resampled rasters in separate lists for each resolution as a test
rasters_5m <- list()
rasters_10m <- list()
rasters_15m <- list()
rasters_30m <- list()
rasters_50m <- list()
rasters_100m <- list()

for (loc in locations) {
  for (res in target_resolutions) {
    resampled_dir <- file.path("processed_tifs", loc, paste0("resampled_", res, "m"))
    resampled_files <- list.files(path = resampled_dir, pattern = "\\.tif$", full.names = TRUE)
    for (tif in resampled_files) {
      r <- rast(tif)
      if (res == 5) {
        rasters_5m[[paste(loc, basename(tif), sep="_")]] <- r
      } else if (res == 10) {
        rasters_10m[[paste(loc, basename(tif), sep="_")]] <- r
      } else if (res == 15) {
        rasters_15m[[paste(loc, basename(tif), sep="_")]] <- r
      } else if (res == 30) {
        rasters_30m[[paste(loc, basename(tif), sep="_")]] <- r
      } else if (res == 50) {
        rasters_50m[[paste(loc, basename(tif), sep="_")]] <- r
      } else if (res == 100) {
        rasters_100m[[paste(loc, basename(tif), sep="_")]] <- r
      }
    }
  }
}

#define list of functions to run from geodiv and spatialEco packages
functions_list <- (c("sa", "sq", "s10z", "sdq", "sdq6", #geodiv functions
                     "sdr", "sbi","sci","ssk","sku","sds", "sfd","srw", "std", 
                     "svi","stxr","ssc","sv","sph","sk","smean","spk","svk", "scl", "sdc", 
                     "curvature", "tpi", "tri", "vrm",  # SpatialEco functions"
                     "sar", "raster.entropy", 
                     "AdjSD", "RIE")) #MultiscaleDTM functions

# Define tasks: All combinations of files and functions
tasks_5m <- expand.grid(file = rasters_5m, func = functions_list, stringsAsFactors = FALSE)
tasks_10m <- expand.grid(file = rasters_10m, func = functions_list, stringsAsFactors = FALSE)
tasks_15m <- expand.grid(file = rasters_15m, func = functions_list, stringsAsFactors = FALSE)
tasks_30m <- expand.grid(file = rasters_30m, func = functions_list, stringsAsFactors = FALSE)
tasks_50m <- expand.grid(file = rasters_50m, func = functions_list, stringsAsFactors = FALSE)
tasks_100m <- expand.grid(file = rasters_100m, func = functions_list, stringsAsFactors = FALSE)




