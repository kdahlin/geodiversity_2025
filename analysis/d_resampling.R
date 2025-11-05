# Resampling rasters to multiple resolutions
# Origninal NEON Data are 1x1 Meter resolution
# This script loads the mosaicked rasters and resamples to various resolutions

# Author: LAB
# Date: 2025-11-05

#Load necessary libraries
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

#--------------------------
# DATA AND TASK PREPARATION 
#--------------------------

#iterate over each NEON site code as a subdirectory of processed_tifs to get tif files
tifs <- c()
for (loc in locations) {
  dir_path <- file.path("processed_tifs", loc)
  tif_files <- list.files(path = dir_path, pattern = "\\.tif$", full.names = TRUE)
  tifs <- c(tifs, tif_files)
}

#--------------------------
# RESAMPLING RASTERS
#--------------------------

# Function to resample and save rasters at multiple resolutions
resample_and_save <- function(input_tif, site_code, resolutions, base_dir = "processed_tifs") {
  # Load input raster
  r <- rast(input_tif)
  
  for (res in resolutions) {
    # Compute new resolution in x and y 
    new_res <- c(res, res)
    
    # Resample raster to new resolution using bilinear interpolation (3x3 window)
    r_resampled <- resample(r, rast(ext(r), resolution=new_res, crs=crs(r)), method = "bilinear")
    
    # Create output directory if it doesn't exist
    out_dir <- file.path(base_dir, site_code, paste0("resampled_", res, "m"))
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    # Define output file path
    file_base <- sub("\\.tif$", "", basename(input_tif), ignore.case = TRUE)
    out_file <- file.path(out_dir, paste0(file_base, "_", res, "m.tif"))
    
    # Write resampled raster
    writeRaster(r_resampled, filename = out_file, overwrite = TRUE)
    
    message("Saved resampled raster to: ", out_file)
  }
}

# Loop over all files and sites
for (loc in locations) {
  dir_path <- file.path("processed_tifs", loc)
  tif_files <- list.files(path = dir_path, pattern = "\\.tif$", full.names = TRUE)
  
  for (tif in tif_files) {
    resample_and_save(tif, loc, target_resolutions)
  }
}


#load the first resampled raster from each site and each resolution to verify
for (loc in locations) {
  for (res in target_resolutions) {
    resampled_dir <- file.path("processed_tifs", loc, paste0("resampled_", res, "m"))
    resampled_files <- list.files(path = resampled_dir, pattern = "\\.tif$", full.names = TRUE)
    
    if (length(resampled_files) > 0) {
      r_test <- rast(resampled_files[1])
      print(paste("Site:", loc, "Resolution:", res, "m"))
      print(r_test)
    }
  }
}

#plot example raster - panel of different resolutions for one site
example_site <- locations[1]
par(mfrow=c(2,3))
for (res in target_resolutions) {
  resampled_dir <- file.path("processed_tifs", example_site, paste0("resampled_", res, "m"))
  resampled_files <- list.files(path = resampled_dir, pattern = "\\.tif$", full.names = TRUE)
  if (length(resampled_files) > 0) {
    r_plot <- rast(resampled_files[1])
    plot(r_plot, main=paste(example_site, "-", res, "m"))
  }
}

# End of script
