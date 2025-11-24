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
#setwd("/Users/leobaldiga/Documents/GitHub/geodiversity_2025/") 

#Define NEON site location codes (processed_tifs/subdirectory names)
locations <- c("ORNL", "RMNP", "CPER", "WOOD", 
               "CLBJ", "MLBS", "OAES", "OSBS", 
               "TEAK", "UNDE", "WREF")

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
  tif_files <- list.files(path = dir_path, pattern = "\\_1m.tif$", full.names = TRUE)
  tifs <- c(tifs, tif_files)
}

# Function to aggregate raster to each target resolution
resample_and_save <- function(input_tif, resolutions) {
  r <- rast(input_tif)
  input_dir <- dirname(input_tif)
  file_base <- sub("_1m$", "", sub("\\.tif$", "", basename(input_tif), ignore.case = TRUE))
  
  # Get original resolution (assumes square pixels)
  orig_res <- res(r)[1]
  
  for (res in resolutions) {
    # Calculate aggregation factor
    agg_factor <- res / orig_res
    
    # Only use integer factors; otherwise, skip or warn about rounding
    if (abs(round(agg_factor) - agg_factor) > .Machine$double.eps^0.5) {
      warning(sprintf("Target res %s not an integer multiple of original (%s). Skipping.", res, orig_res))
      next
    }
    agg_factor <- round(agg_factor)
    
    # Aggregate raster
    r_agg <- aggregate(r, fact = agg_factor, fun = mean, na.rm = TRUE)
    
    # Output filename:
    out_file <- file.path(input_dir, paste0(file_base, "_", res, "m.tif"))
    writeRaster(r_agg, filename = out_file, overwrite = F)
    message("Saved aggregated raster to: ", out_file)
  }
}

#--------------------------
# RESAMPLING RASTERS
#--------------------------

# Loop over files in all locations
for (loc in locations) {
  dir_path <- file.path("processed_tifs", loc)
  tif_files <- list.files(path = dir_path, pattern = "_1m\\.tif$", full.names = TRUE)
  for (tif in tif_files) {
    resample_and_save(tif, target_resolutions)
  }
}

#--------------------------
# PLOT ONE EXAMPLE RASTER AND CHECK RESOLUTIONS
#--------------------------

# Example: Plot original and resampled rasters for the first file
example_tif <- tif_files[1]
print(example_tif)

r_list <- list()

# Load original raster (1m)
r_list[["1m"]] <- rast(example_tif)

# Load each resampled raster
for (res in target_resolutions) {
  resampled_file <- sub("_1m\\.tif$", paste0("_", res, "m.tif"), example_tif)
  r_list[[paste0(res, "m")]] <- rast(resampled_file)
}

# Plot all: original first, then resampled
par(mfrow = c(2, 4))
plot(r_list[["1m"]], main = "Original 1m")
for (res in target_resolutions) {
  plot(r_list[[paste0(res, "m")]], main = paste0(res, "m"))
}

# Print info for each raster
print(r_list[["1m"]])
for (res in target_resolutions) {
  print(r_list[[paste0(res, "m")]])
}

# End of script!