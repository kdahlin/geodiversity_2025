library(terra)


#Set a working directory
setwd("/Users/hallerdi/Documents/Dillon_Classes/geodiversity_2025")

#Define NEON site location codes (subdirectory names)
locations <- c("ORNL",
               "RMNP",
               "CPER",
               "WOOD")

#iterate over each NEON site  as a subdir of processed_tifs to get tif files
tifs <- c()
for (loc in locations) {
  dir_path <- file.path("processed_tifs", loc)
  tif_files <- list.files(path = dir_path, pattern = "\\.tif$", full.names = TRUE)
  tifs <- c(tifs, tif_files)
}
for (tif in tifs) {
  r1 = rast(tif)
  if (any(is.na(values(r1)))){
    print(tif)
  }
}
