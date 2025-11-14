# Title: Get NEON Lidar DEM
# Date: 09/24/2025
# Author: CNR adapted from KMD
# Enter a geographic location to identify a set number of nearby
# NEON tiles, visualize them in GEE, and download to a folder

# load packages 
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(cowplot)
library(terra)
library(neonUtilities)

# prevent R from reporting values in scientific notation
options(scipen=999)

# --------------------------------------------------
# USER-DEFINED VARIABLES- should change before runs
# --------------------------------------------------

# today's date
date <- "20251114"

# directory to save "raw" neon data to (working in github)
save.directory <- "./NEON_data/"

# site code and year
site <- "CPER" 
year <- "2020"  

# make a directory for the data you eventually download
dir.create(paste0(save.directory, site))
save.directory <- paste0(save.directory, site, "/")


# Set buffer (KILOMETERS) to determine how many adjoining tiles to add around
# the sampled points
# ex. buffer = 0 fills gaps in the input data but doesn't add a buffer
# ex. buffer = 3 creates a 7x7km grid and pulls a total of 49 tiles
buffer = 3

# define EPSG code of your spatial data UTM zone
# four NEON sites are listed below

# WOOD
#epsg <- 32614

# CPER
epsg <- 32613

# RMNP
#epsg <- 32613

# ORNL
#epsg <- 32616

# what is the approx centroid of where you want data from (in lat/lon)
# four NEON sites are listed below

# WOOD
#lon <- -99.26000
#lat <- 47.12000

# CPER
lon <- -104.74559
lat <- 40.81554

# RMNP
#lon <- -105.5171737
#lat <- 40.2627492

# ORNL
#lon <- -84.3261184
#lat <- 35.9337824

# --------------------------------
# CREATE data frame OF PT LOCATION
# --------------------------------

# turn that lat/lon into a data frame
centroid <- as.data.frame(matrix(data = c(lon, lat), nrow = 1, ncol = 2))
names(centroid) <- c("lon", "lat")

#-----------------------------------------------
# RUN LIST AOP TILES FUNCTION TO IDENTIFY TILES
#-----------------------------------------------
# note, this function is written assuming you have a list of lat/lon points,
# but here we are only doing one 

list_AOP_Tiles <- function(coords, input_crs = 4326) {
  # coords: matrix or data.frame with two columns (X, Y), or optionally three 
  # with an ID column
  # input_crs: EPSG code (numeric or character) or proj4 string
  
  if (!is.matrix(coords) && !is.data.frame(coords)) {
    stop("coords must be a matrix or data.frame with two columns (X, Y) or
         optionally three with an ID")
  }
  if (ncol(coords) < 2) {
    stop("coords must have at least two columns")
  }
  
  # Check for ID column
  has_id <- is.data.frame(coords) && "ID" %in% colnames(coords)
  ids <- if (has_id) coords$ID else seq_len(nrow(coords))
  
  # Warn if default input_crs was assumed
  if (identical(input_crs, 4326) && !("input_crs" %in% names(match.call()))) {
    warning(
      "input_crs not specified; assuming EPSG:4326 (WGS84 lon/lat).\n",
      "If your coordinates are in a different CRS, please specify input_crs 
      explicitly.\n",
      "You can obtain the CRS of your spatial data in R using 
      sf::st_crs(your_data)."
    )
  }
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' required but not installed.")
  }
  
  # Create sf object from X and Y only
  pts_sf <- sf::st_as_sf(
    data.frame(x = coords[,1], y = coords[,2]),
    coords = c("x", "y"),
    crs = input_crs
  )
  
  pts_wgs84 <- sf::st_transform(pts_sf, crs = 4326)
  lonlat <- sf::st_coordinates(pts_wgs84)
  
  results <- vector("list", nrow(coords))
  
  for (i in seq_len(nrow(coords))) {
    lon <- lonlat[i, 1]
    lat <- lonlat[i, 2]
    
    utm_zone <- floor((lon + 180) / 6) + 1
    crs_utm <- if (lat >= 0) 32600 + utm_zone else 32700 + utm_zone
    
    pt_utm <- sf::st_transform(pts_sf[i, ], crs = crs_utm)
    coords_utm <- sf::st_coordinates(pt_utm)
    
    easting <- round(coords_utm[1], -3)
    northing <- round(coords_utm[2], -3)
    
    results[[i]] <- data.frame(
      ID = ids[i],
      easting = coords_utm[1],
      northing = coords_utm[2],
      tile_easting = easting,
      tile_northing = northing,
      EPSG = crs_utm
    )
  }
  
  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}

#----------------------------------------------------------------
# GENERATE UTM COORDINATES, ESPG CODE, AND NEON TILE COORDINATES
#----------------------------------------------------------------
# call the function to list the UTM coordinates and coordinates of all the tiles 
UTM_coords_df <- list_AOP_Tiles(centroid) 

#----------------
# GET NEON TILES
#----------------
# create a data frame of unique tile eastings and northings 
easting <- UTM_coords_df$tile_easting

northing <- UTM_coords_df$tile_northing

tile_coords <- as.data.frame(unique(cbind(easting,northing)))

# Convert the data frame to an sf object
tile_points <- st_as_sf(tile_coords, 
                        coords = c("easting", "northing"), 
                        crs = epsg) 

#----------
# PLOT Point
#----------

# load basemap 
#world <- ne_countries(scale = "large", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE, crs = epsg))

# Reproject basemaps into UTM
# world_utm  <- st_transform(world, epsg)
states_utm <- st_transform(states, epsg)

# specify the tile map bounding box coordinates 
eastmax <- max(tile_coords$easting) + 10000
eastmin <- min(tile_coords$easting) - 10000
northmax <- max(tile_coords$northing) + 10000
northmin <- min(tile_coords$northing) - 10000

# expand bbox limits for inset map
eastrefmax <- eastmax + 500000
eastrefmin <- eastmin - 500000
northrefmax <- northmax + 500000
northrefmin <- northmin - 500000

# create polygon of the bbox for inset map
bbox <- as.data.frame(cbind(c(eastmax, eastmin), c(northmax, northmin))) %>% 
  st_as_sf(coords = c("V1", "V2"), crs = epsg) %>% 
  st_bbox() %>%
  st_as_sfc

# create inset map
inset <- ggplot(data = states_utm) +
  geom_sf(fill = "lightblue") +
  geom_sf(data = states_utm, fill = "grey")+
  geom_sf(data = bbox, fill = "red")+
  coord_sf(xlim = c(eastrefmin, eastrefmax), ylim = c(northrefmin, northrefmax),
           crs = epsg) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# create map of tiles 
tiles_map <- ggplot(data = states_utm) +
  geom_sf() +
  geom_sf(data = states_utm, fill = NA)+
  geom_sf(data = tile_points, size = 3, shape = 24, fill = "lightgreen")+
  coord_sf(xlim = c(eastmin, eastmax), ylim = c(northmin, northmax), 
           crs = epsg, datum = epsg) +
  theme_bw()

# combine inset and tile maps together and plot 
ggdraw(tiles_map) +
  draw_plot(inset, width = 0.3, height = 0.3, x = 0.15, y = 0.05)

#------------------------
# SELECT ADJOINING TILES
#------------------------
# run function
## need a list of 6 northings and eastings
## each 1000m apart of eachother

adjoin_neon_tiles <- function(coords, kmbuffer = 1) {
  # coords: data.frame with columns "easting" and "northing"
  # kmbuffer: number of km outward to extend in all directions
  
  if(!all(c("easting", "northing") %in% names(coords))) {
    stop("coords must have columns 'easting' and 'northing'")
  }
  
  # Convert buffer to meters
  buf_m <- kmbuffer * 1000
  
  # Find bounding box with buffer
  east_min <- min(coords$easting) - buf_m
  east_max <- max(coords$easting) + buf_m
  north_min <- min(coords$northing) - buf_m
  north_max <- max(coords$northing) + buf_m
  
  # Generate full grid of tiles
  east_seq <- seq(east_min, east_max, by = 1000)
  north_seq <- seq(north_min, north_max, by = 1000)
  
  grid <- expand.grid(easting = east_seq,
                      northing = north_seq)
  
  # Return full set
  return(grid[order(grid$northing, grid$easting), ])
}

# This step is for users who want to download a continuous bounding box 
# of tiles or add a buffer to the tiles. 

tile_coords_new <- adjoin_neon_tiles(tile_coords, kmbuffer = buffer)

# convert tile coordinates to sf object 
tile_points2 <- st_as_sf(tile_coords_new, 
                         coords = c("easting", "northing"), 
                         crs = epsg) 

# create map of adjoining tiles for verification
tiles_map2 <- ggplot(data = states_utm) +
  geom_sf() +
  #geom_sf(data = states_utm, fill = NA)+
  geom_sf(data = tile_points2, size = 3, shape = 24, fill = "blue")+
  geom_sf(data = tile_points, size = 3, shape = 24, fill = "lightgreen")+
  coord_sf(xlim = c(eastmin, eastmax), ylim = c(northmin, northmax), 
           crs = epsg, datum = epsg) +
  theme_bw()

# combine inset and new tile map together and plot 
ggdraw(tiles_map2) +
  draw_plot(inset, width = 0.3, height = 0.3, x = 0.15, y = 0.05)

#-----------------------------------------------
# DOWNLOAD PTS TO VISUALIZE IN GEE (optional)
#-----------------------------------------------

# st_write(tile_points2, paste0(save.directory, "/", site,  "_49points.shp"), 
         # driver = "ESRI Shapefile")
# open in GEE

#-----------------------------------------------
# DOWNLOAD DATA FROM NEON AND SAVE TO DIRECTORY
#-----------------------------------------------
for (i in 1:nrow(tile_coords_new)) {
  neonUtilities::byTileAOP(dpID = "DP3.30024.001",
                           site = site,
                           year = year,
                           easting = tile_coords_new$easting[i],
                           northing = tile_coords_new$northing[i],
                           check.size = FALSE,
                           savepath = save.directory)
  print(i)
}
