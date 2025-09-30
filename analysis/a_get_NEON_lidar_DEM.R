# Title: Get NEON Lidar DEM
# Date: 09/24/2025
# Author: KMD adapted from TRG
# For entering a geographic location and getting the names of the four closest
# NEON tiles, then downloading to this rproject folder

#load packages 
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(cowplot)
library(terra)
library(neonUtilities)

# -----------------------------------
# USER-DEFINED VARIABLES (note there are other hard coded things too, so don't 
# run this without checking! especially at step 5 and beyond)
# -----------------------------------

# today's date
date <- "20250925"

# directory to save "raw" neon data to (working in github)
save.directory <- "./NEON_data/"

# Site Code and Year
site <- "ORNL" 
year <- "2018"  
siteyear <- paste0(site, "/", year, "/")

# make a directory for the data you eventually download
dir.create(paste0(save.directory, site))
save.directory <- paste0(save.directory, site, "/")

# define EPSG code of your spatial data UTM zone (change for new location!)
epsg <- 32616

# what is the approx centroid of where you want data from (in lat/lon)
lon <- -84.33937
lat <- 35.95373

# turn that lat/lon into a sf object for R (with lat/lon epsg)
centroid <- as.data.frame(matrix(data = c(lat, lon), nrow = 1, ncol = 2))
names(centroid) <- c("lat", "lon")

centroid.pt <- st_as_sf(centroid, coords = c("lon", "lat"))

st_crs(centroid.pt) <- 4326 # this is the GEE epsg

#-------------
# RUN LIST AOP TILES FUNCTION
#-------------
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


#-------------
# STEP 1: generate UTM coordinates, ESPG code, and NEON tile coordinates
#-------------

# project data into target UTM Zone 
centroid.rpj <- st_transform(centroid.pt, crs = epsg)

# plot to take a look (would be nice to make a better map here...)
plot(centroid.rpj)

raw_coords_df <- as.data.frame(st_coordinates(centroid.rpj))

# call the function to list the UTM coordinates and coordinates of all the tiles 
UTM_coords_df <- list_AOP_Tiles(raw_coords_df, input_crs = epsg) 

#---------------
# Step 2: get tiles
#---------------

# Set buffer to determine how many adjoining tiles to add around the sampled
# points
# buffer = 0 fills gaps in the input data but doesn't add a buffer
buffer = 1 # in kilometers

# create a data frame of unique tile eastings and northings 
easting <- UTM_coords_df$tile_easting

northing <- UTM_coords_df$tile_northing

tile_coords <- as.data.frame(unique(cbind(easting,northing)))

# Convert the data frame to an sf object
tile_points <- st_as_sf(tile_coords, 
                        coords = c("easting", "northing"), 
                        crs = epsg) 

#-------
# Step 3: plot your points!
#-------

# load basemap 
#world <- ne_countries(scale = "large", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE, crs = epsg))

# Reproject basemaps into UTM
#world_utm  <- st_transform(world, epsg)
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
  coord_sf(xlim = c(eastrefmin, eastrefmax), ylim = c(northrefmin, northrefmax), crs = epsg) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# create map of tiles 
tiles_map <- ggplot(data = states_utm) +
  geom_sf() +
  geom_sf(data = states_utm, fill = NA)+
  geom_sf(data = tile_points, size = 3, shape = 24, fill = "lightgreen")+
  coord_sf(xlim = c(eastmin, eastmax), ylim = c(northmin, northmax), crs = epsg, datum = epsg) +
  theme_bw()

# combine inset and tile maps together and plot 
ggdraw(tiles_map) +
  draw_plot(inset, width = 0.3, height = 0.3, x = 0.15, y = 0.05)

#------
# Step 4 (Optional): Select adjoining tiles 
#------
# run function
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
# If you are satisfied with the tiles plotted in step 2, proceed to step 4. 

tile_coords_new <- adjoin_neon_tiles(tile_coords, kmbuffer = buffer)

# convert tile coordinates to sf object 
tile_points2 <- st_as_sf(tile_coords_new, 
                         coords = c("easting", "northing"), 
                         crs = epsg) 

# create map of adjoining tiles for verification
tiles_map2 <- ggplot(data = states_utm) +
  geom_sf() +
  geom_sf(data = states_utm, fill = NA)+
  geom_sf(data = tile_points2, size = 3, shape = 24, fill = "blue")+
  geom_sf(data = tile_points, size = 3, shape = 24, fill = "lightgreen")+
  coord_sf(xlim = c(eastmin, eastmax), ylim = c(northmin, northmax), crs = epsg, datum = epsg) +
  theme_bw()

# combine inset and new tile map together and plot 
ggdraw(tiles_map2) +
  draw_plot(inset, width = 0.3, height = 0.3, x = 0.15, y = 0.05)

#------
# Step 5: Take a look in GEE to pick the 4 points you want
#------

st_write(tile_points2, paste0(save.directory, "/", site,  "_9points.kml"), 
         driver = "KML")
# open in GEE

# select which coordinates you want to keep (remember you only want the lower
# left corner of each tile)
tile_coords <- tile_coords_new[c(1,2,4,5), ]

#------
# Step 6: Get data from NEON!
#------
for (i in 1:nrow(tile_coords)) {
  neonUtilities::byTileAOP(dpID = "DP3.30024.001",
                         site = site,
                         year = year,
                         easting = tile_coords$easting[i],
                         northing = tile_coords$northing[i],
                         check.size = FALSE,
                         savepath = save.directory)
  print(i)
}

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






