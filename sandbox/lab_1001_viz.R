#Leo Baldiga
#geodiversity testing methods

library(terra)
library(geodiv)
library(sf)
library(tidyverse)
library(viridis)
library(parallel)

mosaic_path <- ("~/Documents/GitHub/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif")

#ORNL ANALYSIS (OAK RIDGE NATIONAL LAB) ----------------------------------------

r1 <- rast(mosaic_path)

#Normalize the raster
stats <- global(r1,fun = c("mean", "sd"))
mu <- stats$mean
sigma <- stats$sd
normr <- (r1-mu)/sigma

#remove plane
r1_rem <- remove_plane(r1)

#remove plane from normalized raster
normr_rem <- remove_plane(normr)

plot(r1)
title("Original Raster")
plot(r1_rem)
title("Plane Removed Raster")
plot(normr)
title("Normalized Raster")
plot(normr_rem)
title("Normalized Plane Removed Raster")

#aspect
aspect <- terrain(r1, v="aspect")
aspect_rem <- terrain(r1_rem, v="aspect")
aspect_norm <- terrain(normr, v="aspect")
aspect_norm_rem <- terrain(normr_rem, v="aspect")

plot(aspect) 
title("Aspect on original raster")
plot(aspect_rem)
title("Aspect on plane removed raster")
plot(aspect_norm)
title("Aspect on normalized raster")
plot(aspect_norm_rem)
title("Aspect on normalized plane removed raster")

#slope
slope <- terrain(r1, v="slope")
slope_rem <- terrain(r1_rem, v="slope")
slope_norm <- terrain(normr, v="slope")
slope_norm_rem <- terrain(normr_rem, v="slope")

plot(slope)
title("Slope on original raster")
plot(slope_rem)
title("Slope on plane removed raster")
plot(slope_norm)
title("Slope on normalized raster")
plot(slope_norm_rem)
title("Slope on normalized plane removed raster")

#roughness
roughness <- terrain(r1, v="roughness")
roughness_rem <- terrain(r1_rem, v="roughness")
roughness_norm <- terrain(normr, v="roughness")
roughness_norm_rem <- terrain(normr_rem, v="roughness")

plot(roughness)
title("Roughness on original raster")
plot(roughness_rem)
title("Roughness on plane removed raster")
plot(roughness_norm)
title("Roughness on normalized raster")
plot(roughness_norm_rem)
title("Roughness on normalized plane removed raster")

#roughness on plane removed raster
sa <- sa(r1)
sa_rem <- sa(r1_rem)
sa_norm <- sa(normr)
sa_norm_rem <- sa(normr_rem)

print(sa)
print(sa_rem)
print(sa_norm)
print(sa_norm_rem)

#surface bearing index on plane removed raster
sbi <- sbi(r1)
sbi_rem <- sbi(r1_rem)
sbi_norm <- sbi(normr)
sbi_norm_rem <- sbi(normr_rem)

print(sbi)
print(sbi_rem)
print(sbi_norm)
print(sbi_norm_rem)

#Root Mean Square Roughness on plane removed raster
sq <- sq(r1)
sq_rem <- sq(r1_rem)
sq_norm <- sq(normr)
sq_norm_rem <- sq(normr_rem)

print(sq)
print(sq_rem)
print(sq_norm)
print(sq_norm_rem)

#Reduced Peak Height
spk <- spk(r1)
spk_rem <- spk(r1_rem)
spk_norm <- spk(normr)
spk_norm_rem <- spk(normr_rem)

print(spk)
print(spk_rem)
print(spk_norm)
print(spk_norm_rem)

#Ten Point Height
s10z <- s10z(r1)
s10z_rem <- s10z(r1_rem)
s10z_norm <- s10z(normr)
s10z_norm_rem <- s10z(normr_rem)

print(s10z)
print(s10z_rem)
print(s10z_norm)
print(s10z_norm_rem)

#Core Fluid Retention Index
sci <- sci(r1)
sci_rem <- sci(r1_rem)
sci_norm <- sci(normr)
sci_norm_rem <- sci(normr_rem)

print(sci)
print(sci_rem)
print(sci_norm)
print(sci_norm_rem)

#texture direction
std <- std(r1, create_plot = TRUE)
std_rem <- std(r1_rem, create_plot = TRUE)
std_norm <- std(normr, create_plot = TRUE)
std_norm_rem <- std(normr_rem, create_plot = TRUE)

print(std)
print(std_rem)
print(std_norm)
print(std_norm_rem)

#fractal dimension
sfd <- sfd(r1)
sfd_rem <- sfd(r1_rem)
sfd_norm <- sfd(normr)
sfd_norm_rem <- sfd(normr_rem)

print(sfd)
print(sfd_rem)
print(sfd_norm)
print(sfd_norm_rem)


#RMNP ANALYSIS (ROCKY MOUNTAIN NATIONAL PARK) ----------------------------------

#Get path to RMNP mosaic data
mosaic_path_rmnp <- ("~/Documents/GitHub/geodiversity_2025/processed_tifs/RMNP_2020_DEM_mosaic_20251005.tif")

r2 <- rast(mosaic_path_rmnp)
functions_list <- (c("sa", "sq", "s10z", "sdq", "sdq6", 
                   "sdr", "sbi","sci","ssk","sku","sds","sfd","srw", "std", 
                   "svi","stxr","ssc","sv","sph","sk",
                   "smean","spk","svk", "scl", "sdc"))
  
results_list <- list()
results_list <- mclapply(functions_list, function(func) {
  metric_fun <- get(func, envir = asNamespace("geodiv"))
  if (func == "sdc") {
    value <- metric_fun(r2, low = 0, high = 0.05)
  } else {
    value <- metric_fun(r2)
  }
  list(func = func, value = value)
}, mc.cores = 10)

# turn results list into a dataframe 
results_df <- data.frame(
  func = sapply(results_list, function(x) x$func),
  value = I(lapply(results_list, function(x) x$value))
)

print(results_df)

#create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

#save results table to csv
write.csv(results_df, "results/rmnp_geodiversity_metrics.csv", row.names = F)

#do the same for r1 - ORNL
results_list_ORNL <- list()

results_list_ORNL <- mclapply(functions_list, function(func) {
  metric_fun <- get(func, envir = asNamespace("geodiv"))
  if (func == "sdc") {
    value <- metric_fun(r1, low = 0, high = 0.05)
  } else {
    value <- metric_fun(r1)
  }
  list(func = func, value = value)
}, mc.cores = 10)

# turn results list into a dataframe 
results_df_ORNL <- data.frame(
  func = sapply(results_list_ORNL, function(x) x$func),
  value = I(lapply(results_list_ORNL, function(x) x$value))
)

print(results_df_ORNL)

#save results table to csv
write.csv(results_df_ORNL, "results/ornl_geodiversity_metrics.csv", row.names = F)


#HEXAGONAL GRIDS BELOW - ENTER AT YOUR OWN RISK --------------------------------

# Create hexagonal grid over the raster
hexes <- st_make_grid(r1, cellsize = 100, square = F)
# Assign unique hex_id to each hexagon
hexes_sf <- st_sf(geometry = hexes) %>%
  mutate(hex_id = row_number())
plot(r1)
plot(st_geometry(hexes_sf), add = TRUE, border = 'red')
title("Hexagonal Grid over Raster")

# Add hex_id to raster cells based on which hexagon they fall into
r1_hex <- rasterize(vect(hexes_sf), r1, field = "hex_id", touches = TRUE)

#calculate mean elevation for each hexagon
hex_means <- zonal(r1, r1_hex, fun = 'mean', na.rm = TRUE)
head(hex_means)

#plot the hex means on a map with viridis color scale
hex_means <- as.data.frame(hex_means)
colnames(hex_means) <- c("hex_id", "mean_elevation")
hex_means <- left_join(hexes_sf, hex_means, by = "hex_id")

#plot the hex means
plot_sf <- function(sf_data) {
  ggplot() +
    geom_sf(data = sf_data, aes(fill = mean_elevation), color = NA) +
    scale_fill_viridis(option = "C", na.value = "transparent") +
    theme_minimal() +
    labs(title = "Mean Elevation in 100m Hexes ORNL", fill = "Mean Elevation")
}
plot_sf(hex_means)

#calculate mean roughness for each hexagon
hex_roughness <- zonal(roughness, r1_hex, fun = 'mean', na.rm = TRUE)
head(hex_roughness)

#plot the hex roughness on a map with viridis color scale
hex_roughness <- as.data.frame(hex_roughness)
colnames(hex_roughness) <- c("hex_id", "mean_roughness")
hex_roughness <- left_join(hexes_sf, hex_roughness, by = "hex_id")

#plot the hex roughness
plot_sf <- function(sf_data) {
  ggplot() +
    geom_sf(data = sf_data, aes(fill = mean_roughness), color = NA) +
    scale_fill_viridis(option = "C", na.value = "transparent") +
    theme_minimal() +
    labs(title = "Mean Roughness in 100m Hexes ORNL", fill = "Mean Roughness")
}
plot_sf(hex_roughness)

#calculate mean slope for each hexagon
hex_slope <- zonal(slope, r1_hex, fun = 'mean', na.rm = TRUE)
head(hex_slope)
#plot the hex slope on a map with viridis color scale
hex_slope <- as.data.frame(hex_slope)
colnames(hex_slope) <- c("hex_id", "mean_slope")
hex_slope <- left_join(hexes_sf, hex_slope, by = "hex_id")
#plot the hex slope
plot_sf <- function(sf_data) {
  ggplot() +
    geom_sf(data = sf_data, aes(fill = mean_slope), color = NA) +
    scale_fill_viridis(option = "C", na.value = "transparent") +
    theme_minimal() +
    labs(title = "Mean Slope in 100m Hexes ORNL", fill = "Mean Slope")
}
#save the plot to figures directory
plot_sf(hex_slope)
ggsave("figures/mean_slope_hexes.png", width = 8, height = 6)


#calculate mean aspect for each hexagon
hex_aspect <- zonal(aspect, r1_hex, fun = 'mean', na.rm = TRUE)
head(hex_aspect)
#plot the hex aspect on a map with viridis color scale
hex_aspect <- as.data.frame(hex_aspect)
colnames(hex_aspect) <- c("hex_id", "mean_aspect")
hex_aspect <- left_join(hexes_sf, hex_aspect, by = "hex_id")
#plot the hex aspect
plot_sf <- function(sf_data) {
  ggplot() +
    geom_sf(data = sf_data, aes(fill = mean_aspect), color = NA) +
    scale_fill_viridis(option = "C", na.value = "transparent") +
    theme_minimal() +
    labs(title = "Mean Aspect in 100m Hexes ORNL", fill = "Mean Aspect")
}
plot_sf(hex_aspect)

#RMNP Plotting --------------------------------

plot(r2)
#get min and max of r2
r2_min <- global(r2, fun = 'min', na.rm = TRUE)
r2_max <- global(r2, fun = 'max', na.rm = TRUE)
r2_min
r2_max

#get mean of r2
r2_mean <- global(r2, fun = 'mean', na.rm = TRUE)
r2_mean
#get the extent of r2
r2_ext <- ext(r2)
r2_ext
#convert UTM extent to Decimal Degrees
library(sf)
library(sp)
library(rgdal)
# Create a SpatialPoints object with the UTM coordinates
utm_coords <- data.frame(
  x = c(r2_ext[1], r2_ext[2]),
  y = c(r2_ext[3], r2_ext[4])
)
coordinates(utm_coords) <- ~x+y
proj4string(utm_coords) <- CRS("+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs")
# Transform to latitude/longitude
latlong_coords <- spTransform(utm_coords, CRS("+proj=longlat +datum=WGS84"))
latlong_coords
#get the centroid of the extent
centroid <- c((r2_ext[1] + r2_ext[2]) / 2, (r2_ext[3] + r2_ext[4]) / 2)
centroid
#convert centroid to lat long
centroid_sp <- data.frame(x = centroid[1], y = centroid[2])
coordinates(centroid_sp) <- ~x+y
proj4string(centroid_sp) <- CRS("+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs")
centroid_latlong <- spTransform(centroid_sp, CRS("+proj=longlat +datum=WGS84"))
centroid_latlong




