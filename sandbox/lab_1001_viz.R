#Leo Baldiga
#geodiversity testing methods

library(terra)
library(geodiv)

mosaic_path <- ("~/Documents/GitHub/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif")

r1 <- rast(mosaic_path)
plot(r1)
aspect <- terrain(r1, v="aspect")
plot(aspect) 
slope <- terrain(r1, v="slope")
plot(slope)
roughness <- terrain(r1, v="roughness")
plot(roughness)

#Normalize the raster
stats <- global(r1,fun = c("mean", "sd"))
mu <- stats$mean
sigma <- stats$sd
normr <- (r1-mu)/sigma
plot(normr)

#remove plane
r1_rem <- remove_plane(r1)
plot(r1_rem)

#remove plane from normalized raster
normr_rem <- remove_plane(normr)
plot(normr_rem)

#roughness on plane removed raster
sa <- sa(r1)
print(sa)
sa_rem <- sa(r1_rem)
print(sa_rem)
sa_norm <- sa(normr)
print(sa_norm)
sa_norm_rem <- sa(normr_rem)
print(sa_norm_rem)

#surface bearing index on plane removed raster
sbi <- sbi(r1_rem)
print(sbi)

#Root Mean Square Roughness on plane removed raster
sq <- sq(r1_rem)
print(sq)

#Reduced Peak Height
spk <- spk(r1_rem)
print(spk)

#Ten Point Height
s10z <- s10z(r1_rem, create_plot = TRUE)
print(s10z)

#Core Fluid Retention Index
sci <- sci(r1_rem)
print(sci_rem)

#texture direction
std <- std(r1, create_plot = TRUE)
print(std)
std_rem <- std(r1_rem, create_plot = TRUE)
print(std_rem)
std_norm <- std(normr, create_plot = TRUE)
print(std_norm)
std_norm_rem <- std(normr_rem, create_plot = TRUE)
print(std_norm_rem)

#fractal dimension
sfd <- sfd(r1)
print(sfd)
sfd_rem <- sfd(r1_rem)
print(sfd_rem)
sfd_norm <- sfd(normr)
print(sfd_norm)
sfd_norm_rem <- sfd(normr_rem)
print(sfd_norm_rem)


