#install packages
library(terra)
library(geodiv)
library(raster)
library(mapdata)
library(ggplot2)
library(cowplot)
library(sf)
library(rasterVis)

setwd("C://Users/rache/Documents/geodiversity_2025/geodiversity_2025/")

ORNL <- "C://Users/rache/Documents/geodiversity_2025/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif"
CPER <- "C://Users/rache/Documents/geodiversity_2025/geodiversity_2025/processed_tifs/CPER_2020_DEM_mosaic_20251005.tif"
RMNP <- "C://Users/rache/Documents/geodiversity_2025/geodiversity_2025/processed_tifs/RMNP_2020_DEM_mosaic_20251005.tif"
WOOD <- "C://Users/rache/Documents/geodiversity_2025/geodiversity_2025/processed_tifs/WOOD_2020_DEM_mosaic_20251005.tif"

r1<-rast(ORNL)
plot(r1)

ORNL_slope <- terrain(r1, v="slope")
ORNL_aspect <- terrain(r1, v="aspect")
plot(ORNL_slope)
plot(ORNL_aspect)

r2<-rast(RMNP)
plot(RMNP_r1)

###Calulate global surface gradient metrics
#Average roughness
ORNL_roughness <- sa(r1)
ORNL_roughness

#root mean square roughness
ORNL_rms_roughness <- sq(r1)
ORNL_rms_roughness

#ten-point height
ORNL_10ptht <- s10z(r1)
ORNL_10ptht

#root mean square slope of surface, 2-point method
ORNL_rms_slope_2pt <- sdq(r1)
ORNL_rms_slope_2pt

#root mean square slope of surface, 7-point method
ORNL_rms_slope_7pt <- sdq6(r1)
ORNL_rms_slope_7pt

#surface area ratio
ORNL_area <- sdr(r1)
ORNL_area

#Surface bearing index
ORNL_surfacebearing <- sbi(r1)
ORNL_surfacebearing

#core fluid retention index
ORNL_fluid <- sci(r1)
ORNL_fluid

#skewness
ORNL_skewness <- ssk(r1)
ORNL_skewness

#kurtosis
ORNL_kurtosis <- sku(r1)
ORNL_kurtosis

#summit density
ORNL_summit <- sds(r1)
ORNL_summit

#3d fractal dimension
ORNL_fractal <- sfd(r1)
ORNL_fractal

#dominant radial wavelength, radial wavelength index, mean half wavelength
ORNL_radial <- srw(r1)
ORNL_radial

#angle of dominating texture, texture direction index
ORNL_texturedirection <- std(r1)
ORNL_texturedirection

#valley fluid retention index
ORNL_valleyfluid <- svi(r1)
ORNL_valleyfluid

#texture aspect ratio
ORNL_textureaspect <- stxr(r1)
ORNL_textureaspect

#mean summit curvature
ORNL_curvature <- ssc(r1)
ORNL_curvature

#maximum valley depth
ORNL_valleydepth <- sv(r1)
ORNL_valleydepth

#maximum peak height
ORNL_peak<- sph(r1)
ORNL_peak

#bearing area curve height interval
ORNL_bearing<- sdc(r1, low=0, high=0.05)
ORNL_bearing

19. 'sph': maximum peak height
20. 'sk': core roughness depth
21. 'smean': mean peak height
22. 'svk': reduced valley depth
23. 'spk': reduced peak height
24. 'scl': correlation length



#texture direction metrics
texturedirection <- std(r1, create_plot = FALSE, option=1)
texturedirection

###Texture image creation
window <- matrix(1, nrow=7, ncol=7)
system.time(
  output_raster<-focal_metrics(r1, window, metrics=list('sa'), progress=TRUE))
print(output_raster)
rasterVis::levelplot(output_raster[[1]], margin=F, par.settings=eviTheme, 
          ylab=NULL, xlab=NULL, main='Sa')
