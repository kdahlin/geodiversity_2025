library(terra)
library(geodiv)
install.packages("devtools")
devtools::install_github("bioXgeo/geodiv")
library(raster)
library(rasterVis)
library(mapdata)
library(maptools)
library(rgeos)
library(ggplot2)
library(tidyverse)
library(parallel)
library(sf)
library(rasterVis)
library(ggmap)
library(corrplot)
library(gridExtra)
library(cowplot)
library(factoextra)
library(cluster)

#Practice using the Geodiv Package
data(orforest)
orforest
terra::unwrap(orforst)

eviCols <- colorRampPalette(c('lightyellow1', 'darkgreen'))(100)
eviTheme <- rasterTheme(region=eviCols)
(orig_ndvi<-rasterVis::levelplot(orforest, margin=F, 
                                 par.settings = eviTheme, xlab='Longitude', 
                                 ylab='Latitude', main='orforest original'))

#get path to mosaic data
mosaic_dem<- "C:/Users/Emma/871/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif"
cper_dem<- "C:/Users/Emma/871/geodiversity_2025/processed_tifs/CPER_2020_DEM_mosaic_20251005.tif" #data that courtney loaded
wood_dem<- "C:/Users/Emma/871/geodiversity_2025/processed_tifs/WOOD_2020_DEM_mosaic_20251005.tif"
ornl_dem<- "C:/Users/Emma/871/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif"


r1<-rast(mosaic_dem)
r2<-rast(cper_dem) #data from courtney
r3<-rast(wood_dem)
r4<-rast(ornl_dem)

plot(r1)
plot(r2)
plot(r3)
plot(r4)

#average roughness of the surface
roughness <- sa(r4)
roughness

rtmean <- sq(r4)
rtmean

tenpoint <- s10z(r4)
tenpoint

rms2pt <- sdq(r4)
rms2pt

rms7pt <- sdq6(r4)
rms7pt

#surface area ratio
surfarearatio <- sdr(r4)
surfarearatio

#surface bearing index
surfbi <- sbi(r4)
surfbi

#height bearing area curve
heightbac <- sdc(r4, low=0, high=0.05)
heightbac

#correlation length
corrlength <- scl(r4)
corrlength

#reduced valley depth
redvaldepth <- svk(r4)
redvaldepth

redpk <- spk(r4)
redpk

meanpkht <- smean(r4)
meanpkht

coreroughnessdepth <- sk(r4)
coreroughnessdepth

maxpeakh <- sph(r4)
maxpeakh

valleydep




metrics<-terrain(r1)
plot(metrics)
help(terrain)

aspect<-terrain(r1, v="aspect")
plot(aspect)

library(geodiv)
