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


gitrm\\cached

#get path to mosaic data
mosaic_dem<- "C:/Users/Emma/871/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif"

r1<-rast(mosaic_dem)
plot(r1)

metrics<-terrain(r1)
plot(metrics)
help(terrain)

aspect<-terrain(r1, v="aspect")
plot(aspect)

#copy/paste run this if you haven't gotten HTTPS push to work yet:
create_github_token()

#once you get the token, run this and paste it in
gitcreds::gitcreds_set()

library(geodiv)
