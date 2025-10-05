#install packages
library(terra)
library(geodiv)
library(raster)
library(mapdata)
library(ggplot2)
library(cowplot)
library(sf)
library(rasterVis)

setwd("C://Users/rache/Documents/geodiversity_2025/geodiversity_2025/sandbox/")

mosaic <- "C://Users/rache/Documents/geodiversity_2025/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif"

r1<-rast(mosaic)
plot(r1)

metrics_slope <- terrain(r1, v="slope")
metrics_aspect <- terrain(r1, v="aspect")
plot(metrics_slope)
plot(metrics_aspect)

###Calulate global surface gradient metrics
#Average roughness
roughness <- sa(r1)
roughness

#Surface bearing index
surfacebearing <- sbi(r1)
surfacebearing

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
