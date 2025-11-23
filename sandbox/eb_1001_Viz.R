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

setwd("C:/Users/Emma/871/geodiversity_2025")

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
ornl_dem<- "C:/Users/Emma/871/geodiversity_2025/processed_tifs/ORNL/ORNL_2018_DEM_mosaic_20250925.tif"
rmnp_dem<- "C:/Users/Emma/871/geodiversity_2025/processed_tifs/RMNP_2020_DEM_mosaic_20251005.tif"

r1<-rast(mosaic_dem)
r2<-rast(cper_dem) #data from courtney
r3<-rast(wood_dem)
r4<-rast(ornl_dem)
r5<-rast(rmnp_dem)

plot(r1)
plot(r2)
plot(r3)
plot(r4)
plot(r5)

ORNL_slope <- terrain(r4, v="slope")
ORNL_aspect <- terrain(r4, v="aspect")
plot(ORNL_slope)
plot(ORNL_aspect)
ORNL_aspect
mean_aspect <- global(ORNL_aspect, fun="mean", na.rm=TRUE)
mean_aspect

##Correct northness aspect
ORNL_aspect2<-cos(terrain(r4, v = "aspect", neighbors = 8, unit = "radians"))
plot(ORNL_aspect2)
ORNL_aspect2
mean_aspect2 <- global(ORNL_aspect2, fun="mean", na.rm=TRUE)
mean_aspect2

##Correct eastness aspect
ORNL_aspect2e<-sin(terrain(r4, v = "aspect", neighbors = 8, unit = "radians"))
plot(ORNL_aspect2e)
ORNL_aspect2e
mean_aspect2e <- global(ORNL_aspect2e, fun="mean", na.rm=TRUE)
mean_aspect2e

################ spatialEco ##############

#NOT geodiv package- trying to measure curvature using spatialEco package
# step 1 install package
install.packages("spatialEco")
library(spatialEco)
library(terra) #already did it above but just to be sure

# step 2 load DEM
elev<-rast("C:/Users/Emma/871/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif")

# step 3 calculate curvature types

profile_curv<-curvature(elev,type="profile")  #profile curvature

planform_curv<-curvature(elev,type="planform") # plan curvature

total_curv<-curvature(elev,type="total")  #total/sum curvature

# step 4 plot

plot(profile_curv,main="ORNL Profile curvature")
plot(planform_curv,main="ORNL Plan curvature")
plot(total_curv,main="ORNL total curvature")

print(total_curv)

########### change scale of the curvature ? ##########3

class(r4)

library(raster)
w <- matrix(1, 3, 3)
ornl <- focal(elev, w = w, fun = mean, na.rm = TRUE)

# calculate curvature w specific scale
curvature_result <- curvature(ornl, scale = 3)  # Example: scale = 3?????

# step 2 load DEM
elev<-rast("C:/Users/Emma/871/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif")

# step 3 calculate curvature types

profile_curv<-curvature(elev,type="profile")  #profile curvature

planform_curv<-curvature(elev,type="planform") # plan curvature

total_curv<-curvature(elev,type="total")  #total/sum curvature

# step 4 plot

plot(profile_curv,main="ORNL Profile curvature")
plot(planform_curv,main="ORNL Plan curvature")
plot(total_curv,main="ORNL total curvature")

print(total_curv)



############### Repeat! ########### Rocky mountain for maybe some more variation...
# step 2 load DEM
elevrm<-rast("C:/Users/Emma/871/geodiversity_2025/processed_tifs/RMNP_2020_DEM_mosaic_20251005.tif")

# step 3 calculate curvature types

profile_curvrm<-curvature(elevrm,type="profile")  #profile curvature

planform_curvrm<-curvature(elevrm,type="planform") # plan curvature

total_curvrm<-curvature(elevrm,type="total")  #total/sum curvature

# step 4 plot

plot(profile_curvrm,main="RMNP Profile curvature")
plot(planform_curvrm,main="RMNP Plan curvature")
plot(total_curvrm,main="RMNP total curvature")

print(total_curvrm)



######### geodiv metrics ###############

#average roughness of the surface
roughness <- sa(r4)
roughness

rtmean <- sq(r4)
rtmean

tenpoint <- s10z(r4)
tenpoint

#2 point slope
rms2pt <- sdq(r4)
rms2pt

#7 point slope
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

#reduced peak height
redpk <- spk(r4)
redpk

#mean peak height
meanpkht <- smean(r4)
meanpkht

#core roughness depth
coreroughnessdepth <- sk(r4)
coreroughnessdepth

#maximum peak height
maxpeakh <- sph(r4)
maxpeakh


metrics<-terrain(r1)
plot(metrics)
help(terrain)

aspect<-terrain(r1, v="aspect")
plot(aspect)


################ New Comparisons - SpatialEco vs Geodiv ###############
library(spatialEco)

#Raster entropy
rEnt <- raster.entropy(r1, d=3, categorical = FALSE, global = TRUE)
opar <- par(no.readonly=TRUE)
par(mfcol=c(2,1))

plot(rEnt, limits=c(2.1962, 2.1974))
par(opar)summary(rEnt)
rEnt
mean_rEnt <- global(rEnt, fun="mean", na.rm=TRUE)
mean_rEnt

