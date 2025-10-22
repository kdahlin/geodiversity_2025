## Get NEON data

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