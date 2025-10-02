library(terra)
library(geodiv)

#Get path to mosaic data
mosaic_dem <- "C:\\Users\\dillo\\Documents\\GitHub\\geodiversity_2025\\processed_tifs\\ORNL_2018_DEM_mosaic_20250925.tif"

r1 <- rast(mosaic_dem)
plot(r1)

#Normalize the raster
stats <- global(r1,fun = c("mean", "sd"))
mu <- stats$mean
sigma <- stats$sd
normr <- (r1-mu)/sigma

#Fractal dimension is a complex metric that measures the degree to which an image looks similar at all scales
#This returns 2 regardless of if the raster is normalized
frac <- sfd(r1)
frac_n <- sfd(normr)

#Radial wavelength, returns three values: the dominant wavelength, the dominance of that wavelength, and the wavelength to encompass half of values
rad_dom <- srw(r1, create_plot = TRUE)
rad_dom_n <-srw(normr)


#Dominant texture direction? This may be an interesting one. Does it give us more than aspect?
dir <- std(r1, create_plot = TRUE)
dir_n <- std(normr)
