library(terra)
library(geodiv)

#Get path to mosaic data
mosaic_dem <- "C:\\Users\\hallerdi\\Documents\\Dillon_Classes\\geodiversity_2025\\processed_tifs\\ORNL_2018_DEM_mosaic_20250925.tif"

r1 <- rast(mosaic_dem)
plot(r1)

#Normalize the raster
stats <- global(r1,fun = c("mean", "sd"))
mu <- stats$mean
sigma <- stats$sd
normr <- (r1-mu)/sigma

##In this section, I'm just going to test to see if normalizing the raster changes the results for any of these.

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
#Texture aspect ratio (not the same as aspect)
##THIS ONE DOES NOT STAY THE SAME AFTER NORMALIZATION
asp_ratio <- stxr(r1)
asp_ratio_n <- stxr(normr)
#Correlation length. This is the minimum and maximum lengths to threshold autocorrelation
#Thus these two numbers form the aspect ratio. This maximum seems not to be affected by normalization
#But the minimum is strongly affected, which explains why the aspect ratio is.
corrlen <- scl(r1, threshold = c(0.2))
corrlen_n <- scl(normr, threshold = c(0.2))
