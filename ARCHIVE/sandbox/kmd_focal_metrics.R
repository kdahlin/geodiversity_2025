# making maps of geodiv metrics and (maybe) troubleshooting
library(geodiv)

in.data <- rast("./processed_tifs/ORNL/ORNL_2018_DEM_mosaic_2025-10-22_1.tif")

test.sa <- texture_image(in.data, 
                         window_type = "square",
                         size = 10,
                         in_meters = TRUE, 
                         metric = "sa",
                         parallel = TRUE, 
                         nclumps = 100)

focal.sds <- focal_metrics(in.data,
                          window = matrix(1, nrow = 11, ncol = 11),
                          metrics = 'sds',
                          progress = TRUE)
