#load libraries
library(terra)
library(geodiv)
library(parallel)

#Set a working directory
setwd("C:\\Users\\hallerdi\\Documents\\Dillon_Classes\\geodiversity_2025") #Change this to the location of your own geodiv folder

#create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

#define list of functions to run from geodiv package
functions_list <- (c("sa", "sq", "s10z", "sdq", "sdq6", 
                     "sdr", "sbi","sci","ssk","sku","sds","sfd","srw", "std", 
                     "svi","stxr","ssc","sv","sph","sk",
                     "smean","spk","svk", "scl", "sdc"))

#Define location codes
locations <- c("ORNL",
               "RMNP",
               "CPER",
               "WOOD")

#Define folder locations of tifs
tifs <- c("processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif",
          "processed_tifs/RMNP_2020_DEM_mosaic_20251005.tif",
          "processed_tifs/CPER_2020_DEM_mosaic_20251005.tif",
          "processed_tifs/WOOD_2020_DEM_mosaic_20251005.tif")

#Use mclapply to do all four sites at once
listed_full_results <- mclapply(tifs, function(tif) {
  r1 <- rast(tif)
  #Just use a for loop to traverse the list of functions
  vals <- list()
  for (func in functions_list){
    metric_fun <- get(func, envir = asNamespace("geodiv"))
    if (func == "sdc") { 
      vals <- append(vals, metric_fun(r1, low = 0, high = 0.05))
    }
    else if (func == "stxr" | func == "scl") {
      vals <- append(vals, metric_fun(r1, threshold = 0.2))
    }
    else {
      vals <- append(vals, metric_fun(r1))
    }
  }
  return(vals)
}, mc.cores = 1)#You can change this, but apparently not on Windows

# turn results list into a dataframe 
indices_list <- (c("sa", "sq", "s10z", "sdq", "sdq6", 
                   "sdr", "sbi","sci","ssk","sku","sds","sfd","srw",
                   "srwi","shw","std","stdi",
                   "svi","stxr","ssc","sv","sph","sk",
                   "smean","spk","svk","scllow","sclhigh","sdc"))

results_df_wood <- data.frame(
  func = indices_list,
  ORNL = listed_full_results[[1]],
  RMNP = listed_full_results[[2]],
  CPER = listed_full_results[[3]],
  WOOD = listed_full_results[[4]]
)

print(results_df_wood)

#save results table to csv
write.csv(results_df_wood, "results/wood_geodiversity_metrics.csv", row.names = F)
# DONE -----
