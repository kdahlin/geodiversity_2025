#load libraries
library(terra)
library(geodiv)
library(parallel)
library(parallelly)
library(snow)

#create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

#define list of functions to run from geodiv package
functions_list <- (c("sa", "sq", "s10z", "sdq", "sdq6", 
                     "sdr", "sbi","sci","ssk","sku","sds","sfd","srw", "std", 
                     "svi","stxr","ssc","sv","sph","sk",
                     "smean","spk","svk", "scl", "sdc"))

# Cross-platform parallel processing function
cross_platform_parallel <- function(func_list, raster_data, num_cores) {
  
  # Detect OS and use appropriate parallel method
  if (.Platform$OS.type == "windows") {
    
    # Windows: Use parLapply with PSOCK cluster
    cl <- makePSOCKcluster(num_cores)
    
    # Export necessary objects to cluster
    clusterExport(cl, c("raster_data"), envir = environment())
    
    # Load required packages on each worker
    clusterEvalQ(cl, {
      library(terra)
      library(geodiv)
    })
    
    # Run parallel computation
    results <- parLapply(cl, func_list, function(func) {
      metric_fun <- get(func, envir = asNamespace("geodiv"))
      if (func == "sdc") {
        value <- metric_fun(raster_data, low = 0, high = 0.05)
      } else {
        value <- metric_fun(raster_data)
      }
      list(func = func, value = value)
    })
    
    # Stop cluster
    stopCluster(cl)
    
  } else {
    
    # Unix systems (Linux/macOS): Use mclapply
    results <- mclapply(func_list, function(func) {
      metric_fun <- get(func, envir = asNamespace("geodiv"))
      if (func == "sdc") {
        value <- metric_fun(raster_data, low = 0, high = 0.05)
      } else {
        value <- metric_fun(raster_data)
      }
      list(func = func, value = value)
    }, mc.cores = num_cores)
    
  }
  return(results)
}


#ORNL -----
mosaic_path_ornl <- ("~/Documents/GitHub/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif")
r1 <- rast(mosaic_path_ornl)

# Run the cross-platform parallel computation
results_list_ornl <- cross_platform_parallel(functions_list, r1, num_cores = parallelly::availableCores(omit = 2))

# turn results list into a dataframe 
results_df_ornl <- data.frame(
  func = sapply(results_list_ornl, function(x) x$func),
  value = I(lapply(results_list_ornl, function(x) x$value))
)
print(results_df_ornl)

# turn results list into a dataframe 
results_df_ornl <- data.frame(
  func = sapply(results_list_ornl, function(x) x$func),
  value = I(lapply(results_list_ornl, function(x) x$value))
)
print(results_df_ornl)

#save results table to csv
write.csv(results_df_ornl, "results/ornl_geodiversity_metrics.csv", row.names = F)

#RMNP -----
mosaic_path_rmnp <- ("~/Documents/GitHub/geodiversity_2025/processed_tifs/RMNP_2020_DEM_mosaic_20251005.tif")
r2 <- rast(mosaic_path_rmnp)

results_list_rmnp <- list()
results_list_rmnp <- mclapply(functions_list, function(func) {
  metric_fun <- get(func, envir = asNamespace("geodiv"))
  if (func == "sdc") {
    value <- metric_fun(r2, low = 0, high = 0.05)
  } else {
    value <- metric_fun(r2)
  }
  list(func = func, value = value)
}, mc.cores = 10)

# turn results list into a dataframe 
results_df_rnmp <- data.frame(
  func = sapply(results_list_rmnp, function(x) x$func),
  value = I(lapply(results_list_rmnp, function(x) x$value))
)

print(results_df_rnmp)

#save results table to csv
write.csv(results_df_rnmp, "results/rmnp_geodiversity_metrics.csv", row.names = F)

#CPER -----
mosaic_path_cper <- ("~/Documents/GitHub/geodiversity_2025/processed_tifs/CPER_2020_DEM_mosaic_20251005.tif")
r3 <- rast(mosaic_path_cper)

results_list_cper <- list()
results_list_cper <- mclapply(functions_list, function(func) {
  metric_fun <- get(func, envir = asNamespace("geodiv"))
  if (func == "sdc") {
    value <- metric_fun(r3, low = 0, high = 0.05)
  } else {
    value <- metric_fun(r3)
  }
  list(func = func, value = value)
}, mc.cores = 10)

# turn results list into a dataframe 
results_df_cper <- data.frame(
  func = sapply(results_list_cper, function(x) x$func),
  value = I(lapply(results_list_cper, function(x) x$value))
)

print(results_df_cper)

#save results table to csv
write.csv(results_df_cper, "results/cper_geodiversity_metrics.csv", row.names = F)

# WOOD -----
mosaic_path_wood <- ("~/Documents/GitHub/geodiversity_2025/processed_tifs/WOOD_2020_DEM_mosaic_20251005.tif")
r4 <- rast(mosaic_path_wood)

results_list_wood <- list()
results_list_wood <- mclapply(functions_list, function(func) {
  metric_fun <- get(func, envir = asNamespace("geodiv"))
  if (func == "sdc") {
    value <- metric_fun(r4, low = 0, high = 0.05)
  } else {
    value <- metric_fun(r4)
  }
  list(func = func, value = value)
}, mc.cores = 10)

# turn results list into a dataframe 
results_df_wood <- data.frame(
  func = sapply(results_list_wood, function(x) x$func),
  value = I(lapply(results_list_wood, function(x) x$value))
)

print(results_df_wood)

#save results table to csv
write.csv(results_df_wood, "results/wood_geodiversity_metrics.csv", row.names = F)
# DONE -----
