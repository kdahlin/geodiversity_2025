#load libraries
library(terra)
library(geodiv)
library(parallel)
library(parallelly)
library(snow)

#Set a working directory
setwd("C:/Users/leobaldiga/Documents/GitHub/geodiversity_2025") #Change this to the location of your own geodiv folder

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
tifs <- c("processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif")

#"processed_tifs/RMNP_2020_DEM_mosaic_20251005.tif",
#"processed_tifs/CPER_2020_DEM_mosaic_20251005.tif",
#"processed_tifs/WOOD_2020_DEM_mosaic_20251005.tif"

#robust detect cores available, won't result in 0 or negative
num_cores <- parallelly::availableCores(omit = 2)

# Define worker task: compute all metrics for one TIFF
compute_metrics_for_tif <- function(filepath, funcs) {
  library(terra)
  library(geodiv)
  r <- rast(filepath)
  
  res_list <- lapply(funcs, function(f) {
    metric_fun <- get(f, envir = asNamespace("geodiv"))
    if (f == "sdc") {
      val <- metric_fun(r, low = 0, high = 0.05)
    } else {
      val <- metric_fun(r)
    }
    list(func = f, value = val)
  })
  
  # Combine metrics for this raster
  data.frame(
    file = basename(filepath),
    func = sapply(res_list, `[[`, "func"),
    value = I(lapply(res_list, `[[`, "value"))
  )
}

# --- Cross‑platform run ---
if (.Platform$OS.type == "windows") {
  cl <- makePSOCKcluster(num_cores)
  clusterExport(cl, c("tifs", "functions_list", "compute_metrics_for_tif"))
  clusterEvalQ(cl, { library(terra); library(geodiv) })
  
  results_all <- parLapply(cl, tifs, compute_metrics_for_tif, funcs = functions_list)
  stopCluster(cl)
  
} else {
  # Unix uses forking (faster)
  results_all <- mclapply(tifs, compute_metrics_for_tif, funcs = functions_list,
                          mc.cores = num_cores)
}

# Bind all results together
results_df <- do.call(rbind, results_all)

# Save combined results
write.csv(results_df, "results/all_geodiversity_metrics.csv", row.names = FALSE)

print(head(results_df))

#DONE -----


#----- OLD CODE BELOW -----
# Choose parallel method based on OS
if (.Platform$OS.type == "windows") {
  cl <- makePSOCKcluster(num_cores)
  clusterExport(cl, c("mosaic_path_ornl", "functions_list"))
  clusterEvalQ(cl, {
    library(terra)
    library(geodiv)
  })
  results_list_ornl <- parLapply(cl, functions_list, function(func) {
    r1 <- terra::rast(mosaic_path_ornl)  # Recreate within each worker function
    metric_fun <- get(func, envir = asNamespace("geodiv"))
    if (func == "sdc") {
      value <- metric_fun(r1, low = 0, high = 0.05)
    } else {
      value <- metric_fun(r1)
    }
    list(func = func, value = value)
  })
  stopCluster(cl)
  
} else {
  # Unix version (original)
  results_list_ornl <- mclapply(functions_list, function(func) {
    metric_fun <- get(func, envir = asNamespace("geodiv"))
    if (func == "sdc") {
      value <- metric_fun(r1, low = 0, high = 0.05)
    } else {
      value <- metric_fun(r1)
    }
    list(func = func, value = value)
  }, mc.cores = num_cores)
}



#ORNL -----
mosaic_path_ornl <- ("C:/Users/leobaldiga/Documents/GitHub/geodiversity_2025/processed_tifs/ORNL_2018_DEM_mosaic_20250925.tif")
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
