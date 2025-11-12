# Batch compute geodiversity metrics across multiple NEON sites and DEM tiles
# This script processes all DEM tiles in the 'processed_tifs' directory for 
# specified NEON sites, computing a suite of geodiversity metrics using
# functions from the 'geodiv', 'spatialEco', 'MultiscaleDTM', and 'terra' packages.
# Results are saved in both long and wide CSV formats for further analysis.

#Author: LAB
#Date: 11/05/2025

#load libraries
library(terra)
library(geodiv)
library(parallel)
library(parallelly)
library(snow)
library(dplyr)
library(tidyr)
library(moments)
library(spatialEco)
library(MultiscaleDTM)

### ONLY MODIFY THE FOLLOWING SECTION AS NEEDED ### ------

#Set a working directory
setwd("/Users/leobaldiga/Documents/GitHub/geodiversity_2025/") 
#Change this to the location of your own geodiv folder

#Define NEON site location codes (processed_tifs/subdirectory names)
locations <- c("ORNL", "RMNP", "CPER", "WOOD")
################################################################## 

#--------------------------
# DATA AND TASK PREPARATION 
#--------------------------

#create a results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

#iterate over each NEON site code as a subdirectory of processed_tifs to get tif files
tifs <- c()
for (loc in locations) {
  dir_path <- file.path("processed_tifs", loc)
  tif_files <- list.files(path = dir_path, pattern = "\\.tif$", full.names = TRUE)
  tifs <- c(tifs, tif_files)
}

#define list of functions to run from geodiv and spatialEco packages
functions_list <- (c("sa", "sq", "s10z", "sdq", "sdq6", #geodiv functions
                     "sdr", "sbi","sci","ssk","sku","sds", "sfd","srw", "std", 
                     "svi","stxr","ssc","sv","sph","sk","smean","spk","svk", "scl", "sdc", 
                     "curvature", "tpi", "tri", "vrm",  # SpatialEco functions"
                     "sar", "raster.entropy", 
                     "AdjSD", "RIE")) #MultiscaleDTM functions

# Define tasks: All combinations of files and functions
tasks <- expand.grid(file = tifs, func = functions_list, stringsAsFactors = FALSE)

#------------------------------
# DEFINE A TASK RUNNER FUNCTION
#------------------------------

run_one_task <- function(task) {
  r <- rast(task$file)
  f <- task$func
  
  # Dynamic function lookup:
  # This checks the geodiv namespace first, then spatialEco, then MultiscaleDTM
  if (exists(f, envir = asNamespace("geodiv"), inherits = FALSE)) {
    metric_fun <- get(f, envir = asNamespace("geodiv"))
  } else if (exists(f, envir = asNamespace("spatialEco"), inherits = FALSE)) {
    metric_fun <- get(f, envir = asNamespace("spatialEco"))
  } else if (exists(f, envir = asNamespace("MultiscaleDTM"), inherits = FALSE)) {
    metric_fun <- get(f, envir = asNamespace("MultiscaleDTM"))
  } else {
    stop(paste("Function", f, "not found in geodiv or spatialEco"))
  }
  
  # Handle special cases for function arguments that require default values
  val <- tryCatch({
    if (f == "sdc") {
      metric_fun(r, low = 0, high = 0.05)
    } else if (f == "scl") {
      metric_fun(r, threshold = 0.2)
    } else {
      metric_fun(r)
    }
  }, error = function(e) {
    warning(sprintf("Error in '%s' on file '%s': %s", f, task$file, e$message))
    NA
  })
  
  # If the result is a raster, aggregate to a single value by taking global mean
  if (inherits(val, "SpatRaster")) {
    val <- tryCatch({
      global(val, fun = "mean", na.rm = TRUE)[1,1]
    }, error = function(e) {
      warning(sprintf("Error aggregating raster from '%s' on file '%s': %s", f, task$file, e$message))
      NA
    })
  }
  
  data.frame(file = basename(task$file), func = f, value = I(list(val)))
}

#--------------------------
# PARALLEL PROCESSING SETUP
#--------------------------

# Determine number of cores to use, leaving some free for system processes
# parallelly::availableCores is better than parallel:: detectCores() for this purpose, 
# it won't result in 0 on a single-core machine
num_cores <- parallelly::availableCores(omit = 2) 
print(paste("Using", num_cores, "cores for parallel processing."))

#--------------------------
# RUN TASKS IN PARALLEL
#--------------------------

# Parallel execution - this may take a while, depending on how many files! 
# Go get lunch, or take a nap while it runs. 

# The execution method differs between Windows and Unix-like systems.

# The Windows version uses PSOCKcluster, and parLapply without forking processes.
# Each worker is a separate R session initialized from scratch, 
# so you must explicitly export libraries, functions, and variables to the workers.

if (.Platform$OS.type == "windows") { 
  cl <- makePSOCKcluster(num_cores)
  clusterExport(cl, c("tasks", "run_one_task")) #Export task list and runner function to each worker
  clusterEvalQ(cl, { #Load libraries on each worker
    library(terra)
    library(geodiv)
    library(spatialEco)
    library(moments) #A dependency of spatialEco
    library(MultiscaleDTM)
  })
  
  # parLapply funcs over rows of tasks
  results_list <- parLapply(cl, seq_len(nrow(tasks)), function(i) {
    run_one_task(tasks[i, ])
  })
  stopCluster(cl)
  
#The Unix version uses mclapply, forking internally. 
#This eliminates the need to export variables to and load libraries indivdually 
#on each worker, because they share a memory environment.
#This could get weird with clustered systems sharing resources (i.e. a HPC), 
#but should be fine for local multicore use (one machine with multiple cores).
  
} else { 
  results_list <- mclapply(seq_len(nrow(tasks)), function(i) {
    library(terra) #load libraries inside shared memory space
    library(geodiv)
    library(spatialEco)
    library(moments)
    library(MultiscaleDTM)
    run_one_task(tasks[i, ])
  }, mc.cores = num_cores)
}

# Print results list structure
str(results_list, max.level = 3)

#convert the list to a data frame
results_all <- lapply(results_list, function(df) {
  data.frame(
    file = df$file,
    func = df$func,
    value = unlist(df$value)  # expands vectors with multiple values to rows 
  )
})
print(results_all)

# Bind all results together
results_df <- do.call(rbind, results_all)
print(results_df)

#--------------------------
# TERRAIN METRICS - TERRA PACKAGE
#--------------------------

# Define terrain metrics to compute
terrain_metrics <- c("slope", "aspect", "TPI", "TRI", 
                    "TRIriley", "TRIrmsd", "roughness") #for terrain function

# Create all combinations of files and terrain metrics
terrain_tasks <- expand.grid(file = tifs, metric = terrain_metrics, stringsAsFactors = FALSE)

# Define a task runner function for terrain metrics
run_terrain_task <- function(task) {
  r <- rast(task$file)
  metric <- task$metric
  
  # Special handling for aspect to compute northness and eastness (sin/cos)
  if (metric == "aspect") {
    val_northness <- tryCatch({
      aspect_rad <- terrain(r, v = "aspect", neighbors = 8, unit = "radians")
      northness <- cos(aspect_rad)
      global(northness, fun = "mean", na.rm = TRUE)[1,1]
    }, error = function(e) {
      warning(sprintf("Error in northness on file '%s': %s", task$file, e$message))
      NA
    })
    val_eastness <- tryCatch({
      aspect_rad <- terrain(r, v = "aspect", neighbors = 8, unit = "radians")
      eastness <- sin(aspect_rad)
      global(eastness, fun = "mean", na.rm = TRUE)[1,1]
    }, error = function(e) {
      warning(sprintf("Error in eastness on file '%s': %s", task$file, e$message))
      NA
    })
    data.frame(
      file = basename(task$file),
      func = c("northness", "eastness"),
      value = c(val_northness, val_eastness)
    )
  } else {
    # For other terrain metrics, compute normally
    val <- tryCatch({
      terrain_raster <- terrain(r, v = metric, unit = "degrees", neighbors = 8)
      global(terrain_raster, fun = "mean", na.rm = TRUE)[1,1]
    }, error = function(e) {
      warning(sprintf("Error in terrain '%s' on file '%s': %s", metric, task$file, e$message))
      NA
    })
    data.frame(file = basename(task$file), func = metric, value = val)
  }
}

#--------------------------
# RUN TERRAIN TASKS IN PARALLEL
#--------------------------

# Parallel execution for terrain metrics
if (.Platform$OS.type == "windows") { 
  cl <- makePSOCKcluster(num_cores)
  clusterExport(cl, c("terrain_tasks", "run_terrain_task"))
  clusterEvalQ(cl, {
    library(terra)
  })
  
  terrain_results_list <- parLapply(cl, seq_len(nrow(terrain_tasks)), function(i) {
    run_terrain_task(terrain_tasks[i, ])
  })
  stopCluster(cl)
  
} else {
  terrain_results_list <- mclapply(seq_len(nrow(terrain_tasks)), function(i) {
    library(terra)
    run_terrain_task(terrain_tasks[i, ])
  }, mc.cores = num_cores)
}

# Bind terrain results together
terrain_results_df <- do.call(rbind, terrain_results_list)

#------------------------
# Run Moments metrics
#------------------------

# Compute and print moments metrics
# for each raster file using parallel processing, + range 
compute_moments <- function(file) {
  r <- rast(file)
  vals <- values(r, na.rm = TRUE)
  mom <- moments(vals)
  data.frame(
    file = basename(file),
    stdv = mom["stdv"],
    mad = mom["mad"],
    nmodes = mom["nmodes"],
    range = max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE)
  )
}

# Parallel execution for moments metrics
if (.Platform$OS.type == "windows") { 
  cl <- makePSOCKcluster(num_cores)
  clusterExport(cl, c("tifs", "compute_moments"))
  clusterEvalQ(cl, {
    library(terra)
    library(moments)
  })
  
  moments_results <- parLapply(cl, tifs, function(file) {
    compute_moments(file)
  })
  stopCluster(cl)
  } else {
  moments_results <- mclapply(tifs, function(file) {
    library(terra)
    library(moments)
    compute_moments(file)
  }, mc.cores = num_cores)
}

# Combine moments results into a single data frame
moments_results_df <- do.call(rbind, moments_results)
print(moments_results_df)

#reshape to be file, func, value
moments_results_df <- moments_results_df %>%
  pivot_longer(
    cols = -file,
    names_to = "func",
    values_to = "value"
  )


#--------------------------
# COMBINE AND RESHAPE RESULTS
#--------------------------

# Combine terrain results with previous results
all_results_df <- rbind(results_df, terrain_results_df, moments_results_df)

# Reshape data to LONG format first
long_df <- terrain_results_df %>%
  unnest(value) %>%   
  mutate(value = as.numeric(value)) %>%          
  group_by(func, file) %>%
  mutate(
    row_id = row_number(),
    func = paste0(func, "_", row_id),
    site = sub("_.*", "", file),
    # Extract tile: last number before _Xxm.tif
    tile = sub(".*_(\\d+)_\\d+m\\.tif$", "\\1", file),
    # Extract resolution: number before m.tif
    res  = sub(".*_(\\d+)m\\.tif$", "\\1", file)
  ) %>%
  select(-row_id) %>%
  ungroup()

#make sure value is unlisted (no list cols)
long_df$value <- unlist(long_df$value)

# Pivot WIDER!
wide_df <- long_df %>%
  pivot_wider(
    id_cols = func,
    names_from = c(file),
    values_from = value
  )

#--------------------------
# SAVE RESULTS
#--------------------------

# Write to CSV
write.csv(long_df, "sandbox/results/all_metrics_long.csv", row.names = FALSE)
view(long_df)

write.csv(wide_df, "sandbox/results/all_metrics_wide.csv", row.names = FALSE)
view(wide_df)

#DONE! -----
