#load libraries
library(terra)
library(geodiv)
library(parallel)
library(parallelly)
library(snow)
library(dplyr)
library(tidyr)
library(purrr)
library(moments)
library(spatialEco)

#Set a working directory
setwd("/Users/leobaldiga/Documents/GitHub/geodiversity_2025") #Change this to the location of your own geodiv folder

#create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

#define list of functions to run from geodiv and spatialEco packages
functions_list <- (c("sa", "sq", "s10z", "sdq", "sdq6", #geodiv
                     "sdr", "sbi","sci","ssk","sku","sds","sfd","srw", "std", 
                     "svi","stxr","ssc","sv","sph","sk",
                     "smean","spk","svk", "scl", "sdc", 
                     "curvature", "tpi", "tri", "vrm",  # Spatial Eco fns "
                     "sar", "raster.entropy"))

#Define location codes (subdirectory names)
locations <- c("ORNL",
               "RMNP",
               "CPER",
               "WOOD")

#iterate over each location as a subdirectory to get tifs
tifs <- c()
for (loc in locations) {
  dir_path <- file.path("processed_tifs", loc)
  tif_files <- list.files(path = dir_path, pattern = "\\.tif$", full.names = TRUE)
  tifs <- c(tifs, tif_files)
}

# All combinations of files and functions
tasks <- expand.grid(file = tifs, func = functions_list, stringsAsFactors = FALSE)

# Define a task-runner function
run_one_task <- function(task) {
  r <- rast(task$file)
  f <- task$func
  # dynamic function lookup
  if (exists(f, envir = asNamespace("geodiv"), inherits = FALSE)) {
    metric_fun <- get(f, envir = asNamespace("geodiv"))
  } else if (exists(f, envir = asNamespace("spatialEco"), inherits = FALSE)) {
    metric_fun <- get(f, envir = asNamespace("spatialEco"))
  } else {
    stop(paste("Function", f, "not found in geodiv or spatialEco"))
  }
  
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

num_cores <- parallelly::availableCores(omit = 2)

# Wrap the parallel call with required libraries loaded
if (.Platform$OS.type == "windows") {
  cl <- makePSOCKcluster(num_cores)
  clusterExport(cl, c("tasks", "run_one_task"))
  clusterEvalQ(cl, {
    library(terra)
    library(geodiv)
    library(spatialEco)
    library(moments)
  })
  
  # parLapply over rows of tasks
  library(pbapply)  # optional progress bar wrapper
  results_list <- parLapply(cl, seq_len(nrow(tasks)), function(i) {
    run_one_task(tasks[i, ])
  })
  stopCluster(cl)
} else { #unix version
  results_list <- mclapply(seq_len(nrow(tasks)), function(i) {
    library(terra)
    library(geodiv)
    library(spatialEco)
    library(moments)
    run_one_task(tasks[i, ])
  }, mc.cores = num_cores)
}

print(results_list)
#print only values from the results list
results_all <- lapply(results_list, function(df) {
  data.frame(file = df$file, func = df$func, value = df$value[[1]])
})

# Bind all results together
results_df <- do.call(rbind, results_all)

print(results_df)

# Expand list-column 'value' into multiple columns
expanded_df <- results_df %>%
  mutate(value = map(value, as.numeric)) %>%  # ensure numeric vectors
  unnest_wider(value, names_sep = "_val_")   # split multi-values into separate columns

# Pivot wider by 'func' so each metric is row, columns are file + value part
wide_df <- expanded_df %>%
  pivot_wider(
    id_cols = func,
    names_from = file,
    values_from = starts_with("value"),
    names_glue = "{file}_{.value}"
  )

# Save to CSV
write.csv(wide_df, "sandbox/results/all_geodiversity_metrics_expanded.csv", row.names = FALSE)
print(wide_df)

#DONE -----
