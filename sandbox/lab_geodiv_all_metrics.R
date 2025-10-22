#load libraries
library(terra)
library(geodiv)
library(parallel)
library(parallelly)
library(snow)
library(dplyr)
library(tidyr)
library(purrr)

#Set a working directory
setwd("/Users/leobaldiga/Documents/GitHub/geodiversity_2025") #Change this to the location of your own geodiv folder

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

#Define folder locations of tifs and get all tifs in those folders
tifs <- c("processed_tifs/ORNL",
          "processed_tifs/RMNP",
          "processed_tifs/CPER",
          "processed_tifs/WOOD")

# Get all TIFF files from the directories
tifs <- unlist(lapply(tifs, function(dir) {
  list.files(path = dir, pattern = "\\.tif$", full.names = TRUE)
}))

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
    } else if ( f == "scl") {
      val <- metric_fun(r, threshold = 0.2)
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
