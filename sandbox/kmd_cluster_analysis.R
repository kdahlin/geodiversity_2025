

library(dbscan)

# see https://cran.r-project.org/web/packages/dbscan/vignettes/hdbscan.html for
# info and links re: hdbscan

#------------------------------------ data org ---------------------------------
# make sure R displays enough digits (some metrics have very small values)
options(digits = 10)

# read in the wide data
in.data <- read.csv("./results/all_metrics_wide.csv", row.names = 1)

# transpose data so each row is a subsite and each column is a metric
in.data.t <- as.data.frame(t(in.data))

# get row names
samples <- rownames(in.data.t)

# split out resolutions
sample.res <- strsplit(samples, "[_\\.]")
res <- unlist(lapply(sample.res, function(x) {x[9]}))

# rename rows so they aren't so clunky
rnames <- row.names(in.data.t)
split.rnames <- strsplit(rnames, "_")

rnames.new <- paste0(split.rnames[[1]][1], "_", split.rnames[[1]][6],
                     "_", res[1])

for (i in 2:length(rnames)) {
  rnames.new <- rbind(rnames.new, 
                      paste0(split.rnames[[i]][1], "_", split.rnames[[i]][6],
                             "_", res[i]))
}

# attach the better names
row.names(in.data.t) <- rnames.new[,1]

################################################################################
### data exploration at 1 m resolution ### need to do this on 36 variables!

# get just 1m resolution
in.data.1m <- subset(in.data.t, res == "1m")

################################################################################
# clusters on subsites as rows
cl <- hdbscan(in.data.1m, minPts = 4)

print(cl)

plot(cl, show_flat = TRUE)

out.clusters <- cbind(rownames(in.data.1m), cl$cluster)

# clusters on metrics as rows (all data?)

cl.metrics <- hdbscan(in.data, minPts = 4)




