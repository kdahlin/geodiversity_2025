

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
### data exploration at 1 m resolution 

# get just 1m resolution
in.data.1m <- subset(in.data.t, res == "1m")

# look at column variances (to remove columns with no variance)
variances <- sapply(in.data.1m, var)

# get the column names with no variance
which.var <- which(variances == 0)

# remove columns with no variance
in.data.1m <- in.data.1m[,-(which.var)]

### OPTIONAL - REMOVE HIGHLY CORRELATED VARIABLES AND WEIRD ONES
remove.some <- c("tri_1", "tpi_1", "sq_1", "sa_1", "s10z_1", "sdq6_1", "sph_1", 
                 "svi_1", "scl_1", "TRI_1", "TRIriley_1", "TRIrmsd_1", 
                 "roughness_1", "std_1")

remove.cols <- which(names(in.data.1m) %in% remove.some)

small.data.1m <- in.data.1m[,-c(remove.cols)]

################################################################################
# clusters on subsites as rows
cl <- hdbscan(small.data.1m, minPts = 3)

print(cl)

plot(cl, show_flat = TRUE)

out.clusters <- cbind(rownames(in.data.1m), cl$cluster)
plot(small.data.1m$range_1, small.data.1m$curvature_total_1,
     col = cl$cluster+1)


# clusters on metrics as rows - transpose (again!)
small.data.t <- as.data.frame(t(small.data.1m))

cl.metrics <- hdbscan(small.data.t, minPts = 2)
print(cl.metrics)
plot(cl.metrics, show_flat = TRUE)
metric.clusters <- cbind(rownames(small.data.t), cl.metrics$cluster)
plot(small.data.t$range_1, small.data.t$curvature_total_1,
     col = cl.metrics$cluster+1)

######################### 
# cluster on PCA metrics - not sure if this is bananas.

pca <- prcomp(small.data.t, center = TRUE, scale = TRUE)
# get the scores/PCs
scores <- as.data.frame(summary(pca)$x)
row.names(scores) <- row.names(small.data.t)

cl.pcs <- hdbscan(scores, minPts = 3)
print(cl.pcs)
pc.clusters <- cbind(rownames(scores), cl.pcs$cluster)
plot(scores[,c(1,2)], type = "n")
text(scores[,c(1,2)], label = row.names(scores), col = cl.pcs$cluster+1)


