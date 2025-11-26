# ------------------------------------------------------------------------------
# Title: Principal Components analysis (PCA) across scales
# Date Created: 10/20/2025
# Last Update: 11/26/2025
# Authors: KMD
# Purpose: calculate PCAs for different scales of data
# ------------------------------------------------------------------------------

# load packages 
library(corrplot)
library(factoextra)
library(FactoMineR)

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
### data exploration at 1 m resolution 

# get just 1m resolution
in.data.1m <- subset(in.data.t, res == "1m")

# look at column variances (to remove columns with no variance)
variances <- sapply(in.data.1m, var)

# get the column names with no variance
which.var <- which(variances == 0)

# remove columns with no variance
in.data.1m <- in.data.1m[,-(which.var)]

# look at correlations
cors <- cor(in.data.1m)

write.csv(cors, "./sandbox/results/correlation_matrix_1m.csv", row.names = TRUE)

# visualize with corrplot
x11()
corrplot(cors, order = "FPC", method = "circle")

### OPTIONAL - REMOVE HIGHLY CORRELATED VARIABLES AND WEIRD ONES
remove.some <- c("tri_1", "tpi_1", "sq_1", "sa_1", "s10z_1", "sdq6_1", "sph_1", 
                 "svi_1", "scl_1", "TRI_1", "TRIriley_1", "TRIrmsd_1", 
                 "roughness_1", "std_1")

remove.cols <- which(names(in.data.1m) %in% remove.some)

small.data.1m <- in.data.1m[,-c(remove.cols)]
# look at correlations
cors <- cor(small.data.1m)

# visualize with corrplot
x11()
corrplot(cors, order = "FPC", method = "circle")

# now let's remove those columns from the entire data set (all resolutions) so
# we can do PCA on each resolution
remove.cols <- which(names(in.data.t) %in% remove.some)
in.data.small <- in.data.t[,-c(remove.cols)]

# now let's loop through each resolution
################################################################################
sizes <- unique(res)

for (i in 1:length(sizes)) {
  # get just 1m resolution
  in.data <- subset(in.data.small, res == sizes[i])
  
  # look at column variances (to remove columns with no variance)
  variances <- sapply(in.data, var)
  
  # get the column names with no variance
  which.var <- which(variances == 0)
  which.na <- which(is.na(colSums(in.data)))
  
  # remove columns with no variance
  in.data <- in.data[,-c(which.var, which.na)]

  # do pca
  pca <- prcomp(in.data, center = TRUE, scale = TRUE)
  
  # write out variances explained
  write.csv(summary(pca)$importance, paste0("./results/PCA/", "PCA_var_exp_", 
                                 sizes[i], ".csv"), row.names = TRUE)
  # get the scores/PCs
  scores <- as.data.frame(summary(pca)$x)
  
  # let's look at how the PCs correlate with the variables
  pc.corrs <- cor(in.data, scores)
  

  png(filename = paste0("./results/PCA/", "PCA_corrplot_", 
                        sizes[i], ".png"), width = 6.5, height = 5.5, units = "in",
      res = 200)
  corrplot(pc.corrs, tl.cex = 0.6, title = paste(sizes[i], "Resolution"),
           mar = c(0.3,1,1.5,1))
  dev.off()
  
  p <- dim(pc.corrs)[1]
  
  # let's look at which variables are most strongly correlated with which scores
  cors.sorted <- as.data.frame(matrix(NA, nrow = p, ncol = p))
  cors.Rvals.sorted <- as.data.frame(matrix(NA, nrow = p, ncol = p))
  
  for (j in 1:p) {
    names(cors.sorted)[j] <- colnames(pc.corrs)[j]
    get.sort <- sort(abs(pc.corrs[,j]), decreasing = TRUE)
    cors.sorted[,j] <- names(get.sort)
    names(cors.Rvals.sorted)[j] <- paste0(colnames(pc.corrs)[j], "_Rval")
    cors.Rvals.sorted[,j] <- sort(abs(pc.corrs[,j]))
  }
  
  write.csv(cors.sorted, paste0("./results/PCA/", "PCA_cors_sorted_", 
                                sizes[i], ".csv"), row.names = TRUE)
  write.csv(cors.Rvals.sorted, paste0("./results/PCA/", 
                                      "PCA_corsRvals_sorted_", 
                                sizes[i], ".csv"), row.names = TRUE)
  print(paste("done with", i))
}



###### some data visualization steps to be used on an individual scale if wanted
# look at variances
x11()
fviz_eig(pca, addlabels = TRUE)

# biplot
fviz_pca_biplot(pca, repel = TRUE,                 
                col.var = "#2e9fdf", # Variables color
                col.ind = "#696969")  # Individuals color

