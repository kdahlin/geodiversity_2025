#CLUSTER ANALYSIS
#Written by: Leo Baldiga (borrowing some of Kyla's Code for data preparation)

#load libraries
library(factoextra)
library(Hmisc)
library(plotly)
library(dendextend)
library(ggplot2)
library(ggdendro)

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


### OPTIONAL - REMOVE HIGHLY CORRELATED VARIABLES AND WEIRD ONES
remove.metrics <- c("Sa", "Sq", "s10z", "Sdq6", "Sph", "Svi", "Scl1",
                    "TRIv1", "TRIriley", "TRIrmsd", "rough",
                    "Std1", "TPIv2", "TRIv2", "Srw1", "Scl2")

remove.cols <- which(names(in.data.t) %in% remove.metrics)

small.data.t <- in.data.t[,-c(remove.cols)]

# get 1m data
small.data.1m <- subset(small.data.t, res == "1m")

#remove the res from the row names
rownames(small.data.1m) <- gsub("_1m", "", rownames(small.data.1m))

#transform all variables to a z-score
small.data.1m.z <- as.data.frame(scale(small.data.1m))

# K-means clustering of sites ----
set.seed(123) # for reproducibility
png(filename = "./results/kmeans_screeplot_1_40_centers.png", 
    width = 6.5, height = 4, units = "in", res = 300)
fviz_nbclust(small.data.1m.z, kmeans,
             method = "wss",     # elbow: total within‑cluster SS
             k.max = 40) +       # adjust upper bound as needed
  labs(title = "Scree plot of k-means within cluster sum of squares",
       ) +
  theme(axis.text.x = element_text(size = 8))
dev.off()

k_range <- 1:40  
km_list <- lapply(k_range, function(k) kmeans(small.data.1m.z, centers = k,
                                              nstart = 25))

var_explained <- sapply(km_list, function(km) km$betweenss / km$totss)

var_explained_df<- data.frame(
  k = k_range,
  var_explained = sprintf("%.4f", var_explained),          # 4 decimal places
  pct_explained = sprintf("%.2f", 100 * var_explained)     # e.g. 87.34
)

#it takes 20 clusters to explain 80% of the variance, 37 to explain 90%

k_opt <- 6  # replace with chosen k from elbow
kmeans_result <- kmeans(small.data.1m.z,
                        centers = k_opt,
                        iter.max = 100,
                        nstart = 25)
#print results
kmeans_result$centers

# Ward's D minimum variance hierarchical clustering ----
dissimilarity <- dist(small.data.1m.z, method = "euclidean")
hc_km <- hclust(dissimilarity, method = "ward.D2")

# cut into k clusters and color branches
k <- 7
dend <- as.dendrogram(hc_km)
dend_col <- color_branches(dend, k = k) 
#color labels to match branches
dend_col <- color_labels(dend_col, k = k)

png("./results/hierarchical_clustering_1m_wards_colored.png",
    width = 5, height = 9, units = "in", res = 300)
  par(cex = 0.5, mar = c(4.5,1,0,3))
  plot(rev(dend_col), horiz = TRUE, cex.axis = 1.5, 
       cex.lab = 2, xlab = "distance")
dev.off()


#END OF SCRIPT
