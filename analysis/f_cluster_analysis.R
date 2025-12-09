#CLUSTER ANALYSIS
#Written by: Leo Baldiga (borrowing some of Kyla's Code for data preparation)

#load libraries
library(corrplot)
library(factoextra)
library(Hmisc)
library(plotly)
library(dendextend)

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

#remove the res from the row names
rownames(in.data.1m) <- gsub("_1m", "", rownames(in.data.1m))

# look at column variances (to remove columns with no variance)
variances <- sapply(in.data.1m, var)

# get the column names with no variance
which.var <- which(variances == 0)

# remove columns with no variance
in.data.1m <- in.data.1m[,-(which.var)]

# look at correlations
cors <- cor(in.data.1m)

#write.csv(cors, "./results/correlation_matrix_1m.csv", row.names = TRUE)

# visualize with corrplot
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
corrplot(cors, order = "FPC", method = "circle")

#transform all variables to a z-score
small.data.1m.z <- as.data.frame(scale(small.data.1m))

# look at correlations
cors <- cor(small.data.1m.z)
# visualize with corrplot (Should be no change)
corrplot(cors, order = "FPC", method = "circle")

# Variable Cluster Analyses -----

#using spearmans rank and Hmisc::varclus
vc_spearmans <- varclus(~ ., data = small.data.1m.z , similarity = "spearman") 

#Plot dendrogram (axes show 1 - R^2 with own cluster)
png(filename = "./results/varclus_spearmans_1m.png",
    width = 12, height = 8, units = "in", res = 300)
plot(vc_spearmans, main = "Hmisc::varclus variable clustering")
dev.off()

#Extract the hclust object and cut the tree to see clusters
hc_spearmans <- vc_spearmans$hclust
hc_spearmans

# Example: cut at height corresponding roughly to |r| >= 0.9
# varclus uses 1 - R^2 as height, so R^2 = 0.9^2 = 0.81 -> height ~ 1 - 0.81 = 0.19
clust_vc_spear <- cutree(hc_spearmans, k = 5)
clusters_vc_spear <- split(names(clust_vc_spear), clust_vc_spear)
clusters_vc_spear

#Using pearson and Hmisc::varclus
vc_pearson <- varclus(~ ., data = small.data.1m.z , similarity = "pearson")

#Plot dendrogram (axes show 1 - R^2 with own cluster)
png(filename = "./results/varclus_pearson_1m.png",
    width = 12, height = 8, units = "in", res = 300)
plot(vc_pearson, main = "Hmisc::varclus variable clustering (Pearson)")
dev.off()

#Extract the hclust object and cut the tree
hc_pearson <- vc_pearson$hclust
hc_pearson
# Example: cut at height corresponding roughly to |r| >= 0.9
# varclus uses 1 - R^2 as height, so R^2 = 0.9^2 = 0.81 -> height ~ 1 - 0.81 = 0.19
clust_vc_pear <- cutree(hc_pearson, k = 5)
clusters_vc_pear <- split(names(clust_vc_pear), clust_vc_pear)
clusters_vc_pear

# K-means clustering of sites ----
set.seed(123) # for reproducibility
png(filename = "./results/kmeans_screeplot_1_20_centers.png", 
    width = 10, height = 8, units = "in", res = 300)
fviz_nbclust(small.data.1m.z, kmeans,
             method = "wss",     # elbow: total within‑cluster SS
             k.max = 20) +       # adjust upper bound as needed
  labs(title = "Elbow method for k-means")
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

k_opt <- 11  # replace with chosen k from elbow - I chose 11 because there are 11 sites
kmeans_result <- kmeans(small.data.1m.z,
                        centers = k_opt,
                        iter.max = 100,
                        nstart = 25)
#print results
kmeans_result$centers

# Assign variables to clusters based on highest loading
variable_clusters_km <- split(names(kmeans_result$cluster), kmeans_result$cluster)
variable_clusters_km

#turn the variable clusters into a table for easier viewing
max_cluster_size <- max(sapply(variable_clusters_km, length))
variable_clusters_km_df <- as.data.frame(matrix(NA, nrow = max_cluster_size,
                                                ncol = length(variable_clusters_km)))
colnames(variable_clusters_km_df) <- paste0("Cluster_", seq_along(variable_clusters_km))
for (i in seq_along(variable_clusters_km)) {
  cluster_vars <- variable_clusters_km[[i]]
  variable_clusters_km_df[1:length(cluster_vars), i] <- cluster_vars
}

#export to excel for table viewing
write.csv(variable_clusters_km_df, "./results/kmeans_clusters_k11_table.csv",
          row.names = FALSE)

# Visualize k-means clustering results of Principal Components
png(filename = "./results/kmeans_clustering_1m_PC1_2.png", width = 10, height = 8, units = "in", res = 300)
fviz_cluster(kmeans_result, data = small.data.1m.z,
             palette = "jco",
             title = "K-means Clustering of sites (PC 1,2) (k=11, 100 iterations, 25 starts)",
             ggtheme = theme_minimal())
dev.off()

png(filename = "./results/kmeans_clustering_1m_PC1_3.png", width = 10, height = 8, units = "in", res = 300)
fviz_cluster(kmeans_result, data = small.data.1m.z,
             palette = "jco",
             title = "K-means Clustering of sites (PC 1,3) (k=11, 100 iterations, 25 starts)",
             axes = c(1, 3),
             ggtheme = theme_minimal())
dev.off()

png(filename = "./results/kmeans_clustering_1m_PC2_3.png", width = 10, height = 8, units = "in", res = 300)
fviz_cluster(kmeans_result, data = small.data.1m.z,
             palette = "jco",
             title = "K-means Clustering of sites (PC 2,3) (k=11, 100 iterations, 25 starts)",
             axes = c(2, 3),
             ggtheme = theme_minimal())
dev.off()

# PCA on the same data used for k-means
pca_res <- prcomp(small.data.1m.z, scale. = FALSE)

# First 3 PCs as a data frame
pc_scores <- as.data.frame(pca_res$x[, 1:3])
colnames(pc_scores) <- c("PC1", "PC2", "PC3")

# Add cluster labels from k-means
pc_scores$cluster <- factor(kmeans_result$cluster)

plot_ly(pc_scores,
        x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~cluster,
        stroke = "Set3",
        type = "scatter3d",
        mode = "markers") %>%
  layout(title = "K-means (k=11) clusters on first 3 principal components")

# Ward's D minimum variance hierarchical clustering ----
dissimilarity <- dist(small.data.1m.z, method = "euclidean")
hc_km <- hclust(dissimilarity, method = "ward.D2")

# cut into k clusters and color branches
k <- 11
dend <- as.dendrogram(hc_km)
dend_col <- color_branches(dend, k = k) 
#color labels to match branches
dend_col <- color_labels(dend_col, k = k)

png("./results/hierarchical_clustering_1m_wardsD2_colored.png",
    width = 20, height = 10, units = "in", res = 300)
plot(dend_col, main = "Ward's D minimum variance (k=11)",
     xlab = "Sites", ylab = "Height")
dev.off()

#get dendrogram clusters as a table
hc_clusters <- cutree(hc_km, k = k)
variable_clusters_hc <- split(names(hc_clusters), hc_clusters)
variable_clusters_hc

#END OF SCRIPT
