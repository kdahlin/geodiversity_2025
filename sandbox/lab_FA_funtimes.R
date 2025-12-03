#LAB CFA FUNTIMES

library(corrplot)
library(factoextra)
library(FactoMineR)
library(lavaan)
library(semPlot)
library(psych)
library(Hmisc)

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

#write.csv(cors, "./sandbox/results/correlation_matrix_1m.csv", row.names = TRUE)

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

#transform all variables to a z-score
small.data.1m.z <- as.data.frame(scale(small.data.1m))

# look at column variances (to remove columns with no variance)
variances <- sapply(small.data.1m.z, var)
which.var <- which(variances == 0)

# remove columns with no variance
in.data.1m.z <- in.data.1m.z[,-(which.var)]

# look at correlations
cors <- cor(small.data.1m.z)

# visualize with corrplot
corrplot(cors, order = "FPC", method = "circle")

#Run Factor Analysis (using psych package) -----


#Minres factoring with varimax rotation
fa.parallel(small.data.1m.z, fa = "fa", fm ="minres", n.iter = 100, show.legend = FALSE)
fa.results.minres <- fa(small.data.1m.z, nfactors = 5, rotate = "varimax", fm = "minres")
print(fa.results.minres)


# FA loadings from your 5-factor minres solution
fa_loadings <- as.data.frame(unclass(fa.results.minres$loadings))[, 1:5]

# Make sure rownames match variable names
fa_loadings <- fa_loadings[rownames(pca_loadings), ]

# Compare loadings
comparison <- data.frame(Variable = rownames(pca_loadings),
                         PCA_PC1 = pca_loadings$PC1,
                         FA_Factor1 = fa_loadings$MR1,
                         PCA_PC2 = pca_loadings$PC2,
                         FA_Factor2 = fa_loadings$MR2,
                         PCA_PC3 = pca_loadings$PC3,
                         FA_Factor3 = fa_loadings$MR3,
                         PCA_PC4 = pca_loadings$PC4,
                         FA_Factor4 = fa_loadings$MR4,
                         PCA_PC5 = pca_loadings$PC5,
                         FA_Factor5 = fa_loadings$MR5)
head(comparison)


# Variable clustering using spearmans rank -----
vc <- varclus(~ ., data = small.data.1m.z , similarity = "spearman") 

#Plot dendrogram (axes show 1 - R^2 with own cluster)
plot(vc, main = "Hmisc::varclus variable clustering")

#Extract the hclust object and cut the tree
hc <- vc$hclust

#change hc to data frame for easier viewing
hc_df <- as.data.frame(hc$labels)
hc_df$height <- hc$height
hc_df$order <- seq(1, nrow(hc_df))
rownames(hc_df) <- hc$labels

# Example: cut at height corresponding roughly to |r| >= 0.9
# varclus uses 1 - R^2 as height, so R^2 = 0.9^2 = 0.81 -> height ~ 1 - 0.81 = 0.19
clust_vc <- cutree(hc, k = 5)
clusters_vc <- split(names(clust_vc), clust_vc)
clusters_vc

# varclus uses 1 - R^2 as height, so R^2 = 0.9^2 = 0.81 -> height ~ 1 - 0.81 = 0.19
clust_vc_7 <- cutree(hc, k = 7)
clusters_vc <- split(names(clust_vc_7), clust_vc_7)
clusters_vc

# kmeans clustering of sites ----
set.seed(123) # for reproducibility
kmeans_result <- kmeans(small.data.1m.z, centers = 11, iter.max = 100, nstart = 25)
kmeans_result$centers
# Assign variables to clusters based on highest loading
variable_clusters_km <- split(names(kmeans_result$cluster), kmeans_result$cluster)
variable_clusters_km

# Visualize k-means clustering results
png(filename = "./sandbox/results/kmeans_clustering_1m_PC1_2.png", width = 10, height = 8, units = "in", res = 300)
fviz_cluster(kmeans_result, data = small.data.1m.z,
             palette = "jco",
             title = "K-means Clustering of sites (PC 1,2) (k=11, 100 iterations, 25 starts)",
             ggtheme = theme_minimal())
dev.off()

png(filename = "./sandbox/results/kmeans_clustering_1m_PC1_3.png", width = 10, height = 8, units = "in", res = 300)
fviz_cluster(kmeans_result, data = small.data.1m.z,
             palette = "jco",
             title = "K-means Clustering of sites (PC 1,3) (k=11, 100 iterations, 25 starts)",
             axes = c(1, 3),
             ggtheme = theme_minimal())
dev.off()

png(filename = "./sandbox/results/kmeans_clustering_1m_PC2_3.png", width = 10, height = 8, units = "in", res = 300)
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

# Add cluster labels from your k-means
pc_scores$cluster <- factor(kmeans_result$cluster)

library(plotly)
x11()
plot_ly(pc_scores,
        x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~cluster,
        colors = "Set3",
        type = "scatter3d",
        mode = "markers") %>%
  layout(title = "K-means clusters on first 3 principal components")


# Ward's D minres hierarchical clustering
dissimilarity <- dist(small.data.1m.z, method = "euclidean")
hc_km <- hclust(dissimilarity, method = "ward.D2")
png(filename = "./sandbox/results/hierarchical_clustering_1m_wardsD2.png", width = 20, height = 10, units = "in", res = 300)
plot(hc_km, main = "Dendrogram of Hierarchical Clustering (Ward's Method)", xlab = "Sites", ylab = "Height")
dev.off()

# cut into k clusters and color branches
k <- 11
dend <- as.dendrogram(hc_km)
dend_col <- color_branches(dend, k = k)  # or specify col = c("red","blue",...)
# optionally color labels to match branches
dend_col <- color_labels(dend_col, k = k)

png("./sandbox/results/hierarchical_clustering_1m_wardsD2_colored.png",
    width = 20, height = 10, units = "in", res = 300)
plot(dend_col, main = "Ward D2 dendrogram with colored clusters",
     xlab = "Sites", ylab = "Height")
dev.off()

