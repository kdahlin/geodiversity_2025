#LAB CFA and CLUSTER ANALYSIS FUNTIMES

library(corrplot)
library(factoextra)
library(FactoMineR)
library(lavaan)
library(semPlot)
library(psych)
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

#write.csv(cors, "./sandbox/results/correlation_matrix_1m.csv", row.names = TRUE)

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
plot(vc_spearmans, main = "Hmisc::varclus variable clustering")

#Using pearson and Hmisc::varclus
vc_pearson <- varclus(~ ., data = small.data.1m.z , similarity = "pearson")

#Plot dendrogram (axes show 1 - R^2 with own cluster)
plot(vc_pearson, main = "Hmisc::varclus variable clustering (Pearson)")


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

# K-means clustering of sites ----
set.seed(123) # for reproducibility
png(filename = "./sandbox/results/kmeans_screeplot_1_20_centers.png", 
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
write.csv(variable_clusters_km_df, "./sandbox/results/kmeans_clusters_k11_table.csv",
          row.names = FALSE)

# Visualize k-means clustering results of Principal Components
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

png("./sandbox/results/hierarchical_clustering_1m_wardsD2_colored.png",
    width = 20, height = 10, units = "in", res = 300)
plot(dend_col, main = "Ward's D minimum variance (k=11)",
     xlab = "Sites", ylab = "Height")
dev.off()

#get dendrogram clusters as a table
hc_clusters <- cutree(hc_km, k = k)
variable_clusters_hc <- split(names(hc_clusters), hc_clusters)
variable_clusters_hc

#Run Factor Analysis (using psych package) -----

#Minres factoring with varimax rotation
fa.parallel(small.data.1m.z, fa = "fa", fm ="minres", n.iter = 100, show.legend = FALSE)
fa.results.minres <- fa(small.data.1m.z, nfactors = 5, rotate = "varimax", fm = "minres")
print(fa.results.minres)
fa.diagram(fa.results.minres, main = "Factor Analysis (Minres, Varimax Rotation)")

#get factors and loadings as a table
fa_loadings <- as.data.frame(fa.results.minres$loadings[1:ncol(small.data.1m.z),])
fa_loadings <- fa_loadings[, c(ncol(fa_loadings), 1:(ncol(fa_loadings)-1))]

#nmodes, sku, and eastness don't load well on any factors - remove from SEM specification
vars_for_sem <- c("sdq_1", "slope_1", "sar_1", "range_1", "stdv_1",
                  "mad_1", "sk_1", "spk_1", "sdc_1", "vrm_1", "svk_1",
                  "srw_3", "sdr_1", "std_2", "ssc_1", "sfd_1", "sds_1",
                  "srw_2", "curvature_planform_1", "curvature_profile_1",
                  "curvature_total_1", "raster.entropy_1", "sci_1",
                  "sv_1", "smean_1", "sbi_1", "northness_1",
                  "stxr_2", "stxr_1", "TPI_1", "ssk_1")
sem_data <- small.data.1m.z[, vars_for_sem]

library(Matrix)
S <- cov(sem_data, use = "pairwise.complete.obs")
S_pd <- as.matrix(nearPD(S, keepDiag = TRUE)$mat)

#specify SEM Using MRs from FA
sem_model <- '
MR1 =~ sdq_1 + slope_1 + sar_1 + range_1 + stdv_1 + mad_1 + sk_1 + spk_1 + sdc_1 + vrm_1 + svk_1
MR2 =~ srw_3 + sdr_1 + std_2 + ssc_1 + sfd_1 + sds_1 + srw_2
MR3 =~ curvature_planform_1 + curvature_profile_1 + curvature_total_1 + raster.entropy_1 
MR4 =~ sci_1 + sv_1 + smean_1 + sbi_1 + northness_1 
MR5 =~ stxr_2 + stxr_1 + TPI_1 + ssk_1
'
fit <- sem(sem_model, sample.cov = S_pd,
           sample.nobs = nrow(sem_data),
           meanstructure = FALSE)
summary(fit, fit.measures = TRUE, standardized = TRUE)

est <- parameterEstimates(fit)
est[est$op == "~~" & est$lhs == est$rhs & est$est < 0, ]

semPaths(
  fit,
  whatLabels = "std",
  layout = "tree",
  edge.label.cex = 0.8,
  sizeMan = 5,
  sizeLat = 7,
  nCharNodes = 0,
  title = FALSE
)

sem_model_reduced <- '
  MR1 =~ sdq_1 + slope_1 + sar_1 + range_1 + stdv_1 + mad_1 + sk_1 + spk_1 + sdc_1 + vrm_1 + svk_1
  MR2 =~ srw_3 + sdr_1 + std_2 + ssc_1 + sfd_1 + sds_1 + srw_2
  MR3 =~ curvature_planform_1 + curvature_profile_1 + raster.entropy_1
  MR4 =~ sci_1 + sv_1 + smean_1 + sbi_1 + northness_1
  MR5 =~ stxr_1 + TPI_1 + ssk_1
'
fit2 <- sem(sem_model_reduced, data = sem_data)  # or with S_pd if still needed
summary(fit2, fit.measures = TRUE, standardized = TRUE)
est2 <- parameterEstimates(fit2)
est2[est2$op == "~~" & est2$lhs == est2 $rhs & est2$est < 0, ]
semPaths(
  fit2,
  whatLabels = "std",
  layout = "tree",
  edge.label.cex = 0.8,
  sizeMan = 5,
  sizeLat = 7,
  nCharNodes = 0,
  title = FALSE
)