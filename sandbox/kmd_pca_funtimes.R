
library(corrplot)
library(factoextra)
library(FactoMineR)

options(digits = 10)

# read in the wide data
in.data <- read.csv("./results/all_metrics_wide.csv", row.names = 1)

# transpose data so each row is a subsite
in.data.t <- as.data.frame(t(in.data))

# get row names
samples <- rownames(in.data.t)

# split out resolutions
sample.res <- strsplit(samples, "[_\\.]")
res <- unlist(lapply(sample.res, function(x) {x[9]}))

# get just 1m resolution
in.data.1m <- subset(in.data.t, res == "1m")

# rename rows so they aren't so clunky
rnames <- row.names(in.data.1m)
split.rnames <- strsplit(rnames, "_")

rnames.new <- paste0(split.rnames[[1]][1], "_", split.rnames[[1]][6])

for (i in 2:length(rnames)) {
  rnames.new <- rbind(rnames.new, 
                      paste0(split.rnames[[i]][1], "_", split.rnames[[i]][6]))
}

# attach the better names
row.names(in.data.1m) <- rnames.new[,1]

# look at column variances (to remove columns with no variance)
variances <- sapply(in.data.1m, var)

# get the column names with no variance
which.var <- which(variances == 0)

# remove columns with no variance
in.data.1m <- in.data.1m[,-(which.var)]

# look at correlations
cors <- cor(in.data.1m)

# visualize with corrplot
x11()
corrplot(cors, order = "FPC", method = "circle")

# ### OPTIONAL - REMOVE HIGHLY CORRELATED VARIABLES AND WEIRD ONES
# remove.some <- c("stdv_1", "sq_1", "mad_1", "tri_1", "TRIriley_1", "TRIrmsd_1",
#                  "raster.entropy_1", "std_1", "sdc_1", "sdq6_1", "smean_1", 
#                  "tpi_1")
# 
# remove.cols <- which(names(in.data.1m) %in% remove.some)
# 
# small.data.1m <- in.data.1m[,-c(remove.cols)]
##########

# do pca
pca <- prcomp(in.data.1m, center = TRUE, scale = TRUE)

scale.data.1m <- scale(in.data.1m, center = TRUE, scale = TRUE)

pca.test <- princomp(scale.data.1m)
# note some plotting from here on out using the factoextra package comes from 
# scripts provided here: https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

# look at variances
fviz_eig(pca, addlabels = TRUE)

# look at variance explained
summary(pca)

# get the PCA loadings (aka rotations)
loadings <- pca$rotation

# get the scores/PCs
scores <- as.data.frame(summary(pca)$x)
scores <- pca.test$scores

# biplot
fviz_pca_biplot(pca.test, repel = TRUE,                 
                col.var = "#2e9fdf", # Variables color
                col.ind = "#696969")  # Individuals color

# let's look at how the PCs correlate with the variables
pca.var <- get_pca_var(pca)

x11()
corrplot(pca.var$cos2, is.corr = FALSE)

# I (Kyla) likes this better because I understand it, but weird about PC36!?
pc.corrs <- cor(scores, in.data.1m)
corrplot(pc.corrs)

# let's take a look
plot(scores$PC36, in.data.1m$stxr_1, type = "n")
text(scores$PC36, in.data.1m$stxr_1, labels = row.names(in.data.1m))
# weird. not sure about this.

# let's look at which variables are most strongly correlated with which scores
cors.sorted <- as.data.frame(matrix(NA, nrow = 34, ncol = 34))
cors.Rvals.sorted <- as.data.frame(matrix(NA, nrow = 34, ncol = 34))

for (i in 1:34) {
  names(cors.sorted)[i] <- colnames(pc.corrs)[i]
  get.sort <- sort(abs(pc.corrs[,i]), decreasing = TRUE)
  cors.sorted[,i] <- names(get.sort)
  names(cors.Rvals.sorted)[i] <- paste0(colnames(pc.corrs)[i], "_Rval")
  cors.Rvals.sorted[,i] <- sort(abs(pc.corrs[,i]))
}



