#LAB CFA FUNTIMES

library(corrplot)
library(factoextra)
library(FactoMineR)
library(lavaan)
library(semPlot)
library(psych)

options(digits = 10)

#Prepare the data -----

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

#transform all variables to a z-score
in.data.1m.z <- as.data.frame(scale(in.data.1m))

# look at column variances (to remove columns with no variance)
variances <- sapply(in.data.1m.z, var)
which.var <- which(variances == 0)

# remove columns with no variance
in.data.1m.z <- in.data.1m.z[,-(which.var)]

# look at correlations
cors <- cor(in.data.1m.z)

# visualize with corrplot
corrplot(cors, order = "FPC", method = "circle")

#check VIF and tolerance to look for multicollinearity
library(car)
vif_values <- vif(lm(as.matrix(in.data.1m.z) ~ 1))
tolerance_values <- 1 / vif_values

#Run Factor Analysis (using psych package) -----

#Principal axis factoring with varimax rotation. 
#This method doesn't assume a normal distribution.

fa.parallel(in.data.1m.z, fa = "fa", fm ="pa", n.iter = 100, show.legend = FALSE)
fa.results <- fa(in.data.1m.z, nfactors = 5, rotate = "varimax", fm = "pa")
print(fa.results)

#Visualize factor loadings
fa.diagram(fa.results)

#Minres factoring with varimax rotation
fa.parallel(in.data.1m.z, fa = "fa", fm ="minres", n.iter = 100, show.legend = FALSE)
fa.results.minres <- fa(in.data.1m.z, nfactors = 5, rotate = "varimax", fm = "minres")
print(fa.results.minres)

#need to remove some variables to make the matrix invertible for SEM 
#remove the variables with low loadings (<.3) across all factors

#get the loadings matrix
loadings <- as.data.frame(fa.results.minres$loadings[1:ncol(in.data.1m.z), 1:5])

#find variables with all loadings <|.3|
low.loadings <- rownames(loadings[which(abs(loadings$MR1) < 0.3 & 
                                           abs(loadings$MR2) < 0.3 & 
                                           abs(loadings$MR3) < 0.3 &
                                          abs(loadings$MR4) < 0.3 &
                                          abs(loadings$MR5) < 0.3), ])
#remove those variables from the data
in.data.1m.z.reduced <- in.data.1m.z[, !(names(in.data.1m.z) %in% low.loadings)]

#Remove cross-loading variables (>|.4| on more than one factor)
cross.loadings <- rownames(loadings[which((abs(loadings$MR1) > 0.4 & abs(loadings$MR2) > 0.4) |
                                            (abs(loadings$MR1) > 0.4 & abs(loadings$MR3) > 0.4) |
                                            (abs(loadings$MR2) > 0.4 & abs(loadings$MR3) > 0.4) |
                                            (abs(loadings$MR1) > 0.4 & abs(loadings$MR4) > 0.4) |
                                            (abs(loadings$MR2) > 0.4 & abs(loadings$MR4) > 0.4) |
                                            (abs(loadings$MR3) > 0.4 & abs(loadings$MR4) > 0.4) |
                                            (abs(loadings$MR1) > 0.4 & abs(loadings$MR5) > 0.4) |
                                            (abs(loadings$MR2) > 0.4 & abs(loadings$MR5) > 0.4) |
                                            (abs(loadings$MR3) > 0.4 & abs(loadings$MR5) > 0.4) |
                                            (abs(loadings$MR4) > 0.4 & abs(loadings$MR5) > 0.4)
                                             ), ])
#remove those variables from the data
in.data.1m.z.reduced <- in.data.1m.z.reduced[, !(names(in.data.1m.z.reduced) %in% cross.loadings)]

# rerun FA on reduced dataset
fa.parallel(in.data.1m.z.reduced, fa = "fa", fm ="minres", n.iter = 100, show.legend = FALSE)
fa.results.reduced <- fa(in.data.1m.z.reduced, nfactors = 5, rotate = "varimax", fm = "minres")
print(fa.results.reduced)

#MORE var REDUCTION! -----
library(Hmisc)

corrplot(cor(in.data.1m.z.reduced), order = "FPC", method = "circle")

# 1. Variable clustering (by default uses squared Pearson correlations)
vc <- varclus(~ ., data = in.data.1m.z.reduced, similarity = "spearman") 

# 2. Plot dendrogram (axes show 1 - R^2 with own cluster)
plot(vc, main = "Hmisc::varclus variable clustering")

# 3. Extract the hclust object and cut the tree
hc <- vc$hclust

# Example: cut at height corresponding roughly to |r| >= 0.9
# varclus uses 1 - R^2 as height, so R^2 = 0.9^2 = 0.81 -> height ~ 1 - 0.81 = 0.19
clust_vc <- cutree(hc, k = 5)
clusters_vc <- split(names(clust_vc), clust_vc)
clusters_vc

