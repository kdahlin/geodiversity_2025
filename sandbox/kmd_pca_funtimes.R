
library(corrplot)
library(factoextra)

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

for (i in 1:length(rnames)) {
  rnames.new <- rbind(rnames.new, 
                      paste0(split.rnames[[i]][1], "_", split.rnames[[i]][6]))
}

# attach the better names
row.names(in.data.1m) <- as.data.frame(rnames.new[,1])

# look at column summaries
summary(in.data.1m)

# columns with no variance
no.var <- c("srw_1", "scl_2")

# remove columns with no variance
no.var.cols <- which((names(in.data.1m) %in% no.var))
in.data.1m <- in.data.1m[,-(no.var.cols)]

# look at correlations
cors <- cor(in.data.1m)

# visualize with corrplot
x11()
corrplot(cors, order = "FPC", method = "circle")

# do pca
pca <- prcomp(in.data.1m, center = TRUE, scale = TRUE)

# look at variances
plot(pca)

# look at variance explained
summary(pca)

# get the PCA 









