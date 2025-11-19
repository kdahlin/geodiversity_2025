#LAB CFA FUNTIMES

library(corrplot)
library(factoextra)
library(FactoMineR)
library(lavaan)
library(semPlot)
library(psych)

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

# remove columns with no variance
which.var <- which(variances == 0)
in.data.1m <- in.data.1m[,-(which.var)]

#scale data
scale.in.data.1m <- scale(in.data.1m)

#run efa on in.data.1m
fa.parallel(scale.in.data.1m, fa = "fa", fm = "minres")
fa.parallel(scale.in.data.1m, fa = "fa", fm = "pa")
fa.parallel(scale.in.data.1m, fa = "fa", fm = "ols")
fa.parallel(scale.in.data.1m, fa = "fa", fm = "gls")
#fa.parallel(scale.in.data.1m, fa = "fa", fm = "ml")
#fa.parallel(scale.in.data.1m, fa = "fa", fm = "wls")

fit_minres <- fa(scale.in.data.1m, nfactors = 3, rotate = "varimax", fm = "minres")
print(fit_minres, cut = 0.3, sort = TRUE)
fa.diagram(fit_minres)
loadings_minres <- as.data.frame(fit_minres$loadings)

fit_pa <- fa(scale.in.data.1m, nfactors = 4, rotate = "varimax", fm = "pa")
print(fit_pa, cut = 0.3, sort = TRUE)
fa.diagram(fit_pa)
loadings_pa <- as.data.frame(fit_pa$loadings)

fit_ols <- fa(scale.in.data.1m, nfactors = 4, rotate = "varimax", fm = "ols")
print(fit_ols, cut = 0.3, sort = TRUE)
fa.diagram(fit_ols)
loadings_ols <- as.data.frame(fit_ols$loadings)





#compare fit1 and fit2 factor groupings


# look at correlations
cors <- cor(in.data.1m)
# visualize with corrplot
corrplot(cors, order = "FPC", method = "circle")

### OPTIONAL - REMOVE HIGHLY CORRELATED VARIABLES AND WEIRD ONES
remove.some <- c("stdv_1", "sq_1", "mad_1", "tri_1", "TRIriley_1", "TRIrmsd_1",
                 "raster.entropy_1", "std_1", "sdc_1", "sdq6_1", "smean_1", 
                 "tpi_1")

remove.cols <- which(names(in.data.1m) %in% remove.some)

small.data.1m <- in.data.1m[,-c(remove.cols)]

scale.small.data.1m <- scale(small.data.1m)

fa.parallel(scale.small.data.1m, fa = "fa")

fit2 <- fa(scale.small.data.1m, nfactors = 4, rotate = "varimax", fm = "pa")
print(fit2, cut = 0.3, sort = TRUE)
fa.diagram(fit2)

#okay this is really wacky, gotta reduce the multicollinearity (HACK AND SLASH METHOD)

# Set correlation threshold
corr_thresh <- 0.9

# Create the correlation matrix of your data
cor_matrix <- cor(small.data.1m, use = "pairwise.complete.obs")

# Find highly correlated pairs (upper triangle only to avoid duplicates)
high_corr_pairs <- which(abs(cor_matrix) > corr_thresh, arr.ind = TRUE)
high_corr_pairs <- high_corr_pairs[high_corr_pairs[,1] < high_corr_pairs[,2], ]

# Identify columns to remove
vars_to_remove <- unique(colnames(small.data.1m)[high_corr_pairs[,2]])

# Remove highly correlated variables
reduced.data <- small.data.1m[, !(colnames(small.data.1m) %in% vars_to_remove)]

# Check result
dim(reduced.data)
colnames(reduced.data)
corrplot(cor(reduced.data), order = "FPC", method = "circle")

scale.reduced.data <- scale(reduced.data)

# --------------------------------
# Confirmatory Factor Analysis
# --------------------------------

#Specify the CFA model with the metric groupings we agreed on in the last meeting
#(from Table 1 in maintext)
geodiv_model <- '
  Amplitude =~ s10z_1 + spk_1 + svk_1 
  Slope =~ northness_1 + eastness_1 + sdq_1 + sdr_1 
  Distribution =~ sbi_1 + ssk_1 + sku_1 + sds_1 + nmodes_1 
  Roughness =~ TPI_1 + sa_1 + curvature_1
  Angular_Rad =~ std_2 + stxr_1 + srw_2 + srw_3 
  
  #regressions among factors 
  Amplitude ~~ Slope + Distribution + Roughness + Angular_Rad
  Slope ~~ Distribution + Roughness + Angular_Rad
  Distribution ~~ Roughness + Angular_Rad
  Roughness ~~ Angular_Rad
  '
# Fit the CFA model
fit <- cfa(geodiv_model, data = scale.reduced.data)
summary(fit, fit.measures = TRUE, standardized = TRUE)
semPaths(fit, "std", layout = "tree", whatLabels = "std", 
         edge.label.cex = 0.8, sizeMan = 5, sizeLat = 7)

#a very bad model fit. try efa ?

KMO(scale.reduced.data)
cortest.bartlett(scale.reduced.data)

fa.parallel(scale.reduced.data, fa = "fa")
efa_results_4 <- fa(scale.reduced.data, nfactors = 4, rotate = "varimax", fm = "ml")
print(efa_results_4, cut = 0.3, sort = TRUE)
fa.diagram(efa_results_4)

efa_results_3 <- fa(scale.reduced.data, nfactors = 3, rotate = "varimax", fm = "ml")
print(efa_results_3, cut = 0.3, sort = TRUE)
fa.diagram(efa_results_3)

#use EFA 4 factor results to specify a new CFA model
geodiv_model_efa <- '
  Factor1 =~ s10z_1 + sdr_1 + sa_1 + ssc_1 + svk_1 + std_2 + curvature_1 + northness_1
  Factor2 =~ srw_2 + srw_3 + sbi_1 + sds_1 + nmodes_1
  Factor3 =~ ssk_1 + spk_1 + sdq_1 + sku_1
  Factor4 =~ eastness_1 + TPI_1 

  '
# Fit the CFA model
fit_cfa <- cfa(geodiv_model_efa, data = scale.reduced.data)
summary(fit_cfa, fit.measures = TRUE, standardized = TRUE)
semPaths(fit_cfa, "std", layout = "tree", whatLabels = "std", 
         edge.label.cex = 0.8, sizeMan = 5, sizeLat = 7)
