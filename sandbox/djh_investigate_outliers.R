library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
Metrics <- read_csv("sandbox/results/all_metrics_wide.csv")
View(Metrics)
n <- Metrics$func
tpose <- as.data.frame(t(Metrics[,-1]))
colnames(tpose) <- n
View(tpose)

#Calculate IQR per column
IQRs <- sapply(tpose, IQR)
Q1s <- sapply(tpose, quantile, 0.25)
Q3s <- sapply(tpose, quantile, 0.75)

#Outlier threshold = Q3 + 1.5*IQR
upper_threshold <- Q3s + 1.5*IQRs
lower_threshold <- Q1s - 1.5*IQRs

#Get names of outliers
out_of_bounds <- Map(function(vals, lo, hi) {
  rownames(tpose)[vals < lo | vals > hi]
}, tpose, lower_threshold, upper_threshold)
out_of_bounds

counts <- table(unlist(out_of_bounds))
