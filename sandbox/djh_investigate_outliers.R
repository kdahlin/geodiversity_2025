library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

setwd("C:/Users/hallerdi/Documents/Dillon_Classes/geodiversity_2025") 

Metrics <- read_csv("sandbox/results/all_metrics_wide.csv")
n <- Metrics$func
tpose <- as.data.frame(t(Metrics[,-1]))
colnames(tpose) <- n

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

#Make new dataframe which merely holds a 1 or a zero indicating whether the value is an outlier
flags <- as.data.frame(
  Map(function(vals, lo, hi) {
    as.integer(vals < lo | vals > hi)
  }, tpose, lower_threshold, upper_threshold)
)
rownames(flags) <- rownames(tpose)
flags

#Save out csvs
write.csv(tpose, "sandbox/results/all_metrics_transposed.csv", row.names = TRUE)
write.csv(flags, "sandbox/results/outlier_flags.csv", row.names = TRUE)
