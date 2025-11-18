######Plotting######################
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
Metrics <- read_csv("./results/all_metrics_long.csv")
View(Metrics)

# subset to just 1m data for plotting
metrics.1m <- subset(Metrics, Metrics$res == 1)

#Boxplot for all
ggplot(metrics.1m, aes(x = site, y = value, color = site)) +
  geom_boxplot(outlier.colour = "black", outlier.size = 1) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.8) +
  facet_wrap(~ func, scales = "free_y") +  # one panel per metric
  scale_color_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3"))+
  theme_minimal() +
  labs(
    x = "Raster ID",
    y = "Value",
    title = "Metrics Across 4 different surfaces")

### pulling out just a few of the metrics (based on PCA)
pca.imp <- metrics.1m %>%
  filter(func %in% c("range_1", "sbi_1", "ssk_1", "TPI_1", "stxr_1", "sku_1",
                     "raster.entropy_1", "std_1", "eastness_1", "nmodes_1"))

#Beeswarm plot for subset
library(ggbeeswarm)
ggplot(pca.imp, aes(x = site, y = value, color = site)) +
  geom_beeswarm() +
  facet_wrap(~ func, scales = "free_y") +  # one panel per metric
  theme_bw() +
  scale_color_manual(values = c("CPER" = "aquamarine4", "ORNL" = "royalblue4", "RMNP"="purple3", "WOOD"="darkorange3"))+
  labs(
    x = "Raster ID",
    y = "Value",
    title = "Metrics Across 4 different surfaces")



