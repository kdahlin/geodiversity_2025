#Script to plot metrics against increasing scale
library(ggplot2)
library(stringr) #Allows python-style string formatting
library(dplyr)

# read in the long data
in.data <- read.csv("./results/all_metrics_long_clean.csv")

remove.metrics <- c("Sa", "Sq", "s10z", "Sdq6", "Sph", "Svi", "Sv", "Scl1",
                    "TRIv1", "TRIriley", "TRIrmsd", "rough",
                    "Std1", "TPIv2", "TRIv2", "Srw1", "Scl2")

check <- !(in.data$func %in% remove.metrics)
Metrics.sub <- filter(in.data, check)


grouped.data <- Metrics.sub %>% group_by(func, site, res) %>% 
  summarise(mean = mean(value, na.rm = TRUE)) %>%
  mutate(res = as.numeric(res),
         func = factor(func,
                       levels = c("range", "MAD", "Smean", "Spk", "Svk",
                                  "Slope", "nor", "east", "SAR", "Sdq", "Sdr",
                                  "nmodes", "SBI", "SCI", "Ssk", "Sku", "Sds", 
                                  "Sdc", "CuPro", "CuPl", "CuTot", "Ssc", 
                                  "TPIv1", "stdv", "VRM", "rEnt", "Sk", "Sfd", 
                                  "Srw2", "Srw3", "Std2", "Stxr1", "Stxr2"))
         )


ggplot(grouped.data, aes(x = res, y = mean, color = site, group=site)) +
  geom_point() +
  geom_line() +
  #geom_errorbar(aes(x = res, ymin = min, ymax = max), width = .2)+
  facet_wrap(~ func, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = sort(unique(grouped.data$res)))+
  theme_bw() +
  scale_color_manual(values = c("CPER" = "#a6cee3", "OAES" = "#a6cee3", "CLBJ" = "#a6cee3", 
                                "ORNL" = "#b2df8a", "MLBS" = "#b2df8a", 
                                "RMNP"="#33a02c", "TEAK"="#33a02c", "WREF"="#33a02c", 
                                "WOOD"="#1f78b4", "OSBS"="#1f78b4", "UNDE"="#1f78b4"))+
  labs(x = "", y = "", title = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 5)
  ) +
  theme(legend.position = "bottom", 
        legend.key.height = unit(0.1, "in"))
ggsave("./results/all_metrics_resolution.tif", width = 6.5 , height = 9.5, 
       units = "in", dpi = 300)

# look at VRM only for discussion section writing
vrm.only <- subset(grouped.data, grouped.data$func == "VRM")
x11()
ggplot(vrm.only, aes(x = res, y = mean, color = site, group=site)) +
  geom_point() +
  geom_line() +
  #geom_errorbar(aes(x = res, ymin = min, ymax = max), width = .2)+
  scale_x_continuous(breaks = sort(unique(grouped.data$res)))+
  theme_bw() +
  scale_color_manual(values = c("CPER" = "#a6cee3", "OAES" = "#a6cee3", "CLBJ" = "#a6cee3", 
                                "ORNL" = "#b2df8a", "MLBS" = "#b2df8a", 
                                "RMNP"="#33a02c", "TEAK"="#33a02c", "WREF"="#33a02c", 
                                "WOOD"="#1f78b4", "OSBS"="#1f78b4", "UNDE"="#1f78b4"))+
  labs(x = "", y = "", title = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 5)
  ) +
  theme(legend.position = "bottom", 
        legend.key.height = unit(0.1, "in"))
