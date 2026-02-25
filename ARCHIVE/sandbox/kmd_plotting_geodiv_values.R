## looking at metrics

library(readxl)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

in.file <- read_excel("./analysis/table_o_metrics.xlsx", sheet = 1)

simple.df <- in.file[c(5:32), c(1,2,9:12)]

simple.df$CPER <- as.numeric(simple.df$CPER)
simple.df$ORNL <- as.numeric(simple.df$ORNL)
simple.df$RMNP <- as.numeric(simple.df$RMNP)
simple.df$WOOD <- as.numeric(simple.df$WOOD)


df.piv <- pivot_longer(data = simple.df,
                       cols = c("CPER", "ORNL", "RMNP", "WOOD"), 
                       values_to = "Values")

ggplot(df.piv, aes(x = 1, y = Values, color = name)) +
  geom_beeswarm(size = 2) +
  facet_wrap(~`Full name`, scales = "free_y") + 
  xlim(0.995, 1.005)
