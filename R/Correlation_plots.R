library(data.table)
library(ggplot2)
library(psych)
library(readxl)
library(grDevices)
library(dplyr)
library(broom)
library(ggpubr)
library(tidyverse)
library(knitr)

ChainLength_KF <- ggplot(data = summary_stats_single) +
  geom_point(mapping = aes(x = nr_CF2, y = K_F, group = biochar), color = "grey45", size = 1) + 
  geom_smooth(mapping = aes(x = nr_CF2, y = K_F, group = biochar), color = "black", formula = y ~ x, method=lm, se = TRUE, fullrange = TRUE) +
  geom_errorbar(aes(nr_CF2, ymin=K_F-K_F_std_error, ymax=K_F+K_F_std_error), width=.2,
                position=position_dodge(.9), color = "grey45") +
  labs(x = expression(number~of~CF[2]~moieties), y = expression(log~K[F])) +
  facet_wrap(~ biochar) +
  geom_label(data = summary_stats_single, size = 2, inherit.aes = T, 
    aes(x = 6.5, y = 2, label = paste("slope =",round(n, digits = 2),","," ","R^2",round(r_squared, digits = 2)))) +
  theme_bw() +
  guides(size = "none", fill = "none")
ChainLength_KF
ggsave(filename = "figs/chainlength_KF.png")
