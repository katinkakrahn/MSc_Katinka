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
  geom_point(mapping = aes(x = nr_CF2, y = K_F, color = biochar)) + 
  geom_smooth(mapping = aes(x = nr_CF2, y = K_F, color = biochar), formula = y ~ x, method=lm, se = TRUE, fullrange = TRUE) + 
  labs(x = expression(number~of~CF[2]~moieties), y = expression(K[F]), col = "Biochar", title = "Chain length vs KF") + 
  stat_regline_equation(
    aes(x = nr_CF2 -1, y = K_F, color = biochar, label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x) +
  facet_wrap(~ biochar) +
  theme_bw() +
  guides(size = "none", fill = "none")
ChainLength_KF
