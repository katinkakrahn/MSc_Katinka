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
library(scales)
library(jcolors)

Biochar <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/250322_biochar_parameters.xlsx")
as.data.table(Biochar)
Biochar <- as.data.table(Biochar)
setnames(Biochar, "Biochar", "biochar")



#Hans Peter: choose a low concentration Kd value to get from KF expression to make correlation plots on


Ca <- subset(Biochar, Parameter == "Ca")
Ca_biochar <- merge(Ca, summary_stats_single, by = "biochar")

Ca_biochar_corr <- ggplot(data = Ca_biochar, aes(x = log(Mean), y = K_F), group = compound) +
  geom_point(size = 2) + 
  geom_smooth(formula = y ~ x, method=lm, 
              se = FALSE, fullrange = FALSE) +
  labs(x = expression("log [Ca] g/kg"), y = expression(log~K[F])) +
  facet_wrap(~ compound) +
  # geom_richtext(
  #   data = summary_stats_CLandKF_label,
  #   aes(label = label, x = nr_CF2, y = K_F),
  #   hjust = 0
  # ) +
  theme_bw() +
  theme(panel.grid = element_blank()) #+
  # guides(size = "none", fill = "none")
Ca_biochar_corr
ggsave(filename = "R/figs/Ca_biochar_corr.pdf")

Carbon <- subset(Biochar, Parameter == "C")
Carbon_biochar <- merge(Carbon, summary_stats_single, by = "biochar")

Carbon_biochar_corr <- ggplot(data = Carbon_biochar, aes(x = Mean, y = K_F), group = compound) +
  geom_point(size = 2) + 
  geom_smooth(formula = y ~ x, method=lm, 
              se = FALSE, fullrange = FALSE) +
  labs(x = expression("% total C"), y = expression(log~K[F])) +
  facet_wrap(~ compound) +
  # geom_richtext(
  #   data = summary_stats_CLandKF_label,
  #   aes(label = label, x = nr_CF2, y = K_F),
  #   hjust = 0
  # ) +
  theme_bw() +
  theme(panel.grid = element_blank()) #+
# guides(size = "none", fill = "none")
Carbon_biochar_corr
ggsave(filename = "R/figs/Ca_biochar_corr.pdf")
