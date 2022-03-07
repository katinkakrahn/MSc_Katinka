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

Soil <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/040322_soil_rawdata.xlsx")
as.data.table(Soil)
Soil <- as.data.table(Soil)


CECions <- filter(Soil, Type == "Exchangeable ions")
Total_element <- filter(Soil, Type == "Total element conc")

# Element concentrations from triplicate soil samples
SoilsummaryCEC <- CECions[, .(mean_conc = mean(Concentration), 
                            sd_conc = sd(Concentration)
),
keyby = .(Compound)]

SoilsummaryTotal <- Total_element[, .(mean_conc = mean(Concentration),
                                      sd_conc = sd(Concentration)
                                      ),
                                  keyby = .(Compound)]

#CEC ions plot
CECions_plot <- ggplot(data = SoilsummaryCEC) + 
  geom_pointrange(mapping = aes(x = reorder(Compound, mean_conc), y = mean_conc, 
                                ymin = mean_conc-sd_conc, ymax = mean_conc+sd_conc, group = Compound)) + 
  labs(x = "Exchangeable ion", y = "Total concentration (meqv/100g)") + 
  theme_bw() +
  guides(size = "none")
CECions_plot
ggsave(filename = "figs/CECions_plot.png")

Totalelement_plot <- ggplot(data = SoilsummaryTotal) + 
  geom_pointrange(aes(x = reorder(Compound, mean_conc),y= mean_conc, ymin=mean_conc-sd_conc, 
                      ymax=mean_conc+sd_conc, group = Compound)) + 
  labs(x = "Element", y = "Total concentration (mg/kg dw)") + 
  scale_y_log10() +
  theme_bw() +
  guides(size = "none")
Totalelement_plot
ggsave(filename = "figs/totalelement_plot.png")
