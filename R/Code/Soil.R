# Library ----
library(data.table)
library(readxl)
library(latex2exp)
library(ggtext)
library(scales)
library(writexl)
library(tidyverse)

Soil <- read_excel("R/data_raw/040322_soil_rawdata.xlsx")


CECions <- filter(Soil, Type == "Exchangeable ions")
Total_element <- filter(Soil, Type == "Total element conc")

# Element concentrations from triplicate soil samples
SoilsummaryCEC <- CECions %>% 
  group_by(Compound) %>% 
  summarise(mean_conc = mean(Concentration), 
            sd_conc = sd(Concentration))


SoilsummaryTotal <- Total_element %>% 
  group_by(Compound) %>% 
  summarise(mean_conc = mean(Concentration),
            sd_conc = sd(Concentration)
  )

#CEC ions plot
CECions_plot <- ggplot(data = SoilsummaryCEC) + 
  geom_point(mapping = aes(x = reorder(Compound, mean_conc), y = mean_conc, 
                               group = Compound), size = 4) + 
  labs(x = "Exchangeable ion", y = "Total concentration (meqv/100g)") + 
  theme_bw() +
  scale_y_log10() +
  guides(size = "none")
CECions_plot
ggsave(filename = "R/figs/CECions_plot.pdf")

Totalelement_plot <- ggplot(data = SoilsummaryTotal) + 
  geom_point(aes(x = reorder(Compound, mean_conc), y = mean_conc, 
                 group = Compound), size = 4) + 
  labs(x = "Element", y = "Total concentration (mg/kg dw)") + 
  scale_y_log10() +
  theme_bw() +
  guides(size = "none")
Totalelement_plot
ggsave(filename = "R/figs/totalelement_plot.pdf")

