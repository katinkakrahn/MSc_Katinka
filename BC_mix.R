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
library(plotrix)

Sorption_BC_mix_summary <- Sorption_BC_mix[, .(mean_logCw = mean(log_Cw), 
                                               mean_logCs = mean(log_Cs),
                                               se_logCw = std.error(log_Cw),
                                               se_logCs = std.error(log_Cs),
                                               log_Kd = mean(log_Cs/log_Cw),
                                               se_logKd = std.error(log_Cs/log_Cw),
                                               mixLogic = mixLogic
                                               ),
                                           keyby = .(Compound, Biochar)]

Sorption_BC_mix_summary$Compound <- factor(Sorption_BC_mix_summary$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                      "PFOA", "PFNA", "PFDA"))

BC_mix_Kd <- ggplot(data = Sorption_BC_mix_summary, aes(x = Compound, y = log_Kd, color = Biochar)) + 
  geom_point(size = 2)+ 
  geom_errorbar(aes(ymin=log_Kd-se_logKd, ymax=log_Kd+se_logKd), width = .05) + 
  labs(x = "", y = expression(log~K[d]), col = "") + 
  # scale_color_manual(
  #   values = c('black','grey60'),
  #   breaks = c("TRUE", "FALSE"),
  #   labels = c("cocktail", "single compound")
  # ) +
  theme_bw()
BC_mix_Kd
set_palette(BC_mix_Kd, "uchicago")
ggsave(filename="R/figs/BC_mix_Kd.pdf")


#Sorption attenuation
Sorption_BC_single_C10 <- subset(Sorption_BC_single, Conc_point == 10)
Sorption_BC_single_C10_PFPeA <- subset(Sorption_BC_single_C10, Compound == "PFPeA" & Biochar == "CWC")
Sorption_BC_single_C10_PFHxA <- subset(Sorption_BC_single_C10, Compound == "PFHxA"& Biochar == "ULS")
Sorption_BC_single_C10_PFHpA <- subset(Sorption_BC_single_C10, Compound == "PFHpA"& Biochar == "ULS")
Sorption_BC_single_C10_PFOA <- subset(Sorption_BC_single_C10, Compound == "PFOA")
Sorption_BC_single_C10_PFNA <- subset(Sorption_BC_single_C10, Compound == "PFNA")
Sorption_BC_single_C10_PFDA <- subset(Sorption_BC_single_C10, Compound == "PFDA")
Sorption_BC_single_C10_common <- merge(Sorption_BC_single_C10_PFPeA, Sorption_BC_single_C10_PFHxA, all=T)
Sorption_BC_single_C10_common <- merge(Sorption_BC_single_C10_common, Sorption_BC_single_C10_PFHpA, all=T)
Sorption_BC_single_C10_common <- merge(Sorption_BC_single_C10_common, Sorption_BC_single_C10_PFOA, all=T)
Sorption_BC_single_C10_common <- merge(Sorption_BC_single_C10_common, Sorption_BC_single_C10_PFNA, all=T)
Sorption_BC_single_C10_common <- merge(Sorption_BC_single_C10_common, Sorption_BC_single_C10_PFDA, all=T)

Sorption_BC_single_C10_common$Compound <- factor(Sorption_BC_single_C10_common$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))
Sorption_BC_mix_summary$Compound <- factor(Sorption_BC_mix_summary$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                          "PFOA", "PFNA", "PFDA"))

Sorption_attenuation_BC <- ggplot() +
  geom_point(data = Sorption_BC_single_C10_common, mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar), shape = mixLogic, size = 2), 
             ) + 
  geom_point(data = Sorption_BC_mix_summary, mapping = aes(x = mean_logCw, y = mean_logCs, color = factor(Biochar), shape = mixLogic, size = 2), 
             ) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), shape = "Cocktail", color = "Biochar") +
  guides(size = "none") +
  facet_wrap(~Compound) +
  theme_bw() +
  theme(panel.grid = element_blank())
Sorption_attenuation_BC
set_palette(Sorption_attenuation_BC, "uchicago")
ggsave(filename="R/figs/Sorption_attenuation_BC.pdf")

#table with attenuation factors


