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

PFPeA_sorption_single <- Sorption_BC_single[Compound == "PFPeA"]
PFHxA_sorption_single <- Sorption_BC_single[Compound == "PFHxA"]
PFHpA_sorption_single <- Sorption_BC_single[Compound == "PFHpA"]
PFOA_sorption_single <- Sorption_BC_single[Compound == "PFOA"]
PFNA_sorption_single <- Sorption_BC_single[Compound == "PFNA"]
PFDA_sorption_single <- Sorption_BC_single[Compound == "PFDA"]

#PFPeA
PFPeA_isotherm <- ggplot(data = PFPeA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(method = "lm", mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFPeA") + 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFPeA_isotherm

#PFHxA
PFHxA_isotherm <- ggplot(data = PFHxA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFHxA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFHxA_isotherm

#PFHpA
PFHpA_isotherm <- ggplot(data = PFHpA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFHpA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFHpA_isotherm

#PFOA
PFOA_isotherm <- ggplot(data = PFOA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFOA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFOA_isotherm

#PFNA
PFNA_isotherm <- ggplot(data = PFNA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFNA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFNA_isotherm

#PFDA
PFDA_isotherm <- ggplot(data = PFDA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFDA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFDA_isotherm



