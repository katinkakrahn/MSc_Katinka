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

#Outputs isotherms by biochar:
CWC_isotherm
ULS_isotherm
BRL_isotherm

#Outputs isotherms by compound:
PFPeA_isotherm
PFHxA_isotherm
PFHpA_isotherm
PFOA_isotherm
PFNA_isotherm
PFDA_isotherm

#Compare sorption of biochars across char types
summary_stats_CWC_single[, nr_CF2 := 4:9]
summary_stats_ULS_single[, nr_CF2 := 4:9]
summary_stats_DSL_single[, nr_CF2 := 4:9]

summary_stats_single <- merge(summary_stats_CWC_single, summary_stats_ULS_single, all = TRUE)
summary_stats_single <- merge(summary_stats_single, summary_stats_DSL_single, all = TRUE)
summary_stats_single$compound <- factor(summary_stats_single$compound, levels = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA"))

#Summary stats of each compound
summary_stats_PFPeA <- filter(summary_stats_single, compound == "PFPeA")
summary_stats_PFHxA <- filter(summary_stats_single, compound == "PFHxA")
summary_stats_PFHpA <- filter(summary_stats_single, compound == "PFHpA")
summary_stats_PFOA <- filter(summary_stats_single, compound == "PFOA")
summary_stats_PFNA <- filter(summary_stats_single, compound == "PFNA")
summary_stats_PFDA <- filter(summary_stats_single, compound == "PFDA")