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

Standards <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/180222_Standards_rawdata.xlsx")
as.data.table(Standards)
Standards <- as.data.table(Standards)
Standards_single <- filter(Standards, Type == "Std"| Type == "FB")
Standard_mix <- filter(Standards, Type == "MIX-STD")

#Rawdata table of standard concentrations
#print(Standards, digits = 2)
Standards_table <- kable(Standards, "latex", booktabs = TRUE)


#Mean and standard deviation
Standards_single_summary <- Standards_single[, .(mean_conc = mean(C13), 
                                                 sd_conc = sd(C13)
                                                 ),
                                             keyby = .(Compound)]
Standards_single_summary

#Plot to check for outliers
Standards_single_summary$Compound <- factor(Standards_single_summary$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA"))
Standard_conc <- ggplot(data = Standards_single_summary, aes(x = Compound, y = mean_conc)) +
  geom_point() + 
  geom_errorbar(aes(ymin=mean_conc-sd_conc, ymax=mean_conc+sd_conc), width=.2, position=position_dodge(.9)) +
  labs(x = "Compound", y = "concentration (\u03bcg/L)") +
  theme_bw()
Standard_conc
