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

Standards <- read_excel("/Users/katinkakrahn/Documents/Skole/VOW/Lab/180222_Standards_rawdata.xlsx")
as.data.table(Standards)
Standards <- as.data.table(Standards)
Standards_single <- filter(Standards, Type == "Std"| Type == "FB")
Standard_mix <- filter(Standards, Type == "MIX-STD")

#Rawdata table of standard concentrations
print()
formatSignif(
  Standards,
  digits = 2,
  interval = 3,
  mark = ",",
  dec.mark = getOption("OutDec")
)

print(Standards, digits = 2)
Standards_table <- kable(Standards, "latex", booktabs = TRUE)


#Mean and standard deviation
Standards_single_Summary <- Standards_single[, .(mean_conc = mean(C13), 
                                                 sd_conc = sd(C13),
                                                 ),
                                             keyby = .(Compound)]
Standards_single_Summary

#Plot to check for outliers
Standards$Compound <- factor(Standards$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA"))
Standard_conc <- ggplot(Standards, aes(x = Compound, y = C13)) +
  geom_point()
Standard_conc <- Standard_conc + geom_errorbar(aes(ymin=mean_ph-sd_ph, ymax=mean_ph+sd_ph), width=.2,
                                               position=position_dodge(.9))
Standard_conc