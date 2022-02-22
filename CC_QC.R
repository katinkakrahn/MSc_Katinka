library(data.table)
library(ggplot2)
library(psych)
library(readxl)
library(grDevices)
library(dplyr)
library(broom)
library(ggpubr)
library(tidyverse)
library(magrittr)
library(sandwich)
library(tidyverse)

CCQC <- read_excel("/Users/katinkakrahn/Documents/Skole/VOW/Lab/200222_Calibration_data.xlsx")
CCQC <- as.data.table(CCQC)
CCQC

#Mean for each compound keyed by concentration, CC/QC and Area/Area/IS
CCQCsummary <- CCQC[, .(mean_PFPeA = mean(PFPeA),
                        mean_PFHxA = mean(PFHxA),
                        mean_HFHpA = mean(PFHpA), 
                        mean_PFOA = mean(PFOA),
                        mean_PFNA = mean(PFNA),
                        mean_PFDA = mean(PFDA),
                        mean_PFOA_IS = mean(`PFOA-13C8`)
),
keyby = .(Concentration_ppb, Cal_type, Correction)]
CCQCsummary

#Subset CC
CCsummary <- CCQCsummary[Cal_type == "CC"]
CCsummary <- CCsummary[order(Correction)]
CCsummary

#Subset QC
QCsummary <- CCQCsummary[Cal_type == "QC"]
QCsummary <- QCsummary[order(Correction)]
QCsummary

CC <- ggplot(data = pHcondSummary, aes(x = reorder(Sample, mean_ph), y = mean_ph))
pH <- pH + geom_point()
pH <- pH + geom_errorbar(aes(ymin=mean_ph-sd_ph, ymax=mean_ph+sd_ph), width=.2,
                         position=position_dodge(.9))
pH <- pH + labs(x = "", y = "pH")
pH + theme_bw()
pH
