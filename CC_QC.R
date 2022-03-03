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

CCQC <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/010322_Calibration_data.xlsx")
as.data.table(CCQC)
CCQC <- as.data.table(CCQC)

#Mean for each compound keyed by concentration, CC/QC and Area/Area/IS
CCQCsummary <- CCQC[, .(mean_signal = mean(Signal)
),
keyby = .(Compound, Cal_type, Correction, Concentration_ppb)]

#Subset CC and Area/area IS
CCsummary <- CCQCsummary[Cal_type == "CC"]
CCsummary_IS <- CCsummary[Correction == "Area/Area IS"]

PFPeA_CCsummary <- CCsummary_IS[Compound == "PFPeA"]
PFHxA_CCsummary <- CCsummary_IS[Compound == "PFHxA"]
PFHpA_CCsummary <- CCsummary_IS[Compound == "PFHpA"]
PFOA_CCsummary <- CCsummary_IS[Compound == "PFOA"]
PFNA_CCsummary <- CCsummary_IS[Compound == "PFNA"]
PFDA_CCsummary <- CCsummary_IS[Compound == "PFDA"]


lm_PFPeA <- lm(mean_signal~Concentration_ppb, data = PFPeA_CCsummary)
summary(lm_PFPeA)

CC_PFPeA <- ggplot(data = PFPeA_CCsummary) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal)) + 
  geom_smooth(method = "lm", mapping = aes(x = Concentration_ppb, y = mean_signal), formula = y ~ x, se=FALSE, fullrange = TRUE) +
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Calibraton curve PFPeA") + 
  stat_regline_equation(
    aes(x = Concentration_ppb, y = mean_signal, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
    formula = y ~ x
    ) +
  theme_bw()
CC_PFPeA

CC_PFHxA <- ggplot(data = PFHxA_CCsummary) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal)) + 
  geom_smooth(method = "lm", mapping = aes(x = Concentration_ppb, y = mean_signal), formula = y ~ x, se=FALSE, fullrange = TRUE) +
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Calibraton curve PFHxA") + 
  stat_regline_equation(
    aes(x = Concentration_ppb, y = mean_signal, label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x
  ) +
  theme_bw()
CC_PFHxA


CC_PFHpA <- ggplot(data = PFHpA_CCsummary) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal)) + 
  geom_smooth(method = "lm", mapping = aes(x = Concentration_ppb, y = mean_signal), formula = y ~ x, se=FALSE, fullrange = TRUE) +
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Calibraton curve PFHpA") + 
  stat_regline_equation(
    aes(x = Concentration_ppb, y = mean_signal, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
    formula = y ~ x
  ) +
  theme_bw()
CC_PFHpA

CC_PFOA <- ggplot(data = PFOA_CCsummary) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal)) + 
  geom_smooth(method = "lm", mapping = aes(x = Concentration_ppb, y = mean_signal), formula = y ~ x, se=FALSE, fullrange = TRUE) +
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Calibraton curve PFOA") + 
  stat_regline_equation(
    aes(x = Concentration_ppb, y = mean_signal, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
    formula = y ~ x
  ) +
  theme_bw()
CC_PFOA

CC_PFNA <- ggplot(data = PFNA_CCsummary) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal)) + 
  geom_smooth(method = "lm", mapping = aes(x = Concentration_ppb, y = mean_signal), formula = y ~ x, se=FALSE, fullrange = TRUE) +
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Calibraton curve PFNA") + 
  stat_regline_equation(
    aes(x = Concentration_ppb, y = mean_signal, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
    formula = y ~ x
  ) +
  theme_bw()
CC_PFNA

CC_PFDA <- ggplot(data = PFDA_CCsummary) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal)) + 
  geom_smooth(method = "lm", mapping = aes(x = Concentration_ppb, y = mean_signal), formula = y ~ x, se=FALSE, fullrange = TRUE) +
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Calibraton curve PFDA") + 
  stat_regline_equation(
    aes(x = Concentration_ppb, y = mean_signal, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
    formula = y ~ x
  ) +
  theme_bw()
CC_PFDA

CCsummary_IS$Compound <- factor(CCsummary_IS$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))

CC_all <- ggplot(data = CCsummary_IS) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Calibraton curve", col = "Compound") + 
  stat_regline_equation(
    aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound), label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
    formula = y ~ x
  ) +
  theme_bw() +
  theme(legend.position = c(0.3, 0.87))
CC_all

#Subset QC and Area/Area IS
QCsummary <- CCQCsummary[Cal_type == "QC"]
QCsummary_IS <- QCsummary[Correction == "Area/Area IS"]
QCsummary_IS <- na.omit(QCsummary_IS)

QCsummary_IS$Compound <- factor(QCsummary_IS$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                  "PFOA", "PFNA", "PFDA"))

QC_all <- ggplot(data = QCsummary_IS) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Matrix matched calibraton curve", col = "Compound") + 
  stat_regline_equation(
    aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound), label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
    formula = y ~ x
  ) +
  theme_bw() +
  theme(legend.position = c(0.3, 0.87))
QC_all



