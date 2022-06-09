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
library(RColorBrewer)

CCQC <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/010322_Calibration_data.xlsx")
as.data.table(CCQC)
CCQC <- as.data.table(CCQC)

#Mean for each compound keyed by concentration, CC/QC and Area/Area/IS
CCQCsummary <- CCQC[, .(mean_signal = mean(Signal),
                        sd_signal = sd(Signal)
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

CC_allinone <- ggplot(data = CCsummary_IS) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound))) +
  geom_smooth(mapping = aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound)),
              formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  scale_colour_brewer(palette = "Set2") +
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Calibraton curve", col = "Compound") +
  facet_wrap(~Compound) +
  theme_bw()
CC_allinone
ggsave(filename = "R/figs/CC_allinone.pdf")

facetCC <- CCsummary_IS |>
  transform(pre_compound = Compound)

facetCC <- rbind(
  transform(facetCC, Compound = unique(CCsummary_IS$Compound)[1]),
  transform(facetCC, Compound = unique(CCsummary_IS$Compound)[2]),
  transform(facetCC, Compound = unique(CCsummary_IS$Compound)[3]),
  transform(facetCC, Compound = unique(CCsummary_IS$Compound)[4]),
  transform(facetCC, Compound = unique(CCsummary_IS$Compound)[5]),
  transform(facetCC, Compound = unique(CCsummary_IS$Compound)[6])
)



CC_all <- ggplot(data = CCsummary_IS) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal, group = factor(Compound)), 
             color = "gray45", size = 1) + 
  # #geom_smooth(mapping = aes(x = Concentration_ppb, y = mean_signal, group = pre_compound), 
  #             formula = y ~ x, 
  #             method=lm, 
  #             se=FALSE, 
  #             colour = "grey", 
  #             size = 0.5,
  #             data = facetCC
  # ) +
  geom_smooth(mapping = aes(x = Concentration_ppb, y = mean_signal, group = factor(Compound)), 
              color = "black", 
              formula = y ~ x, 
              method=lm, 
              se=FALSE, 
              fullrange = FALSE) + 
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)") + 
  #ggtitle("Calibration curve") +
  facet_wrap(~Compound) +
  theme_bw() #+
  #theme(panel.grid = element_blank())
CC_all
ggsave(filename="R/figs/CC_all.pdf")

#Subset QC and Area/Area IS
QCsummary <- CCQCsummary[Cal_type == "QC"]
QCsummary_IS <- QCsummary[Correction == "Area/Area IS"]
QCsummary_IS <- na.omit(QCsummary_IS)

QCsummary_IS$Compound <- factor(QCsummary_IS$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                  "PFOA", "PFNA", "PFDA"))

QC_allinone <- ggplot(data = QCsummary_IS) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = Concentration_ppb, y = mean_signal, color = factor(Compound)), 
              formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  scale_colour_brewer(palette = "Set2") +
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", title = "Matrix matched calibraton curve", 
       col = "Compound") +
  theme_bw()
QC_allinone
ggsave(filename = "R/figs/QC_allinone.pdf")


facetQC <- QCsummary_IS |>
  transform(pre_compound = Compound)

facetQC <- rbind(
  transform(facetQC, Compound = unique(QCsummary_IS$Compound)[1]),
  transform(facetQC, Compound = unique(QCsummary_IS$Compound)[2]),
  transform(facetQC, Compound = unique(QCsummary_IS$Compound)[3]),
  transform(facetQC, Compound = unique(QCsummary_IS$Compound)[4]),
  transform(facetQC, Compound = unique(QCsummary_IS$Compound)[5]),
  transform(facetQC, Compound = unique(QCsummary_IS$Compound)[6])
)



QC_all <- ggplot(data = QCsummary_IS) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal, group = factor(Compound)), 
             color = "gray45", size = 1) + 
  # #geom_smooth(mapping = aes(x = Concentration_ppb, y = mean_signal, group = pre_compound), 
  #             formula = y ~ x, 
  #             method=lm, 
  #             se=FALSE, 
  #             colour = "grey", 
  #             size = 0.5,
  #             data = facetCC
  # ) +
  geom_smooth(mapping = aes(x = Concentration_ppb, y = mean_signal, group = factor(Compound)), 
              color = "black", 
              formula = y ~ x, 
              method=lm, 
              se=FALSE, 
              fullrange = FALSE) + 
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)") + 
  #ggtitle("Calibration curve") +
  facet_wrap(~Compound) +
  theme_bw() #+
#theme(panel.grid = element_blank())
QC_all
ggsave(filename="R/figs/QC_all.pdf")

#matrix effect
CCQCsummary_IS <- filter(CCQCsummary, Correction == "Area/Area IS")
CCQCsummary_IS <- filter(CCQCsummary_IS, Compound != "PFOA_13C8")

CCQCsummary_IS$Compound <- factor(CCQCsummary_IS$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                  "PFOA", "PFNA", "PFDA"))
CCQC_matrixeffect <- ggplot(data = CCQCsummary_IS) +
  geom_point(mapping = aes(x = Concentration_ppb, y = mean_signal, group = interaction(Compound, Cal_type), color = Cal_type), size = 1) + 
  geom_smooth(mapping = aes(x = Concentration_ppb, y = mean_signal, group = interaction(Compound, Cal_type), color = Cal_type), 
              formula = y ~ x, 
              method=lm, 
              se=FALSE, 
              fullrange = FALSE) + 
  labs(x = "Concentration (ppb)", y = "Signal (area/area IS)", col = "Calibration") + 
  #guides(color = "none") +
  #ggtitle("Calibration curve") +
  facet_wrap(~Compound) +
  theme_bw() #+
#theme(panel.grid = element_blank())
CCQC_matrixeffect
ggsave(filename="R/figs/CCQC_matrixeffect.pdf")
