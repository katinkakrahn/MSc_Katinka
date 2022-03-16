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

Sorption_soil <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/010322_sorption_rawdata_soil.xlsx")
as.data.table(Sorption_soil)
Sorption_soil <- as.data.table(Sorption_soil)

Sorption_soil$SoilLogic <- as.logical(Sorption_soil$Soil_binary)
Sorption_soil$mixLogic <- as.logical(Sorption_soil$mix_binary)
Sorption_soil <- subset(Sorption_soil,select = -c(Soil_binary,mix_binary))
Sorption_BC_BS <- kable(Sorption_soil, "latex", booktabs = TRUE, digits = 2)

# Subset biochar and cocktail/single compound, change from Sorption_soil to Sorption_soil_NAomit when data is updated
# Sorption_soil_NAomit <- na.omit(Sorption_soil) #this one must be turned on when actual data is there
# Sorption_soil_NA_C1omit <- Sorption_NAomit %>% slice(-c(column numbers))
Sorption_soil_blank <- subset(Sorption_soil, Biochar == "no") 
Sorption_soil_blank_mix <- subset(Sorption_soil_blank, mixLogic == TRUE)
Sorption_soil_blank_PFOA <- subset(Sorption_soil_blank, mixLogic == FALSE)
Sorption_soil_BC <- Sorption_soil[Biochar != 'no'] 
Sorption_soil_BC_PFOA <- subset(Sorption_soil_BC, mixLogic == FALSE) 
Sorption_soil_BC_mix <- subset(Sorption_soil_BC, mixLogic == TRUE)

CWC_soil_mix <- filter(Sorption_soil_BC_mix, Biochar == "CWC")
ULS_soil_mix <- filter(Sorption_soil_BC_mix, Biochar == "ULS")
DSL_soil_mix <- filter(Sorption_soil_BC_mix, Biochar == "DSL")

CWC_soil_PFOA <- filter(Sorption_soil_BC_PFOA, Biochar == "CWC")
ULS_soil_PFOA <- filter(Sorption_soil_BC_PFOA, Biochar == "ULS")
DSL_soil_PFOA <- filter(Sorption_soil_BC_PFOA, Biochar == "DSL")

Sorption_soil_blank_mix$Compound <- factor(Sorption_soil_blank_mix$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                  "PFOA", "PFNA", "PFDA"))

Sorption_soil_blank_mix_points <- ggplot(data = Sorption_soil_blank_mix) + 
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)))+ 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "Soil blank cocktail triplicates") + 
  theme_bw()
Sorption_soil_blank_mix_points

Sorption_soil_blank_PFOA$Compound <- factor(Sorption_soil_blank_PFOA$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                        "PFOA", "PFNA", "PFDA"))

Sorption_soil_blank_PFOA_points <- ggplot(data = Sorption_soil_blank_PFOA) + 
  geom_point(mapping = aes(x = log_Cw, y = log_Cs), shape = 21, size = 3, fill = "#077DAA")+ 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), title = "Soil blank PFOA triplicates") + 
  theme_bw()
Sorption_soil_blank_PFOA_points

# Soil cocktail isotherms
CWC_soil_mix$Compound <- factor(CWC_soil_mix$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))

CWC_isotherm_soil_mix <- ggplot(data = CWC_soil_mix) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "CWC soil cocktail isotherm") + 
  theme_bw()
CWC_isotherm_soil_mix

ULS_soil_mix$Compound <- factor(ULS_soil_mix$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                  "PFOA", "PFNA", "PFDA"))
ULS_isotherm_soil_mix <- ggplot(data = ULS_soil_mix) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "ULS soil cocktail isotherm") + 
  theme_bw()
ULS_isotherm_soil_mix

DSL_soil_mix$Compound <- factor(DSL_soil_mix$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                  "PFOA", "PFNA", "PFDA"))

DSL_isotherm_soil_mix <- ggplot(data = DSL_soil_mix) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "DSL soil cocktail isotherm") + 
  theme_bw()
DSL_isotherm_soil_mix

# Soil PFOA isotherms
CWC_isotherm_soil_PFOA <- ggplot(data = CWC_soil_PFOA) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs)) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), title = "CWC soil PFOA isotherm") + 
  theme_bw()
CWC_isotherm_soil_PFOA

ULS_isotherm_soil_PFOA <- ggplot(data = ULS_soil_PFOA) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs)) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), title = "ULS soil PFOA isotherm") + 
  theme_bw()
ULS_isotherm_soil_PFOA

DSL_isotherm_soil_PFOA <- ggplot(data = DSL_soil_PFOA) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs)) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), title = "DSL soil PFOA isotherm") + 
  theme_bw()
DSL_isotherm_soil_PFOA

PFOA_soil_isotherm <- ggplot(data = Sorption_soil_BC_PFOA) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFOA biochar soil isotherm") + 
  theme_bw()
PFOA_soil_isotherm
