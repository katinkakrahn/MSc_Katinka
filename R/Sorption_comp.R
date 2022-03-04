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

#Personal access token
#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

Sorption <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/160222_sorption_rawdata.xlsx")
Sorption <- as.data.table(Sorption)

#Convert 1 and 0 to TRUE and FALSE and delete integer columns
Sorption$SoilLogic <- as.logical(Sorption$Soil_binary)
Sorption$mixLogic <- as.logical(Sorption$mix_binary)
Sorption <- subset(Sorption,select = -c(Soil_binary,mix_binary))
Sorption_BC <- kable(Sorption, "latex", booktabs = TRUE, digits = 2)

# Subset biochar and cocktail/single compound
Sorption_NAomit <- na.omit(Sorption)
Sorption_NA_C1omit <- Sorption_NAomit %>% slice(-c(20, 30, 40, 50, 79, 88, 98, 108, 118, 157, 167, 177))
Sorption_BC_single <- subset(Sorption_NA_C1omit, mixLogic == FALSE)
Sorption_BC_mix <- subset(Sorption_NA_C1omit, mixLogic == TRUE)

CWC_single <- filter(Sorption_BC_single, Biochar == "CWC")
ULS_single <- filter(Sorption_BC_single, Biochar == "ULS")
DSL_single <- filter(Sorption_BC_single, Biochar == "DSL")

#CWC Freundlich isotherm plot
CWC_single$Compound <- factor(CWC_single$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))

CWC_isotherm <- ggplot(data = CWC_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  # geom_smooth(mapping = aes(x = log_Cw, y = log_Cs), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE, colour = "grey",
  #             data = CWC_single) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("CWC isotherm") +
  facet_wrap(~ Compound) +
  geom_label(data = summary_stats_CWC_single, inherit.aes = FALSE, aes(x = 0, y = 0.2, label = paste("K_F =",K_F,","," ","n =",n,","," ","R^2",r_squared))) +
  theme_bw() +
  guides(color = "none")
  #theme(legend.position = c(0.9, 0.15))
CWC_isotherm

nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

summary_stats_CWC_single <- data.table(K_F = rep(0, nr_compounds), 
                                       K_F_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "CWC")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = CWC_single[Compound == compounds[i]])
  summary_stats_CWC_single[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_CWC_single[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_CWC_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_CWC_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_CWC_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_CWC_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_CWC_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                       summary(fit)$fstatistic[3],lower.tail=F)]
}

#ULS Freundlich isotherm plot
ULS_single$Compound <- factor(ULS_single$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))

ULS_isotherm <- ggplot(data = ULS_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "ULS isotherm") + 
  theme_bw() #+
  #theme(legend.position = c(0.9, 0.15))
ULS_isotherm

summary_stats_ULS_single <- data.table(K_F = rep(0, nr_compounds), 
                                       K_F_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "ULS")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = ULS_single[Compound == compounds[i]])
  summary_stats_ULS_single[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_ULS_single[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_ULS_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_ULS_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_ULS_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_ULS_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_ULS_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

#DSL Freundlich isotherm plot
DSL_single$Compound <- factor(DSL_single$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))
DSL_isotherm <- ggplot(data = DSL_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "DSL isotherm") + 
  theme_bw() #+
  #theme(legend.position = c(0.9, 0.15))
DSL_isotherm

summary_stats_DSL_single <- data.table(K_F = rep(0, nr_compounds), 
                                       K_F_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "DSL")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = DSL_single[Compound == compounds[i]])
  summary_stats_DSL_single[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_DSL_single[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_DSL_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_DSL_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_DSL_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_DSL_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_DSL_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}
