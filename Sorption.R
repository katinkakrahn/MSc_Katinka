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

Sorption <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Lab/160222_sorption_rawdata.xlsx")
as.data.table(Sorption)
Sorption <- as.data.table(Sorption)
order <- data.table(order = c(1:6))
Sorption[order, order := "PFPeA" == 1, "PFHxA" == 2, "PFHpA" == 3, "PFOA" == 4, "PFNA" == 5, "PFDA" == 6]

#Convert 1 and 0 to TRUE and FALSE and delete integer columns
Sorption$SoilLogic <- as.logical(Sorption$Soil_binary)
Sorption$mixLogic <- as.logical(Sorption$mix_binary)
Sorption <- subset(Sorption,select = -c(Soil_binary,mix_binary))
Sorption_BC <- kable(Sorption, "latex", booktabs = TRUE, digits = 2)

# Subset biochar and cocktail/single compound
Sorption_NAomit <- na.omit(Sorption)
Sorption_C1omit <- filter(Sorption, between(Conc_point, 2, 10))
Sorption_C1NAomit <- na.omit(Sorption_C1omit)
Sorption_BCsingleComp_all <- subset(Sorption_NAomit, mixLogic == FALSE)
Sorption_BCsingleComp_C1omit <- subset(Sorption_C1NAomit, mixLogic == FALSE)
Sorption_BCmixComp <- subset(Sorption_NAomit, mixLogic == TRUE)

CWC_single_all <- filter(Sorption_BCsingleComp_all, Biochar == "CWC")
ULS_single_all <- filter(Sorption_BCsingleComp_all, Biochar == "ULS")
BRL_single_all <- filter(Sorption_BCsingleComp_all, Biochar == "BRL")

CWC_single_C1omit <- filter(Sorption_BCsingleComp_C1omit, Biochar == "CWC")
ULS_single_C1omit <- filter(Sorption_BCsingleComp_C1omit, Biochar == "ULS")
BRL_single_C1omit <- filter(Sorption_BCsingleComp_C1omit, Biochar == "BRL")

#CWC Freundlich isotherm plot all points
CWC_single_all$Compound <- factor(CWC_single_all$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                      "PFOA", "PFNA", "PFDA"))
CWC_isotherm_all <- ggplot(data = CWC_single_all) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE)
CWC_isotherm_all <- CWC_isotherm_all + labs(x = "log C_w", y = "log C_s")
CWC_isotherm_all

nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)

# her initsialiseres et skjelett til en data table som du kan fylle inn i for løkken
summary_stats_CWC_single_all <- data.table(K_F = rep(0, nr_compounds), 
                                           K_F_std_error = rep(0, nr_compounds),
                                           K_F_t_value = rep (0, nr_compounds),
                                           K_F_p_value_t = rep(0, nr_compounds),
                                           n = rep(0, nr_compounds),
                                           n_std_error = rep(0, nr_compounds),
                                           n_t_value = rep(0, nr_compounds),
                                           n_F_p_value_t = rep(0, nr_compounds),
                                           r_squared = rep(0, nr_compounds),
                                           r_squared_adj = rep(0, nr_compounds),
                                           residual_std_error = rep(0, nr_compounds),
                                           F_statistic = rep(0, nr_compounds),
                                           p_value = rep(0, nr_compounds),
                                           nr_points = rep("all", nr_compounds),
                                           compound = compounds)

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = CWC_single_all[Compound == compounds[i]])
  #print(summary(fit))
  #print(compounds[i])
  summary_stats_CWC_single_all[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_CWC_single_all[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_CWC_single_all[compound == compounds[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_CWC_single_all[compound == compounds[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_CWC_single_all[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_CWC_single_all[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_CWC_single_all[compound == compounds[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_CWC_single_all[compound == compounds[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_CWC_single_all[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_CWC_single_all[compound == compounds[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_CWC_single_all[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_CWC_single_all[compound == compounds[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_CWC_single_all[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                       summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_CWC_single_all

#CWC Freundlich isotherm plot omit C1
CWC_single_C1omit$Compound <- factor(CWC_single_C1omit$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                            "PFOA", "PFNA", "PFDA"))
CWC_isotherm_C1omit <- ggplot(data = CWC_single_C1omit) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE)
CWC_isotherm_C1omit <- CWC_isotherm_C1omit + labs(x = "log C_w", y = "log C_s")
CWC_isotherm_C1omit

summary_stats_CWC_single_C1omit <- data.table(K_F = rep(0, nr_compounds), 
                                              K_F_std_error = rep(0, nr_compounds),
                                              K_F_t_value = rep (0, nr_compounds),
                                              K_F_p_value_t = rep(0, nr_compounds),
                                              n = rep(0, nr_compounds),
                                              n_std_error = rep(0, nr_compounds),
                                              n_t_value = rep(0, nr_compounds),
                                              n_F_p_value_t = rep(0, nr_compounds),
                                              r_squared = rep(0, nr_compounds),
                                              r_squared_adj = rep(0, nr_compounds),
                                              residual_std_error = rep(0, nr_compounds),
                                              F_statistic = rep(0, nr_compounds),
                                              p_value = rep(0, nr_compounds),
                                              nr_points = rep("C1omit", nr_compounds),
                                              compound = compounds)

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = CWC_single_C1omit[Compound == compounds[i]])
  #print(summary(fit))
  #print(compounds[i])
  summary_stats_CWC_single_C1omit[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_CWC_single_C1omit[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_CWC_single_C1omit[compound == compounds[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_CWC_single_C1omit[compound == compounds[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_CWC_single_C1omit[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_CWC_single_C1omit[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_CWC_single_C1omit[compound == compounds[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_CWC_single_C1omit[compound == compounds[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_CWC_single_C1omit[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_CWC_single_C1omit[compound == compounds[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_CWC_single_C1omit[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_CWC_single_C1omit[compound == compounds[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_CWC_single_C1omit[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                       summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_CWC_single_C1omit

#Compare C1omit points to C1 omitted CWC 
compare_fit_CWC_singleComp <- merge(summary_stats_CWC_single_all,summary_stats_CWC_single_C1omit, all = TRUE)
compare_fit_CWC_singleComp$compound <- factor(compare_fit_CWC_singleComp$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                              "PFOA", "PFNA", "PFDA"))


###############################################################################################################################################################################
###############################################################################################################################################################################

#ULS Freundlich isotherm plot all points
ULS_single_all$Compound <- factor(ULS_single_all$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                      "PFOA", "PFNA", "PFDA"))
ULS_isotherm_all <- ggplot(data = ULS_single_all) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE)
ULS_isotherm_all <- ULS_isotherm_all + labs(x = "log C_w", y = "log C_s")
ULS_isotherm_all

summary_stats_ULS_single_all <- data.table(K_F = rep(0, nr_compounds), 
                                           K_F_std_error = rep(0, nr_compounds),
                                           K_F_t_value = rep (0, nr_compounds),
                                           K_F_p_value_t = rep(0, nr_compounds),
                                           n = rep(0, nr_compounds),
                                           n_std_error = rep(0, nr_compounds),
                                           n_t_value = rep(0, nr_compounds),
                                           n_F_p_value_t = rep(0, nr_compounds),
                                           r_squared = rep(0, nr_compounds),
                                           r_squared_adj = rep(0, nr_compounds),
                                           residual_std_error = rep(0, nr_compounds),
                                           F_statistic = rep(0, nr_compounds),
                                           p_value = rep(0, nr_compounds),
                                           nr_points = rep("all", nr_compounds),
                                           compound = compounds)

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = ULS_single_all[Compound == compounds[i]])
  #print(summary(fit))
  #print(compounds[i])
  summary_stats_ULS_single_all[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_ULS_single_all[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_ULS_single_all[compound == compounds[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_ULS_single_all[compound == compounds[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_ULS_single_all[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_ULS_single_all[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_ULS_single_all[compound == compounds[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_ULS_single_all[compound == compounds[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_ULS_single_all[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_ULS_single_all[compound == compounds[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_ULS_single_all[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_ULS_single_all[compound == compounds[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_ULS_single_all[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                       summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_ULS_single_all

#ULS Freundlich isotherm plot omit C1
ULS_single_C1omit$Compound <- factor(ULS_single_C1omit$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                            "PFOA", "PFNA", "PFDA"))
ULS_isotherm_C1omit <- ggplot(data = ULS_single_C1omit) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE)
ULS_isotherm_C1omit <- ULS_isotherm_C1omit + labs(x = "log C_w", y = "log C_s")
ULS_isotherm_C1omit

summary_stats_ULS_single_C1omitomit <- data.table(K_F = rep(0, nr_compounds), 
                                           K_F_std_error = rep(0, nr_compounds),
                                           K_F_t_value = rep (0, nr_compounds),
                                           K_F_p_value_t = rep(0, nr_compounds),
                                           n = rep(0, nr_compounds),
                                           n_std_error = rep(0, nr_compounds),
                                           n_t_value = rep(0, nr_compounds),
                                           n_F_p_value_t = rep(0, nr_compounds),
                                           r_squared = rep(0, nr_compounds),
                                           r_squared_adj = rep(0, nr_compounds),
                                           residual_std_error = rep(0, nr_compounds),
                                           F_statistic = rep(0, nr_compounds),
                                           p_value = rep(0, nr_compounds),
                                           nr_points = rep("C1omitomit", nr_compounds),
                                           compound = compounds)

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = ULS_single_C1omitomit[Compound == compounds[i]])
  #print(summary(fit))
  #print(compounds[i])
  summary_stats_ULS_single_C1omit[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_ULS_single_C1omit[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_ULS_single_C1omit[compound == compounds[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_ULS_single_C1omit[compound == compounds[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_ULS_single_C1omit[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_ULS_single_C1omit[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_ULS_single_C1omit[compound == compounds[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_ULS_single_C1omit[compound == compounds[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_ULS_single_C1omit[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_ULS_single_C1omit[compound == compounds[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_ULS_single_C1omit[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_ULS_single_C1omit[compound == compounds[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_ULS_single_C1omit[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                       summary(fit)$fstatistic[3],lower.tail=F)]
}

#Compare C1omit points to C1omit omitted ULS
compare_fit_ULS_singleComp <- merge(summary_stats_ULS_single_all,summary_stats_ULS_single_C1omit, all = TRUE)
compare_fit_ULS_singleComp$compound <- factor(compare_fit_ULS_singleComp$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                              "PFOA", "PFNA", "PFDA"))

###############################################################################################################################################################################
###############################################################################################################################################################################


#BRL Freundlich isotherm plot all points
BRL_single_all$Compound <- factor(BRL_single_all$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                      "PFOA", "PFNA", "PFDA"))
BRL_isotherm_all <- ggplot(data = BRL_single_all) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE)
BRL_isotherm_all <- BRL_isotherm_all + labs(x = "log C_w", y = "log C_s")
BRL_isotherm_all

summary_stats_BRL_single_all <- data.table(K_F = rep(0, nr_compounds), 
                                           K_F_std_error = rep(0, nr_compounds),
                                           K_F_t_value = rep (0, nr_compounds),
                                           K_F_p_value_t = rep(0, nr_compounds),
                                           n = rep(0, nr_compounds),
                                           n_std_error = rep(0, nr_compounds),
                                           n_t_value = rep(0, nr_compounds),
                                           n_F_p_value_t = rep(0, nr_compounds),
                                           r_squared = rep(0, nr_compounds),
                                           r_squared_adj = rep(0, nr_compounds),
                                           residual_std_error = rep(0, nr_compounds),
                                           F_statistic = rep(0, nr_compounds),
                                           p_value = rep(0, nr_compounds),
                                           nr_points = rep("all", nr_compounds),
                                           compound = compounds)

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = BRL_single_all[Compound == compounds[i]])
  #print(summary(fit))
  #print(compounds[i])
  summary_stats_BRL_single_all[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_BRL_single_all[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_BRL_single_all[compound == compounds[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_BRL_single_all[compound == compounds[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_BRL_single_all[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_BRL_single_all[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_BRL_single_all[compound == compounds[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_BRL_single_all[compound == compounds[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_BRL_single_all[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_BRL_single_all[compound == compounds[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_BRL_single_all[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_BRL_single_all[compound == compounds[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_BRL_single_all[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                       summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_BRL_single_all


#BRL Freundlich isotherm plot omit C1
BRL_single_C1omit$Compound <- factor(BRL_single_C1omit$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                            "PFOA", "PFNA", "PFDA"))
BRL_isotherm_C1omit <- ggplot(data = BRL_single_C1omit) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE)
BRL_isotherm_C1omit <- BRL_isotherm_C1omit + labs(x = "log C_w", y = "log C_s")
BRL_isotherm_C1omit

summary_stats_BRL_single_C1omit <- data.table(K_F = rep(0, nr_compounds), 
                                           K_F_std_error = rep(0, nr_compounds),
                                           K_F_t_value = rep (0, nr_compounds),
                                           K_F_p_value_t = rep(0, nr_compounds),
                                           n = rep(0, nr_compounds),
                                           n_std_error = rep(0, nr_compounds),
                                           n_t_value = rep(0, nr_compounds),
                                           n_F_p_value_t = rep(0, nr_compounds),
                                           r_squared = rep(0, nr_compounds),
                                           r_squared_adj = rep(0, nr_compounds),
                                           residual_std_error = rep(0, nr_compounds),
                                           F_statistic = rep(0, nr_compounds),
                                           p_value = rep(0, nr_compounds),
                                           nr_points = rep("C1omit", nr_compounds),
                                           compound = compounds)
for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = BRL_single_C1omit[Compound == compounds[i]])
  #print(summary(fit))
  #print(compounds[i])
  summary_stats_BRL_single_all[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_BRL_single_all[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_BRL_single_all[compound == compounds[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_BRL_single_all[compound == compounds[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_BRL_single_all[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_BRL_single_all[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_BRL_single_all[compound == compounds[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_BRL_single_all[compound == compounds[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_BRL_single_all[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_BRL_single_all[compound == compounds[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_BRL_single_all[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_BRL_single_all[compound == compounds[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_BRL_single_all[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                       summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_BRL_single_all

#Compare C1omit points to C1 omitted BRL
compare_fit_BRL_singleComp <- merge(summary_stats_BRL_single_all,summary_stats_BRL_single_C1omit, all = TRUE)
compare_fit_BRL_singleComp$compound <- factor(compare_fit_BRL_singleComp$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                              "PFOA", "PFNA", "PFDA"))
