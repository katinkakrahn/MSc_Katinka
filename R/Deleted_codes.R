#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

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
                                           #resid_min = rep(0, nr_compounds),
                                           #resid_1Q = rep(0, nr_compounds),
                                           #resid_median = rep(0, nr_compounds),
                                           #resid_3Q = rep(0, nr_compounds),
                                           #resid_max = rep(0, nr_compounds),
                                           residual_std_error = rep(0, nr_compounds),
                                           F_statistic = rep(0, nr_compounds),
                                           p_value = rep(0, nr_compounds),
                                           nr_points = rep("all", nr_compounds),
                                           compound = compounds)

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = CWC_single_all[Compound == compounds[i]])
  print(summary(fit))
  print(compounds[i])
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
  #summary_stats_CWC_single_all[compound == compounds[i], resid_min := summary(fit$residuals)[,1]]
  #summary_stats_CWC_single_all[compound == compounds[i], resid_1Q := summary(fit$residuals)[,2]]
  #summary_stats_CWC_single_all[compound == compounds[i], resid_median := summary(fit$residuals)[,3]]
  #summary_stats_CWC_single_all[compound == compounds[i], resid_3Q := summary(fit$residuals)[,4]]
  #summary_stats_CWC_single_all[compound == compounds[i], resid_max := summary(fit$residuals)[,5]]
  summary_stats_CWC_single_all[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_CWC_single_all[compound == compounds[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_CWC_single_all[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_CWC_single_all

order <- data.table(order = c(1:6))
Sorption[order, order := "PFPeA" == 1, "PFHxA" == 2, "PFHpA" == 3, "PFOA" == 4, "PFNA" == 5, "PFDA" == 6]

summary_stats_compounds <- merge(summary_stats_PFPeA, summary_stats_PFHxA, all = TRUE)
summary_stats_compounds <- merge(summary_stats_compounds, summary_stats_PFHpA, all = TRUE)
summary_stats_compounds <- merge(summary_stats_compounds, summary_stats_PFOA, all = TRUE)
summary_stats_compounds <- merge(summary_stats_compounds, summary_stats_PFNA, all = TRUE)
summary_stats_compounds <- merge(summary_stats_compounds, summary_stats_PFDA, all = TRUE)

#Compare C1omit points to C1 omitted CWC 
compare_fit_CWC_singleComp <- merge(summary_stats_CWC_single_all,summary_stats_CWC_single_C1omit, all = TRUE)
compare_fit_CWC_singleComp$compound <- factor(compare_fit_CWC_singleComp$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                              "PFOA", "PFNA", "PFDA"))

#Compare C1omit points to C1omit omitted ULS
compare_fit_ULS_singleComp <- merge(summary_stats_ULS_single_all,summary_stats_ULS_single_C1omit, all = TRUE)
compare_fit_ULS_singleComp$compound <- factor(compare_fit_ULS_singleComp$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                              "PFOA", "PFNA", "PFDA"))

summary_stats_PFPeA <- data.table(K_F = rep(0, nr_biochars),
                                  K_F_std_error = rep(0, nr_biochars),
                                  K_F_t_value = rep (0, nr_biochars),
                                  K_F_p_value_t = rep(0, nr_biochars),
                                  n = rep(0, nr_biochars),
                                  n_std_error = rep(0, nr_biochars),
                                  n_t_value = rep(0, nr_biochars),
                                  n_F_p_value_t = rep(0, nr_biochars),
                                  r_squared = rep(0, nr_biochars),
                                  r_squared_adj = rep(0, nr_biochars),
                                  residual_std_error = rep(0, nr_biochars),
                                  F_statistic = rep(0, nr_biochars),
                                  p_value = rep(0, nr_biochars),
                                  nr_points = rep("all", nr_biochars),
                                  Biochar = biochars)

for(i in 1:nr_biochars){
  fit <- lm(log_Cs ~ log_Cw, data = PFPeA[Biochar == biochars[i]])
  #print(summary(fit))
  #print(biochars[i])
  summary_stats_PFPeA[Biochar == biochars[i], K_F := fit$coefficients[1]]
  summary_stats_PFPeA[Biochar == biochars[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_PFPeA[Biochar == biochars[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_PFPeA[Biochar == biochars[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_PFPeA[Biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_PFPeA[Biochar == biochars[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_PFPeA[Biochar == biochars[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_PFPeA[Biochar == biochars[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_PFPeA[Biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_PFPeA[Biochar == biochars[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_PFPeA[Biochar == biochars[i], residual_std_error := summary(fit)$sigma]
  summary_stats_PFPeA[Biochar == biochars[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_PFPeA[Biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                            summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_PFPeA

summary_stats_PFHxA <- data.table(K_F = rep(0, nr_biochars),
                                  K_F_std_error = rep(0, nr_biochars),
                                  K_F_t_value = rep (0, nr_biochars),
                                  K_F_p_value_t = rep(0, nr_biochars),
                                  n = rep(0, nr_biochars),
                                  n_std_error = rep(0, nr_biochars),
                                  n_t_value = rep(0, nr_biochars),
                                  n_F_p_value_t = rep(0, nr_biochars),
                                  r_squared = rep(0, nr_biochars),
                                  r_squared_adj = rep(0, nr_biochars),
                                  residual_std_error = rep(0, nr_biochars),
                                  F_statistic = rep(0, nr_biochars),
                                  p_value = rep(0, nr_biochars),
                                  nr_points = rep("all", nr_biochars),
                                  Biochar = biochars)

for(i in 1:nr_biochars){
  fit <- lm(log_Cs ~ log_Cw, data = PFHxA[Biochar == biochars[i]])
  #print(summary(fit))
  #print(biochars[i])
  summary_stats_PFHxA[Biochar == biochars[i], K_F := fit$coefficients[1]]
  summary_stats_PFHxA[Biochar == biochars[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_PFHxA[Biochar == biochars[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_PFHxA[Biochar == biochars[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_PFHxA[Biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_PFHxA[Biochar == biochars[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_PFHxA[Biochar == biochars[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_PFHxA[Biochar == biochars[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_PFHxA[Biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_PFHxA[Biochar == biochars[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_PFHxA[Biochar == biochars[i], residual_std_error := summary(fit)$sigma]
  summary_stats_PFHxA[Biochar == biochars[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_PFHxA[Biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                            summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_PFHxA

summary_stats_PFHpA <- data.table(K_F = rep(0, nr_biochars),
                                  K_F_std_error = rep(0, nr_biochars),
                                  K_F_t_value = rep (0, nr_biochars),
                                  K_F_p_value_t = rep(0, nr_biochars),
                                  n = rep(0, nr_biochars),
                                  n_std_error = rep(0, nr_biochars),
                                  n_t_value = rep(0, nr_biochars),
                                  n_F_p_value_t = rep(0, nr_biochars),
                                  r_squared = rep(0, nr_biochars),
                                  r_squared_adj = rep(0, nr_biochars),
                                  residual_std_error = rep(0, nr_biochars),
                                  F_statistic = rep(0, nr_biochars),
                                  p_value = rep(0, nr_biochars),
                                  nr_points = rep("all", nr_biochars),
                                  Biochar = biochars)

for(i in 1:nr_biochars){
  fit <- lm(log_Cs ~ log_Cw, data = PFHpA[Biochar == biochars[i]])
  #print(summary(fit))
  #print(biochars[i])
  summary_stats_PFHpA[Biochar == biochars[i], K_F := fit$coefficients[1]]
  summary_stats_PFHpA[Biochar == biochars[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_PFHpA[Biochar == biochars[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_PFHpA[Biochar == biochars[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_PFHpA[Biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_PFHpA[Biochar == biochars[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_PFHpA[Biochar == biochars[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_PFHpA[Biochar == biochars[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_PFHpA[Biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_PFHpA[Biochar == biochars[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_PFHpA[Biochar == biochars[i], residual_std_error := summary(fit)$sigma]
  summary_stats_PFHpA[Biochar == biochars[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_PFHpA[Biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                            summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_PFHpA

summary_stats_PFOA <- data.table(K_F = rep(0, nr_biochars),
                                 K_F_std_error = rep(0, nr_biochars),
                                 K_F_t_value = rep (0, nr_biochars),
                                 K_F_p_value_t = rep(0, nr_biochars),
                                 n = rep(0, nr_biochars),
                                 n_std_error = rep(0, nr_biochars),
                                 n_t_value = rep(0, nr_biochars),
                                 n_F_p_value_t = rep(0, nr_biochars),
                                 r_squared = rep(0, nr_biochars),
                                 r_squared_adj = rep(0, nr_biochars),
                                 residual_std_error = rep(0, nr_biochars),
                                 F_statistic = rep(0, nr_biochars),
                                 p_value = rep(0, nr_biochars),
                                 nr_points = rep("all", nr_biochars),
                                 Biochar = biochars)

for(i in 1:nr_biochars){
  fit <- lm(log_Cs ~ log_Cw, data = PFOA[Biochar == biochars[i]])
  #print(summary(fit))
  #print(biochars[i])
  summary_stats_PFOA[Biochar == biochars[i], K_F := fit$coefficients[1]]
  summary_stats_PFOA[Biochar == biochars[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_PFOA[Biochar == biochars[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_PFOA[Biochar == biochars[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_PFOA[Biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_PFOA[Biochar == biochars[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_PFOA[Biochar == biochars[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_PFOA[Biochar == biochars[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_PFOA[Biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_PFOA[Biochar == biochars[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_PFOA[Biochar == biochars[i], residual_std_error := summary(fit)$sigma]
  summary_stats_PFOA[Biochar == biochars[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_PFOA[Biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                           summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_PFOA

summary_stats_PFNA <- data.table(K_F = rep(0, nr_biochars),
                                 K_F_std_error = rep(0, nr_biochars),
                                 K_F_t_value = rep (0, nr_biochars),
                                 K_F_p_value_t = rep(0, nr_biochars),
                                 n = rep(0, nr_biochars),
                                 n_std_error = rep(0, nr_biochars),
                                 n_t_value = rep(0, nr_biochars),
                                 n_F_p_value_t = rep(0, nr_biochars),
                                 r_squared = rep(0, nr_biochars),
                                 r_squared_adj = rep(0, nr_biochars),
                                 residual_std_error = rep(0, nr_biochars),
                                 F_statistic = rep(0, nr_biochars),
                                 p_value = rep(0, nr_biochars),
                                 nr_points = rep("all", nr_biochars),
                                 Biochar = biochars)

for(i in 1:nr_biochars){
  fit <- lm(log_Cs ~ log_Cw, data = PFNA[Biochar == biochars[i]])
  #print(summary(fit))
  #print(biochars[i])
  summary_stats_PFNA[Biochar == biochars[i], K_F := fit$coefficients[1]]
  summary_stats_PFNA[Biochar == biochars[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_PFNA[Biochar == biochars[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_PFNA[Biochar == biochars[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_PFNA[Biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_PFNA[Biochar == biochars[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_PFNA[Biochar == biochars[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_PFNA[Biochar == biochars[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_PFNA[Biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_PFNA[Biochar == biochars[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_PFNA[Biochar == biochars[i], residual_std_error := summary(fit)$sigma]
  summary_stats_PFNA[Biochar == biochars[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_PFNA[Biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                           summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_PFNA



summary_stats_PFDA <- data.table(K_F = rep(0, nr_biochars),
                                 K_F_std_error = rep(0, nr_biochars),
                                 K_F_t_value = rep (0, nr_biochars),
                                 K_F_p_value_t = rep(0, nr_biochars),
                                 n = rep(0, nr_biochars),
                                 n_std_error = rep(0, nr_biochars),
                                 n_t_value = rep(0, nr_biochars),
                                 n_F_p_value_t = rep(0, nr_biochars),
                                 r_squared = rep(0, nr_biochars),
                                 r_squared_adj = rep(0, nr_biochars),
                                 residual_std_error = rep(0, nr_biochars),
                                 F_statistic = rep(0, nr_biochars),
                                 p_value = rep(0, nr_biochars),
                                 nr_points = rep("all", nr_biochars),
                                 Biochar = biochars)

for(i in 1:nr_biochars){
  fit <- lm(log_Cs ~ log_Cw, data = PFDA[Biochar == biochars[i]])
  #print(summary(fit))
  #print(biochars[i])
  summary_stats_PFDA[Biochar == biochars[i], K_F := fit$coefficients[1]]
  summary_stats_PFDA[Biochar == biochars[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_PFDA[Biochar == biochars[i], K_F_t_value := summary(fit)$coefficients[1,3]]
  summary_stats_PFDA[Biochar == biochars[i], K_F_p_value_t := summary(fit)$coefficients[1,4]]
  summary_stats_PFDA[Biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_PFDA[Biochar == biochars[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_PFDA[Biochar == biochars[i], n_t_value := summary(fit)$coefficients[2,3]]
  summary_stats_PFDA[Biochar == biochars[i], n_F_p_value_t := summary(fit)$coefficients[2,4]]
  summary_stats_PFDA[Biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_PFDA[Biochar == biochars[i], r_squared_adj := summary(fit)$adj.r.squared]
  summary_stats_PFDA[Biochar == biochars[i], residual_std_error := summary(fit)$sigma]
  summary_stats_PFDA[Biochar == biochars[i], F_statistic := summary(fit)$fstatistic[1]]
  summary_stats_PFDA[Biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                           summary(fit)$fstatistic[3],lower.tail=F)]
}
summary_stats_PFDA

#Compare C1omit points to C1 omitted BRL
compare_fit_BRL_singleComp <- merge(summary_stats_BRL_single_all,summary_stats_BRL_single_C1omit, all = TRUE)
compare_fit_BRL_singleComp$compound <- factor(compare_fit_BRL_singleComp$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                              "PFOA", "PFNA", "PFDA"))
#C1 omit
Sorption_C1omit <- filter(Sorption_NAomit, between(Conc_point, 2, 10))
Sorption_BCsingleComp_C1omit <- subset(Sorption_C1omit, mixLogic == FALSE)
CWC_single_C1omit <- filter(Sorption_BCsingleComp_C1omit, Biochar == "CWC")
ULS_single_C1omit <- filter(Sorption_BCsingleComp_C1omit, Biochar == "ULS")
BRL_single_C1omit <- filter(Sorption_BCsingleComp_C1omit, Biochar == "BRL")

#CWC Freundlich isotherm plot omit C1
CWC_single_C1omit$Compound <- factor(CWC_single_C1omit$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                            "PFOA", "PFNA", "PFDA"))
CWC_isotherm_C1omit <- ggplot(data = CWC_single_C1omit) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "CWC isotherm C1 omitted") + 
  theme_bw() +
  theme(legend.position = c(0.9, 0.15))
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
                                              compound = compounds,
                                              biochar = "CWC")

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

#Compare C1omit points to C1 omitted
compare_fit_CWC_singleComp <- filter(summary_stats_single_allC1omit, biochar == "CWC")
compare_fit_ULS_singleComp <- filter(summary_stats_single_allC1omit, biochar == "ULS")
compare_fit_BRL_singleComp <- filter(summary_stats_single_allC1omit, biochar == "BRL")

pHcondSummarySummary <- pHcondSummary[, .(mean_ph = mean(mean_ph), 
                                          mean_cond = mean(mean_cond),
                                          sd_cond = sd(sd_cond),
                                          sd_ph = sd(sd_ph)
)]

pHcondSummary_S <- pHcondSummary %>% slice(-c(1))
pHcondSummarySummary_S <- pHcondSummary_S[, .(mean_ph = mean(mean_ph), 
                                              mean_cond = mean(mean_cond),
                                              sd_cond = sd(sd_cond),
                                              sd_ph = sd(sd_ph))]

#TukeyHSD(pH_ANOVA_all, conf.level = 0.95)
pH_TKHSD <- TukeyHSD(pH_ANOVA_all, "Sample", ordered = TRUE)
pH_TKHSD <- as.data.frame(pH_TKHSD$Sample)

pH_ANOVA_BC <- aov(pH ~ Sample, data = pHcond_BC)
summary(pH_ANOVA_BC)
pH_BC_TKHSD <- TukeyHSD(pH_ANOVA_BC, "Sample", ordered = TRUE)
pH_BC_TKHSD <- as.data.frame(pH_BC_TKHSD$Sample)

pH_ANOVA_S <- aov(pH ~ Sample, data = pHcond_BC_S)
summary(pH_ANOVA_S)
pH_S_TKHSD <- TukeyHSD(pH_ANOVA_S, "Sample", ordered = TRUE)
pH_S_TKHSD <- as.data.frame(pH_S_TKHSD$Sample)

#ANOVA and Tukey HSD cond
cond_ANOVA_all <- aov(Conductivity ~ Sample, data = pHcond)
summary(cond_ANOVA_all)
cond_TKHSD <- TukeyHSD(cond_ANOVA_all, "Sample", ordered = TRUE)
cond_TKHSD <- as.data.frame(cond_TKHSD$Sample)

cond_ANOVA_BC <- aov(Conductivity ~ Sample, data = pHcond_BC)
summary(cond_ANOVA_BC)
cond_BC_TKHSD <- TukeyHSD(cond_ANOVA_BC, "Sample", ordered = TRUE)
cond_BC_TKHSD <- as.data.frame(cond_BC_TKHSD$Sample)

cond_ANOVA_S <- aov(Conductivity ~ Sample, data = pHcond_BC_S)
summary(cond_ANOVA_S)
cond_S_TKHSD <- TukeyHSD(cond_ANOVA_S, "Sample", ordered = TRUE)
cond_S_TKHSD <- as.data.frame(cond_S_TKHSD$Sample)

# ANOVA and Tukey HSD pH
pHcond_BC <- filter(pHcond, Sample == "CWC" | Sample == "ULS" | Sample == "BRL")
pHcond_BC_S <- filter (pHcond, Sample == "CWC+S" | Sample == "ULS+S" | Sample == "BRL+S")

pH_ANOVA_all <- aov(pH ~ Sample, data = pHcond)
summary(pH_ANOVA_all)

#copy of CWC isotherm facet before changing labels
CWC_facet_isotherm <- ggplot(data = CWC_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_compound), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetCWC) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("CWC isotherm") +
  geom_label(data = transform(summary_stats_CWC_single, Compound = compound), size = 3, inherit.aes = T, 
             aes(x = 0.25, y = 0.5, label = (paste("K_F =",round(K_F, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2))))) +
  facet_wrap(~Compound) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none")

#DSL Freundlich isotherm plot
DSL_single$Compound <- factor(DSL_single$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))

facetDSL <- DSL_single |>
  transform(pre_compound = Compound)

facetDSL <- rbind(
  transform(facetDSL, Compound = unique(DSL_single$Compound)[1]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[2]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[3]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[4]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[5]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[6])
)

DSL_isotherm <- ggplot(data = DSL_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_compound), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetDSL) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("DSL isotherm") +
  geom_label(data = transform(summary_stats_CWC_single, Compound = compound), size = 2, inherit.aes = T, aes(x = 0, y = 0.3, label = paste("K_F=",round(K_F, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~ Compound) +
  scale_colour_brewer(palette = "Set2") +
  theme_bw() +
  guides(color = "none")
DSL_isotherm
ggsave(filename="figs/DSL_isotherm.pdf")

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

#DSL Freundlich isotherm plot
DSL_single$Compound <- factor(DSL_single$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))
facetDSL$Compound <- factor(facetDSL$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                          "PFOA", "PFNA", "PFDA"))
summary_stats_DSL_single$compound <- factor(summary_stats_DSL_single$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                          "PFOA", "PFNA", "PFDA"))

facetDSL <- DSL_single |>
  transform(pre_compound = Compound)

facetDSL <- rbind(
  transform(facetDSL, Compound = unique(DSL_single$Compound)[1]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[2]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[3]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[4]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[5]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[6])
)

DSL_isotherm <- ggplot(data = DSL_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_compound), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetDSL) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("DSL isotherm") +
  geom_label(data = transform(summary_stats_DSL_single, Compound = compound), size = 2, inherit.aes = FALSE, aes(x = -0.5, y = 3.5, label = paste("K_F =",round(K_F, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =",round(r_squared, digits = 2)))) +
  facet_wrap(~ Compound) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(group = "none")
DSL_isotherm

ggsave(filename="R/figs/DSL_isotherm.pdf")

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

Sorption_attenuation_BC <- ggplot() +
  geom_point(data = Sorption_BC_single_C10_common, mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar), shape = mixLogic, size = 2), 
  ) + 
  geom_point(data = Sorption_BC_mix_summary, mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar), shape = mixLogic, size = 2), 
  ) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), shape = "Cocktail", color = "Biochar") +
  guides(size = "none") +
  facet_wrap(~Compound) +
  theme_bw() +
  theme(panel.grid = element_blank())
Sorption_attenuation_BC
set_palette(Sorption_attenuation_BC, "uchicago")
ggsave(filename="R/figs/Sorption_attenuation_BC.pdf")

SA <- ggplot(data = subset(SA_PV, Porosity %in% "SA"), 
             aes(x = Gas, y = Mean_diffunit, color = Biochar)) +
  geom_point() +
  labs(x = "", y = "m2/g") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
SA

PV <- ggplot(data = subset(SA_PV, Porosity %in% "PV"), 
             aes(x = Gas, y = Mean_diffunit, color = biochar)) +
  geom_point() +
  labs(x = "", y = "cm3/g") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
PV

C_Ca_ratios_plot <- ggplot(data = subset(Elements_ratios, Parameter %in% "C_Ca"), 
                           aes(x = Biochar, y = Ratio, group = Biochar)) +
  geom_point() +
  # geom_smooth(formula = y ~ x, 
  #             method=lm, 
  #             se = FALSE, 
  #             fullrange = FALSE,
  #             color = "grey") +
  labs(x = expression(""), y = "C/Ca ratio") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
C_Ca_ratios_plot


C_Fe_ratios_plot <- ggplot(data = subset(Elements_ratios, Parameter %in% "C_Fe"), 
                           aes(x = Biochar, y = Ratio, group = Biochar)) +
  geom_point() +
  # geom_smooth(formula = y ~ x, 
  #             method=lm, 
  #             se = FALSE, 
  #             fullrange = FALSE,
  #             color = "grey") +
  labs(x = expression(""), y = "C/Fe ratio") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
C_Fe_ratios_plot

Elements_CWC <- subset(Elements_biochar, biochar == "CWC")
Elements_ULS <- subset(Elements_biochar, biochar == "ULS")
Elements_DSL <- subset(Elements_biochar, biochar == "DSL")

Elements_biochar_join <- inner_join(Elements_CWC, Elements_ULS, by = "Parameter")
Elements_biochar_join <- inner_join(Elements_biochar_join, Elements_DSL, by = "Parameter")
Elements_biochar_join <- Elements_biochar_join[, .(Parameter, Mean_diffunit.x, Mean_diffunit.y, Mean_diffunit)]
setnames(Elements_biochar_join, c("Mean_diffunit.x", "Mean_diffunit.y", "Mean_diffunit"), c("CWC", "ULS", "DSL"))
Elements_biochar_join_t <- setNames(data.frame(t(Elements_biochar_join[, - 1])), Elements_biochar_join[, 1])
colnames(Elements_biochar_join_t) <- t(names)


Elements_ratios <- Elements_biochar_join_t
Elements_ratios$C_Ca <- Elements_ratios$C/Elements_ratios$Ca
Elements_ratios$C_Fe <- Elements_ratios$C/Elements_ratios$Fe
Elements_ratios <- as.data.table(Elements_ratios, keep.rownames = T)
setnames(Elements_ratios, "rn", "Biochar")