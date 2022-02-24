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