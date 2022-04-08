Kd_1ugL <- summary_stats_single[, .(log_Cs = log_KF + n*log_Cw,
                                    log_Kd = log_KF + n*log_Cw,
                                    Kd = 10^(log_KF + n*log_Cw)
  
), keyby = .(compound, biochar)]

Kd_1ugL_select <- Kd_1ugL[compound != "PFHpA" | biochar != "ULS"]
Kd_1ugL_select <- Kd_1ugL_select[compound != "PFPeA" | biochar != "DSL"]
Kd_1ugL_select <- Kd_1ugL_select[compound != "PFHxA" | biochar != "DSL"]
Kd_1ugL_select <- Kd_1ugL_select[compound != "PFPeA" | biochar != "CWC"]
Kd_1ugL_select <- Kd_1ugL_select[compound != "PFHxA" | biochar != "CWC"]

Biochar_use <- subset(Biochar, Key == 1)
Biochar_CWC <- subset(Biochar_use, biochar == "CWC")
Biochar_ULS <- subset(Biochar_use, biochar == "ULS")
Biochar_DSL <- subset(Biochar_use, biochar == "DSL")

Biochar1 <- subset(Biochar_use, biochar == "CWC")
names <- Biochar1[, Parameter]
names <- as.data.table(names)
  
Biochar_join <- inner_join(Biochar_CWC, Biochar_ULS, by = "Parameter")
Biochar_join <- inner_join(Biochar_join, Biochar_DSL, by = "Parameter")
Biochar_join <- Biochar_join[, .(Parameter, Mean_diffunit.x, Mean_diffunit.y, Mean_diffunit)]
setnames(Biochar_join, c("Mean_diffunit.x", "Mean_diffunit.y", "Mean_diffunit"), c("CWC", "ULS", "DSL"))
Biochar_join_t <- setNames(data.frame(t(Biochar_join[, - 1])), Biochar_join[, 1])
colnames(Biochar_join_t) <- t(names)

Biochar_ratios <- Biochar_join_t
Biochar_ratios$C_Ca <- Biochar_ratios$C/Biochar_ratios$Ca
Biochar_ratios$C_Fe <- Biochar_ratios$C/Biochar_ratios$Fe
Biochar_ratios$SAN2_Fe <- Biochar_ratios$N2_BET_SA/Biochar_ratios$Fe
Biochar_ratios$SAN2_Ca <- Biochar_ratios$N2_BET_SA/Biochar_ratios$Ca
Biochar_ratios$SACO2_Fe <- Biochar_ratios$CO2_DFT_SA/Biochar_ratios$Fe
Biochar_ratios$SACO2_Ca <- Biochar_ratios$CO2_DFT_SA/Biochar_ratios$Ca
Biochar_ratios$PVN2_Fe <- Biochar_ratios$N2_BJH_PV/Biochar_ratios$Fe
Biochar_ratios$PVN2_Ca <- Biochar_ratios$N2_BJH_PV/Biochar_ratios$Ca
Biochar_ratios$PVCO2_Fe <- Biochar_ratios$CO2_DFT_PV/Biochar_ratios$Fe
Biochar_ratios$PVCO2_Ca <- Biochar_ratios$CO2_DFT_PV/Biochar_ratios$Ca
Biochar_ratios$SA_PV <- Biochar_ratios$N2_SA/Biochar_ratios$N2_PV
Biochar_ratios$SA_PV_Ca <- (Biochar_ratios$N2_SA/Biochar_ratios$N2_PV)/Biochar_ratios$Ca
Biochar_ratios$SA_PV_C <- (Biochar_ratios$N2_SA/Biochar_ratios$N2_PV)/Biochar_ratios$C
Biochar_ratios$PV_C <- Biochar_ratios$N2_PV/Biochar_ratios$C
Biochar_ratios$SA_C <- Biochar_ratios$N2_SA/Biochar_ratios$C


Biochar_ratios <- as.data.table(Biochar_ratios, keep.rownames = T)
setnames(Biochar_ratios, "rn", "biochar")
write_xlsx(Biochar_ratios, "/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/280322_Biochar_ratios.xlsx")
setnames(Biochar_ratios, c("C",	"Ca",	"Fe",	"N2_BET_SA",	"N2_BJH_PV",	"CO2_DFT_SA",	"CO2_DFT_PV",	"C_Ca",	"C_Fe",	"SAN2_Fe",	"SAN2_Ca",	"SACO2_Fe",	"SACO2_Ca",	"PVN2_Fe",	"PVN2_Ca",	"PVCO2_Fe",	"PVCO2_Ca"), 
         c("C",	"Ca",	"Fe",	"N2_SA2",	"N2_PV",	"CO2_SA",	"CO2_PV",	"C_Ca",	"C_Fe",	"SAN2_Fe",	"SAN2_Ca",	"SACO2_Fe",	"SACO2_Ca",	"PVN2_Fe",	"PVN2_Ca",	"PVCO2_Fe",	"PVCO2_Ca"))

Biochar_ratios_1ugL_select <- full_join(Biochar_ratios, Kd_1ugL_select, by = "biochar")
Biochar_ratios_1ugL_all <- full_join(Biochar_ratios, Kd_1ugL, by = "biochar")

nr_compounds2 <- length(unique(Biochar_ratios_1ugL_select$compound))
compounds2 <- unique(Biochar_ratios_1ugL_select$compound)

Corr_table <- data.table(r_squared = rep(0, nr_compounds2),
                         p_value = rep(0, nr_compounds2),
                         compound = compounds2)

for(i in 1:nr_compounds2){
  fit <- lm(log_Kd ~ SA_PV_C, data = Biochar_ratios_1ugL_select[compound == compounds2[i]])
  Corr_table[compound == compounds2[i], r_squared := summary(fit)$r.squared]
  Corr_table[compound == compounds2[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                           summary(fit)$fstatistic[3],lower.tail=F)]
}




