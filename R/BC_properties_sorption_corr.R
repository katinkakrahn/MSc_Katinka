# Biochar ratios ---- 
Biochar <- read_excel("R/data_raw/250322_biochar_parameters.xlsx")
as.data.table(Biochar)
Biochar <- setnames(Biochar, "Biochar", "biochar")

Biochar_CWC <- subset(Biochar, biochar == "CWC")
Biochar_ULS <- subset(Biochar, biochar == "ULS")
Biochar_DSL <- subset(Biochar, biochar == "DSL")

Biochar1 <- subset(Biochar, biochar == "CWC")
names <- unique(Biochar1$Parameter)
names <- as.data.table(names)

Biochar_join <- inner_join(Biochar_CWC, Biochar_ULS, by = "Parameter")
Biochar_join <- inner_join(Biochar_join, Biochar_DSL, by = "Parameter") 
Biochar_join <- select(Biochar_join, Parameter, Mean_diffunit.x, Mean_diffunit.y, Mean_diffunit)
Biochar_join <- setnames(Biochar_join, c("Mean_diffunit.x", "Mean_diffunit.y", "Mean_diffunit"), c("CWC", "ULS", "DSL"))
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
Biochar_ratios$SA_PV <- Biochar_ratios$N2_BET_SA/Biochar_ratios$N2_BJH_PV
Biochar_ratios$SA_PV_Ca <- (Biochar_ratios$N2_BET_SA/Biochar_ratios$N2_BJH_PV)/Biochar_ratios$Ca
Biochar_ratios$SA_PV_C <- (Biochar_ratios$N2_BET_SA/Biochar_ratios$N2_BJH_PV)/Biochar_ratios$C
Biochar_ratios$PV_C <- Biochar_ratios$N2_BJH_PV/Biochar_ratios$C
Biochar_ratios$SA_C <- Biochar_ratios$N2_BET_SA/Biochar_ratios$C

Biochar_ratios <- as.data.table(Biochar_ratios, keep.rownames = T)
setnames(Biochar_ratios, "rn", "biochar")
write_xlsx(Biochar_ratios, "R/data_manipulated/280322_Biochar_ratios.xlsx")
setnames(Biochar_ratios, c("C",	"Ca",	"Fe",	"N2_BET_SA",	"N2_BJH_PV",	"CO2_DFT_SA",	"CO2_DFT_PV",	"C_Ca",	"C_Fe",	"SAN2_Fe",	"SAN2_Ca",	"SACO2_Fe",	"SACO2_Ca",	"PVN2_Fe",	"PVN2_Ca",	"PVCO2_Fe",	"PVCO2_Ca"), 
         c("C",	"Ca",	"Fe",	"N2_SA2",	"N2_PV",	"CO2_SA",	"CO2_PV",	"C_Ca",	"C_Fe",	"SAN2_Fe",	"SAN2_Ca",	"SACO2_Fe",	"SACO2_Ca",	"PVN2_Fe",	"PVN2_Ca",	"PVCO2_Fe",	"PVCO2_Ca"))


# Summary stats single compound isotherms ----
Sorption <- read_excel("R/data_raw/160222_sorption_rawdata.xlsx")
as.data.table(Sorption)
Sorption <- as.data.table(Sorption)

#Convert 1 and 0 to TRUE and FALSE and delete integer columns
Sorption$SoilLogic <- as.logical(Sorption$Soil_binary)
Sorption$mixLogic <- as.logical(Sorption$mix_binary)
Sorption <- select(Sorption, -Soil_binary, -mix_binary)

# Subset biochar and cocktail/single compound
Sorption_NAomit <- na.omit(Sorption)
Sorption_NA_C1omit <- Sorption_NAomit %>% 
  filter(Conc_point != 1) %>% 
  mutate(Kd = Cs/Cw, log_Kd = log10(Cs/Cw))
Sorption_BC_single <- subset(Sorption_NA_C1omit, mixLogic == FALSE)
Sorption_BC_mix <- subset(Sorption_NA_C1omit, mixLogic == TRUE)

CWC_single <- filter(Sorption_BC_single, Biochar == "CWC")
ULS_single <- filter(Sorption_BC_single, Biochar == "ULS")
DSL_single <- filter(Sorption_BC_single, Biochar == "DSL")

nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

summary_stats_CWC_single <- data.table(log_KF = rep(0, nr_compounds), 
                                       log_KF_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       log_Cw = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "CWC")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = CWC_single[Compound == compounds[i]])
  summary_stats_CWC_single[compound == compounds[i], log_KF := fit$coefficients[1]]
  summary_stats_CWC_single[compound == compounds[i], log_KF_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_CWC_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_CWC_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_CWC_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_CWC_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_CWC_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}


summary_stats_ULS_single <- data.table(log_KF = rep(0, nr_compounds), 
                                       log_KF_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       log_Cw = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "ULS")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = ULS_single[Compound == compounds[i]])
  summary_stats_ULS_single[compound == compounds[i], log_KF := fit$coefficients[1]]
  summary_stats_ULS_single[compound == compounds[i], log_KF_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_ULS_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_ULS_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_ULS_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_ULS_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_ULS_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_DSL_single <- data.table(log_KF = rep(0, nr_compounds), 
                                       log_KF_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       log_Cw = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "DSL")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = DSL_single[Compound == compounds[i]])
  summary_stats_DSL_single[compound == compounds[i], log_KF := fit$coefficients[1]]
  summary_stats_DSL_single[compound == compounds[i], log_KF_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_DSL_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_DSL_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_DSL_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_DSL_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_DSL_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_CWC_single[, nr_CF2 := 4:9]
summary_stats_ULS_single[, nr_CF2 := 4:9]
summary_stats_DSL_single[, nr_CF2 := 4:9]

summary_stats_single <- merge(summary_stats_CWC_single, summary_stats_ULS_single, all = TRUE)
summary_stats_single <- merge(summary_stats_single, summary_stats_DSL_single, all = TRUE)


# Kd at 1 ug/L ----
Kd_1ugL <- summary_stats_single[, .(log_Cs = log_KF + n*log_Cw,
                                    log_Kd = log_KF + n*log_Cw,
                                    Kd = 10^(log_KF + n*log_Cw),
                                    logKd_error = sqrt(log_KF_std_error^2+n_std_error^2)
                                    
), keyby = .(compound, biochar)]

# Kd for the compounds with isotherm within the 1 ug/L range
Kd_1ugL_select <- Kd_1ugL %>% 
  filter(compound != "PFPeA", compound != "PFHxA", compound !="PFHpA")

Biochar_ratios_1ugL_select <- merge(Biochar_ratios, Kd_1ugL_select, by = "biochar")

Kd_1ugL_select <- merge(Kd_1ugL_select, summary_stats_single, all=FALSE, by = c("compound", "biochar")) 


nr_compounds2 <- length(unique(Kd_1ugL_select$compound))
compounds2 <- unique(Kd_1ugL_select$compound)

Corr_table <- data.table(r_squared = rep(0, nr_compounds2),
                         p_value = rep(0, nr_compounds2),
                         compound = compounds2)

for(i in 1:nr_compounds2){
  fit <- lm(log_Kd ~ SA_PV_C, 
            data = Biochar_ratios_1ugL_select[compound == compounds2[i]])
  Corr_table[compound == compounds2[i], r_squared := summary(fit)$r.squared]
  Corr_table[compound == compounds2[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                      summary(fit)$fstatistic[3],lower.tail=F)]
}


#Chain length vs Kd 1ug/L plot ----
Kd_1ugL[compound=="PFPeA", nr_CF2:=4]
Kd_1ugL[compound=="PFHxA", nr_CF2:=5]
Kd_1ugL[compound=="PFHpA", nr_CF2:=6]
Kd_1ugL[compound=="PFOA", nr_CF2:=7]
Kd_1ugL[compound=="PFNA", nr_CF2:=8]
Kd_1ugL[compound=="PFDA", nr_CF2:=9]

Kd_1ugL_chain_length <- Kd_1ugL[compound != "PFPeA"]
Kd_1ugL_chain_length <- Kd_1ugL_chain_length[compound != "PFHxA" | biochar != "DSL"]
Kd_1ugL_chain_length <- Kd_1ugL_chain_length[compound != "PFHxA" | biochar != "CWC"]
Kd_1ugL_chain_length <- Kd_1ugL_chain_length[compound != "PFHpA" | biochar != "ULS"]

Kd_1ugL_select$biochar <- factor(Kd_1ugL_select$biochar, levels = c("ULS", "DSL", "CWC"))

Kd_1ugL_plot <- ggplot(data = Kd_1ugL_chain_length, 
                       aes(x = nr_CF2, y = log_Kd, color = biochar)) + 
  geom_point(size = 4) + 
  geom_line(size = 1) + 
  geom_errorbar(aes(ymin = log_Kd-logKd_error, 
                    ymax = log_Kd+logKd_error), 
                width = 0.05) +
  labs(x = TeX(r'($CF_2~chain~length$)'), 
       y = TeX(r'($log~K_d~(at~C_w~1 \mu g/L)$)'), 
       color = "", 
       shape = "") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values= c("#767676FF","#800000FF","#FFB547FF"))+
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom", 
        text = element_text(size = 16))
Kd_1ugL_plot
ggsave(filename="R/figs/Kd_1ugL_plot.pdf")

summary_stats_1ugL_CWC <- lm(log_Kd ~ nr_CF2, data = subset(Kd_1ugL_chain_length, biochar == "CWC"))
summary(summary_stats_1ugL_CWC)
summary_stats_1ugL_ULS <- lm(log_Kd ~ nr_CF2, data = subset(Kd_1ugL_chain_length, biochar == "ULS"))
summary(summary_stats_1ugL_ULS)
summary_stats_1ugL_DSL <- lm(log_Kd ~ nr_CF2, data = subset(Kd_1ugL_chain_length, biochar == "DSL"))
summary(summary_stats_1ugL_DSL)

# Kd 1 ng/L plot (not good) ----
Kd_1ugL_n <- merge(Kd_1ugL_chain_length, summary_stats_single, all=FALSE, by = c("compound", "biochar"))
Kd_1ngL <- Kd_1ugL_n[, Kd_1ngL:=log_Kd+3*(1-n)]
Kd_1ngL <- Kd_1ngL[, logKd_error := sqrt(log_KF_std_error^2+n_std_error^2)]
Kd_1ugL <- Kd_1ugL_n[, logKd_error := sqrt(log_KF_std_error^2+n_std_error^2)]


Kd_1ngL_plot <- ggplot(data = Kd_1ngL, aes(x = nr_CF2.x, y = Kd_1ngL, color = biochar)) + 
  geom_point(size = 4) + 
  geom_line(size = 1) + 
  labs(x = "", y = expression(log~K[d]~(1~ng~L^-1))) + 
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Kd_1ngL_plot
set_palette(Kd_1ngL_plot, "uchicago")
ggsave(filename="R/figs/Kd_1ngL_plot.pdf")

# Trace and main elements plots together ----
Trace_elements_biochar_plot <- ggplot(data = subset(Biochar, Unit %in% "mg/kg"), 
                                      aes(x = Parameter, y = Mean_diffunit, color = biochar),
) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values= c("#767676FF","#800000FF","#FFB547FF"))+
  labs(x = "", 
       y = "mg/kg", 
       color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom",
        text = element_text(size = 16))
Trace_elements_biochar_plot
ggsave("R/figs/Trace_elements_biochar.pdf")

CHON_biochar_plot <- ggplot(data = subset(Biochar, Unit %in% "%"), 
                            aes(x = Parameter, y = Mean_diffunit, color = biochar),
) +
  geom_point() +
  labs(x = "", y = "%") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
CHON_biochar_plot

Main_elements_biochar_plot <- ggplot(data = subset(Biochar, Unit %in% "g/kg"), 
                                      aes(x = Parameter, y = Mean_diffunit, color = biochar),
) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values= c("#767676FF","#800000FF","#FFB547FF"))+
  labs(x = "", 
       y = "g/kg", 
       color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom",
        text = element_text(size = 16))
Main_elements_biochar_plot
ggsave("R/figs/Main_elements_biochar.pdf")


# Main elements separately ----
#need to fix error bars in Biochar_ratios_1ugL_select
Biochar_ratios_1ugL_select$compound <- factor(Biochar_ratios_1ugL_select$compound, levels = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA"))

Kd_1ugL_Ca <- ggplot(data = Biochar_ratios_1ugL_select,
                     aes(x = log10(Ca), y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point(size = 8) +
  # geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.01)+
  # geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) + 
  labs(x = TeX(r'(log Ca $(g~kg^{-1})$)'), y = TeX(r'($log~K_d~(at~C_w~1 \mu g~L^{-1})$)'), color = "", shape = "") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme_bw() +
  #guides(shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_Ca
ggsave(filename="R/figs/Kd_1ugL_Ca.pdf")

Kd_1ugL_Fe <- ggplot(data = subset(Elements_ratios_1ugL, Parameter %in% "Fe"),
                     aes(x = Mean_sameunit, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) + 
  labs(x = "Fe (g/kg)", y = expression(log~K[d]), color = "", shape = "") +
  guides(color = "none", shape = "none") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_Fe
set_palette(Kd_1ugL_Fe, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_Fe.pdf")

Kd_1ugL_C <- ggplot(data = subset(Elements_ratios_1ugL, Parameter %in% "C"),
                    aes(x = Mean_sameunit, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) + 
  labs(x = "C (g/kg)", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_C
set_palette(Kd_1ugL_C, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_C.pdf")

Kd_1ugL_SACO2 <- ggplot(data = SA_PV,
                        aes(x = CO2_SA, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "SA CO2", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  #guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SACO2
set_palette(Kd_1ugL_SACO2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SACO2.pdf")

Kd_1ugL_PVCO2 <- ggplot(data = SA_PV,
                        aes(x = CO2_PV, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.001)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "PV CO2", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVCO2
set_palette(Kd_1ugL_PVCO2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PVCO2.pdf")

Kd_1ugL_SAN2 <- ggplot(data = SA_PV,
                       aes(x = N2_SA, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "SA N2", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  #guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SAN2
set_palette(Kd_1ugL_SAN2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SAN2.pdf")

Kd_1ugL_PVN2 <- ggplot(data = SA_PV,
                       aes(x = log10(N2_PV), y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point(size = 8) +
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.001)+
  labs(x = TeX(r'(log PV $(cm^{3}~g^{-1})$)'), y = TeX(r'($log~K_d~(at~C_w~1 \mu g~L^{-1})$)'), color = "", shape = "") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  theme_bw() +
  guides(shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 30))
Kd_1ugL_PVN2
ggsave(filename="R/figs/Kd_1ugL_PVN2.pdf")


Kd_1ugL_PVN2_Ca <- ggplot(data = SA_PV,
                          aes(x = log10(PVN2_Ca), y = log_Kd, shape = compound, color = biochar), 
) +
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.005)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 8) +
  labs(x = TeX(r'(log PV/Ca)'), y = TeX(r'($log~K_d~(at~C_w~1 \mu g~L^{-1})$)'), color = "", shape = "") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme_bw() +
  #guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVN2_Ca
ggsave(filename="R/figs/Kd_1ugL_PVN2_Ca.pdf")

Kd_1ugL_SA_PV_Ca <- ggplot(data = SA_PV,
                           aes(x = log10(SA_PV_Ca), y = log_Kd, shape = compound, color = biochar), 
) +
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.03)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 8) +
  labs(x = TeX(r'(log (SA/PV)/Ca)'), y = TeX(r'($log~K_d~(at~C_w~1 \mu g~L^{-1})$)'), color = "", shape = "") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SA_PV_Ca
ggsave(filename="R/figs/Kd_1ugL_SA_PV_Ca.pdf")




Biochar_elements <- Biochar_SAPV_elements %>% 
  filter(PV_SA == 0) %>% 
  filter(Key == 1)
write_xlsx(Biochar_elements, "R/data_manipulated/100422_Biochar_elements.xlsx")
Biochar_SAPV <- filter(Biochar, PV_SA == 1)

#Correlation plots C3 ----
SA_PV_C3 <- filter(Biochar_sorption_soil_C3_PFOA, PV_SA == "1")

# Random ----
Elements_ratios <- read_excel("R/data_raw/270322_Elements_ratios.xlsx")
Elements_ratios <- as.data.table(Elements_ratios)
Sorption_BC_single_C3 <- filter(Sorption_BC_single, Conc_point == 3)
Biochar_sorption_soil_C3 <- full_join(x = Sorption_BC_single_C3, y = Biochar, by = c("Biochar" = "biochar"))

Biochar_sorption_soil_C3_PFOA <- filter(Biochar_sorption_soil_C3, Compound == "PFOA")
Biochar_sorption_soil_C3$Compound <- factor(Biochar_sorption_soil_C3$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                "PFOA", "PFNA", "PFDA"))
#Ca
Ca_Kd_corr <- ggplot(data = subset(Biochar_sorption_soil_C3, Parameter %in% "Ca"), aes(x = log(Mean_diffunit), y = log_Kd), group = compound) +
  geom_point(size = 2, aes(shape = Biochar), color = "black") + 
  # geom_smooth(formula = y ~ x, 
  #             method=lm, 
  #             se = FALSE, 
  #             fullrange = FALSE,
  #             color = "grey") +
  labs(x = expression("log [Ca] g/kg"), y = expression(log~K[d]), shape = "") +
  facet_wrap(~ Compound) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Ca_Kd_corr
ggsave(filename = "R/figs/Ca_Kd_corr.pdf")

#Carbon
Carbon_Kd_corr <- ggplot(data = subset(Biochar_sorption_soil_C3, Parameter %in% "C"), aes(x = log(Mean_diffunit), y = log_Kd), group = compound) +
  geom_point(size = 2, aes(shape = Biochar), color = "black") + 
  # geom_smooth(formula = y ~ x, 
  #             method=lm, 
  #             se = FALSE, 
  #             fullrange = FALSE,
  #             color = "grey") +
  labs(x = expression("% C"), y = expression(log~K[d]), shape = "") +
  facet_wrap(~ Compound) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Carbon_Kd_corr
ggsave(filename = "R/figs/Carbon_Kd_corr.pdf")

#Fe
Fe_Kd_corr <- ggplot(data = subset(Biochar_sorption_soil_C3, Parameter %in% "Fe"), aes(x = log(Mean_diffunit), y = log_Kd), group = compound) +
  geom_point(size = 2, aes(shape = Biochar), color = "black") + 
  # geom_smooth(formula = y ~ x, 
  #             method=lm, 
  #             se = FALSE, 
  #             fullrange = FALSE,
  #             color = "grey") +
  labs(x = expression("log [Fe] g/kg"), y = expression(log~K[d]), shape = "") +
  facet_wrap(~ Compound) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Fe_Kd_corr
ggsave(filename = "R/figs/Fe_Kd_corr.pdf")

# SA PV plots ----
#SA PV plots
SA_PV <- filter(Elements_ratios_1ugL, Use == 1)
SA_PV <- filter(Biochar, PV_SA == 1)

SA_PV <- filter(SA_PV, Use == 1)

SA <- ggplot(data = subset(SA_PV, Porosity %in% "SA"),
             aes(x = Mean_diffunit, y = log_Kd, group = Gas, color = Gas, shape = Biochar), 
) +
  geom_point(size = 3) +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  labs(x = expression(m^2~g^-1), y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
SA

PV <- ggplot(data = subset(SA_PV, Porosity %in% "PV"),
             aes(x = Mean_diffunit, y = log_Kd, group = Gas, color = Gas, shape = Biochar), 
) +
  geom_point(size = 3) +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  labs(x = expression(cm^3~g^-1), y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
PV
#So far all parameters follow the trend DSL<ULS<CWC except PV N2

PVN2_Ca <- ggplot(data = Biochar_ratios_1ugL_all,
             aes(x = PVN2_Ca, y = log_Kd, color = compound) 
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm",se = F,
              formula = y ~ x)+
  labs(x = "PVN2/Ca ratio", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
PVN2_Ca

