Biochar <- read_excel("/Users/katinkakrahn/Library/CloudStorage/OneDrive-NGI/VOW/Data/250322_biochar_parameters.xlsx")
as.data.table(Biochar)
Biochar <- as.data.table(Biochar)
setnames(Biochar, "Biochar", "biochar")

#Main and trace elements plot
Elements_biochar <- filter(Biochar, PV_SA == 0)
SA_PV <- filter(Biochar, PV_SA == 1)

Main_elements_1ugL <- Kd_1ugL_select %>% filter(Parameter %in% Elements_biochar)
Main_elements_biochar_plot <- ggplot(data = subset(Elements_biochar, Unit %in% "g/kg"), 
                                aes(x = Parameter, y = Mean_diffunit, color = biochar),
                                ) +
  geom_point(position = position_jitter(w=0.1, h=0)) +
  labs(x = "", y = "g/kg") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
set_palette(Main_elements_biochar_plot, "uchicago")
Main_elements_biochar_plot


Trace_elements_biochar_plot <- ggplot(data = subset(Elements_biochar, Unit %in% "mg/kg"), 
                                     aes(x = Parameter, y = Mean_diffunit, color = biochar),
) +
  geom_point() +
  labs(x = "", y = "mg/kg") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Trace_elements_biochar_plot

CHON_biochar_plot <- ggplot(data = subset(Elements_biochar, Unit %in% "%"), 
                                     aes(x = Parameter, y = Mean_diffunit, color = biochar),
) +
  geom_point() +
  labs(x = "", y = "%") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
CHON_biochar_plot

#Correlation plots
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


#SA PV plots
SA_PV <- filter(Biochar_sorption_soil_C3_PFOA, PV_SA == "1")
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
