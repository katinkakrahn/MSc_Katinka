Biochar <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/250322_biochar_parameters.xlsx")
as.data.table(Biochar)
Biochar <- as.data.table(Biochar)
setnames(Biochar, "Biochar", "biochar")

#Main and trace elements plot
Elements_biochar <- filter(Biochar, PV_SA == 0)
SA_PV <- filter(Biochar, PV_SA == 1)

Main_elements_biochar_plot <- ggplot(data = subset(Elements_biochar, Unit %in% "g/kg"), 
                                aes(x = Parameter, y = Mean_diffunit, color = biochar),
                                ) +
  geom_point() +
  labs(x = "", y = "g/kg") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
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

#SA PV plots
SA_PV <- filter(Biochar, Model != "0")
SA_PV <- filter(SA_PV, Use == 1)

SA <- ggplot(data = subset(SA_PV, Porosity %in% "SA"), 
                            aes(x = Gas, y = Mean_diffunit, color = biochar)) +
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

#Correlation plots
Sorption_BC_single_C3 <- filter(Sorption_BC_single, Conc_point == 3)
Biochar_sorption_soil_C3 <- full_join(x = Sorption_BC_single_C3, y = Biochar, by = c("Biochar" = "biochar"))


Biochar_sorption_soil_C3$Compound <- factor(Biochar_sorption_soil_C3$Compound, 
                                            levels = c("PFPeA", "PFHxA", "PFHpA","PFOA", "PFNA", "PFDA"))

#Ca
Ca_biochar_corr <- ggplot(data = subset(Biochar_sorption_soil_C3, Parameter %in% "Ca"), aes(x = log(Mean_diffunit), y = log_Kd), group = compound) +
  geom_point(size = 2, aes(shape = Biochar), color = "black") + 
  geom_smooth(formula = y ~ x, 
              method=lm, 
              se = FALSE, 
              fullrange = FALSE,
              color = "grey") +
  labs(x = expression("log [Ca] g/kg"), y = expression(log~K[d]), shape = "") +
  facet_wrap(~ Compound) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Ca_biochar_corr
ggsave(filename = "R/figs/Ca_biochar_corr.pdf")

#Carbon
Carbon_biochar_corr <- ggplot(data = subset(Biochar_sorption_soil_C3, Parameter %in% "C"), aes(x = log(Mean_diffunit), y = log_Kd), group = compound) +
  geom_point(size = 2, aes(shape = Biochar), color = "black") + 
  geom_smooth(formula = y ~ x, 
              method=lm, 
              se = FALSE, 
              fullrange = FALSE,
              color = "grey") +
  labs(x = expression("% C"), y = expression(log~K[d]), shape = "") +
  facet_wrap(~ Compound) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Carbon_biochar_corr
ggsave(filename = "R/figs/Carbon_biochar_corr.pdf")

#Fe
Fe_biochar_corr <- ggplot(data = subset(Biochar_sorption_soil_C3, Parameter %in% "Fe"), aes(x = log(Mean_diffunit), y = log_Kd), group = compound) +
  geom_point(size = 2, aes(shape = Biochar), color = "black") + 
  geom_smooth(formula = y ~ x, 
              method=lm, 
              se = FALSE, 
              fullrange = FALSE,
              color = "grey") +
  labs(x = expression("log [Fe] g/kg"), y = expression(log~K[d]), shape = "") +
  facet_wrap(~ Compound) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Fe_biochar_corr
ggsave(filename = "R/figs/Fe_biochar_corr.pdf")


