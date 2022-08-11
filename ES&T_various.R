Sorption_BC <- read_excel("R/data_raw/160222_sorption_rawdata.xlsx") %>% 
  na.omit() %>% 
  filter(Conc_point != 1) %>% 
  mutate(SoilLogic = as.logical(Soil_binary),
         mixLogic = as.logical(mix_binary),
         isothermLogic = as.logical(isotherm_binary),
         Kd = Cs / Cw) %>% 
  filter(mixLogic == FALSE) %>% 
  select(Conc_point, Compound, Biochar, type, SoilLogic, mixLogic, isothermLogic, 
         Cs, Cw, Kd, Ci) 

Sorption_soil <- read_excel("R/data_raw/010322_sorption_rawdata_soil.xlsx") %>% 
  drop_na(log_Cs) %>% 
  mutate(BClogic = if_else(Biochar == "no",
                           TRUE,
                           FALSE),
         SoilLogic = as.logical(Soil_binary),
         mixLogic = as.logical(mix_binary),
         isothermLogic = as.logical(isotherm_binary),
         Kd = Cs / Cw) %>% 
  filter(Biochar != "no",
        # Compound %in% c("PFOA", "PFNA", "PFDA")
         ) %>% 
  select(Conc_point, Compound, Biochar, type, SoilLogic, mixLogic, isothermLogic, 
         BClogic, Cs, Cw, K_ds, Kd)

Sorption_BC_soil <- full_join(Sorption_BC, Sorption_soil)

# Summary statistics ----
regression_glance <- Sorption_BC_soil %>% 
  group_by(Compound, Biochar, type) %>%
  do(fit_isotherms = glance(lm(log10(Cs) ~ log10(Cw), data = .))) %>% 
  unnest(fit_isotherms)

regression_tidy <- Sorption_BC_soil %>% 
  group_by(Compound, Biochar, type) %>%
  do(fit_isotherms = tidy(lm(log10(Cs) ~ log10(Cw), data = .))) %>% 
  unnest(fit_isotherms) %>% 
  pivot_wider(names_from = term,
              values_from = c(estimate, std.error, statistic, p.value))

regression_statistics_BC_soil <- full_join(regression_glance, 
                                          regression_tidy) %>% 
  rename(K_F = "estimate_(Intercept)",
         n_F = "estimate_log10(Cw)",
         K_F_se = "std.error_(Intercept)",
         n_F_se = "std.error_log10(Cw)"
  ) %>% 
  select(Compound, Biochar, type, r.squared, p.value, K_F, n_F, K_F_se, n_F_se)

write_xlsx(regression_statistics_BC_soil, "R/data_manipulated/150622_Freundlich_coefficients.xlsx")

# Sorption isotherm plots ----
Sorption_isotherms_all <- Sorption_BC_soil %>% 
  filter(Compound != "PFPeA" | Biochar != "DSL",
         Compound != "PFPeA" | Biochar != "CWC",
         Compound != "PFHxA" | Biochar != "CWC") %>% 
  mutate(factor(Compound, 
                levels = c("PFDA", "PFNA", "PFOA", 
                           "PFHpA", "PFHxA","PFPeA"))) %>% 
  group_by(Biochar, Compound, type) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
  geom_point(mapping = aes(x = log10(Cw), 
                           y = log10(Cs), 
                           color = factor(Biochar),
                           shape = type
  ), 
  #color = "grey",
  size = 1) + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            color = factor(Biochar),
                            linetype = type), 
              formula = y ~ x, 
              method=lm, 
              se=F, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~Compound,
             scales = "free",
             ncol = 2) +
  theme_bw() +
  theme(
    #text = element_text(size = 16),
    panel.spacing = unit(0.2, "cm"),
    legend.position = "bottom",
    panel.grid = element_blank(),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,0,-10)) +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values=c("#FFB547FF","#4E9C81","#40E0CF"))
Sorption_isotherms_all




Attenuation <- Sorption_BC_soil %>% 
  filter(Compound != "PFPeA" | Biochar != "DSL",
         Compound != "PFPeA" | Biochar != "CWC",
         Compound != "PFHxA" | Biochar != "CWC") %>% 
  mutate(Compound = factor(Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                "PFOA", "PFNA", "PFDA")),
         type = factor(type, levels = c("BC_sing", "BC_S_sing", 
                                        "BC_S_mix"))) %>% 
  ggplot(mapping = aes(x = log10(Cw), 
                       y = log10(Cs),
                       color = type)) +
  geom_point(size = 0.5,
             color = "grey") + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            color = type), 
              formula = y ~ x, 
              method=lm, 
              se=F, fullrange = FALSE,
              alpha = 0.5) + 
  facet_grid(rows = vars(Biochar),
             cols = vars(Compound)) +
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'),
       color = "") +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil mixed",
                                "BC mixed (n=3)")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 12)
  )
Attenuation

              






