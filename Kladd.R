Soil_BC_join_triplicate_mean <- Soil_BC_join %>%
  filter(isothermLogic == FALSE) %>% 
  group_by(Biochar, Compound, SoilLogic, mixLogic, isothermLogic, BClogic, Conc_point, type) %>% 
  summarise(Cw_sd = sd(Cw),
            Cs_sd = sd(Cs),
            Kd_sd = sqrt(Cw_sd^2+Cs_sd^2),
            Kd = mean(Kd),
            Cw = mean(Cw),
            Cs = mean(Cs),
            log_Kd = log10(mean(Kd)),
            log_Kd_sd = log10(sqrt(Cw_sd^2+Cs_sd^2)),
            n = n())

Soil_BC_join_isotherm <- Soil_BC_join %>%
  filter(isothermLogic == TRUE)

Soil_BC_join <- full_join(Soil_BC_join_triplicate_mean, Soil_BC_join_isotherm)

Soil_dummy <- Soil_BC_join %>% 
  filter(Biochar == "no") |>
  transform(pre_compound = Compound,
            pre_biochar = Biochar)

Soil_dummy <- rbind(
  transform(Soil_dummy, Compound = unique(Soil_BC_join$Compound)[1]),
  transform(Soil_dummy, Compound = unique(Soil_BC_join$Compound)[2]),
  transform(Soil_dummy, Compound = unique(Soil_BC_join$Compound)[3]),
  transform(Soil_dummy, Compound = unique(Soil_BC_join$Compound)[4]),
  transform(Soil_dummy, Compound = unique(Soil_BC_join$Compound)[5]),
  transform(Soil_dummy, Compound = unique(Soil_BC_join$Compound)[6])
)

Soil_BC_join_biochars <- Soil_BC_join %>% 
  filter(Biochar == "no")

Soil_dummy2 <- Soil_BC_join %>% 
  filter(Biochar == "no") |>
  transform(pre_biochar = Biochar)

Soil_dummy3 <- rbind(
  transform(Soil_dummy2, Biochar = unique(Soil_BC_join_biochars$Biochar)[1]),
  transform(Soil_dummy2, Biochar = unique(Soil_BC_join_biochars$Biochar)[2]),
  transform(Soil_dummy2, Biochar = unique(Soil_BC_join_biochars$Biochar)[3])
)

Soil_BC_join$Compound <- 
  factor(Soil_BC_join$Compound, 
         levels = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA"))

Soil_BC_join$type <- 
  factor(Soil_BC_join$type, 
         levels = c("BC_sing", "BC_S_sing", "BC_S_mix", "BC_mix"))

Attenuation <- Soil_BC_join %>% 
  filter(Biochar != "no") %>% 
  ggplot() +
  # geom_errorbar(data = subset(Soil_BC_join,
  #                             !is.na(Kd_sd)),
  #               xmin = log10(Cw) - log10(Cw_sd),
  #               xmax = log10(Cw) + log10(Cw_sd),
  #               ymin = log10(Cs) - log10(Cs_sd),
  #               ymax = log10(Cs) + log10(Cs_sd)) +
  geom_point(mapping = aes(x = log10(Cw), 
                           y = log10(Cs),
                           color = type)) + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            color = type), 
              formula = y ~ x, 
              method=lm, 
              se=F, fullrange = FALSE) + 
  # geom_point(data = Soil_dummy,
  #            aes(x = log10(Cw),
  #                y = log10(Cs),
  #                group = pre_compound)) + 
  #facet_wrap(~ Compound + Biochar) +
  facet_grid(rows = vars(Biochar),
             cols = vars(Compound)) +
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/g)$)'),
       color = "",
       ) +
  scale_color_discrete(labels = c("BC single", 
                                 "BC soil single", 
                                 "BC soil cocktail",
                                 "BC cocktail")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")
Attenuation
# Can choose to have triplicate points as triplicates or as mean point (but then error bars become an issue)

C10 <- Soil_BC_join %>% 
  filter(Conc_point == 10,
         Biochar != "no") %>%
  ggplot(aes(x = log10(Cw),
             y = log10(Cs),
             color = type
  )) +
  geom_point() +
  facet_wrap(~ Compound)
C10



# filter(Compound != "PFPeA" | Biochar == "DSL",
#        Compound != "PFHxA" | Biochar == "DSL",
#        Compound != "PFHpA" | Biochar == "DSL"
# ) %>%


