# Library ----
library(data.table)
library(readxl)
library(latex2exp)
library(ggtext)
library(scales)
library(writexl)
library(broom)
library(tidyverse)
library(knitr)
library(moderndive)

# Data manipulation ----
Sorption_soil <- read_excel("R/data_raw/010322_sorption_rawdata_soil.xlsx") %>% 
  drop_na(log_Cs) %>% 
  mutate(BClogic = if_else(Biochar == "no",
                 TRUE,
                 FALSE),
         SoilLogic = as.logical(Soil_binary),
         mixLogic = as.logical(mix_binary),
         isothermLogic = as.logical(isotherm_binary),
         Kd = Cs / Cw) %>% 
  select(Conc_point, Compound, Biochar, type, SoilLogic, mixLogic, isothermLogic, 
         BClogic, Cs, Cw, K_ds, Kd)
  
write_xlsx(Sorption_soil_summary,"R/data_manipulated/010322_sorption_soil_summary.xlsx")

Sorption_soil_summary$Compound <- 
  factor(Sorption_soil_summary$Compound, 
         levels = c("PFPeA", 
                    "PFHxA", 
                    "PFHpA", 
                    "PFOA", 
                    "PFNA", 
                    "PFDA")
         )

Sorption_BC <- read_excel("R/data_raw/160222_sorption_rawdata.xlsx") %>% 
  na.omit() %>% 
  filter(Conc_point != 1) %>% 
  mutate(SoilLogic = as.logical(Soil_binary),
         mixLogic = as.logical(mix_binary),
         isothermLogic = as.logical(isotherm_binary),
         BClogic = FALSE,
         Kd = Cs / Cw) %>% 
  select(Conc_point, Compound, Biochar, type, SoilLogic, mixLogic, isothermLogic, 
         BClogic, Cs, Cw, Kd)

Soil_BC_join <- full_join(Sorption_BC, Sorption_soil)

# Summary statistics ----
# https://stackoverflow.com/questions/22713325/fitting-several-regression-models-with-dplyr

Soil_BC_join_isotherm <- Soil_BC_join %>%
  filter(isothermLogic == TRUE) %>% 
  group_by(Biochar, Compound, type, Conc_point) %>% 
  mutate(n = n())

summary_statistics_isotherms <- Soil_BC_join_isotherm %>% 
  group_by(Compound, Biochar, type, SoilLogic, mixLogic) %>%
  do(fit_isotherms = glance(lm(Cs ~ Cw, data = .))) %>% 
  unnest(fit_isotherms)

Soil_BC_join_triplicate_mean <- Soil_BC_join %>%
  filter(isothermLogic == FALSE) %>% 
  group_by(Biochar, Compound, type, Conc_point) %>% 
  summarise(n = n(),
            Cw_se = sd(Cw) / sqrt(n),
            Cs_se = sd(Cs) / sqrt(n),
            Kd_se = sqrt(Cw_se^2+Cs_se^2),
            Kd = mean(Kd),
            Cw = mean(Cw),
            Cs = mean(Cs))

Soil_BC_join_mean_and_isotherm <- full_join(Soil_BC_join_triplicate_mean, 
                                            Soil_BC_join_isotherm)

Attenuation_C10 <- Soil_BC_join_mean_and_isotherm %>% 
  filter(Conc_point == 10) %>% 
  mutate(log_Kd = log10(Kd),
         log_Kd_se = log10(Kd_se)) %>% 
  select(Compound, Biochar, type, log_Kd, log_Kd_se, n) %>% 
  group_by(Compound, Biochar, type) %>% 
  arrange(Biochar = factor(Biochar, levels = c("ULS", "DSL", "CWC", "no")),
          type = factor(type, levels = c("BC_sing", "BC_S_sing", 
                                         "BC_S_mix", "BC_mix")),
          Compound = factor(Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                 "PFOA", "PFNA", "PFDA"))
          ) %>% 
  write_xlsx("R/data_manipulated/190422_Attenuation_factors.xlsx")

# Attenuation all combinations ----
Attenuation <- Soil_BC_join_mean_and_isotherm %>% 
  filter(Biochar != "no") %>% 
  mutate(Compound = factor(Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                      "PFOA", "PFNA", "PFDA")),
         type = factor(type, levels = c("BC_sing", "BC_S_sing", 
                                  "BC_S_mix", "BC_mix"))) %>% 
  ggplot(mapping = aes(x = log10(Cw), 
                       y = log10(Cs),
                       color = type)) +
  # geom_errorbar(aes(xmin = log10(Cw) - log10(Cw_sd),
  #                   xmax = log10(Cw) + log10(Cw_sd),
  #                   )) +
  # geom_errorbar(aes(ymin = log10(Cs) - log10(Cs_sd),
  #                   ymax = log10(Cs) + log10(Cs_sd))) +
  geom_point(alpha = 0.7) + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            color = type), 
              formula = y ~ x, 
              method=lm, 
              se=F, fullrange = FALSE) + 
  facet_grid(rows = vars(Biochar),
             cols = vars(Compound)) +
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'),
       color = "") +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil cocktail",
                                "BC cocktail (n=3)")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 12)
        )
Attenuation
ggsave(filename = "R/figs/Attenuation.pdf")
# Can choose to have triplicate points as triplicates or as mean point (but then error bars become an issue)

# Attenuation at C10 ----
C10 <- Soil_BC_join_mean_and_isotherm %>% 
  filter(Conc_point == 10) %>%
  mutate(Compound = factor(Compound, 
                           levels = c("PFPeA", "PFHxA", "PFHpA", 
                                      "PFOA", "PFNA", "PFDA")),
         Biochar = factor(Biochar,
                          levels = c("no", "CWC", "DSL", "ULS")),
         type = factor(type,
                       levels = c("BC_sing", "BC_S_sing", "BC_S_mix", 
                                  "BC_mix", "S_sing", "S_mix"))) %>% 
  ggplot(aes(x = Biochar,
             y = log10(Kd),
             color = type,
             shape = type
  )) +
  geom_point(size = 4,
             alpha = 1) +
  facet_wrap(~ Compound) +
  labs(x = "Biochar", 
       y = TeX(r'($log~K_{d}~(L/kg)$)'),
       color = "",
       shape = ""
  ) +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil cocktail",
                                "BC cocktail (n=3)",
                                "Soil single (n=3)",
                                "Soil cocktail (n=3)")) +
  scale_shape_manual(name = "",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil cocktail",
                                "BC cocktail (n=3)",
                                "Soil single (n=3)",
                                "Soil cocktail (n=3)"),
                     values = c(16, 17, 17, 16, 17, 17)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20))
C10
ggsave(filename = "R/figs/C10.pdf")

Attenuation_factors <- read_xlsx("R/data_manipulated/190422_Attenuation_factors_manual.xlsx") %>% 
  mutate(Attenuation_percent = Attenuation * 100)
  
Attenuation_C10 <- Attenuation_factors %>% 
  drop_na(Attenuation) %>% 
  mutate(Compound = factor(Compound, 
                           levels = c("PFPeA", "PFHxA", "PFHpA", 
                                      "PFOA", "PFNA", "PFDA")),
         Biochar = factor(Biochar,
                          levels = c("CWC", "DSL", "ULS", "no_CWC", "no_DSL", "no_ULS"))
         ) %>% 
  ggplot(aes(x = Biochar,
             y = Attenuation_percent,
             color = type,
             shape = type
  )) +
  geom_point(size = 4,
             alpha = 1) +
  facet_wrap(~ Compound) +
  labs(x = "Biochar",
       y = "Attenuation factor (% sorption reduction from BC single)",
       color = "",
       shape = ""
  ) +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC cocktail (n=3)", 
                                "BC soil cocktail",
                                "BC soil single",
                                "Soil cocktail (n=3)",
                                "Soil single (n=3)")) +
  scale_shape_manual(name = "",
                     labels = c("BC cocktail (n=3)", 
                                "BC soil cocktail",
                                "BC soil single",
                                "Soil cocktail (n=3)",
                                "Soil single (n=3)"),
                     values = c(16, 17, 17, 17, 17)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 12))
Attenuation_C10
ggsave(filename = "R/figs/Attenuation_factors_C10.pdf")

# PFOA soil isotherm and BC isotherm ----
PFOA_isotherm_attenuation <- Soil_BC_join_isotherm %>% 
  filter(Compound == "PFOA") %>% 
  mutate(type = factor(type, levels = c("BC_sing", "BC_S_sing", 
                                        "BC_S_mix"))) %>% 
  group_by(type) %>% 
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw), y = log10(Cs), color = type),
                           size = 1) + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            color = type), 
              formula = y ~ x, 
              method=lm, 
              se=F, fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'),
       color = "") + 
  facet_grid(rows = vars(Biochar),
             cols = vars(Compound)) +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil cocktail")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20)
  )
PFOA_isotherm_attenuation
ggsave(filename = "R/figs/Attenuation_isotherms_PFOA.pdf")
