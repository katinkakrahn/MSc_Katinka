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
  do(fit_isotherms = glance(lm(log10(Cs) ~ log10(Cw), data = .))) %>% 
  unnest(fit_isotherms)

Freundlich_statistics_isotherms <- Soil_BC_join_isotherm %>% 
  group_by(Compound, Biochar, type, SoilLogic, mixLogic) %>%
  do(fit_isotherms = tidy(lm(log10(Cs) ~ log10(Cw), data = .))) %>% 
  unnest(fit_isotherms) %>% 
  pivot_wider(names_from = "term",
              values_from = "estimate") %>% 
  rename(log_KF = "(Intercept)",
         n = "log10(Cw)")

Freundlich_KF <- Freundlich_statistics_isotherms %>% 
  drop_na(log_KF) %>% 
  select(-n)

Freundlich_n <- Freundlich_statistics_isotherms %>% 
  drop_na(n) %>% 
  select(-log_KF)

Freundlich_soil <- full_join(Freundlich_KF, Freundlich_n, by = c("type", "Compound", "Biochar"))
r_squared <- summary_statistics_isotherms %>% 
  select(Compound, Biochar, type, r.squared,p.value)
Freundlich_soil <- merge(Freundlich_soil, r_squared, by = c("Compound", "Biochar", "type")) %>% 
  filter(Compound %in% c("PFOA"),
         type %in% c("BC_S_mix", "BC_S_sing", "BC_sing")) %>% 
  write_xlsx("R/data_manipulated/280422_PFOA_Freundlich.xlsx")

  
PFOA_linear <- Soil_BC_join_isotherm %>% 
  filter(Compound == "PFOA") %>% 
  mutate(type = factor(type, levels = c("BC_sing", "BC_S_sing", 
                                          "BC_S_mix")),
         Cw = Cw/1000,
         Cs = Cs/1000) %>% 
  group_by(type) %>% 
  ggplot(mapping = aes(x = Cw, y = Cs, color = type)) +
  geom_point(size = 3) + 
  geom_smooth(method = "lm", 
              formula = y ~ poly(log(x), 2), 
              se = FALSE) +
  labs(x = TeX(r'($C_{w}~(mg/L)$)'), 
         y = TeX(r'($C_{s}~(mg/kg)$)'),
         color = "") + 
  facet_wrap(~Biochar,
               scales = "free_x") +
  scale_y_continuous(labels = comma) +
  scale_color_brewer(palette = "Paired",
                       labels = c("BC single", 
                                  "BC soil single", 
                                  "BC soil mixed")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20),
        panel.spacing = unit(2, "lines"))
PFOA_linear
ggsave(filename = "R/figs/PFOA_linear.pdf")  


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

TEST <- merge(summary_statistics_BC_single,
              Soil_BC_join_isotherm)
Soil_BC_join_mean_and_isotherm <- full_join(Soil_BC_join_triplicate_mean, 
                                            Soil_BC_join_isotherm)

Attenuation_C10 <- Soil_BC_join_mean_and_isotherm %>% 
  filter(Conc_point == 10) %>% 
  mutate(log_Kd = log10(Kd),
         log_Kd_se = log10(Kd_se),
         Kd = Kd) %>% 
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
  geom_point(alpha = 1) + 
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
ggsave(filename = "R/figs/Attenuation.pdf")
# Can choose to have triplicate points as triplicates or as mean point (but then error bars become an issue)

# Attenuation at C10 ----
C10 <- Soil_BC_join_mean_and_isotherm %>% 
  filter(Conc_point == 10,
         Compound %in% c("PFOA", "PFNA", "PFDA")) %>%
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
  geom_point(size = 6,
             alpha = 1) +
  facet_wrap(~ Compound) +
  labs(x = "", 
       y = TeX(r'($log~K_{d}~(L/kg)$)'),
       color = "",
       shape = ""
  ) +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil mixed",
                                "BC mixed (n=3)",
                                "Soil single (n=3)",
                                "Soil mixed (n=3)")) +
  scale_shape_manual(name = "",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil mixed",
                                "BC mixed (n=3)",
                                "Soil single (n=3)",
                                "Soil mixed (n=3)"),
                     values = c(16, 17, 17, 16, 17, 17)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 25),
        axis.text.x = element_text(angle = 0, 
                                   vjust = 0.5, 
                                   hjust=0.5))
C10
ggsave(filename = "R/figs/C10.pdf")

Attenuation_factors <- read_xlsx("R/data_manipulated/190422_Attenuation_factors_manual.xlsx") %>% 
  mutate(percent_drop = (1-Attenuation)*100) %>% 
  filter(Compound %in% c("PFOA", "PFNA", "PFDA"))

Attenuation_C10_OND <- Attenuation_factors %>% 
  drop_na(Attenuation) %>% 
  filter(Biochar %in% c("CWC", "DSL", "ULS"),
         Compound %in% c("PFOA", "PFNA", "PFDA"),
         ) %>% 
  mutate(Compound = factor(Compound, 
                           levels = c("PFOA", "PFNA", "PFDA")),
         Biochar = factor(Biochar,
                          levels = c("CWC", "DSL", "ULS")),
         type = factor(type,
                       levels = c("BC_sing", "BC_S_sing", 
                                  "BC_S_mix", "BC_mix",
                                  "S_sing", "S_mix"
                                  )) 
         ) %>% 
  ggplot(aes(x = Biochar,
             y = Attenuation,
             color = type,
             shape = type
  )) +
  geom_point(size = 6,
             alpha = 0.7) +
  facet_wrap(~ Compound) +
  labs(x = "",
       y = "AF",
       color = "",
       shape = ""
  ) +
  #scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC single (ref.)",
                                "BC soil single",
                                "BC soil mixed",
                                "BC mixed (n=3)",
                                "Soil single (n=3)",
                                "Soil mixed (n=3)"
                                )) +
  scale_shape_manual(name = "",
                     labels = c("BC single (ref.)",
                                "BC soil single",
                                "BC soil mixed",
                                "BC mixed (n=3)",
                                "Soil single (n=3)",
                                "Soil mixed (n=3)"
                                ),
                     values = c(16, 17, 17, 16)) +
  scale_y_continuous(breaks=c(25,50,75,100,125,150)) +
  guides(color=guide_legend(nrow=2),byrow=TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 25),
        axis.text.x = element_text(angle = 0, 
                                   vjust = 0.5, 
                                   hjust=0.5))
Attenuation_C10_OND
ggsave(filename = "R/figs/Attenuation_factors_C10_OND.pdf")

# Attenuation C10 PFOA + PFNA ----
C10_PFOA_PFNA <- Soil_BC_join_mean_and_isotherm %>% 
  filter(Conc_point == 10,
         Compound %in% c("PFOA", "PFNA")) %>%
  mutate(Compound = factor(Compound, 
                           levels = c("PFOA", "PFNA")),
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
  geom_point(size = 6,
             alpha = 1) +
  facet_wrap(~ Compound) +
  labs(x = "", 
       y = TeX(r'($log~K_{d}~(L/kg)$)'),
       color = "",
       shape = ""
  ) +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil mixed",
                                "BC mixed (n=3)",
                                "Soil single (n=3)",
                                "Soil mixed (n=3)")) +
  scale_shape_manual(name = "",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil mixed",
                                "BC mixed (n=3)",
                                "Soil single (n=3)",
                                "Soil mixed (n=3)"),
                     values = c(16, 17, 17, 16, 17, 17)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        text = element_text(size = 20))
C10_PFOA_PFNA
ggsave(filename = "R/figs/SETAC/C10_PFOA_PFNA.pdf")

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
                                "BC soil mixed")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20)
  )
PFOA_isotherm_attenuation
ggsave(filename = "R/figs/Attenuation_isotherms_PFOA.pdf")
