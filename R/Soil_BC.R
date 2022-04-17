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
  mutate(BClogic = if_else(Biochar == "no",
                 TRUE,
                 FALSE))

Sorption_soil$SoilLogic <- as.logical(Sorption_soil$Soil_binary)
Sorption_soil$mixLogic <- as.logical(Sorption_soil$mix_binary)
Sorption_soil$isothermLogic <- as.logical(Sorption_soil$isotherm_binary)
Sorption_soil <- subset(Sorption_soil,
                        select = -c(Soil_binary,mix_binary,isotherm_binary))
Sorption_soil <- Sorption_soil %>%
  drop_na(log_Cs) %>% 
  select(-c(V_w, m_tot, m_aq, M_s, M_bc, y)) %>% 
  mutate(Kd = Cs/Cw, 
         log_Kd = log10(Cs/Cw))

Sorption_soil_summary <- Sorption_soil %>% 
  drop_na(log_Cs) %>% 
  filter(Compound != "PFPeA") %>% 
  group_by(mixLogic, BClogic, Biochar, Compound) %>% 
  summarise(K_ds_mean = mean(K_ds),
            log_Kds_mean = log10(mean(K_ds)),
            sd_Kds= sd(K_ds),
            Kd_mean = mean(Kd),
            log_Kd_mean = log10(mean(Kd)),
            sd_Kd = sd(Kd),
            n = n(),
            se_Kd = sd_Kd / sqrt(n),
            se_Kds = sd_Kds / sqrt(n)
            )

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

Sorption <- read_excel("R/data_raw/160222_sorption_rawdata.xlsx")

#Convert 1 and 0 to TRUE and FALSE and delete integer columns
Sorption$SoilLogic <- as.logical(Sorption$Soil_binary)
Sorption$mixLogic <- as.logical(Sorption$mix_binary)
Sorption$isothermLogic <- as.logical(Sorption$isotherm_binary)
Sorption <- select(Sorption, -Soil_binary, -mix_binary, -isotherm_binary)

# Subset biochar and cocktail/single compound
Sorption_NAomit <- na.omit(Sorption)
Sorption_NA_C1omit <- Sorption_NAomit %>% 
  filter(Conc_point != 1) %>% 
  mutate(Kd = Cs/Cw, log_Kd = log10(Cs/Cw),
         BClogic = FALSE)

Soil_BC_join <- full_join(Sorption_NA_C1omit, Sorption_soil)

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

Soil_BC_join$Compound <- 
  factor(Soil_BC_join$Compound, 
         levels = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA"))

Soil_BC_join$type <- 
  factor(Soil_BC_join$type, 
         levels = c("BC_sing", "BC_S_sing", "BC_S_mix", "BC_mix"))


# Attenuation all combinations ----
Attenuation <- Soil_BC_join %>% 
  filter(Biochar != "no") %>% 
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw), 
                           y = log10(Cs),
                           color = type),
             alpha = 0.7) + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            color = type), 
              formula = y ~ x, 
              method=lm, 
              se=F, fullrange = FALSE) + 
  facet_grid(rows = vars(Biochar),
             cols = vars(Compound)) +
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/g)$)'),
       color = "",
  ) +
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

# Attenuation at C10
Soil_BC_join$Compound <- 
  factor(Soil_BC_join$Compound, 
         levels = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA"))

Soil_BC_join$Biochar <- 
  factor(Soil_BC_join$Biochar, 
         levels = c("no", "CWC", "DSL", "ULS"))

Soil_BC_join$type <- 
  factor(Soil_BC_join$type, 
         levels = c("BC_sing", "BC_S_sing", "BC_S_mix", "BC_mix", "S_sing", "S_mix"))

C10 <- Soil_BC_join %>% 
  filter(Conc_point == 10) %>%
  ggplot(aes(x = Biochar,
             y = log10(Kd),
             color = type
  )) +
  geom_point(size = 3,
             alpha = 0.7) +
  facet_wrap(~ Compound) +
  labs(x = "", 
       y = TeX(r'($log~K_{d}~(L/kg)$)'),
       color = "",
  ) +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC single", 
                                "BC soil single", 
                                "BC soil cocktail",
                                "BC cocktail (n=3)",
                                "Soil single (n=3)",
                                "Soil cocktail (n=3)")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 12))
C10
ggsave(filename = "R/figs/C10.pdf")

# PFOA soil isotherm and BC isotherm ----
PFOA_isotherm_attenuation <- Soil_BC_join %>% 
  filter(Compound == "PFOA",
         Biochar != "no",
         mixLogic == FALSE) %>% 
  group_by(SoilLogic, Biochar) %>% 
  ggplot() +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = SoilLogic),
                           size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, 
                            y = log_Cs, 
                            color = SoilLogic), 
              formula = y ~ x, 
              method=lm, 
              se=F, fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/g)$)'),
       color = "") + 
  ggtitle("PFOA") +
  facet_grid(rows = vars(Biochar)) +
  theme_bw() +
  theme(panel.grid = element_blank())
PFOA_isotherm_attenuation

# PFOA summary statistics ----
summary_stats_PFOA <- Soil_BC_join %>% 
  filter(Compound == "PFOA",
         Biochar != "no",
         mixLogic == FALSE,
         SoilLogic == TRUE) %>% 
  group_by(Biochar) %>% 
  do(model = lm(log_Cs ~ log_Cw, data = .)) %>% 
  tidy(model)

summary_stats_PFOA_lm <- lm(log_Cs ~ log_Cw, data = summary_stats_PFOA) %>% 
  get_regression_table() %>% 
  mutate(biochar = "CWC") #need to adjust filter above to the right biochar




