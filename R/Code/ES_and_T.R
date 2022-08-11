# Library ----
library(data.table)
library(readxl)
library(latex2exp)
library(ggtext)
library(scales)
library(writexl)
library(tidyverse)
library(ggtext)
library(glue)
library(plotrix)
library(moderndive)
library(broom)
library(knitr)


# Data manipulation sorption isotherms ----
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
         #Compound %in% c("PFOA", "PFNA", "PFDA")
         ) %>% 
  select(Conc_point, Compound, Biochar, type, SoilLogic, mixLogic, isothermLogic, 
         BClogic, Cs, Cw, K_ds, Kd)

Sorption_BC_soil <- full_join(Sorption_BC, Sorption_soil)

Sorption_BC_soil %>% filter(isothermLogic == F) %>% write_xlsx("R/data_manipulated/090822_rawdata_joined.xlsx")

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
  select(Compound, Biochar, type, r.squared, p.value, K_F, n_F, K_F_se, n_F_se, nobs)

write_xlsx(regression_statistics_BC_soil, "R/data_manipulated/150622_Freundlich_coefficients.xlsx")

# Sorption isotherm plots ----
Sorption_isotherms_nonsigremoved <- Sorption_BC_soil %>% 
  filter(SoilLogic == FALSE,
         Compound != "PFPeA" | Biochar != "DSL",
         Compound != "PFPeA" | Biochar != "CWC",
         Compound != "PFHxA" | Biochar != "CWC") %>% 
  mutate(Compound = factor(Compound, 
                           levels = c("PFDA", "PFNA", "PFOA", 
                                      "PFHpA", "PFHxA","PFPeA")
                           )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
  geom_point(mapping = aes(x = log10(Cw), 
                           y = log10(Cs), 
                           color = factor(Biochar)
                           ), 
                           #color = "grey",
             size = 3,
             alpha = 0.6) + 
  geom_smooth(mapping = aes(x = log10(Cw), y = log10(Cs), color = factor(Biochar)), 
              formula = y ~ x, 
              method=lm, 
              se=F, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~Compound,
             #scales = "free_x",
             ncol = 2) +
  theme_bw() +
  theme(text = element_text(size = 30),
    #legend.text=element_text(size=20),
    panel.spacing = unit(0.2, "cm"),
    legend.position = "bottom",
    panel.grid = element_blank(),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,0,-10),
    axis.title.x = element_text(margin = margin(10,0,5,0))) +
  scale_color_manual(labels = c("WCBC", "SSBC1", "SSBC2"),
                     breaks = c("CWC", "ULS", "DSL"),
                     values = c("#FFB547FF","#4E9C81","#40E0CF")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
Sorption_isotherms_nonsigremoved
ggsave(filename="R/figs/article/Sorption_isotherms_single_nonsigremoved.pdf")

# Data cleaning ----
PSD <- read_xlsx("R/data_raw/080422_PZD.xlsx") 
PSD <- as.data.table(PSD)
PSD$SA_PV <- PSD$SA/PSD$PV
PSD$SA_PV_C <- PSD$SA_PV/PSD$C

# Pore size distribution ----
PSD_pivot <- PSD %>% 
  mutate(id = row_number()) %>% 
  select(SA, PV, Pore_size, Biochar, Gas) %>%
  pivot_longer(c(where(is.numeric), -Biochar, -Pore_size))

PSD_pivot$name <- factor(PSD_pivot$name,
                         levels = c("SA","PV"),
                         labels = c("SA~(m^2/g)", "PV~(cm^3/g)"))
PSD_pivot$name <- paste(PSD_pivot$Gas, "-", PSD_pivot$name)
PSD_pivot$name <- factor(PSD_pivot$name,
                         levels = c("N2 - SA~(m^2/g)", "N2 - PV~(cm^3/g)", "CO2 - SA~(m^2/g)", "CO2 - PV~(cm^3/g)"),
                         labels = c("N[2]~SA~(m^2/g)", "N[2]~PV~(cm^3/g)", "CO[2]~SA~(m^2/g)", "CO[2]~PV~(cm^3/g)"))

PSD_plot <- PSD_pivot %>%
  mutate(name = factor(name,
                       levels = c("CO[2]~SA~(m^2/g)", "CO[2]~PV~(cm^3/g)", "N[2]~SA~(m^2/g)", "N[2]~PV~(cm^3/g)"
                       )))%>%
  ggplot(aes(
    y = value,
    x = Pore_size,
    color = Biochar
  )) +
  geom_vline(xintercept = c(0.96, 1.08, 1.19, 1.36, 1.54, 1.42),
             linetype = 2,
             color = "gray45") +
  labs(x = "Pore diameter (nm)", y = NULL, color = "", shape = "") +
  geom_point(size = 3) +
  scale_x_log10() +
  facet_wrap(.~ name,
             scales = "free",
             labeller = label_parsed,
             strip.position = "left") +
  scale_color_manual(labels = c("WCBC", "SSBC1", "SSBC2"),
                     breaks = c("CWC", "ULS", "DSL"),
                     values=c("#FFB547FF","#4E9C81","#40E0CF")) +
  #scale_x_continuous(breaks=c(1,3,5,10,20,30)) +
  scale_y_continuous() +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 40),
        strip.placement = "outside",
        strip.background = element_blank(),
        legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10, 10, -10),
        panel.grid = element_blank(),
        panel.spacing = unit(0.5, "cm"),
        axis.title.x = element_text(margin = margin(10,0,5,0))
  ) +
  guides(colour = guide_legend(override.aes = list(size=7)))
PSD_plot
ggsave("R/figs/article/PSD_plot.pdf")

# Attenuation ----
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

Soil_BC_join_isotherm <- Soil_BC_join %>%
  filter(isothermLogic == TRUE) %>% 
  group_by(Biochar, Compound, type, Conc_point) %>% 
  mutate(n = n(),
         log_Kd = log10(Kd))

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

Freundlich_coefficients <- full_join(Freundlich_KF, Freundlich_n, by = c("type", "Compound", "Biochar"))
r_squared <- summary_statistics_isotherms %>% 
  select(Compound, Biochar, type, r.squared,p.value)

Soil_BC_join_triplicate_mean <- Soil_BC_join %>%
  filter(isothermLogic == FALSE) %>% 
  group_by(Biochar, Compound, type, Conc_point, SoilLogic, isothermLogic, BClogic, mixLogic) %>% 
  summarise(Kd_mean = mean(Kd),
            log_Kd = mean(log10(Kd)),
            log_Kd_sd = sd(log10(Kd)),
            n = n())

  write_xlsx(Soil_BC_join_triplicate_mean, "R/data_manipulated/090822_triplicate_tests.xlsx")

Soil_BC_join_mean_and_isotherm <- full_join(Soil_BC_join_triplicate_mean, 
                                            Soil_BC_join_isotherm)

C10 <- Soil_BC_join_mean_and_isotherm %>% 
  filter(Conc_point == 10,
         Compound %in% c("PFOA", "PFNA", "PFDA"),
         Biochar != "no"#,
         #type != "BC_mix"
  ) %>%
  mutate(Compound = factor(Compound, 
                           levels = c("PFPeA", "PFHxA", "PFHpA", 
                                      "PFOA", "PFNA", "PFDA")),
         type = factor(type,
                       levels = c("BC_sing", "BC_S_sing", "BC_S_mix", 
                                  "BC_mix", "S_sing", "S_mix"))) %>% 
  ggplot(aes(x = Biochar,
             y = log_Kd,
             color = type,
             fill = type,
             shape = type
  )) +
  geom_point(size = 7,
             alpha = 1) +
  geom_errorbar(aes(ymin=log_Kd-log_Kd_sd, 
                    ymax=log_Kd+log_Kd_sd),
                width = 0.3,
                show.legend = F) +
  facet_wrap(~ Compound) +
  labs(x = "", 
       y = TeX(r'($log~K_{d}~(L/kg)$)'),
       color = "",
       shape = "",
       fill = ""
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("BC-single",
                                "BC-soil-single",
                                "BC-soil-mix",
                                "BC-mix (n=3)")) +
  scale_fill_manual(values = c("#262239","#6bb120", "#009df7", "grey45"),
                     labels = c("BC-single",
                                "BC-soil-single",
                                "BC-soil-mix",
                                "BC-mix (n=3)")) +
  scale_color_manual(values = c("black", "black", "black", "black"),
                    labels = c("BC-single",
                               "BC-soil-single",
                               "BC-soil-mix",
                               "BC-mix (n=3)")) +
  scale_x_discrete(limit = c("CWC", "DSL", "ULS"),
                   labels = c("WCBC","SSBC2","SSBC1")) +
  theme_bw() +
  theme(text = element_text(size = 30),
        panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0, 
                                   vjust = 0.5, 
                                   hjust=0.5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,10,-10),
        legend.text=element_text(size=20)) #+
  #guides(colour = guide_legend(override.aes = list(size=4)))
C10
ggsave(filename = "R/figs/article/C10_grey_triplicateremoved.pdf")
