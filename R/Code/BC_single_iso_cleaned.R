# Library ----
library(data.table)
library(readxl)
library(latex2exp)
library(ggtext)
library(scales)
library(writexl)
library(tidyverse)
library(glue)
library(plotrix)
library(moderndive)
library(broom)

# Data manipulation ----
Sorption <- read_excel("R/data_raw/160222_sorption_rawdata.xlsx") %>% 
  na.omit() %>% 
  filter(Conc_point != 1) %>% 
  mutate(SoilLogic = as.logical(Soil_binary),
         mixLogic = as.logical(mix_binary),
         isothermLogic = as.logical(isotherm_binary),
         Kd = Cs / Cw) %>% 
  select(Conc_point, Compound, Biochar, type, SoilLogic, mixLogic, isothermLogic, 
         Cs, Cw, Kd, Ci) 
  
Sorption_BC_single <- subset(Sorption, mixLogic == FALSE)
Sorption_BC_mix <- subset(Sorption, mixLogic == TRUE)

write_xlsx(Sorption, "R/data_manipulated/010422_Sorption_BC.xlsx")
write_xlsx(Sorption_BC_single, "R/data_manipulated/010422_Sorption_BC_single.xlsx")
write_xlsx(Sorption_BC_mix, "R/data_manipulated/010422_Sorption_BC_mix.xlsx")

CWC_single <- filter(Sorption_BC_single, Biochar == "CWC")
ULS_single <- filter(Sorption_BC_single, Biochar == "ULS")
DSL_single <- filter(Sorption_BC_single, Biochar == "DSL")

nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

# Summary statistics ----
summary_BC_single <- Sorption_BC_single %>% 
  group_by(Compound, Biochar) %>%
  do(fit_isotherms = glance(lm(Cs ~ Cw, data = .))) %>% 
  unnest(fit_isotherms)

summary_BC_single_tidy <- Sorption_BC_single %>% 
  group_by(Compound, Biochar) %>%
  do(fit_isotherms = tidy(lm(Cs ~ Cw, data = .))) %>% 
  unnest(fit_isotherms) %>% 
  pivot_wider(names_from = term,
              values_from = c(estimate, std.error, statistic, p.value))

summary_statistics_BC_single <- full_join(summary_BC_single, 
                                          summary_BC_single_tidy) %>% 
  rename(K_F = "estimate_(Intercept)",
         n_F = estimate_Cw
         )

# Sorption BC cocktail ----
Sorption_BC_mix_summary <- Sorption_BC_mix %>% 
  group_by(Biochar, Compound) %>% 
  summarise(Conc_point = 10,
            n = n(),
            se_Cw = sd(Cw) / sqrt(n),
            se_Cs = sd(Cs) / sqrt(n),
            Ci = mean(Ci),
            Cw = mean(Cw),
            Cs = mean(Cs),
            Kd = mean(Cs/Cw),
            se_Kd = sqrt(se_Cw^2+se_Cs^2),
            mixLogic = TRUE, 
            SoilLogic = FALSE,
            type = type)


#CWC ---- 
facetCWC <- CWC_single |>
  transform(pre_compound = Compound)

facetCWC <- rbind(
  transform(facetCWC, Compound = unique(CWC_single$Compound)[1]),
  transform(facetCWC, Compound = unique(CWC_single$Compound)[2]),
  transform(facetCWC, Compound = unique(CWC_single$Compound)[3]),
  transform(facetCWC, Compound = unique(CWC_single$Compound)[4]),
  transform(facetCWC, Compound = unique(CWC_single$Compound)[5]),
  transform(facetCWC, Compound = unique(CWC_single$Compound)[6])
)

facetCWC$Compound <- factor(facetCWC$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                          "PFOA", "PFNA", "PFDA"))

CWC_facet_isotherm <- Sorption_BC_single %>% 
  filter(Biochar == "CWC") %>% 
  mutate(Compound = factor(Compound, 
                           levels = c("PFPeA", "PFHxA", "PFHpA", 
                                      "PFOA", "PFNA", "PFDA")),
         Biochar = factor(Biochar,
                          levels = c("CWC", "DSL", "ULS"))) %>% 
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw), 
                           y = log10(Cs), 
                           group = factor(Compound)), 
             color = "gray45", 
             size = 1) + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            group = pre_compound), 
              formula = y ~ x, 
              method=lm, 
              se=FALSE, 
              colour = "grey", 
              size = 0.5,
              data = facetCWC) +
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            group = factor(Compound)), 
              color = "black", 
              formula = y ~ x, 
              method=lm, 
              se=T, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~Compound) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  theme(panel.grid = element_blank())
CWC_facet_isotherm
ggsave(filename="R/figs/CWC_facet_isotherm.pdf")

# ULS ----
facetULS$Compound <- factor(facetULS$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                          "PFOA", "PFNA", "PFDA"))

ULS_facet_isotherm <- Sorption_BC_single %>% 
  filter(Biochar == "ULS") %>% 
  mutate(Compound = factor(Compound, 
                           levels = c("PFPeA", "PFHxA", "PFHpA", 
                                      "PFOA", "PFNA", "PFDA")),
         Biochar = factor(Biochar,
                          levels = c("CWC", "DSL", "ULS"))) %>% 
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw), 
                           y = log10(Cs), 
                           group = factor(Compound)), 
             color = "gray45", 
             size = 1) + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            group = pre_compound), 
              formula = y ~ x, 
              method=lm, 
              se=FALSE, 
              colour = "grey", 
              size = 0.5,
              data = facetULS) +
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            group = factor(Compound)), 
              color = "black", 
              formula = y ~ x, 
              method=lm, 
              se=T, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~Compound) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  theme(panel.grid = element_blank())
ULS_facet_isotherm
ggsave(filename="R/figs/ULS_facet_isotherm.pdf")

facetULS <- ULS_single |>
  transform(pre_compound = Compound)

facetULS <- rbind(
  transform(facetULS, Compound = unique(ULS_single$Compound)[1]),
  transform(facetULS, Compound = unique(ULS_single$Compound)[2]),
  transform(facetULS, Compound = unique(ULS_single$Compound)[3]),
  transform(facetULS, Compound = unique(ULS_single$Compound)[4]),
  transform(facetULS, Compound = unique(ULS_single$Compound)[5]),
  transform(facetULS, Compound = unique(ULS_single$Compound)[6])
)


# DSL ---- 
facetDSL$Compound <- factor(facetDSL$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                          "PFOA", "PFNA", "PFDA"))

DSL_facet_isotherm <- Sorption_BC_single %>% 
  filter(Biochar == "DSL") %>% 
  mutate(Compound = factor(Compound, 
                           levels = c("PFPeA", "PFHxA", "PFHpA", 
                                      "PFOA", "PFNA", "PFDA")),
         Biochar = factor(Biochar,
                          levels = c("CWC", "DSL", "ULS"))) %>% 
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw), 
                           y = log10(Cs), 
                           group = factor(Compound)), 
             color = "gray45", 
             size = 1) + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            group = pre_compound), 
              formula = y ~ x, 
              method=lm, 
              se=FALSE, 
              colour = "grey", 
              size = 0.5,
              data = facetDSL) +
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            group = factor(Compound)), 
              color = "black", 
              formula = y ~ x, 
              method=lm, 
              se=T, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~Compound) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  theme(panel.grid = element_blank())
DSL_facet_isotherm
ggsave(filename="R/figs/DSL_facet_isotherm.pdf")

facetDSL <- DSL_single |>
  transform(pre_compound = Compound)

facetDSL <- rbind(
  transform(facetDSL, Compound = unique(DSL_single$Compound)[1]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[2]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[3]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[4]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[5]),
  transform(facetDSL, Compound = unique(DSL_single$Compound)[6])
)


# Summary stats ----
summary_stats_CWC_single[, nr_CF2 := 4:9]
summary_stats_ULS_single[, nr_CF2 := 4:9]
summary_stats_DSL_single[, nr_CF2 := 4:9]

summary_stats_single <- merge(summary_stats_CWC_single, summary_stats_ULS_single, 
                              all = TRUE)
summary_stats_single <- merge(summary_stats_single, summary_stats_DSL_single, 
                              all = TRUE)
summary_stats_single$compound <- factor(summary_stats_single$compound, 
                                        levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                   "PFOA", "PFNA", "PFDA"))

write_xlsx(summary_stats_single, "R/data_manipulated/310322_summary_stats_single.xlsx")

# Sorption isotherm all chars ----
Sorption_BC_single$Compound <- factor(Sorption_BC_single$Compound, 
                                      levels = c("PFDA", "PFNA", "PFOA", 
                                                 "PFHpA", "PFHxA","PFPeA"))

Sorption_isotherms <- Sorption_BC_single %>% 
  filter(Compound != "PFPeA" | Biochar != "DSL",
         Compound != "PFPeA" | Biochar != "CWC",
         Compound != "PFHxA" | Biochar != "CWC") %>% 
  # mutate(Compound = recode(Compound,
  #                          "PFPeA" = "PFPeA (C5)",
  #                          "PFHxA" = "PFHxA (C6)",
  #                          "PFHpA" = "PFHpA (C7)",
  #                          "PFOA" = "PFOA (C8)",
  #                          "PFNA" = "PFNA (C9)",
  #                          "PFDA" = "PFDA (C10)")) %>%  
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw), y = log10(Cs), color = factor(Biochar)), 
             size = 1) + 
  geom_smooth(mapping = aes(x = log10(Cw), y = log10(Cs), color = factor(Biochar)), 
              formula = y ~ x, 
              method=lm, 
              se=T, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g~L^{-1})$)'), 
       y = TeX(r'($log~C_{s}~(\mu g~kg^{-1})$)'), 
       color = "") + 
  facet_wrap(~Compound,
             scales = "free_x") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.spacing = unit(0.8, "cm")) +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values=c("#FFB547FF","#4E9C81","#40E0CF"))+
  theme(legend.position = "bottom") 
Sorption_isotherms
ggsave(filename="R/figs/Sorption_isotherms_single_BC.pdf")
ggsave(filename="R/figs/article/Sorption_isotherms_single_BC.pdf")

Sorption_isotherms_nolabel <- Sorption_BC_single %>% 
  filter(Compound != "PFPeA" | Biochar != "DSL",
         Compound != "PFHxA" | Biochar != "CWC",
         Compound != "PFPeA" | Biochar != "CWC") %>% 
  # mutate(Compound = recode(Compound,
  #                          "PFPeA" = "PFPeA (C5)",
  #                          "PFHxA" = "PFHxA (C6)",
  #                          "PFHpA" = "PFHpA (C7)",
  #                          "PFOA" = "PFOA (C8)",
  #                          "PFNA" = "PFNA (C9)",
  #                          "PFDA" = "PFDA (C10)")) %>%  
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw), y = log10(Cs), 
                           color = factor(Biochar)),
             alpha = 0.6,
             size = 2) + 
  geom_smooth(mapping = aes(x = log10(Cw), y = log10(Cs), color = factor(Biochar)), 
              formula = y ~ x, 
              method=lm, 
              se=F, 
              fullrange = FALSE,
              size = 2) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g~L^{-1})$)'), 
       y = TeX(r'($log~C_{s}~(\mu g~kg^{-1})$)'), color = "") + 
  facet_wrap(~Compound,
             scales = "free_x") +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values=c("#FFB547FF","#4E9C81","#40E0CF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom")
  #guides(color = "none")
Sorption_isotherms_nolabel
ggsave(filename="R/figs/SETAC/Sorption_isotherms_nolabel.pdf")

# Individual sorption isotherms all ----
facet_all <- full_join(facetCWC, facetDSL)
facet_all <- full_join(facet_all, facetULS)

BC_facet_isotherm <- Sorption_BC_single %>% 
  mutate(Compound = factor(Compound, 
                           levels = c("PFPeA", "PFHxA", "PFHpA", 
                                      "PFOA", "PFNA", "PFDA")),
         Biochar = factor(Biochar,
                          levels = c("CWC", "DSL", "ULS"))) %>% 
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw), 
                           y = log10(Cs), 
                           group = factor(Compound)), 
             color = "gray45", 
             size = 1) + 
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            group = pre_compound), 
              formula = y ~ x, 
              method=lm, 
              se=FALSE, 
              colour = "grey", 
              size = 0.5,
              data = facet_all) +
  geom_smooth(mapping = aes(x = log10(Cw), 
                            y = log10(Cs), 
                            group = factor(Compound)), 
              color = "black", 
              formula = y ~ x, 
              method=lm, 
              se=T, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_grid(rows = vars(Biochar),
             cols = vars(Compound)) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  theme(panel.grid = element_blank())
BC_facet_isotherm
ggsave(filename="R/figs/BC_facet_isotherm.pdf")



