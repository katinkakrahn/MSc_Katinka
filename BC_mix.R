library(data.table)
library(ggplot2)
library(psych)
library(readxl)
library(grDevices)
library(dplyr)
library(broom)
library(ggpubr)
library(tidyverse)
library(knitr)
library(plotrix)

Sorption_BC_mix_summary <- Sorption_BC_mix[, .(mean_logCw = mean(log_Cw), 
                                               mean_logCs = mean(log_Cs),
                                               se_logCw = std.error(log_Cw),
                                               se_logCs = std.error(log_Cs),
                                               log_Kd = mean(log_Cs/log_Cw),
                                               se_logKd = std.error(log_Cs/log_Cw)
                                               ),
                                           keyby = .(Compound, Biochar)]

Sorption_BC_mix_summary$Compound <- factor(Sorption_BC_mix_summary$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                      "PFOA", "PFNA", "PFDA"))

BC_mix_Kd <- ggplot(data = Sorption_soil_blank_Kd, aes(x = Compound, y = log_Kd, color = mixLogic)) + 
  geom_point()+ 
  geom_errorbar(aes(ymin=log_Kd-se_logKd, ymax=log_Kd+se_logKd), width = .1) + 
  labs(x = "", y = expression(log~K[d]), col = "") + 
  scale_color_manual(
    values = c('black','grey45'),
    breaks = c("TRUE", "FALSE"),
    labels = c("cocktail", "single compound")
  ) +
  theme_bw()
SoilBlankKd
ggsave(filename="R/figs/SoilBlankKd.pdf")


#CWC facet soil mix
nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

summary_stats_CWC_soil_mix <- data.table(K_F = rep(0, nr_compounds), 
                                         K_F_std_error = rep(0, nr_compounds),
                                         n = rep(0, nr_compounds),
                                         n_std_error = rep(0, nr_compounds),
                                         r_squared = rep(0, nr_compounds),
                                         residual_std_error = rep(0, nr_compounds),
                                         p_value = rep(0, nr_compounds),
                                         compound = compounds,
                                         biochar = "CWC")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = CWC_soil_mix[Compound == compounds[i]])
  summary_stats_CWC_soil_mix[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_CWC_soil_mix[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_CWC_soil_mix[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_CWC_soil_mix[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_CWC_soil_mix[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_CWC_soil_mix[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_CWC_soil_mix[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                     summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_CWC_soil_mix_label <- summary_stats_CWC_soil_mix %>%
  mutate(
    log_Cw = 0.6, log_Cs = 1,
    label =
      glue("*r<sup>2</sup>* = {round(r_squared, 2)} <br> *log K<sub>F</sub>* = {round(K_F, 2)} <br> *n<sub>F</sub>* = {round(n, 2)}")
  )

summary_stats_CWC_soil_mix_label <- summary_stats_CWC_soil_mix_label |>
  transform(Compound = compound)

#CWC Freundlich isotherm plot
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

CWC_single$Compound <- factor(CWC_single$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))
facetCWC$Compound <- factor(facetCWC$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                          "PFOA", "PFNA", "PFDA"))
summary_stats_CWC_single$compound <- factor(summary_stats_CWC_single$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                          "PFOA", "PFNA", "PFDA"))
summary_stats_CWC_single_label$Compound <- factor(summary_stats_CWC_single_label$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                                      "PFOA", "PFNA", "PFDA"))


CWC_facet_isotherm <- ggplot(data = CWC_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), 
             color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_compound), 
              formula = y ~ x, 
              method=lm, 
              se=FALSE, 
              colour = "grey", 
              size = 0.5,
              data = facetCWC
  ) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), 
              color = "black", 
              formula = y ~ x, 
              method=lm, 
              se=T, 
              fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  facet_wrap(~Compound) +
  #ggtitle("CWC isotherm") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none") +
  geom_richtext(
    data = summary_stats_CWC_single_label,
    aes(label = label, x = log_Cw, y = log_Cs),
    hjust = 0
  )
CWC_facet_isotherm
ggsave(filename="R/figs/CWC_facet_isotherm.pdf")