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
library(latex2exp)
library(glue)
library(ggtext)

Sorption <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/160222_sorption_rawdata.xlsx")
Sorption <- as.data.table(Sorption)
#old <- c("Ci_(ug/L)", "Cw_(ug/L)", "Cs_(ug/kg)")
#new <- c("Ci", "Cw", "Cs")
#setnames(Sorption, old, new, skip_absent = TRUE)

#Convert 1 and 0 to TRUE and FALSE and delete integer columns
Sorption$SoilLogic <- as.logical(Sorption$Soil_binary)
Sorption$mixLogic <- as.logical(Sorption$mix_binary)
Sorption <- subset(Sorption,select = -c(Soil_binary,mix_binary))
Sorption_BC <- kable(Sorption, "latex", booktabs = TRUE, digits = 2)

# Subset biochar and cocktail/single compound
Sorption_NAomit <- na.omit(Sorption)
Sorption_NA_C1omit <- Sorption_NAomit %>% slice(-c(20, 30, 40, 50, 79, 88, 98, 108, 118, 157, 167, 177))
Sorption_BC_single <- subset(Sorption_NA_C1omit, mixLogic == FALSE)
Sorption_BC_mix <- subset(Sorption_NA_C1omit, mixLogic == TRUE)

CWC_single <- filter(Sorption_BC_single, Biochar == "CWC")
ULS_single <- filter(Sorption_BC_single, Biochar == "ULS")
DSL_single <- filter(Sorption_BC_single, Biochar == "DSL")

#Summary statistics CWC
nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

summary_stats_CWC_single <- data.table(K_F = rep(0, nr_compounds), 
                                       K_F_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "CWC")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = CWC_single[Compound == compounds[i]])
  summary_stats_CWC_single[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_CWC_single[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_CWC_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_CWC_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_CWC_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_CWC_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_CWC_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_CWC_single_label <- summary_stats_CWC_single %>%
  mutate(
    log_Cw = 1, log_Cs = 0.8,
    label =
      glue("*r<sup>2</sup>* = {round(r_squared, 2)} <br> *log K<sub>F</sub>* = {round(K_F, 2)} <br> *n<sub>F</sub>* = {round(n, 2)}")
    )

summary_stats_CWC_single_label <- summary_stats_CWC_single_label |>
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



#Summary statistics ULS
nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

summary_stats_ULS_single <- data.table(K_F = rep(0, nr_compounds), 
                                       K_F_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "ULS")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = ULS_single[Compound == compounds[i]])
  summary_stats_ULS_single[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_ULS_single[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_ULS_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_ULS_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_ULS_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_ULS_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_ULS_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_ULS_single_label <- summary_stats_ULS_single %>%
  mutate(
    log_Cw = -2.5, log_Cs = 6.1,
    label =
      glue("*r<sup>2</sup>* = {round(r_squared, 2)} <br> *log K<sub>F</sub>* = {round(K_F, 2)} <br> *n<sub>F</sub>* = {round(n, 2)}")
  )

summary_stats_ULS_single_label <- summary_stats_ULS_single_label |>
  transform(Compound = compound)

#ULS Freundlich isotherm plot
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

ULS_single$Compound <- factor(ULS_single$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))
facetULS$Compound <- factor(facetULS$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                          "PFOA", "PFNA", "PFDA"))
summary_stats_ULS_single$compound <- factor(summary_stats_ULS_single$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                          "PFOA", "PFNA", "PFDA"))
summary_stats_ULS_single_label$Compound <- factor(summary_stats_ULS_single_label$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                                      "PFOA", "PFNA", "PFDA"))

ULS_facet_isotherm <- ggplot(data = ULS_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), 
             color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_compound),
              formula = y ~ x,
              method=lm,
              se=FALSE,
              colour = "grey",
              size = 0.5,
              data = facetULS
  ) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), 
              color = "black", 
              formula = y ~ x, 
              method=lm, 
              se=T, 
              fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  #ggtitle("ULS isotherm") +
  facet_wrap(~Compound) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none") +
  geom_richtext(
    data = summary_stats_ULS_single_label,
    aes(label = label, x = log_Cw, y = log_Cs),
    hjust = 0
  )
ULS_facet_isotherm
ggsave(filename="R/figs/ULS_facet_isotherm.pdf")

#Summary statistics DSL
nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

summary_stats_DSL_single <- data.table(K_F = rep(0, nr_compounds), 
                                       K_F_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "DSL")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = DSL_single[Compound == compounds[i]])
  summary_stats_DSL_single[compound == compounds[i], K_F := fit$coefficients[1]]
  summary_stats_DSL_single[compound == compounds[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_DSL_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_DSL_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_DSL_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_DSL_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_DSL_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_DSL_single_label <- summary_stats_DSL_single %>%
  mutate(
    log_Cw = 0.8, log_Cs = 1.4,
    label =
      glue("*r<sup>2</sup>* = {round(r_squared, 2)} <br> *log K<sub>F</sub>* = {round(K_F, 2)} <br> *n<sub>F</sub>* = {round(n, 2)}")
  )

summary_stats_DSL_single_label <- summary_stats_DSL_single_label |>
  transform(Compound = compound)



#DSL Freundlich isotherm plot
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

DSL_single$Compound <- factor(DSL_single$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))
facetDSL$Compound <- factor(facetDSL$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                          "PFOA", "PFNA", "PFDA"))
summary_stats_DSL_single$compound <- factor(summary_stats_DSL_single$compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                          "PFOA", "PFNA", "PFDA"))
summary_stats_DSL_single_label$Compound <- factor(summary_stats_DSL_single_label$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                                      "PFOA", "PFNA", "PFDA"))

DSL_facet_isotherm <- ggplot(data = DSL_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), 
             color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_compound), 
              formula = y ~ x, 
              method=lm, 
              se=FALSE, 
              colour = "grey", 
              size = 0.5,
              data = facetDSL
  ) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), 
              color = "black", 
              formula = y ~ x, 
              method=lm, 
              se=T, 
              fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  #ggtitle("DSL isotherm") +
  facet_wrap(~Compound) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none") +
  geom_richtext(
    data = summary_stats_DSL_single_label,
    aes(label = label, x = log_Cw, y = log_Cs),
    hjust = 0
  )
DSL_facet_isotherm
ggsave(filename="R/figs/DSL_facet_isotherm.pdf")


