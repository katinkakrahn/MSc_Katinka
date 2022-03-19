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

#Compare sorption of biochars across char types
summary_stats_CWC_single[, nr_CF2 := 4:9]
summary_stats_ULS_single[, nr_CF2 := 4:9]
summary_stats_DSL_single[, nr_CF2 := 4:9]

#Summary stats of entire biochar-water batch tests
summary_stats_single <- merge(summary_stats_CWC_single, summary_stats_ULS_single, all = TRUE)
summary_stats_single <- merge(summary_stats_single, summary_stats_DSL_single, all = TRUE)
summary_stats_single$compound <- factor(summary_stats_single$compound, levels = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA"))

#Summary stats of each compound
summary_stats_PFPeA <- filter(summary_stats_single, compound == "PFPeA")
summary_stats_PFHxA <- filter(summary_stats_single, compound == "PFHxA")
summary_stats_PFHpA <- filter(summary_stats_single, compound == "PFHpA")
summary_stats_PFOA <- filter(summary_stats_single, compound == "PFOA")
summary_stats_PFNA <- filter(summary_stats_single, compound == "PFNA")
summary_stats_PFDA <- filter(summary_stats_single, compound == "PFDA")

#Individual compound sorption
PFPeA_sorption_single <- filter(Sorption_BC_single, Compound == "PFPeA")
PFHxA_sorption_single <- filter(Sorption_BC_single, Compound == "PFHxA")
PFHpA_sorption_single <- filter(Sorption_BC_single, Compound == "PFHpA")
PFOA_sorption_single <- filter(Sorption_BC_single, Compound == "PFOA")
PFNA_sorption_single <- filter(Sorption_BC_single, Compound == "PFNA")
PFDA_sorption_single <- filter(Sorption_BC_single, Compound == "PFDA")


#PFPeA
PFPeA_isotherm <- ggplot(data = PFPeA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar))) + 
  geom_smooth(method = "lm", mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFPeA") + 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFPeA_isotherm

facetPFPeA <- PFPeA_sorption_single |>
  transform(pre_biochar = Biochar)

facetPFPeA <- rbind(
  transform(facetPFPeA, Biochar = unique(PFPeA_sorption_single$Biochar)[1]),
  transform(facetPFPeA, Biochar = unique(PFPeA_sorption_single$Biochar)[2]),
  transform(facetPFPeA, Biochar = unique(PFPeA_sorption_single$Biochar)[3])
)


PFPeA_facet_isotherm <- ggplot(data = PFPeA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_biochar), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetPFPeA) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("PFPeA") +
  geom_label(data = transform(summary_stats_PFPeA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = 1.25, y = 3.5, label = paste("K_F =",round(K_F, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFPeA_facet_isotherm
ggsave(filename="R/figs/PFPeA_facet_isotherm.pdf")

#PFHxA
PFHxA_isotherm <- ggplot(data = PFHxA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFHxA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFHxA_isotherm

facetPFHxA <- PFHxA_sorption_single |>
  transform(pre_biochar = Biochar)

facetPFHxA <- rbind(
  transform(facetPFHxA, Biochar = unique(PFHxA_sorption_single$Biochar)[1]),
  transform(facetPFHxA, Biochar = unique(PFHxA_sorption_single$Biochar)[2]),
  transform(facetPFHxA, Biochar = unique(PFHxA_sorption_single$Biochar)[3])
)


PFHxA_facet_isotherm <- ggplot(data = PFHxA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_biochar), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetPFHxA) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("PFHxA") +
  geom_label(data = transform(summary_stats_PFHxA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = 0.2, y = 0.25, label = paste("K_F =",round(K_F, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFHxA_facet_isotherm
ggsave(filename="R/figs/PFHxA_facet_isotherm.pdf")

#PFHpA
PFHpA_isotherm <- ggplot(data = PFHpA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFHpA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFHpA_isotherm

facetPFHpA <- PFHpA_sorption_single |>
  transform(pre_biochar = Biochar)

facetPFHpA <- rbind(
  transform(facetPFHpA, Biochar = unique(PFHpA_sorption_single$Biochar)[1]),
  transform(facetPFHpA, Biochar = unique(PFHpA_sorption_single$Biochar)[2]),
  transform(facetPFHpA, Biochar = unique(PFHpA_sorption_single$Biochar)[3])
)


PFHpA_facet_isotherm <- ggplot(data = PFHpA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_biochar), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetPFHpA) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("PFHpA") +
  geom_label(data = transform(summary_stats_PFHpA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = -0.75, y = 3.5, label = paste("K_F =",round(K_F, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFHpA_facet_isotherm
ggsave(filename="R/figs/PFHpA_facet_isotherm.pdf")

#PFOA
old <- c("Ci_(ug/L)", "Cw_(ug/L)", "Cs_(ug/kg)")
new <- c("Ci", "Cw", "Cs")
setnames(PFOA_sorption_single, old, new, skip_absent = TRUE)

PFOA_CWC_isotherm_nonlinear <- ggplot(subset(PFOA_sorption_single, Biochar %in% "CWC"), aes(x = Cw, y = Cs)) +
  geom_point(size = 3) + 
  geom_smooth(method = "lm", formula = y ~ poly(log(x), 2), se = FALSE) +
  labs(x = expression(C[w]~ug/L), y = expression(C[s]~ug/kg), title = "Freundlich isoterm CWC-PFOA") + 
  theme_bw()
PFOA_CWC_isotherm_nonlinear
ggsave(filename="R/figs/PFOA_CWC_isotherm_nonlinear.pdf")

PFOA_CWC_isotherm_linear <- ggplot(subset(PFOA_sorption_single, Biochar %in% "CWC"), aes(x = log_Cw, y = log_Cs)) +
  geom_point() + 
  geom_smooth(formula = y ~ x, method = lm, se=FALSE) +
  labs(x = expression(log~C[w]~ug/L), y = expression(log~C[s]~ug/kg), title = "Lineær Freundlich isoterm CWC-PFOA") + 
  theme_bw()
PFOA_CWC_isotherm_linear
ggsave(filename="R/figs/PFOA_CWC_isotherm_linear.pdf")

PFOA_isotherm <- ggplot(data = PFOA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFOA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFOA_isotherm

facetPFOA <- PFOA_sorption_single |>
  transform(pre_biochar = Biochar)

facetPFOA <- rbind(
  transform(facetPFOA, Biochar = unique(PFOA_sorption_single$Biochar)[1]),
  transform(facetPFOA, Biochar = unique(PFOA_sorption_single$Biochar)[2]),
  transform(facetPFOA, Biochar = unique(PFOA_sorption_single$Biochar)[3])
)


PFOA_facet_isotherm <- ggplot(data = PFOA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_biochar), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetPFOA) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("PFOA") +
  geom_label(data = transform(summary_stats_PFOA, Biochar = biochar), 
             size = 2, inherit.aes = T, aes(x = 0.8, y = 4.9, label = paste("K_F =",round(K_F, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFOA_facet_isotherm
ggsave(filename="R/figs/PFOA_facet_isotherm.pdf")

#PFNA
PFNA_isotherm <- ggplot(data = PFNA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFNA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFNA_isotherm

facetPFNA <- PFNA_sorption_single |>
  transform(pre_biochar = Biochar)

facetPFNA <- rbind(
  transform(facetPFNA, Biochar = unique(PFNA_sorption_single$Biochar)[1]),
  transform(facetPFNA, Biochar = unique(PFNA_sorption_single$Biochar)[2]),
  transform(facetPFNA, Biochar = unique(PFNA_sorption_single$Biochar)[3])
)


PFNA_facet_isotherm <- ggplot(data = PFNA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_biochar), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetPFNA) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("PFNA") +
  geom_label(data = transform(summary_stats_PFNA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = 0, y = 4.75, label = paste("K_F =",round(K_F, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFNA_facet_isotherm
ggsave(filename="R/figs/PFNA_facet_isotherm.pdf")

#PFDA
old <- c("Ci_(ug/L)", "Cw_(ug/L)", "Cs_(ug/g)")
new <- c("Ci", "Cw", "Cs")
setnames(PFDA_sorption_single, old, new, skip_absent = TRUE)

PFDA_CWC_isotherm_nonlinear <- ggplot(subset(PFDA_sorption_single, Biochar %in% "CWC"), aes(x = Cw, y = Cs)) +
  geom_point(size = 3) + 
  geom_smooth(method = "lm", formula = y ~ poly(log(x), 2), se = FALSE) +
  labs(x = expression(C[w]~ug/L), y = expression(C[s]~ug/kg), title = "Freundlich isoterm CWC-PFDA") + 
  theme_bw() +
  theme(text = element_text(size = 15), panel.grid = element_blank())
PFDA_CWC_isotherm_nonlinear
ggsave(filename="R/figs/PFD_CWC_isotherm_nonlinear.pdf")

PFDA_CWC_isotherm_linear <- ggplot(subset(PFDA_sorption_single, Biochar %in% "CWC"), aes(x = log_Cw, y = log_Cs)) +
  geom_point(size = 3) + 
  geom_smooth(formula = y ~ x, method = lm, se=FALSE) +
  labs(x = expression(log~C[w]~ug/L), y = expression(log~C[s]~ug/kg), title = "Lineær Freundlich isoterm CWC-PFDA") + 
  theme_bw() +
  theme(text = element_text(size = 15), panel.grid = element_blank())
PFDA_CWC_isotherm_linear
ggsave(filename="R/figs/PFDA_CWC_isotherm_linear.pdf")

PFDA_isotherm <- ggplot(data = PFDA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) +
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFDA")+ 
  theme_bw() +
  theme(legend.position = c(0.9, 0.18))
PFDA_isotherm

facetPFDA <- PFDA_sorption_single |>
  transform(pre_biochar = Biochar)

facetPFDA <- rbind(
  transform(facetPFDA, Biochar = unique(PFDA_sorption_single$Biochar)[1]),
  transform(facetPFDA, Biochar = unique(PFDA_sorption_single$Biochar)[2]),
  transform(facetPFDA, Biochar = unique(PFDA_sorption_single$Biochar)[3])
)


PFDA_facet_isotherm <- ggplot(data = PFDA_sorption_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_biochar), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetPFDA) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("PFDA") +
  geom_label(data = transform(summary_stats_PFDA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = 0, y = 5, label = paste("K_F =",round(K_F, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFDA_facet_isotherm
ggsave(filename="R/figs/PFDA_facet_isotherm.pdf")



