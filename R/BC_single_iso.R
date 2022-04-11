# Sorption biochar water data ---- 
Sorption <- read_excel("R/data_raw/160222_sorption_rawdata.xlsx")
as.data.table(Sorption)
Sorption <- as.data.table(Sorption)

#Convert 1 and 0 to TRUE and FALSE and delete integer columns
Sorption$SoilLogic <- as.logical(Sorption$Soil_binary)
Sorption$mixLogic <- as.logical(Sorption$mix_binary)
Sorption <- select(Sorption, -Soil_binary, -mix_binary)

# Subset biochar and cocktail/single compound
Sorption_NAomit <- na.omit(Sorption)
Sorption_NA_C1omit <- Sorption_NAomit %>% 
  filter(Conc_point != 1) %>% 
  mutate(Kd = Cs/Cw, log_Kd = log10(Cs/Cw))
Sorption_BC_single <- subset(Sorption_NA_C1omit, mixLogic == FALSE)
Sorption_BC_mix <- subset(Sorption_NA_C1omit, mixLogic == TRUE)

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

# Sorption BC cocktail ----

Sorption_BC_mix_summary <- Sorption_BC_mix[, .(Conc_point = 10,
                                               Ci = mean(Ci),
                                               Cw = mean(Cw),
                                               Cs = mean(Cs),
                                               Kd = mean(Cs/Cw),
                                               log_Kd = log10(mean(Cs/Cw)),
                                               log_Cw = mean(log_Cw), 
                                               log_Cs = mean(log_Cs),
                                               se_logCw = std.error(log_Cw),
                                               se_logCs = std.error(log_Cs),
                                               se_logKd = std.error(log10(Cs/Cw)),
                                               mixLogic = TRUE, 
                                               SoilLogic = FALSE
),
keyby = .(Biochar, Compound)]

# CWC ----

#Summary statistics CWC
summary_stats_CWC_single <- data.table(log_KF = rep(0, nr_compounds), 
                                       log_KF_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       log_Cw = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "CWC")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = CWC_single[Compound == compounds[i]])
  summary_stats_CWC_single[compound == compounds[i], log_KF := fit$coefficients[1]]
  summary_stats_CWC_single[compound == compounds[i], log_KF_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_CWC_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_CWC_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_CWC_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_CWC_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_CWC_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_CWC_single_label <- summary_stats_CWC_single %>%
  mutate(
    log_Cw = 0.6, log_Cs = 1,
    label =
      glue("*r<sup>2</sup>* = {round(r_squared, 2)} <br> *log K<sub>F</sub>* = {round(log_KF, 2)} <br> *n<sub>F</sub>* = {round(n, 2)}")
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
Sorption_BC_mix_summary$Compound <- factor(Sorption_BC_mix_summary$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                                      "PFOA", "PFNA", "PFDA"))


 
CWC_facet_isotherm <- ggplot(data = CWC_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), 
             color = "gray45", size = 1) + 
  geom_point(data = subset(Sorption_BC_mix_summary, Biochar %in% "CWC"), 
             mapping = aes(x = log_Cw, y = log_Cs), 
             color = "red",
             group = "Compound",
             size = 1) +
  # geom_errorbar(data = subset(Sorption_BC_mix_summary, Biochar %in% "CWC"),
  #               mapping = aes(x=log_Cw, y=log_Cs, ymin=log_Cs-se_logCs, ymax=log_Cs+se_logCs, xmin=log_Cw-se_logCw, xmax=log_Cw+se_logCw),
  #               width=.1,
  #               position=position_dodge(.9),
  #               color = "red") +
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
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), color = "") + 
  facet_wrap(~Compound) +
  #ggtitle("CWC isotherm") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  theme(panel.grid = element_blank()) #+
  #guides(color = "none") #+
  # geom_richtext(
  #   data = summary_stats_CWC_single_label,
  #   aes(label = label, x = log_Cw, y = log_Cs),
  #   hjust = 0)
CWC_facet_isotherm
ggsave(filename="R/figs/CWC_facet_isotherm.pdf")

# ULS ----
# Summary statistics ULS
summary_stats_ULS_single <- data.table(log_KF = rep(0, nr_compounds), 
                                       log_KF_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       log_Cw = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "ULS")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = ULS_single[Compound == compounds[i]])
  summary_stats_ULS_single[compound == compounds[i], log_KF := fit$coefficients[1]]
  summary_stats_ULS_single[compound == compounds[i], log_KF_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_ULS_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_ULS_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_ULS_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_ULS_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_ULS_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_ULS_single_label <- summary_stats_ULS_single %>%
  mutate(
    log_Cw = -2.5, log_Cs = 6,
    label =
      glue("*r<sup>2</sup>* = {round(r_squared, 2)} <br> *log K<sub>F</sub>* = {round(log_KF, 2)} <br> *n<sub>F</sub>* = {round(n, 2)}")
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
Sorption_BC_mix_summary$Compound <- factor(Sorption_BC_mix_summary$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                        "PFOA", "PFNA", "PFDA"))


ULS_facet_isotherm <- ggplot(data = ULS_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), 
             color = "gray45", size = 1) + 
  geom_point(data = subset(Sorption_BC_mix_summary, Biochar %in% "ULS"), 
             mapping = aes(x = log_Cw, y = log_Cs), 
             color = "red",
             group = "Compound",
             size = 1) +
  # geom_errorbar(data = subset(Sorption_BC_mix_summary, Biochar %in% "ULS"),
  #               mapping = aes(x=log_Cw, y=log_Cs, ymin=log_Cs-se_logCs, ymax=log_Cs+se_logCs, xmin=log_Cw-se_logCw, xmax=log_Cw+se_logCw),
  #               width=.1,
  #               position=position_dodge(.9),
  #               color = "red") +
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
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), color = "") + 
  facet_wrap(~Compound) +
  #ggtitle("ULS isotherm") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  theme(panel.grid = element_blank()) +
  guides(color = "none") #+
  # geom_richtext(
  #   data = summary_stats_ULS_single_label,
  #   aes(label = label, x = log_Cw, y = log_Cs),
  #   hjust = 0
  # )
ULS_facet_isotherm
ggsave(filename="R/figs/ULS_facet_isotherm.pdf")

# DSL ---- 
#Summary statistics DSL
summary_stats_DSL_single <- data.table(log_KF = rep(0, nr_compounds), 
                                       log_KF_std_error = rep(0, nr_compounds),
                                       n = rep(0, nr_compounds),
                                       n_std_error = rep(0, nr_compounds),
                                       r_squared = rep(0, nr_compounds),
                                       residual_std_error = rep(0, nr_compounds),
                                       p_value = rep(0, nr_compounds),
                                       log_Cw = rep(0, nr_compounds),
                                       compound = compounds,
                                       biochar = "DSL")

for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = DSL_single[Compound == compounds[i]])
  summary_stats_DSL_single[compound == compounds[i], log_KF := fit$coefficients[1]]
  summary_stats_DSL_single[compound == compounds[i], log_KF_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_DSL_single[compound == compounds[i], n := fit$coefficients[2]]
  summary_stats_DSL_single[compound == compounds[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_DSL_single[compound == compounds[i], r_squared := summary(fit)$r.squared]
  summary_stats_DSL_single[compound == compounds[i], residual_std_error := summary(fit)$sigma]
  summary_stats_DSL_single[compound == compounds[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_DSL_single_label <- summary_stats_DSL_single %>%
  mutate(
    log_Cw = 0.4, log_Cs = 1.4,
    label =
      glue("*r<sup>2</sup>* = {round(r_squared, 2)} <br> *log K<sub>F</sub>* = {round(log_KF, 2)} <br> *n<sub>F</sub>* = {round(n, 2)}")
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
Sorption_BC_mix_summary$Compound <- factor(Sorption_BC_mix_summary$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                        "PFOA", "PFNA", "PFDA"))


DSL_facet_isotherm <- ggplot(data = DSL_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Compound)), 
             color = "gray45", size = 1) + 
  geom_point(data = subset(Sorption_BC_mix_summary, Biochar %in% "DSL"), 
             mapping = aes(x = log_Cw, y = log_Cs), 
             color = "red",
             group = "Compound",
             size = 1) +
  # geom_errorbar(data = subset(Sorption_BC_mix_summary, Biochar %in% "DSL"),
  #               mapping = aes(x=log_Cw, y=log_Cs, ymin=log_Cs-se_logCs, ymax=log_Cs+se_logCs, xmin=log_Cw-se_logCw, xmax=log_Cw+se_logCw),
  #               width=.1,
  #               position=position_dodge(.9),
  #               color = "red") +
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
  labs(x = TeX(r'($log~C_{w}~(\mu g~L^{-1})$)'), y = TeX(r'($log~C_{s}~(\mu g~kg^{-1})$)'), color = "") +  
  facet_wrap(~Compound) +
  #ggtitle("DSL isotherm") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  theme(panel.grid = element_blank()) +
  guides(color = "none") #+
  # geom_richtext(
  #   data = summary_stats_DSL_single_label,
  #   aes(label = label, x = log_Cw, y = log_Cs),
  #   hjust = 0
  # )
DSL_facet_isotherm
ggsave(filename="R/figs/DSL_facet_isotherm.pdf")

# Summary stats ----
summary_stats_CWC_single[, nr_CF2 := 4:9]
summary_stats_ULS_single[, nr_CF2 := 4:9]
summary_stats_DSL_single[, nr_CF2 := 4:9]

summary_stats_single <- merge(summary_stats_CWC_single, summary_stats_ULS_single, all = TRUE)
summary_stats_single <- merge(summary_stats_single, summary_stats_DSL_single, all = TRUE)
summary_stats_single$compound <- factor(summary_stats_single$compound, levels = c("PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA"))
write_xlsx(summary_stats_single, "R/data_manipulated/310322_summary_stats_single.xlsx")

# Sorption isotherm all chars ----
Sorption_BC_single$Compound <- factor(Sorption_BC_single$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                              "PFOA", "PFNA", "PFDA"))

Sorption_isotherms <- ggplot(data = Sorption_BC_single) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), 
             size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), 
              formula = y ~ x, 
              method=lm, 
              se=T, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g~L^{-1})$)'), y = TeX(r'($log~C_{s}~(\mu g~kg^{-1})$)'), color = "") + 
  facet_wrap(~Compound) +
  #ggtitle("Freundlich linear sorption isotherms") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom") #+
# geom_richtext(
#   data = summary_stats_ULS_single_label,
#   aes(label = label, x = log_Cw, y = log_Cs),
#   hjust = 0
# )
Sorption_isotherms
ggsave(filename="R/figs/Sorption_isotherms_single_BC.pdf")

# Summary stats of each compound ----
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


# PFPeA ----
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
  geom_label(data = transform(summary_stats_PFPeA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = 1.25, y = 3.5, label = paste("log_KF =",round(log_KF, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFPeA_facet_isotherm
ggsave(filename="R/figs/PFPeA_facet_isotherm.pdf")

# PFHxA ----
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
  geom_label(data = transform(summary_stats_PFHxA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = 0.2, y = 0.25, label = paste("log_KF =",round(log_KF, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFHxA_facet_isotherm
ggsave(filename="R/figs/PFHxA_facet_isotherm.pdf")

# PFHpA ----
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
  geom_label(data = transform(summary_stats_PFHpA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = -0.75, y = 3.5, label = paste("log_KF =",round(log_KF, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFHpA_facet_isotherm
ggsave(filename="R/figs/PFHpA_facet_isotherm.pdf")

# PFOA ----
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

summary_stats_PFOA_single_label <- summary_stats_PFOA %>%
  mutate(
    log_Cw = 1.85, log_Cs = 5.3,
    label =
      glue("*r<sup>2</sup>* = {round(r_squared, 2)} <br> *log K<sub>F</sub>* = {round(log_KF, 2)} <br> *n<sub>F</sub>* = {round(n, 2)}")
  )

summary_stats_PFOA_single_label <- summary_stats_PFOA_single_label |>
  transform(Biochar = biochar)

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
  facet_grid(rows = vars(Biochar)) +
  theme_bw() +
  geom_richtext(
    data = summary_stats_PFOA_single_label,
    aes(label = label, x = log_Cw, y = log_Cs),
    hjust = 0
  ) +
  theme(panel.grid = element_blank()) +
  guides(color = "none")
PFOA_facet_isotherm
ggsave(filename="R/figs/PFOA_facet_isotherm.pdf")

# PFNA ----
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
  geom_label(data = transform(summary_stats_PFNA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = 0, y = 4.75, label = paste("log_KF =",round(log_KF, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFNA_facet_isotherm
ggsave(filename="R/figs/PFNA_facet_isotherm.pdf")

# PFDA ----
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
ggsave(filename="R/figs/PFDA_CWC_isotherm_nonlinear.pdf")

PFDA_CWC_isotherm_linear <- ggplot(subset(PFDA_sorption_single, Biochar %in% "CWC"), aes(x = log_Cw, y = log_Cs)) +
  geom_point(size = 3) + 
  geom_smooth(formula = y ~ x, method = lm, se=FALSE) +
  labs(x = expression(log~C[w]~ug/L), y = expression(log~C[s]~ug/kg), title = "Lineær Freundlich isoterm CWC-PFDA") + 
  theme_bw() +
  theme(text = element_text(size = 15), panel.grid = element_blank())
PFDA_CWC_isotherm_linear
ggsave(filename="R/figs/PFDA_CWC_isotherm_linear.pdf")

#ULS PFDA
PFDA_ULS_isotherm_nonlinear <- ggplot(subset(PFDA_sorption_single, Biochar %in% "ULS"), aes(x = Cw, y = Cs)) +
  geom_point(size = 3) + 
  geom_smooth(method = "lm", formula = y ~ poly(log(x), 2), se = FALSE) +
  labs(x = expression(C[w]~ug/L), y = expression(C[s]~ug/kg), title = "Freundlich isoterm ULS-PFDA") + 
  theme_bw() +
  theme(text = element_text(size = 15), panel.grid = element_blank())
PFDA_ULS_isotherm_nonlinear
ggsave(filename="R/figs/PFDA_ULS_isotherm_nonlinear.pdf")

PFDA_ULS_isotherm_linear <- ggplot(subset(PFDA_sorption_single, Biochar %in% "ULS"), aes(x = log_Cw, y = log_Cs)) +
  geom_point(size = 3) + 
  geom_smooth(formula = y ~ x, method = lm, se=FALSE) +
  labs(x = expression(log~C[w]~ug/L), y = expression(log~C[s]~ug/kg), title = "Lineær Freundlich isoterm ULS-PFDA") + 
  theme_bw() +
  theme(text = element_text(size = 15), panel.grid = element_blank())
PFDA_ULS_isotherm_linear
ggsave(filename="R/figs/PFDA_ULS_isotherm_linear.pdf")

#DSL PFDA
PFDA_DSL_isotherm_nonlinear <- ggplot(subset(PFDA_sorption_single, Biochar %in% "DSL"), aes(x = Cw, y = Cs)) +
  geom_point(size = 3) + 
  geom_smooth(method = "lm", formula = y ~ poly(log(x), 2), se = FALSE) +
  labs(x = expression(C[w]~ug/L), y = expression(C[s]~ug/kg), title = "Freundlich isoterm DSL-PFDA") + 
  theme_bw() +
  theme(text = element_text(size = 15), panel.grid = element_blank())
PFDA_DSL_isotherm_nonlinear
ggsave(filename="R/figs/PFDA_DSL_isotherm_nonlinear.pdf")

PFDA_DSL_isotherm_linear <- ggplot(subset(PFDA_sorption_single, Biochar %in% "DSL"), aes(x = log_Cw, y = log_Cs)) +
  geom_point(size = 3) + 
  geom_smooth(formula = y ~ x, method = lm, se=FALSE) +
  labs(x = expression(log~C[w]~ug/L), y = expression(log~C[s]~ug/kg), title = "Lineær Freundlich isoterm DSL-PFDA") + 
  theme_bw() +
  theme(text = element_text(size = 15), panel.grid = element_blank())
PFDA_DSL_isotherm_linear
ggsave(filename="R/figs/PFDA_DSL_isotherm_linear.pdf")

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
  geom_label(data = transform(summary_stats_PFDA, Biochar = biochar), size = 2, inherit.aes = T, aes(x = 0, y = 5, label = paste("log_KF =",round(log_KF, digits = 2),","," ","n =",round(n, digits = 2),","," ","R^2 =", round(r_squared, digits = 2)))) +
  facet_wrap(~Biochar) +
  theme_bw() +
  guides(color = "none")
PFDA_facet_isotherm
ggsave(filename="R/figs/PFDA_facet_isotherm.pdf")





