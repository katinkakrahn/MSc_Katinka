Sorption <- read_excel("C:/Users/KMK/OneDrive - NGI/VOW/Data/160222_sorption_rawdata.xlsx")
Sorption <- as.data.table(Sorption)

#Convert 1 and 0 to TRUE and FALSE and delete integer columns
Sorption$SoilLogic <- as.logical(Sorption$Soil_binary)
Sorption$mixLogic <- as.logical(Sorption$mix_binary)
Sorption <- subset(Sorption,select = -c(Soil_binary,mix_binary))

# Subset biochar and cocktail/single compound
Sorption_NAomit <- na.omit(Sorption)
Sorption_NA_C1omit <- Sorption_NAomit %>% slice(-c(20, 30, 40, 50, 79, 88, 98, 108, 118, 157, 167, 177))
Sorption_NA_C1omit <- Sorption_NA_C1omit %>%
  mutate(Kd = Cs/Cw, log_Kd = log10(Cs/Cw))
Sorption_BC_single <- subset(Sorption_NA_C1omit, mixLogic == FALSE)
Sorption_BC_mix <- subset(Sorption_NA_C1omit, mixLogic == TRUE)
write_xlsx(Sorption, "/Users/katinkakrahn/Library/CloudStorage/OneDrive-NGI/VOW/Data/010422_Sorption_BC.xlsx")
write_xlsx(Sorption_BC_single, "/Users/katinkakrahn/Library/CloudStorage/OneDrive-NGI/VOW/Data/010422_Sorption_BC_single.xlsx")
write_xlsx(Sorption_BC_mix, "/Users/katinkakrahn/Library/CloudStorage/OneDrive-NGI/VOW/Data/010422_Sorption_BC_mix.xlsx")


CWC_single <- filter(Sorption_BC_single, Biochar == "CWC")
ULS_single <- filter(Sorption_BC_single, Biochar == "ULS")
DSL_single <- filter(Sorption_BC_single, Biochar == "DSL")

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

nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

#Sorption isotherm all chars
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
  labs(x = expression(log~C[w]), y = expression(log~C[s]), color = "") + 
  facet_wrap(~Compound) +
  #ggtitle("Freundlich linear sorption isotherms") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom") #+
  # geom_richtext(
  #   data = summary_stats_ULS_single_label,
  #   aes(label = label, x = log_Cw, y = log_Cs),
  #   hjust = 0
  # )
Sorption_isotherms
set_palette(Sorption_isotherms, "uchicago")
ggsave(filename="R/figs/Sorption_isotherms_single_BC.pdf")

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
             color = "#800000FF",
             group = "Compound",
             size = 2) +
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
    hjust = 0)
CWC_facet_isotherm
ggsave(filename="R/figs/CWC_facet_isotherm.pdf")


expression("Diameter of apeture (" * mu ~ "m)")
#Summary statistics ULS
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
             color = "#800000FF",
             group = "Compound",
             size = 2) +
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
  facet_wrap(~Compound) +
  #ggtitle("ULS isotherm") +
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
             color = "#800000FF",
             group = "Compound",
             size = 2) +
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
  facet_wrap(~Compound) +
  #ggtitle("DSL isotherm") +
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

