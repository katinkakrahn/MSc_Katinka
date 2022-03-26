Sorption_soil <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Data/010322_sorption_rawdata_soil.xlsx")
as.data.table(Sorption_soil)
Sorption_soil <- as.data.table(Sorption_soil)

Sorption_soil$SoilLogic <- as.logical(Sorption_soil$Soil_binary)
Sorption_soil$mixLogic <- as.logical(Sorption_soil$mix_binary)
Sorption_soil <- subset(Sorption_soil,select = -c(Soil_binary,mix_binary))
Sorption_soil <- Sorption_soil %>%
  mutate(Kd = Cs/Cw, log_Kd = log10(Cs/Cw))
Sorption_BC_BS <- kable(Sorption_soil, "latex", booktabs = TRUE, digits = 2)


# Subset biochar and cocktail/single compound, change from Sorption_soil to Sorption_soil_NAomit when data is updated
Sorption_soil_NAomit <- na.omit(Sorption_soil)
#Sorption_soil_NA_C1omit <- Sorption_NAomit %>% slice(-c(column numbers))
Sorption_soil_blank <- subset(Sorption_soil_NAomit, Biochar == "no") 
Sorption_soil_blank_mix <- subset(Sorption_soil_blank, mixLogic == TRUE)
Sorption_soil_blank_PFOA <- subset(Sorption_soil_blank, mixLogic == FALSE)
Sorption_soil_BC <- Sorption_soil_NAomit[Biochar != 'no'] 
Sorption_soil_BC_PFOA <- subset(Sorption_soil_BC, mixLogic == FALSE) 
Sorption_soil_BC_mix <- subset(Sorption_soil_BC, mixLogic == TRUE)

CWC_soil_mix <- filter(Sorption_soil_BC_mix, Biochar == "CWC")
ULS_soil_mix <- filter(Sorption_soil_BC_mix, Biochar == "ULS")
DSL_soil_mix <- filter(Sorption_soil_BC_mix, Biochar == "DSL")

CWC_soil_PFOA <- filter(Sorption_soil_BC_PFOA, Biochar == "CWC")
ULS_soil_PFOA <- filter(Sorption_soil_BC_PFOA, Biochar == "ULS")
DSL_soil_PFOA <- filter(Sorption_soil_BC_PFOA, Biochar == "DSL")

Sorption_soil_blank_summary <- Sorption_soil_blank[, .(mean_logCw = mean(log_Cw), 
                            mean_logCs = mean(log_Cs),
                            se_logCw = std.error(log_Cw),
                            se_logCs = std.error(log_Cs),
                            log_Kd = mean(log_Cs/log_Cw),
                            se_logKd = std.error(log_Cs/log_Cw)
),
keyby = .(Compound,mixLogic)]

Sorption_soil_blank_summary$Compound <- factor(Sorption_soil_blank_summary$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))

SoilBlankKd <- ggplot(data = Sorption_soil_blank_summary, aes(x = Compound, y = log_Kd, color = mixLogic)) + 
  geom_point()+ 
  geom_errorbar(aes(ymin=log_Kd-se_logKd, ymax=log_Kd+se_logKd), width = .1) + 
  labs(x = "", y = expression(log~K[d]), col = "") + 
  scale_color_manual(
    values = c('black','grey60'),
    breaks = c("TRUE", "FALSE"),
    labels = c("cocktail", "single compound")
  ) +
  theme_bw()
SoilBlankKd
ggsave(filename="R/figs/SoilBlankKd.pdf")
#for this plot need to normalize for same conc. to get comparable Kd.

#Sorption_soil_BC <- Sorption_soil_BC[, Compound := as.factor(Compound)]
#CWC_soil_mix <- CWC_soil_mix[, Compound := as.factor(Compound)]

#CWC facet soil mix
compounds2 <- unique(CWC_soil_mix$Compound)


summary_stats_CWC_soil_mix <- data.table(K_F = rep(0, nr_compounds-1), 
                                       K_F_std_error = rep(0, nr_compounds-1),
                                       n = rep(0, nr_compounds-1),
                                       n_std_error = rep(0, nr_compounds-1),
                                       r_squared = rep(0, nr_compounds-1),
                                       residual_std_error = rep(0, nr_compounds-1),
                                       p_value = rep(0, nr_compounds-1),
                                       compound = compounds2,
                                       biochar = "CWC")

for(i in 1:(nr_compounds-1)){
  fit <- lm(log_Cs ~ log_Cw, data = CWC_soil_mix[Compound == compounds2[i]])
  summary_stats_CWC_soil_mix[compound == compounds2[i], K_F := fit$coefficients[1]]
  summary_stats_CWC_soil_mix[compound == compounds2[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_CWC_soil_mix[compound == compounds2[i], n := fit$coefficients[2]]
  summary_stats_CWC_soil_mix[compound == compounds2[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_CWC_soil_mix[compound == compounds2[i], r_squared := summary(fit)$r.squared]
  summary_stats_CWC_soil_mix[compound == compounds2[i], residual_std_error := summary(fit)$sigma]
  summary_stats_CWC_soil_mix[compound == compounds2[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                         summary(fit)$fstatistic[3],lower.tail=F)]
}

#ULS facet soil mix
summary_stats_ULS_soil_mix <- data.table(K_F = rep(0, nr_compounds-1), 
                                         K_F_std_error = rep(0, nr_compounds-1),
                                         n = rep(0, nr_compounds-1),
                                         n_std_error = rep(0, nr_compounds-1),
                                         r_squared = rep(0, nr_compounds-1),
                                         residual_std_error = rep(0, nr_compounds-1),
                                         p_value = rep(0, nr_compounds-1),
                                         compound = compounds2,
                                         biochar = "ULS")

for(i in 1:(nr_compounds-1)){
  fit <- lm(log_Cs ~ log_Cw, data = ULS_soil_mix[Compound == compounds2[i]])
  summary_stats_ULS_soil_mix[compound == compounds2[i], K_F := fit$coefficients[1]]
  summary_stats_ULS_soil_mix[compound == compounds2[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_ULS_soil_mix[compound == compounds2[i], n := fit$coefficients[2]]
  summary_stats_ULS_soil_mix[compound == compounds2[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_ULS_soil_mix[compound == compounds2[i], r_squared := summary(fit)$r.squared]
  summary_stats_ULS_soil_mix[compound == compounds2[i], residual_std_error := summary(fit)$sigma]
  summary_stats_ULS_soil_mix[compound == compounds2[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                      summary(fit)$fstatistic[3],lower.tail=F)]
}


#DSL facet soil mix
summary_stats_DSL_soil_mix <- data.table(K_F = rep(0, nr_compounds-1), 
                                         K_F_std_error = rep(0, nr_compounds-1),
                                         n = rep(0, nr_compounds-1),
                                         n_std_error = rep(0, nr_compounds-1),
                                         r_squared = rep(0, nr_compounds-1),
                                         residual_std_error = rep(0, nr_compounds-1),
                                         p_value = rep(0, nr_compounds-1),
                                         compound = compounds2,
                                         biochar = "DSL")

for(i in 1:(nr_compounds-1)){
  fit <- lm(log_Cs ~ log_Cw, data = DSL_soil_mix[Compound == compounds2[i]])
  summary_stats_DSL_soil_mix[compound == compounds2[i], K_F := fit$coefficients[1]]
  summary_stats_DSL_soil_mix[compound == compounds2[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_DSL_soil_mix[compound == compounds2[i], n := fit$coefficients[2]]
  summary_stats_DSL_soil_mix[compound == compounds2[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_DSL_soil_mix[compound == compounds2[i], r_squared := summary(fit)$r.squared]
  summary_stats_DSL_soil_mix[compound == compounds2[i], residual_std_error := summary(fit)$sigma]
  summary_stats_DSL_soil_mix[compound == compounds2[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                      summary(fit)$fstatistic[3],lower.tail=F)]
}


#PFOA soil summary statistics
summary_stats_PFOA_soil <- data.table(K_F = rep(0, nr_biochars), 
                                         K_F_std_error = rep(0, nr_biochars),
                                         n = rep(0, nr_biochars),
                                         n_std_error = rep(0, nr_biochars),
                                         r_squared = rep(0, nr_biochars),
                                         residual_std_error = rep(0, nr_biochars),
                                         p_value = rep(0, nr_biochars),
                                      biochar = biochars,
                                         compound = "PFOA")

for(i in 1:(nr_biochars)){
  fit <- lm(log_Cs ~ log_Cw, data = Sorption_soil_BC_PFOA[Biochar == biochars[i]])
  summary_stats_PFOA_soil[biochar == biochars[i], K_F := fit$coefficients[1]]
  summary_stats_PFOA_soil[biochar == biochars[i], K_F_std_error := summary(fit)$coefficients[1,2]]
  summary_stats_PFOA_soil[biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_PFOA_soil[biochar == biochars[i], n_std_error := summary(fit)$coefficients[2,2]]
  summary_stats_PFOA_soil[biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_PFOA_soil[biochar == biochars[i], residual_std_error := summary(fit)$sigma]
  summary_stats_PFOA_soil[biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                      summary(fit)$fstatistic[3],lower.tail=F)]
}



#PFOA soil isotherm
summary_stats_PFOA_soil_label <- summary_stats_PFOA_soil %>%
  mutate(
    log_Cw = 2.5, log_Cs = 4.9,
    label =
      glue("*r<sup>2</sup>* = {round(r_squared, 2)} <br> *log K<sub>F</sub>* = {round(K_F, 2)} <br> *n<sub>F</sub>* = {round(n, 2)}")
  )

summary_stats_PFOA_soil_label <- summary_stats_PFOA_soil_label |>
  transform(Biochar = biochar)

facetPFOA_soil <- Sorption_soil_BC_PFOA |>
  transform(pre_biochar = Biochar)

facetPFOA_soil <- rbind(
  transform(facetPFOA_soil, Biochar = unique(Sorption_soil_BC_PFOA$Biochar)[1]),
  transform(facetPFOA_soil, Biochar = unique(Sorption_soil_BC_PFOA$Biochar)[2]),
  transform(facetPFOA_soil, Biochar = unique(Sorption_soil_BC_PFOA$Biochar)[3])
)

Sorption_soil_BC_PFOA$Biochar <- factor(Sorption_soil_BC_PFOA$Biochar, levels = c("CWC", "DSL", "ULS"))
facetPFOA_soil$Biochar <- factor(facetPFOA_soil$Biochar, levels = c("CWC", "DSL", "ULS"))
summary_stats_PFOA_soil_label$Biochar <- factor(summary_stats_PFOA_soil_label$Biochar, levels = c("CWC", "DSL", "ULS"))


PFOA_facet_soil_isotherm <- ggplot(data = Sorption_soil_BC_PFOA) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "gray45", size = 1) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = pre_biochar), formula = y ~ x, method=lm, se=FALSE, colour = "grey", size = 0.5,
              data = facetPFOA_soil) +
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, group = factor(Biochar)), color = "black", formula = y ~ x, method=lm, se=T, fullrange = FALSE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s])) + 
  ggtitle("PFOA soil") +
  facet_grid(rows = vars(Biochar)) +
  theme_bw() +
  geom_richtext(
    data = summary_stats_PFOA_soil_label,
    aes(label = label, x = log_Cw, y = log_Cs),
    hjust = 0
  ) +
  theme(panel.grid = element_blank()) +
  guides(color = "none")
PFOA_facet_soil_isotherm
ggsave(filename="R/figs/PFOA_facet_soil_isotherm.pdf")



#Cocktail soil isotherms










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


CWC_facet_isotherm_soil_mix <- ggplot(data = CWC_) +
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












































Sorption_soil_blank_mix$Compound <- factor(Sorption_soil_blank_mix$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                  "PFOA", "PFNA", "PFDA"))

Sorption_soil_blank_mix_points <- ggplot(data = Sorption_soil_blank_mix) + 
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)))+ 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "Soil blank cocktail triplicates") + 
  theme_bw()
Sorption_soil_blank_mix_points

Sorption_soil_blank_PFOA$Compound <- factor(Sorption_soil_blank_PFOA$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                        "PFOA", "PFNA", "PFDA"))

Sorption_soil_blank_PFOA_points <- ggplot(data = Sorption_soil_blank_PFOA) + 
  geom_point(mapping = aes(x = log_Cw, y = log_Cs), shape = 21, size = 3, fill = "#077DAA")+ 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), title = "Soil blank PFOA triplicates") + 
  theme_bw()
Sorption_soil_blank_PFOA_points

# Soil cocktail isotherms
CWC_soil_mix$Compound <- factor(CWC_soil_mix$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                              "PFOA", "PFNA", "PFDA"))

CWC_isotherm_soil_mix <- ggplot(data = CWC_soil_mix) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "CWC soil cocktail isotherm") + 
  theme_bw()
CWC_isotherm_soil_mix

ULS_soil_mix$Compound <- factor(ULS_soil_mix$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                  "PFOA", "PFNA", "PFDA"))
ULS_isotherm_soil_mix <- ggplot(data = ULS_soil_mix) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "ULS soil cocktail isotherm") + 
  theme_bw()
ULS_isotherm_soil_mix

DSL_soil_mix$Compound <- factor(DSL_soil_mix$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                  "PFOA", "PFNA", "PFDA"))

DSL_isotherm_soil_mix <- ggplot(data = DSL_soil_mix) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Compound)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Compound", title = "DSL soil cocktail isotherm") + 
  theme_bw()
DSL_isotherm_soil_mix

# Soil PFOA isotherms
CWC_isotherm_soil_PFOA <- ggplot(data = CWC_soil_PFOA) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs)) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), title = "CWC soil PFOA isotherm") + 
  theme_bw()
CWC_isotherm_soil_PFOA

ULS_isotherm_soil_PFOA <- ggplot(data = ULS_soil_PFOA) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs)) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), title = "ULS soil PFOA isotherm") + 
  theme_bw()
ULS_isotherm_soil_PFOA

DSL_isotherm_soil_PFOA <- ggplot(data = DSL_soil_PFOA) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs)) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), title = "DSL soil PFOA isotherm") + 
  theme_bw()
DSL_isotherm_soil_PFOA

PFOA_soil_isotherm <- ggplot(data = Sorption_soil_BC_PFOA) +
  geom_point(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar))) + 
  geom_smooth(mapping = aes(x = log_Cw, y = log_Cs, color = factor(Biochar)), formula = y ~ x, method=lm, se=FALSE, fullrange = TRUE) + 
  labs(x = expression(log~C[w]), y = expression(log~C[s]), col = "Biochar", title = "PFOA biochar soil isotherm") + 
  theme_bw()
PFOA_soil_isotherm
