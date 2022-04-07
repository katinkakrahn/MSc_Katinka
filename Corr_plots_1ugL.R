Elements_ratios <- read_excel("/Users/katinkakrahn/Library/CloudStorage/OneDrive-NGI/VOW/Data/270322_Elements_ratios.xlsx")
Elements_ratios <- as.data.table(Elements_ratios)
Kd_1ugL_select <- merge(Kd_1ugL_select, summary_stats_single, all=FALSE, by = c("compound", "biochar"))
Kd_1ugL_select <- Kd_1ugL[, logKd_error := sqrt(log_KF_std_error^2+n_std_error^2)]
Elements_ratios_1ugL <- full_join(x = Biochar_ratios, y = Biochar_1ugL_select, by = "biochar")

Kd_1ugL_Ca <- ggplot(data = subset(Elements_ratios_1ugL, Parameter %in% "Ca"),
                                aes(x = Mean_sameunit, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) + 
  labs(x = "Ca (g/kg)", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_Ca
set_palette(Kd_1ugL_Ca, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_Ca.pdf")

Kd_1ugL_Fe <- ggplot(data = subset(Elements_ratios_1ugL, Parameter %in% "Fe"),
                                aes(x = Mean_sameunit, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) + 
  labs(x = "Fe (g/kg)", y = expression(log~K[d]), color = "", shape = "") +
  guides(color = "none", shape = "none") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_Fe
set_palette(Kd_1ugL_Fe, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_Fe.pdf")

Kd_1ugL_C <- ggplot(data = subset(Elements_ratios_1ugL, Parameter %in% "C"),
                               aes(x = Mean_sameunit, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) + 
  labs(x = "C (g/kg)", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_C
set_palette(Kd_1ugL_C, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_C.pdf")


Kd_1ugL_C_Ca <- ggplot(data = Elements_ratios_1ugL,
                                  aes(x = C_Ca, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) + 
  labs(x = "C/Ca", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_C_Ca
set_palette(Kd_1ugL_C_Ca, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_C_Ca.pdf")

Kd_1ugL_C_Fe <- ggplot(data = Elements_ratios_1ugL,
                       aes(x = C_Fe, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "C/Fe", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_C_Fe
set_palette(Kd_1ugL_C_Fe, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_C_Fe.pdf")

#SA PV
SA_PV <- filter(Elements_ratios_1ugL, Use == 1)

Kd_1ugL_SACO2 <- ggplot(data = SA_PV,
                       aes(x = CO2_SA, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "SA CO2", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SACO2
set_palette(Kd_1ugL_SACO2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SACO2.pdf")

Kd_1ugL_PVCO2 <- ggplot(data = SA_PV,
                        aes(x = CO2_PV, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.001)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "PV CO2", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVCO2
set_palette(Kd_1ugL_PVCO2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PVCO2.pdf")

Kd_1ugL_SAN2 <- ggplot(data = SA_PV,
                        aes(x = N2_SA2, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "SA N2", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SAN2
set_palette(Kd_1ugL_SAN2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SAN2.pdf")

Kd_1ugL_PVN2 <- ggplot(data = SA_PV,
                        aes(x = N2_PV, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.001)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "PV N2", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVN2
set_palette(Kd_1ugL_PVN2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PVN2.pdf")

#SA PV Ca
SA_PV <- filter(Elements_ratios_1ugL, Use == 1)

Kd_1ugL_SACO2_Ca <- ggplot(data = SA_PV,
                        aes(x = SACO2_Ca, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "SA CO2/Ca", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SACO2_Ca
set_palette(Kd_1ugL_SACO2_Ca, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SACO2_Ca.pdf")

Kd_1ugL_PVCO2_Ca <- ggplot(data = SA_PV,
                        aes(x = PVCO2_Ca, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.0001)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "PV CO2/Ca", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVCO2_Ca
set_palette(Kd_1ugL_PVCO2_Ca, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PVCO2_Ca.pdf")

Kd_1ugL_SAN2_Ca <- ggplot(data = SA_PV,
                       aes(x = SAN2_Ca, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "SA N2/Ca", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SAN2_Ca
set_palette(Kd_1ugL_SAN2_Ca, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SAN2_Ca.pdf")

Kd_1ugL_PVN2_Ca <- ggplot(data = SA_PV,
                       aes(x = PVN2_Ca, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.00005)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "PV N2/Ca", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVN2_Ca
set_palette(Kd_1ugL_PVN2_Ca, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PVN2_Ca.pdf")

#SA PV Fe
SA_PV <- filter(Elements_ratios_1ugL, Use == 1)

Kd_1ugL_SACO2_Fe <- ggplot(data = SA_PV,
                           aes(x = SACO2_Fe, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "SA CO2/Fe", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SACO2_Fe
set_palette(Kd_1ugL_SACO2_Fe, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SACO2_Fe.pdf")

Kd_1ugL_PVCO2_Fe <- ggplot(data = SA_PV,
                           aes(x = PVCO2_Fe, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.0001)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "PV CO2/Fe", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVCO2_Fe
set_palette(Kd_1ugL_PVCO2_Fe, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PVCO2_Fe.pdf")

Kd_1ugL_SAN2_Fe <- ggplot(data = SA_PV,
                          aes(x = SAN2_Fe, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "SA N2/Fe", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SAN2_Fe
set_palette(Kd_1ugL_SAN2_Fe, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SAN2_Fe.pdf")

Kd_1ugL_PVN2_Fe <- ggplot(data = SA_PV,
                          aes(x = PVN2_Fe, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.00005)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "PV N2/Fe", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVN2_Fe
set_palette(Kd_1ugL_PVN2_Fe, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PVN2_Fe.pdf")
