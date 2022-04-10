Kd_1ugL_Ca_C <- ggplot(data = Biochar_ratios_1ugL_select,
                       aes(x = Ca_C, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.01)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) + 
  scale_x_log10() +
  labs(x = "Ca/C ratio", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  # guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_Ca_C
set_palette(Kd_1ugL_Ca_C, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_Ca_C.pdf")

Kd_1ugL_Fe_C <- ggplot(data = Biochar_ratios_1ugL_select,
                       aes(x = Fe_C, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.01)+
  geom_line(aes(group = compound), color = "black") +
  scale_x_log10() +
  geom_point(size = 2) +
  labs(x = "Fe/C ratio", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_Fe_C
set_palette(Kd_1ugL_Fe_C, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_Fe_C.pdf")

Kd_1ugL_Ca_SACO2 <- ggplot(data = SA_PV,
                           aes(x = Ca_SACO2, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "Ca_SA CO2/", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_Ca_SACO2
set_palette(Kd_1ugL_Ca_SACO2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_Ca_SACO2.pdf")

Kd_1ugL_Ca_PVCO2 <- ggplot(data = SA_PV,
                           aes(x = log10(Ca_PVCO2), y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point(size = 8) +
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.0001)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "Ca_PV CO2", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 30))
Kd_1ugL_PVCO2_Ca
set_palette(Kd_1ugL_Ca_PVCO2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_Ca_PVCO2.pdf")

Kd_1ugL_Ca_SAN2 <- ggplot(data = SA_PV,
                          aes(x = Ca_SAN2, y = log_Kd, shape = compound, color = biochar), 
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
Kd_1ugL_Ca_SAN2
set_palette(Kd_1ugL_Ca_SAN2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_Ca_SAN2.pdf")

Kd_1ugL_SACO2_Fe <- ggplot(data = SA_PV,
                           aes(x = SACO2_Fe, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  scale_x_log10() +
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
  scale_x_log10() +
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
  scale_x_log10() +
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
  scale_x_log10() +
  labs(x = "PV N2/Fe", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  scale_x_log10() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVN2_Fe
set_palette(Kd_1ugL_PVN2_Fe, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PVN2_Fe.pdf")


Kd_1ugL_PV_C <- ggplot(data = SA_PV,
                       aes(x = PV_C, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.00005)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "PV/C", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(Kd_1ugL_PV_C, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PV_C.pdf")

###############################################################################


# Kd_1ugL_C_PV <- ggplot(data = SA_PV,
#                           aes(x = C_PV, y = log_Kd, shape = compound, color = biochar), 
# ) +
#   geom_point() +
#   # geom_smooth(method = "lm",
#   #             formula = y ~ x)+
#   geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.00005)+
#   geom_line(aes(group = compound), color = "black") +
#   geom_point(size = 2) +
#   scale_x_log10() +
#   labs(x = "C/PV", y = expression(log~K[d]), shape = "") +
#   theme_bw() +
#   guides(shape = "none") +
#   theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
# Kd_1ugL_C_PV
# set_palette(Kd_1ugL_C_PV, "uchicago")
# ggsave(filename="R/figs/Kd_1ugL_C_PV.pdf")
