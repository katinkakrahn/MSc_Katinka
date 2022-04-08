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
  guides(shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(Kd_1ugL_Ca, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_Ca.pdf")

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
  #guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(Kd_1ugL_C, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_C.pdf")

Kd_1ugL_SAN2 <- ggplot(data = SA_PV,
                       aes(x = N2_SA, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.05)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  labs(x = "SA N2", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  #guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
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
  guides(shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(Kd_1ugL_PVN2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_PVN2.pdf")


Kd_1ugL_SA_C <- ggplot(data = SA_PV,
                       aes(x = SA_C, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.00005)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "SA/C", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  #guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(Kd_1ugL_SA_C, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SA_C.pdf")

Kd_1ugL_SA_PV <- ggplot(data = SA_PV,
                        aes(x = SA_PV, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.00005)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "SA/PV", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  #guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(Kd_1ugL_SA_PV, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SA_PV.pdf")

Kd_1ugL_SA_PV_C <- ggplot(data = SA_PV,
                          aes(x = SA_PV_C, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.00005)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "(SA/PV)/C", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(Kd_1ugL_SA_PV_C, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SA_PV_C.pdf")

Kd_1ugL_SA_PV_Ca <- ggplot(data = SA_PV,
                           aes(x = SA_PV_Ca, y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.00005)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "(SA/PV)/Ca", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(Kd_1ugL_SA_PV_Ca, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SA_PV_Ca.pdf")

