# Library ----
library(data.table)
library(ggplot2)
library(psych)
library(readxl)
library(grDevices)
library(dplyr)
library(knitr)
library(latex2exp)
library(glue)
library(ggtext)
library(jcolors)
library(plotrix)
library(scales)
library(writexl)
library(ggbreak) 
library(patchwork)
library(scales)
library(mapproj)
library(tidyverse)

# Ca ----
Kd_1ugL_Ca <- ggplot(data = subset(Elements_ratios_1ugL, Parameter %in% "Ca"),
                                aes(x = log10(Mean_sameunit), y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point(size = 8) +
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.01)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) + 
  labs(x = TeX(r'(log Ca $(g~kg^{-1})$)'), y = TeX(r'($log~K_d~(at~C_w~1 \mu g/L)$)'), color = "", shape = "") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme_bw() +
  #guides(shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_Ca
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
  #guides(color = "none", shape = "none") +
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
Kd_1ugL_SAN2
set_palette(Kd_1ugL_SAN2, "uchicago")
ggsave(filename="R/figs/Kd_1ugL_SAN2.pdf")

Kd_1ugL_PVN2 <- ggplot(data = SA_PV,
                        aes(x = log10(N2_PV), y = log_Kd, shape = compound, color = biochar), 
) +
  geom_point(size = 8) +
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.001)+
  labs(x = TeX(r'(log PV $(cm^{3}~g^{-1})$)'), y = TeX(r'($log~K_d~(at~C_w~1 \mu g/L)$)'), color = "", shape = "") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 2) +
  theme_bw() +
  guides(shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 30))
Kd_1ugL_PVN2
ggsave(filename="R/figs/Kd_1ugL_PVN2.pdf")


Kd_1ugL_PVN2_Ca <- ggplot(data = SA_PV,
                       aes(x = log10(PVN2_Ca), y = log_Kd, shape = compound, color = biochar), 
) +
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.005)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 8) +
  labs(x = TeX(r'(log PV/Ca)'), y = TeX(r'($log~K_d~(at~C_w~1 \mu g/L)$)'), color = "", shape = "") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme_bw() +
  #guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_PVN2_Ca
ggsave(filename="R/figs/Kd_1ugL_PVN2_Ca.pdf")

Kd_1ugL_SA_PV_Ca <- ggplot(data = SA_PV,
                          aes(x = log10(SA_PV_Ca), y = log_Kd, shape = compound, color = biochar), 
) +
  geom_errorbar(aes(ymin=log_Kd-logKd_error, ymax=log_Kd+logKd_error), color = "grey", width=.03)+
  geom_line(aes(group = compound), color = "black") +
  geom_point(size = 8) +
  labs(x = TeX(r'(log (SA/PV)/Ca)'), y = TeX(r'($log~K_d~(at~C_w~1 \mu g/L)$)'), color = "", shape = "") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme_bw() +
  guides(color = "none", shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
Kd_1ugL_SA_PV_Ca
ggsave(filename="R/figs/Kd_1ugL_SA_PV_Ca.pdf")


