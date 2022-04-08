PZD <- read_xlsx("/Users/katinkakrahn/Library/CloudStorage/OneDrive-NGI/VOW/Data/080422_PZD.xlsx")
PZD <- as.data.table(PZD)
PZD$SA_PV <- PZD$SA/PZD$PV
PZD$SA_PV_C <- PZD$SA_PV/PZD$C

SA_all <- ggplot(data = PZD, mapping = aes(x = Pore_size, y = SA, color = Biochar, shape = Gas)) +
  geom_point() +
  scale_x_log10() +
  labs(x = "Pore size (nm)", y = TeX(r'(Surface area $(m^{2} g^{-1})$)'), color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(SA_all, "uchicago")

PV_all <- ggplot(data = PZD, mapping = aes(x = Pore_size, y = PV, color = Biochar, shape = Gas)) +
  geom_point() + 
  scale_x_log10() +
  labs(x = "Pore size (nm)", y = TeX(r'(Pore volume $(cm^{3} g^{-1})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(PV_all, "uchicago")

SA_small <- ggplot(data = subset(PZD, Gas %in% "CO2"), mapping = aes(x = Pore_size, y = SA, color = Biochar)) +
  geom_point() +
  scale_x_log10() +
  #scale_y_log10() +
  labs(x = "Pore size (nm)", y = TeX(r'(Surface area $(m^{2} g^{-1})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(SA_small, "uchicago")
ggsave(filename="R/figs/SA_small.pdf")

SA_large <- ggplot(data = subset(PZD, Gas %in% "N2"), mapping = aes(x = Pore_size, y = SA, color = Biochar)) +
  geom_point() +
  #scale_x_log10() +
  #scale_y_log10() +
  labs(x = "Pore size (nm)", y = TeX(r'(Surface area $(m^{2} g^{-1})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(SA_large, "uchicago")
ggsave(filename="R/figs/SA_large.pdf")

PV_small <- ggplot(data = subset(PZD, Gas %in% "CO2"), mapping = aes(x = Pore_size, y = PV, color = Biochar)) +
  geom_point() + 
  scale_x_log10() +
  labs(x = "Pore size (nm)", y = TeX(r'(Pore volume $(cm^{3} g^{-1})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(PV_small, "uchicago")
ggsave(filename="R/figs/PV_small.pdf")

PV_large <- ggplot(data = subset(PZD, Gas %in% "N2"), mapping = aes(x = Pore_size, y = PV, color = Biochar)) +
  geom_point() + 
  scale_x_log10() +
  labs(x = "Pore size (nm)", y = TeX(r'(Pore volume $(cm^{3} g^{-1})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(PV_large, "uchicago")
ggsave(filename="R/figs/PV_large.pdf")

PZD_ULS_DSL <- filter(PZD, Biochar != "CWC")
SAPV_large <- ggplot(data = subset(PZD_ULS_DSL, Gas %in% "N2"), mapping = aes(x = Pore_size, y = SA_PV, color = Biochar)) +
  geom_point() + 
  labs(x = "Pore size (nm)", y = TeX(r'(SA/PV $(m^{2} cm^{-3})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(SAPV_large, "uchicago")
ggsave(filename="R/figsSAPV_large.pdf")

SAPV_C_large <- ggplot(data = subset(PZD_ULS_DSL, Gas %in% "N2"), mapping = aes(x = Pore_size, y = SA_PV_C, color = Biochar)) +
  geom_point() + 
  labs(x = "Pore size (nm)", y = TeX(r'((SA/PV)/C)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
set_palette(SAPV_C_large, "uchicago")
ggsave(filename="R/figs/SAPV_C_large.pdf")
