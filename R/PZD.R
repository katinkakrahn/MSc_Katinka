PZD <- read_xlsx("R/data_raw/080422_PZD.xlsx")
PZD <- as.data.table(PZD)
PZD$SA_PV <- PZD$SA/PZD$PV
PZD$SA_PV_C <- PZD$SA_PV/PZD$C

SA_all <- ggplot(data = PZD, mapping = aes(x = Pore_size, y = SA, color = Biochar, shape = Gas)) +
  geom_point() +
  scale_x_log10() +
  labs(x = "Pore size (nm)", y = TeX(r'(Surface area $(m^{2} g^{-1})$)'), color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))

PV_all <- ggplot(data = PZD, mapping = aes(x = Pore_size, y = PV, color = Biochar, shape = Gas)) +
  geom_point() + 
  scale_x_log10() +
  labs(x = "Pore size (nm)", y = TeX(r'(Pore volume $(cm^{3} g^{-1})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))

SA_small <- ggplot(data = subset(PZD, Gas %in% "CO2"), mapping = aes(x = Pore_size, y = log10(SA), color = Biochar)) +
  geom_point(size = 4) +
  #scale_y_log10() +
  labs(x = "Pore size (nm)", y = TeX(r'(log SA $(m^{2}~g^{-1})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 30))
SA_small
ggsave(filename="R/figs/SA_small.pdf")

SA_large <- ggplot(data = subset(PZD, Gas %in% "N2"), mapping = aes(x = Pore_size, y = log10(SA), color = Biochar)) +
  geom_point() +
  scale_x_continuous(breaks=c(1,2,3,4,5,10,20,30)) +
  labs(x = "Pore size (nm)", y = TeX(r'(log SA $(m^{2}~g^{-1})$)'), color = "") +
  theme_bw() +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20), 
        axis.text.x = element_text(angle=0))
SA_large
ggsave(filename="R/figs/SA_large.pdf")

PV_small <- ggplot(data = subset(PZD, Gas %in% "CO2"), mapping = aes(x = Pore_size, y = log10(PV), color = Biochar)) +
  geom_point(size = 4) + 
  labs(x = "Pore size (nm)", y = TeX(r'(log PV $(cm^{3}~g^{-1})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 30))
PV_small
ggsave(filename="R/figs/PV_small.pdf")

PV_large <- ggplot(data = subset(PZD, Gas %in% "N2"), mapping = aes(x = Pore_size, y = log10(PV), color = Biochar)) +
  geom_point() + 
  labs(x = "Pore size (nm)", y = TeX(r'(log PV $(cm^{3}~g^{-1})$)'), color = "") +
  theme_bw() +
  scale_x_continuous(breaks=c(1,3,5,10,20,30)) +
  # guides(color = "none") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
PV_large
ggsave(filename="R/figs/PV_large.pdf")

SAPV_large <- ggplot(data = subset(PZD, Gas %in% "N2"), mapping = aes(x = Pore_size, y = log10(SA_PV), color = Biochar)) +
  geom_point() + 
  labs(x = "Pore size (nm)", y = TeX(r'(log SA/PV $(m^{2}~cm^{-3})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  scale_x_continuous(breaks=c(1,5,10,20,30)) +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
SAPV_large
ggsave(filename="R/figs/SAPV_large.pdf")

SAPV_small <- ggplot(data = subset(PZD, Gas %in% "CO2"), mapping = aes(x = Pore_size, y = log10(SA_PV), color = Biochar)) +
  geom_point(size = 4) + 
  labs(x = "Pore size (nm)", y = TeX(r'(log SA/PV $(m^{2}~cm^{-3})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 30))
SAPV_small
ggsave(filename="R/figs/SAPV_small.pdf")

PZD_ULS_DSL <- filter(PZD, Biochar != "CWC")
SAPV_large_ULS_DSL <- ggplot(data = subset(PZD_ULS_DSL, Gas %in% "N2"), mapping = aes(x = Pore_size, y = SA_PV, color = Biochar)) +
  geom_point() + 
  labs(x = "Pore size (nm)", y = TeX(r'(log SA/PV $(m^{2} cm^{-3})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
ggsave(filename="R/figs/SAPV_large_ULS_DSL.pdf")

SAPV_C_large <- ggplot(data = subset(PZD_ULS_DSL, Gas %in% "N2"), mapping = aes(x = Pore_size, y = log10(SA_PV_C), color = Biochar)) +
  geom_point() + 
  scale_x_continuous(breaks=c(1,5,10,20,30)) +
  labs(x = "Pore size (nm)", y = TeX(r'(log (SA/PV)/C ($m^2~cm^{-3})$)'), color = "") +
  theme_bw() +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
SAPV_C_large
ggsave(filename="R/figs/SAPV_C_large.pdf")

