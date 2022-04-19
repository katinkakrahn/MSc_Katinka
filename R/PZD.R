# Library ----
library(data.table)
library(readxl)
library(latex2exp)
library(ggtext)
library(scales)
library(writexl)
library(tidyverse)

#Data cleaning ----
PZD <- read_xlsx("R/data_raw/080422_PZD.xlsx")
PZD <- as.data.table(PZD)
PZD$SA_PV <- PZD$SA/PZD$PV
PZD$SA_PV_C <- PZD$SA_PV/PZD$C


# Pivot PZD ----
scaleX <- function(x) sprintf("%.1f", x)
scaleY <- function(y) sprintf("%.0f", y)

PZD_label <- PZD %>% 
  mutate(SA_PV = log10(SA_PV),
         SA_PV_C = log10(SA_PV_C)) %>% 
  filter(Gas == "N2") %>% 
  mutate(id = row_number()) %>% 
  select(SA, PV, SA_PV, SA_PV_C, Pore_size, Biochar) %>%
  pivot_longer(c(where(is.numeric), -Biochar, -Pore_size)) %>% 
  filter(name != "SA_PV_C" | Biochar != "CWC")

PZD_label$name <- factor(PZD_label$name,
                         levels = c("SA","PV","SA_PV", "SA_PV_C"),
                         labels = c("SA~(m^2/g)", "PV~(cm^3/g)", "log~SA/PV", "log~SA/PV/C"))

PZD_SAPV_C_plot <- PZD_label %>% 
  ggplot(aes(
    y = value,
    x = Pore_size,
    color = Biochar
  )) +
  labs(x = "Pore size (mm)", y = NULL, color = "", shape = "") +
  geom_point(size = 2) +
  facet_wrap(.~ name,
             scales = "free_y",
             labeller = label_parsed,
             strip.position = "left") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values=c("#FFB547FF","#4E9C81","#40E0CF")) +
  scale_x_continuous(breaks=c(1,3,5,10,20,30)) +
  scale_y_continuous() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom", 
        text = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_blank())
PZD_SAPV_C_plot
ggsave("R/figs/PZD_SAPV_C_large.pdf")

PZD_SAPV_C_plot_nolabel <- PZD_label %>% 
  ggplot(aes(
    y = value,
    x = Pore_size,
    color = Biochar
  )) +
  labs(x = "Pore size (mm)", y = NULL, color = "", shape = "") +
  geom_point(size = 2) +
  facet_wrap(.~ name,
             scales = "free_y",
             labeller = label_parsed,
             strip.position = "left") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values=c("#FFB547FF","#4E9C81","#40E0CF")) +
  scale_x_continuous(breaks=c(1,3,5,10,20,30)) +
  scale_y_continuous() +
  theme_bw() +
  guides(color = "none") +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom", 
        text = element_text(size = 20),
        strip.placement = "outside",
        strip.background = element_blank())
PZD_SAPV_C_plot_nolabel

# Small pores
PZD_label_small <- PZD %>% 
  mutate(SA_PV = log10(SA_PV)) %>% 
  filter(Gas == "CO2") %>% 
  mutate(id = row_number()) %>% 
  select(SA, PV, SA_PV, Pore_size, Biochar) %>%
  pivot_longer(c(where(is.numeric), -Biochar, -Pore_size)) %>% 
  filter(name != "SA_PV_C" | Biochar != "CWC") 

PZD_label_small$name <- factor(PZD_label_small$name,
                               levels = c("SA","PV","SA_PV"),
                               labels = c("SA~(m^2/g)", "PV~(cm^3/g)", "log~SA/PV"))

PZD_SAPV_C_small_plot <- PZD_label_small %>% 
  ggplot(aes(
    y = value,
    x = Pore_size,
    color = Biochar
  )) +
  labs(x = "Pore size (mm)", y = NULL, color = "", shape = "") +
  geom_point(size = 2) +
  facet_wrap(.~ name,
             scales = "free_y",
             labeller = label_parsed,
             strip.position = "left") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values=c("#FFB547FF","#4E9C81","#40E0CF")) +
  scale_y_continuous() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom", 
        text = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_blank())
PZD_SAPV_C_small_plot
ggsave("R/figs/PZD_SAPV_C_small_plot.pdf")

# Individual plots ----
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

PV_small <- ggplot(data = subset(PZD, Gas %in% "CO2"), mapping = aes(x = Pore_size, y = log10(PV), color = Biochar)) +
  geom_point(size = 4) + 
  labs(x = "Pore size (nm)", y = TeX(r'(log PV $(cm^{3}~g^{-1})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 30))
PV_small

PV_large <- ggplot(data = subset(PZD, Gas %in% "N2"), mapping = aes(x = Pore_size, y = log10(PV), color = Biochar)) +
  geom_point() + 
  labs(x = "Pore size (nm)", y = TeX(r'(log PV $(cm^{3}~g^{-1})$)'), color = "") +
  theme_bw() +
  scale_x_continuous(breaks=c(1,3,5,10,20,30)) +
  # guides(color = "none") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
PV_large

SAPV_large <- ggplot(data = subset(PZD, Gas %in% "N2"), mapping = aes(x = Pore_size, y = log10(SA_PV), color = Biochar)) +
  geom_point() + 
  labs(x = "Pore size (nm)", y = TeX(r'(log SA/PV $(m^{2}~cm^{-3})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  scale_x_continuous(breaks=c(1,5,10,20,30)) +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
SAPV_large

SAPV_small <- ggplot(data = subset(PZD, Gas %in% "CO2"), mapping = aes(x = Pore_size, y = log10(SA_PV), color = Biochar)) +
  geom_point(size = 4) + 
  labs(x = "Pore size (nm)", y = TeX(r'(log SA/PV $(m^{2}~cm^{-3})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 30))
SAPV_small

PZD_ULS_DSL <- filter(PZD, Biochar != "CWC")
SAPV_large_ULS_DSL <- ggplot(data = subset(PZD_ULS_DSL, Gas %in% "N2"), mapping = aes(x = Pore_size, y = SA_PV, color = Biochar)) +
  geom_point() + 
  labs(x = "Pore size (nm)", y = TeX(r'(log SA/PV $(m^{2} cm^{-3})$)'), color = "") +
  theme_bw() +
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
SAPV_large_ULS_DSL

SAPV_C_large <- ggplot(data = subset(PZD_ULS_DSL, Gas %in% "N2"), mapping = aes(x = Pore_size, y = log10(SA_PV_C), color = Biochar)) +
  geom_point() + 
  scale_x_continuous(breaks=c(1,5,10,20,30)) +
  labs(x = "Pore size (nm)", y = TeX(r'(log (SA/PV)/C ($m^2~cm^{-3})$)'), color = "") +
  theme_bw() +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),values=c("#767676FF","#800000FF","#FFB547FF"))+
  # guides(color = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
SAPV_C_large

