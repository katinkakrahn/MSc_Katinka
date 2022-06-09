# Library ----
library(data.table)
library(readxl)
library(latex2exp)
library(ggtext)
library(scales)
library(writexl)
library(tidyverse)

#Data cleaning ----
PSD <- read_xlsx("R/data_raw/080422_PZD.xlsx") 
PSD <- as.data.table(PSD)
PSD$SA_PV <- PSD$SA/PSD$PV
PSD$SA_PV_C <- PSD$SA_PV/PSD$C

#Pore size distribution ----
PSD_pivot <- PSD %>% 
  mutate(id = row_number()) %>% 
  select(SA, PV, Pore_size, Biochar, Gas) %>%
  pivot_longer(c(where(is.numeric), -Biochar, -Pore_size))

PSD_pivot$name <- factor(PSD_pivot$name,
                         levels = c("SA","PV"),
                         labels = c("SA~(m^2/g)", "PV~(cm^3/g)"))
PSD_pivot$name <- paste(PSD_pivot$Gas, "-", PSD_pivot$name)
PSD_pivot$name <- factor(PSD_pivot$name,
                         levels = c("N2 - SA~(m^2/g)", "N2 - PV~(cm^3/g)", "CO2 - SA~(m^2/g)", "CO2 - PV~(cm^3/g)"),
                         labels = c("N[2]~SA~(m^2/g)", "N[2]~PV~(cm^3/g)", "CO[2]~SA~(m^2/g)", "CO[2]~PV~(cm^3/g)"))

PSD_plot <- PSD_pivot %>% 
  mutate(name = factor(name, 
                       levels = c("CO[2]~SA~(m^2/g)", "CO[2]~PV~(cm^3/g)", "N[2]~SA~(m^2/g)", "N[2]~PV~(cm^3/g)" 
                       )))%>% 
  ggplot(aes(
    y = value,
    x = Pore_size,
    color = Biochar
  )) +
  labs(x = "Pore diameter (nm)", y = NULL, color = "", shape = "") +
  geom_point(size = 2) +
  facet_wrap(.~ name,
             scales = "free",
             labeller = label_parsed,
             strip.position = "left") +
  scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
                     values=c("#FFB547FF","#4E9C81","#40E0CF")) +
  #scale_x_continuous(breaks=c(1,3,5,10,20,30)) +
  scale_y_continuous() +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        strip.placement = "outside",
        strip.background = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4)))
PSD_plot
ggsave("R/figs/article/PSD_plot.pdf")
