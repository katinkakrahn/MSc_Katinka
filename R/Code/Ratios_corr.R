# Not cleaned
summary_stats_single2 <- setnames(summary_stats_single, c("compound", "biochar"), 
                                  c("Compound", "Biochar"))
Biochar_Kd_merge <- merge(Sorption_BC_single_C3, summary_stats_single2, by = c("Compound", "Biochar"))

#Ratio with respect to CWC
CWC_C3 <- filter(Biochar_Kd_merge, Biochar == "CWC")
ULS_C3 <- filter(Biochar_Kd_merge, Biochar == "ULS")
DSL_C3 <- filter(Biochar_Kd_merge, Biochar == "DSL")

CWC_ULS <- inner_join(CWC_C3, ULS_C3, by = "Compound")
CWC_DSL <- inner_join(CWC_C3, DSL_C3, by = "Compound")

Biochar_Kd_ratio_merge <- rbind(CWC_ULS, CWC_DSL)
  
Biochar_Kd_ratio <- Biochar_Kd_ratio_merge[, .(log_Kd_ratio = (log_Kd.x/log_Kd.y),
                                                    nr_CF2 = nr_CF2.x
                                                   ),
                                                keyby = .(Biochar.y, Compound)]
Biochar_Kd_ratio <- filter(Biochar_Kd_ratio, Compound != "PFPeA")
Biochar_Kd_ratio[Biochar.y=="DSL", Biochar.y:= "CWC/DSL"]
Biochar_Kd_ratio[Biochar.y=="ULS", Biochar.y:= "CWC/ULS"]

ChainLength_Kd_ratio_plot <- ggplot(data = Biochar_Kd_ratio, 
                                aes(x = nr_CF2, y = log_Kd_ratio, color = Biochar.y)) +
  geom_point() +
  labs(x = expression(nr~CF[2]~moieties), y = expression(log~K[d]~ratio), color = "") +
  geom_smooth(data = Biochar_Kd_ratio,
              aes(x = nr_CF2, y = log_Kd_ratio, group = Biochar.y),
              formula = y ~ x,
              method=lm,
              se=FALSE,
              color = "grey") +
  theme_bw()
ChainLength_Kd_ratio_plot
#Sorption of DSL increases more with chain length than ULS, which is hard to explain...

#Ratio with respect to ULS
ULS_CWC <- inner_join(ULS_C3, CWC_C3, by = "Compound")
ULS_DSL <- inner_join(ULS_C3, DSL_C3, by = "Compound")

Biochar_Kd_ratio_merge2 <- rbind(ULS_CWC, ULS_DSL)

Biochar_Kd_ratio2 <- Biochar_Kd_ratio_merge2[, .(log_Kd_ratio = (log_Kd.x/log_Kd.y),
                                               nr_CF2 = nr_CF2.x
),
keyby = .(Biochar.y, Compound)]
Biochar_Kd_ratio2 <- filter(Biochar_Kd_ratio2, Compound != "PFPeA")
Biochar_Kd_ratio2[Biochar.y=="DSL", Biochar.y:= "ULS/DSL"]
Biochar_Kd_ratio2[Biochar.y=="CWC", Biochar.y:= "ULS/CWC"]

ChainLength_Kd_ratio_plot2 <- ggplot(data = Biochar_Kd_ratio2, 
                                    aes(x = nr_CF2, y = log_Kd_ratio, color = Biochar.y)) +
  geom_point() +
  labs(x = expression(nr~CF[2]~moieties), y = expression(log~K[d]~ratio), color = "") +
  geom_smooth(data = Biochar_Kd_ratio2,
              aes(x = nr_CF2, y = log_Kd_ratio, group = Biochar.y),
              formula = y ~ x,
              method=lm,
              se=FALSE,
              color = "grey") +
  theme_bw()
ChainLength_Kd_ratio_plot2

Elements_ratios <- read_excel("/Users/katinkakrahn/Library/CloudStorage/OneDrive-NGI/VOW/Data/270322_Elements_ratios.xlsx")
Elements_ratios <- as.data.table(Elements_ratios)

Elements_ratios_Kd <- full_join(x = Elements_ratios, y = Sorption_BC_single_C3, by = "Biochar")
Elements_ratios_Kd_PFOA <- filter(Elements_ratios_Kd, Compound == "PFOA")

Elements_ratios_plot_PFOA <- ggplot(data = Elements_ratios_Kd_PFOA,
             aes(x = Ratio, y = log_Kd, color = Parameter, shape = Biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  labs(x = "Element ratio", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Elements_ratios_plot_PFOA

Elements_ratios_plot_PFOA <- ggplot(data = Elements_ratios_Kd,
                                    aes(x = Ratio, y = log_Kd, color = Compound, shape = Biochar), 
) +
  geom_point() +
  # geom_smooth(method = "lm",
  #             formula = y ~ x)+
  labs(x = "Element ratio", y = expression(log~K[d]), color = "", shape = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
Elements_ratios_plot_PFOA

Elements_ratios_Kd_1ugL <- full_join(x = Elements_ratios, y = Sorption_BC_single_C3, by = "Biochar")

