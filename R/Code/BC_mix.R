Sorption_BC_mix_summary <- Sorption_BC_mix[, .(Conc_point = 10,
                                               Ci = mean(Ci),
                                               Cw = mean(Cw),
                                               Cs = mean(Cs),
                                               Kd = mean(Cs/Cw),
                                               se_Kd = std.error(Cs/Cw),
                                               log_Kd = log10(mean(Cs/Cw)),
                                               log_Cw = mean(log_Cw), 
                                               log_Cs = mean(log_Cs),
                                               se_logCw = std.error(log_Cw),
                                               se_logCs = std.error(log_Cs),
                                               se_logKd = std.error(log10(Cs/Cw)),
                                               mixLogic = TRUE, 
                                               SoilLogic = FALSE
                                               ),
                                           keyby = .(Biochar, Compound)]


#Sorption attenuation competition factor
Sorption_BC_single_C10 <- subset(Sorption_BC_single, Conc_point == 10)
Sorption_BC_single_C10_common <- semi_join(x = Sorption_BC_single_C10, y = Sorption_BC_mix_summary, by = c("Compound" = "Compound", "Biochar" = "Biochar"))
C10_mixVsSingle_BC <- rbind(Sorption_BC_single_C10_common, Sorption_BC_mix_summary, fill = TRUE)
C10_mixVsSingle_BC[mixLogic==TRUE, single_mix:="cocktail"]
C10_mixVsSingle_BC[mixLogic==FALSE, single_mix:="single compound"]

C10_mixVsSingle_BC$Compound <- factor(C10_mixVsSingle_BC$Compound, levels = c("PFPeA", "PFHxA", "PFHpA", 
                                                                                        "PFOA", "PFNA", "PFDA"))


  
C10_mixVsSingle_BC_plot <- ggplot(data = C10_mixVsSingle_BC, aes(x = Compound, y = log_Kd, color = single_mix, shape = Biochar)) + 
  geom_point(size = 4) + 
  geom_errorbar(aes(ymin=log_Kd-se_logKd, ymax=log_Kd+se_logKd), width = 0) + 
  labs(x = "", y = expression(log~K[d]), color = "", shape = "") + 
  theme_bw() +
  scale_color_manual(breaks = c("cocktail", "single compound"),values=c("#0F425CFF","#999999"))+
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 20))
C10_mixVsSingle_BC_plot
ggsave(filename="R/figs/C10_mixVsSingle_BC_plot.pdf")

# Each data point from cocktail (red) represents an average of triplicate batch tests.

Sorption_BC_C10_Kdsingle <- setnames(Sorption_BC_single_C10_common, "log_Kd", "log_Kd_s")
Sorption_BC_C10_Kdmix <- setnames(Sorption_BC_mix_summary, "log_Kd", "log_Kd_mix")
Competition_BC <- full_join(x = Sorption_BC_C10_Kdsingle, y = Sorption_BC_C10_Kdmix, by = c("Biochar", "Compound"))

Competition_factor_BC <- Competition_BC[, .(Conc_point = "10",
                                            Competition_factor = Kd.y/Kd.x*100
                                            ),
                                        keyby = .(Biochar, Compound)]
write_xlsx(Competition_factor_BC, "/Users/katinkakrahn/Library/CloudStorage/OneDrive-NGI/VOW/Data/310322_Competition_factor_BC.xlsx")
write_xlsx(Competition_BC, "/Users/katinkakrahn/Library/CloudStorage/OneDrive-NGI/VOW/Data/310322_Competition_Kd_BC.xlsx" )




