Sorption_isotherms_one <- Sorption_BC_single %>% 
  filter(Compound != "PFPeA" | Biochar != "DSL",
         Compound != "PFPeA" | Biochar != "CWC",
         Compound != "PFHxA" | Biochar != "CWC") %>% 
  group_by(Compound,Biochar) %>% 
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw), y = log10(Cs), group = Compound, shape = Biochar), 
             size = 1) + 
  geom_smooth(mapping = aes(x = log10(Cw), y = log10(Cs), linetype = Compound, color = Biochar),
              formula = y ~ x,
              method=lm,
              se=F,
              fullrange = FALSE) +
  labs(x = TeX(r'($log~C_{w}~(\mu g~L^{-1})$)'), 
       y = TeX(r'($log~C_{s}~(\mu g~kg^{-1})$)'), 
       color = "") + 
  theme_bw() +
  theme(panel.grid = element_blank())
  # scale_color_manual(breaks = c("CWC", "ULS", "DSL"),
  #                    values=c("#FFB547FF","#4E9C81","#40E0CF"))+
  # theme(legend.position = "bottom")
Sorption_isotherms_one
ggsave(filename="R/figs/one_isotherm.pdf")

C10_PFOA <- Soil_BC_join_mean_and_isotherm %>% 
  filter(Conc_point == 10,
         Compound %in% c("PFOA")) %>%
  mutate(Compound = factor(Compound, 
                           levels = c("PFPeA", "PFHxA", "PFHpA", 
                                      "PFOA", "PFNA", "PFDA")),
         Biochar = factor(Biochar,
                          levels = c("no", "CWC", "DSL", "ULS")),
         type = factor(type,
                       levels = c("BC_sing", "BC_S_sing", "BC_S_mix", 
                                  "BC_mix", "S_sing", "S_mix"))) %>% 
  ggplot(aes(x = Biochar,
             y = log10(Kd),
             color = type,
             shape = type
  )) +
  geom_point(size = 5,
             alpha = 1) +
  facet_wrap(~ Compound) +
  labs(x = "", 
       y = TeX(r'($log~K_{d}~(L/kg)$)'),
       color = "",
       shape = ""
  ) +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC-single", 
                                "BC-soil-single", 
                                "BC-soil-mix",
                                "BC-mix",
                                "Soil-single",
                                "Soil-mix")) +
  scale_shape_manual(name = "",
                     labels = c("BC-single", 
                                "BC-soil-single", 
                                "BC-soil-mix",
                                "BC-mix",
                                "Soil-single",
                                "Soil-mix"),
                     values = c(16, 17, 17, 16, 17, 17)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, 
                                   vjust = 0.5, 
                                   hjust=0.5)) #+
  #guides(colour = guide_legend(override.aes = list(size=10)))
C10_PFOA
ggsave(filename = "R/figs/C10_PFOA.pdf")


AF_PFOA <- Attenuation_factors %>% 
  drop_na(Attenuation) %>%
  filter(Biochar %in% c("CWC", "DSL", "ULS"),
         Compound %in% c("PFOA")) %>%
  mutate(Compound = factor(Compound, 
                           levels = c("PFOA", "PFNA", "PFDA")),
         Biochar = factor(Biochar,
                          levels = c("CWC", "DSL", "ULS")),
         type = factor(type,
                       levels = c("BC_sing", "BC_S_sing", 
                                  "BC_S_mix", "BC_mix",
                                  "S_sing", "S_mix"
                       ))) %>% 
  ggplot(aes(x = Biochar,
             y = Attenuation,
             color = type,
             shape = type)) +
  geom_point(size = 5,
             alpha = 0.7) +
  facet_wrap(~ Compound) +
  labs(x = "",
       y = "AF",
       color = "",
       shape = "") +
  scale_color_brewer(palette = "Paired",
                     labels = c("BC-single",
                                "BC-soil-single",
                                "BC-soil-mix",
                                "BC-mix",
                                "Soil-single",
                                "Soil-mix")) +
  scale_shape_manual(name = "",
                     labels = c("BC-single",
                                "BC-soil-single",
                                "BC-soil-mix",
                                "BC-mix",
                                "Soil-single",
                                "Soil-mix"),
                     values = c(16, 17, 17, 16)) +
  scale_y_continuous(breaks = c(25,50,75,100,125,150)) +
  #guides(color = guide_legend(nrow=2, 
                              #override.aes = list(size=10)),
         #byrow=TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, 
                                   vjust = 0.5, 
                                   hjust=0.5))
AF_PFOA
ggsave(filename = "R/figs/AF_PFOA.pdf")