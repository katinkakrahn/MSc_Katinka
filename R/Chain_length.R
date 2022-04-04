#Summary statistics chain length and log KF
nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

summary_stats_CLandKd <- data.table(n = rep(0, nr_biochars),
                                       r_squared = rep(0, nr_biochars),
                                       p_value = rep(0, nr_biochars),
                                       biochar = biochars
                                       )

for(i in 1:nr_biochars){
  fit <- lm(log_Kd ~ nr_CF2, data = Sorption_BC_single_C3[Biochar == biochars[i]])
  summary_stats_CLandKd[Biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_CLandKd[Biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_CLandKd[Biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_CLandKd_label <- summary_stats_CLandKF %>%
  mutate(
    nr_CF2 = 6, K_F = 2.7,
    label =
      glue("*r*<sup>2</sup> = {round(r_squared, 2)} <br> n = {round(n, 2)} <br> p-value = {round(p_value, 3)}")
  )

summary_stats_CLandKF <- data.table(n = rep(0, nr_biochars),
                                    r_squared = rep(0, nr_biochars),
                                    p_value = rep(0, nr_biochars),
                                    biochar = biochars
)

for(i in 1:nr_biochars){
  fit <- lm(log_KF ~ nr_CF2, data = summary_stats_single[biochar == biochars[i]])
  summary_stats_CLandKF[biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_CLandKF[biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_CLandKF[biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                              summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_CLandKF_label <- summary_stats_CLandKF %>%
  mutate(
    nr_CF2 = 7, K_F = 3,
    label =
      glue("*r*<sup>2</sup> = {round(r_squared, 2)} <br> p = {round(p_value, 3)}")
  )

summary_stats_single$biochar <- factor(summary_stats_single$biochar, levels = c("ULS", "DSL", "CWC"))
summary_stats_CLandKF_label$biochar <- factor(summary_stats_CLandKF_label$biochar, levels = c("ULS", "DSL", "CWC"))

Sorption_BC_single_C3[Compound=="PFPeA", nr_CF2:=4]
Sorption_BC_single_C3[Compound=="PFHxA", nr_CF2:=5]
Sorption_BC_single_C3[Compound=="PFHpA", nr_CF2:=6]
Sorption_BC_single_C3[Compound=="PFOA", nr_CF2:=7]
Sorption_BC_single_C3[Compound=="PFNA", nr_CF2:=8]
Sorption_BC_single_C3[Compound=="PFDA", nr_CF2:=9]

# ChainLength_Kd <- ggplot(data = Sorption_BC_single_C3) +
#   geom_point(mapping = aes(x = nr_CF2, y = log_Kd, group = Biochar), color = "grey45", size = 1) + 
#   geom_smooth(mapping = aes(x = nr_CF2, y = log_Kd, group = Biochar), color = "black", 
#               formula = y ~ x, method=lm, 
#               se = TRUE, fullrange = TRUE) +
#   labs(x = expression(CF[2]~chain~length), y = expression(log~K[d])) +
#   facet_wrap(~ Biochar) +
#   geom_richtext(
#     data = summary_stats_CLandKF_label,
#     aes(label = label, x = nr_CF2, y = K_F),
#     hjust = 0
#   ) +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   theme(panel.grid = element_blank()) +
#   guides(size = "none", fill = "none")
# ChainLength_KF
# ggsave(filename = "R/figs/chainlength_KF.pdf")

ChainLength_KF <- ggplot(data = summary_stats_single) +
  geom_point(mapping = aes(x = nr_CF2, y = log_KF, group = biochar), color = "grey45", size = 1) + 
  geom_smooth(mapping = aes(x = nr_CF2, y = log_KF, group = biochar), color = "black", 
              formula = y ~ x, method=lm, 
              se = TRUE, fullrange = TRUE) +
  geom_errorbar(aes(nr_CF2, ymin=log_KF-log_KF_std_error, ymax=log_KF+log_KF_std_error), width=.2,
                position=position_dodge(.9), color = "grey45") +
  labs(x = expression(CF[2]~chain~length), y = expression(log~K[F])) +
  facet_wrap(~ biochar) +
  geom_richtext(
    data = summary_stats_CLandKF_label,
    aes(label = label, x = nr_CF2, y = K_F),
    hjust = 0
  ) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  theme(panel.grid = element_blank()) +
  guides(size = "none", fill = "none")
ChainLength_KF
ggsave(filename = "R/figs/chainlength_KF.pdf")

ChainLength_KF_line <- ggplot(data = summary_stats_single, mapping = aes(x = nr_CF2, y = log_KF, color = biochar)) +
  geom_point(size = 1) + 
  geom_errorbar(aes(nr_CF2, ymin=log_KF-log_KF_std_error, ymax=log_KF+log_KF_std_error), width=.05,
                position=position_dodge(.9), color = "grey45") +
  geom_line(size = 1) + 
  labs(x = expression(CF[2]~chain~length), y = expression(log~K[F])) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 25)) +
  guides(size = "none", fill = "none")
ChainLength_KF_line
set_palette(ChainLength_KF_line, "uchicago")
ggsave(filename = "R/figs/chainlength_KF_line.pdf")

n_KF <- ggplot(data = summary_stats_single, mapping = aes(x = nr_CF2, y = n, color = biochar)) +
  geom_point(size = 1) + 
  geom_errorbar(aes(nr_CF2, ymin=n-n_std_error, ymax=n+n_std_error), width=.05,
                position=position_dodge(.9), color = "grey45") +
  geom_line(size = 1) +
  labs(x = expression(CF[2]~chain~length), y = expression(n[F]), color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom", text = element_text(size = 25)) +
  guides(size = "none", fill = "none")
n_KF
set_palette(n_KF, "uchicago")
ggsave(filename = "R/figs/n_KF.pdf")
