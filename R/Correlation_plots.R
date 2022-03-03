ChainLength_KF <- ggplot(data = summary_stats_single) +
  geom_point(mapping = aes(x = nr_CF2, y = K_F, color = biochar)) + 
  geom_smooth(mapping = aes(x = nr_CF2, y = K_F, color = biochar, fill = biochar), formula = y ~ x, method=lm, fullrange = TRUE) + 
  labs(x = expression(number~of~CF[2]~moieties), y = expression(K[F]), col = "Biochar", title = "Chain length vs KF") + 
  stat_regline_equation(
    aes(x = nr_CF2, y = K_F, color = biochar, label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x) +
  theme_bw() +
  guides(size = "none", fill = "none")
ChainLength_KF
