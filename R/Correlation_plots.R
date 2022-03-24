library(data.table)
library(ggplot2)
library(psych)
library(readxl)
library(grDevices)
library(dplyr)
library(broom)
library(ggpubr)
library(tidyverse)
library(knitr)
library(latex2exp)

#Summary statistics chain length and log KF
nr_compounds <- length(unique(Sorption$Compound))
compounds <- unique(Sorption$Compound)
nr_biochars <- length(unique(Sorption$Biochar))
biochars <- unique(Sorption$Biochar)

summary_stats_CLandKF <- data.table(n = rep(0, nr_biochars),
                                       r_squared = rep(0, nr_biochars),
                                       p_value = rep(0, nr_biochars),
                                       biochar = biochars
                                       )

for(i in 1:nr_biochars){
  fit <- lm(K_F ~ nr_CF2, data = summary_stats_single[biochar == biochars[i]])
  summary_stats_CLandKF[biochar == biochars[i], n := fit$coefficients[2]]
  summary_stats_CLandKF[biochar == biochars[i], r_squared := summary(fit)$r.squared]
  summary_stats_CLandKF[biochar == biochars[i], p_value := pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],
                                                                   summary(fit)$fstatistic[3],lower.tail=F)]
}

summary_stats_CLandKF_label <- summary_stats_CLandKF %>%
  mutate(
    nr_CF2 = 6, K_F = 2.7,
    label =
      glue("*r*<sup>2</sup> = {round(r_squared, 2)} <br> n = {round(n, 2)} <br> p-value = {round(p_value, 3)}")
  )

# curve fitting - looks like a linear model fits best, try to use nls function?
#fit2 <- nls(data = summary_stats_CWC_single, K_F~C*(1-exp(k*nr_CF2)))
#summary(fit2)

summary_stats_single$biochar <- factor(summary_stats_single$biochar, levels = c("ULS", "DSL", "CWC"))
summary_stats_CLandKF_label$biochar <- factor(summary_stats_CLandKF_label$biochar, levels = c("ULS", "DSL", "CWC"))

ChainLength_KF <- ggplot(data = summary_stats_single) +
  geom_point(mapping = aes(x = nr_CF2, y = K_F, group = biochar), color = "grey45", size = 1) + 
  geom_smooth(mapping = aes(x = nr_CF2, y = K_F, group = biochar), color = "black", 
              formula = y ~ x, method=lm, 
              se = TRUE, fullrange = TRUE) +
  geom_errorbar(aes(nr_CF2, ymin=K_F-K_F_std_error, ymax=K_F+K_F_std_error), width=.2,
                position=position_dodge(.9), color = "grey45") +
  labs(x = expression(number~of~CF[2]~moieties), y = expression(log~K[F])) +
  facet_wrap(~ biochar) +
  geom_richtext(
    data = summary_stats_CLandKF_label,
    aes(label = label, x = nr_CF2, y = K_F),
    hjust = 0
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(size = "none", fill = "none")
ChainLength_KF
ggsave(filename = "R/figs/chainlength_KF.pdf")
