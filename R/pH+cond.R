library(data.table)
library(ggplot2)
library(psych)
library(readxl)
library(grDevices)
library(latex2exp)

pHcond <- read_excel("/Users/katinkakrahn/Library/Mobile Documents/com~apple~CloudDocs/Documents/Skole/VOW/Lab/pH_Conductivity/pH&conductivity_R.xlsx")
as.data.table(pHcond)
pHcond <- as.data.table(pHcond)

#Mean and standard deviation pH and cond
pHcondSummary <- pHcond[, .(mean_ph = mean(pH), 
                            mean_cond = mean(Conductivity),
                            sd_cond = sd(Conductivity),
                            sd_ph = sd(pH)
                            ),
                        keyby = .(Sample)]
pHcondSummary <- pHcondSummary[order(mean_ph),]


#pH plot
pH <- ggplot(data = pHcondSummary, aes(x = reorder(Sample, mean_ph), y = mean_ph)) + 
  geom_point()+ 
  geom_errorbar(aes(ymin=mean_ph-sd_ph, ymax=mean_ph+sd_ph), width=.2,
                       position=position_dodge(.9))+ 
  labs(x = "", y = "pH") + 
  theme_bw()
pH
ggsave(filename="R/figs/pH.pdf")

#cond plot
cond <- ggplot(data = pHcondSummary, aes(x = reorder(Sample,mean_cond), y = mean_cond)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=mean_cond-sd_cond, ymax=mean_cond+sd_cond), width=.2,
                         position=position_dodge(.9)) + 
  labs(x = "", y = TeX(r'(Conductivity $(\mu S cm^{-1})$)')) + 
  theme_bw()
cond
ggsave(filename = "R/figs/conductivity.pdf")


