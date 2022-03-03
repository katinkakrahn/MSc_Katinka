library(data.table)
library(ggplot2)
library(psych)
library(readxl)
library(grDevices)

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

pH_points <- ggplot(data = pHcond, aes(x = reorder(Sample, pH), y = pH)) + 
  geom_point(shape = 21, size = 3, fill = "#077DAA")+ 
  labs(x = "", y = "pH") + 
  theme_bw()
pH_points

#cond plot
cond <- ggplot(data = pHcondSummary, aes(x = reorder(Sample, mean_cond), y = mean_cond))+ 
  geom_point()+ 
  geom_errorbar(aes(ymin=mean_cond-sd_cond, ymax=mean_cond+sd_cond), width=.2,
                         position=position_dodge(.9)) + 
  labs(x = "", y = "cond (\u03bcS/cm)") + 
  theme_bw()
cond

cond_points <- ggplot(data = pHcond, aes(x = reorder(Sample, Conductivity), y = Conductivity)) + 
  geom_point(position = position_jitter(h=0.2, w=0), shape = 21, size = 3, fill = "#077DAA")+ 
  labs(x = "", y = "cond (\u03bcS/cm)") + 
  theme_bw()
cond_points


# ANOVA and Tukey HSD pH
pHcond_BC <- filter(pHcond, Sample == "CWC" | Sample == "ULS" | Sample == "BRL")
pHcond_BC_S <- filter (pHcond, Sample == "CWC+S" | Sample == "ULS+S" | Sample == "BRL+S")

pH_ANOVA_all <- aov(pH ~ Sample, data = pHcond)
summary(pH_ANOVA_all)
#TukeyHSD(pH_ANOVA_all, conf.level = 0.95)
pH_TKHSD <- TukeyHSD(pH_ANOVA_all, "Sample", ordered = TRUE)
pH_TKHSD <- as.data.frame(pH_TKHSD$Sample)

pH_ANOVA_BC <- aov(pH ~ Sample, data = pHcond_BC)
summary(pH_ANOVA_BC)
pH_BC_TKHSD <- TukeyHSD(pH_ANOVA_BC, "Sample", ordered = TRUE)
pH_BC_TKHSD <- as.data.frame(pH_BC_TKHSD$Sample)

pH_ANOVA_S <- aov(pH ~ Sample, data = pHcond_BC_S)
summary(pH_ANOVA_S)
pH_S_TKHSD <- TukeyHSD(pH_ANOVA_S, "Sample", ordered = TRUE)
pH_S_TKHSD <- as.data.frame(pH_S_TKHSD$Sample)

#ANOVA and Tukey HSD cond
cond_ANOVA_all <- aov(Conductivity ~ Sample, data = pHcond)
summary(cond_ANOVA_all)
cond_TKHSD <- TukeyHSD(cond_ANOVA_all, "Sample", ordered = TRUE)
cond_TKHSD <- as.data.frame(cond_TKHSD$Sample)

cond_ANOVA_BC <- aov(Conductivity ~ Sample, data = pHcond_BC)
summary(cond_ANOVA_BC)
cond_BC_TKHSD <- TukeyHSD(cond_ANOVA_BC, "Sample", ordered = TRUE)
cond_BC_TKHSD <- as.data.frame(cond_BC_TKHSD$Sample)

cond_ANOVA_S <- aov(Conductivity ~ Sample, data = pHcond_BC_S)
summary(cond_ANOVA_S)
cond_S_TKHSD <- TukeyHSD(cond_ANOVA_S, "Sample", ordered = TRUE)
cond_S_TKHSD <- as.data.frame(cond_S_TKHSD$Sample)

