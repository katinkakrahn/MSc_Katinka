library(data.table)
library(ggplot2)
library(psych)
library(readxl)
library(grDevices)
library(dplyr)
library(broom)
library(ggpubr)
library(tidyverse)

Sorption <- read_excel("/Users/aurorahofman/Documents/Diverse/hjelp_katinka/170222_sorption_rawdata_removedLOQ.xlsx")
as.data.table(Sorption)
Sorption <- as.data.table(Sorption)

#TRIKS

Sorption[, SoilLogic:= as.logical(Sorption$Soil_binary)]
Sorption[, mixLogic:= as.logical(Sorption$mix_binary)]


#Convert 1 and 0 to TRUE and FALSE and delete integer columns
# Sorption$SoilLogic <- as.logical(Sorption$Soil_binary)
# Sorption$mixLogic <- as.logical(Sorption$mix_binary)
Sorption <- subset(Sorption,select = -c(Soil_binary,mix_binary))

# Subset biochar and cocktail/single compound
Sorption_NAomit <- na.omit(Sorption)
Sorption_LOQomit <- filter(Sorption, between(Conc_point, 2, 10))
Sorption_LOQNAomit <- na.omit(Sorption_LOQomit)
Sorption_BCsingleComp <- subset(Sorption_NAomit, mixLogic == FALSE)
Sorption_BCmixComp <- subset(Sorption_NAomit, mixLogic == TRUE)

CWC_single_all <- filter(Sorption_NAomit, Biochar == "CWC")
ULS_single_all <- filter(Sorption_NAomit, Biochar == "ULS")
BRL_single_all <- filter(Sorption_NAomit, Biochar == "BRL")

CWC_single_C1omit <- filter(Sorption_LOQNAomit, Biochar == "CWC")
ULS_single_C1omit <- filter(Sorption_LOQNAomit, Biochar == "ULS")
BRL_single_C1omit <- filter(Sorption_LOQNAomit, Biochar == "BRL")

nr_compounds <- length(unique(Sorption$Compound))
for(i in 1:nr_compounds){
  fit <- lm(log_Cs ~ log_Cw, data = Sorption)
  print(summary(fit))
  # her kan du velge hva du vil printe! 
  # Du kan ofså lagre mean og sdt til en tabell for å sammenligne senere
}


