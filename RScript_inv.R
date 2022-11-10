#Inventering####

rm(list=ls())

setwd("~/Documents/Documents/R/Master/Master")

invdat <- read.csv("Inventering data.csv", sep=';')
head(invdat)
str(invdat)

invdat$Bestånd <- as.factor(invdat$Bestånd)
invdat$Behandling <- as.factor(invdat$X)
invdat

invdat$Höjd <- gsub(",",".", invdat$Höjd)
invdat$Höjd <- as.numeric(invdat$Höjd)

invdat$Medel.diameter <- gsub(",",".", invdat$Medel.diameter)
invdat$Medel.diameter <- as.numeric(invdat$Medel.diameter)

invdat$Area <- gsub(",",".", invdat$Area)
invdat$Area <- as.numeric(invdat$Area)

invdat$Medelantal <- gsub(",",".", invdat$Medelantal)
invdat$Medelantal <- as.numeric(invdat$Medelantal)


m = lm(invdat$Höjd~invdat$Behandling)
anova(m)
