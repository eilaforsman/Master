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

summary(m)
qqnorm(invdat$Höjd)

m1 = lm(invdat$Höjd~invdat$Behandling-1)
summary(m1)

m = lm(invdat$Medel.diameter~invdat$Behandling)
anova(m)
summary(m)

qqnorm(invdat$Medel.diameter)

m1 = lm(invdat$Medel.diameter~invdat$Behandling-1)
summary(m1)

qqnorm(invdat$Medelantal)
invdat$logantal <- log(invdat$Medelantal)
qqnorm(invdat$logantal)

m = lm(invdat$logantal~invdat$Behandling) #Ändra ordning så kontroll är intercept
anova(m)
summary(m)

m1 = lm(invdat$logantal~invdat$Behandling-1)
summary(m1)

qqnorm(invdat$Area)
invdat$logarea <- log(invdat$Area)
qqnorm(invdat$logarea)

m = lm(invdat$logarea~invdat$Behandling) #Ändra ordning så kontroll är intercept
anova(m)
summary(m)

m1 = lm(invdat$logarea~invdat$Behandling-1)
summary(m1)



