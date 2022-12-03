#Inventering####

rm(list=ls())

setwd("~/Documents/Documents/R/Master/Master")

invdat <- read.csv("Inventering data.csv", sep=';')
head(invdat)
str(invdat)

#Data sorting####
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

invdat_sub = subset(invdat, select = -c(Diameter..mm., X.1, X.2, X.3, X.4, Antal,
                                        Bestånd.1))

invdat_sub <- invdat_sub[!apply(is.na(invdat_sub) | invdat_sub == "", 1, all),]

str(invdat_sub)

#Data exploration and model fitting####

invdat$Behandling = factor(invdat$Behandling, levels=c("Kontroll", "Klippt", "Heatweed"))
#Change order so Kontroll is intercept

#Height
m = lm(invdat$Höjd~invdat$Behandling)
anova(m)

summary(m)#48.1% of variance described
hist(invdat$Höjd)

m1 = lm(invdat$Höjd~invdat$Behandling-1)
summary(m1)

#Mean diameter
m2 = lm(invdat$Medel.diameter~invdat$Behandling)
anova(m2)
summary(m2)#53.9% of variance described

hist(invdat$Medel.diameter)

m1 = lm(invdat$Medel.diameter~invdat$Behandling-1)
summary(m1)

#Mean number
library(lme4)
library(MASS)

hist(invdat$Medelantal)
invdat$logantal <- log(invdat$Medelantal)
hist(invdat$logantal)

m3 = glm.nb(Medelantal~Behandling, data=invdat_sub) 
anova(m3)
summary(m3) #35.6% of variance is explained 

m1 = lm(invdat$logantal~invdat$Behandling-1) #P only shows difference from 0
summary(m1) #OBS log scale

coefs = summary(m1)$coef
exp(coefs[1,1])
exp(coefs[2,1])
exp(coefs[3,1])

#Area
hist(invdat$Area)
invdat$logarea <- log(invdat$Area)
hist(invdat$logarea)

invdat$Behandling = factor(invdat$Behandling, levels=c("Kontroll", "Klippt", "Heatweed"))
#Change order so kontroll is intercept
m4 = lm(invdat$logarea~invdat$Behandling) 
anova(m4)
summary(m4) #25.2% of variance explained

m1 = lm(invdat$logarea~invdat$Behandling-1)
summary(m1)

#Plotting####

ggplot(invdat_sub, aes(x=Behandling, y=Medel.diameter)) +
  geom_col(width = 0.75, color="black", fill="grey") +
  #geom_errorbar(data=FDAMeans, aes(ymin= mean.antal - se.FDA, ymax=mean.antal + se.FDA),
                #width = 0.13, alpha = 1, position=position_dodge(0.75)) +
  theme_classic() + 
  #scale_fill_manual(values=c("#73ba2e", "#034c77")) +
  #scale_color_manual(values=c("#000000", "#000000")) +
  scale_y_continuous(limits = c(0,150), expand = c(0,0)) +
  labs(y="Medeldiameter (mm)", x="Behandling", 
       title = "Skillnad i medeldiameter mellan behandlingar") +
  theme(plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 20, angle = 60,
                                   hjust = 1, color = "grey1")) +
  theme(axis.ticks.length=unit(.25, "cm"))


