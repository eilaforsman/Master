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
hist(invdat$Höjd)

m1 = lm(invdat$Höjd~invdat$Behandling-1)
summary(m1)

m = lm(invdat$Medel.diameter~invdat$Behandling)
anova(m)
summary(m)

hist(invdat$Medel.diameter)

m1 = lm(invdat$Medel.diameter~invdat$Behandling-1)
summary(m1)

hist(invdat$Medelantal)
invdat$logantal <- log(invdat$Medelantal)
hist(invdat$logantal)


invdat$Behandling = factor(invdat$Behandling, levels=c("Kontroll", "Klippt", "Heatweed"))
#Change order so Kontroll is intercept
m = lm(invdat$logantal~invdat$Behandling) 
anova(m)
summary(m)
SS_T = 9.0346+16.3737 # 25.4083
9.0346/SS_T #0.3556 Meaning 35.56% of variance is explained 


m1 = lm(invdat$logantal~invdat$Behandling-1) #P only shows difference from 0
summary(m1)

hist(invdat$Area)
invdat$logarea <- log(invdat$Area)
hist(invdat$logarea)

invdat$Behandling = factor(invdat$Behandling, levels=c("Kontroll", "Klippt", "Heatweed"))
#Change order so kontroll is intercept
m = lm(invdat$logarea~invdat$Behandling) 
anova(m)
summary(m) #25.2% of variance explained

m1 = lm(invdat$logarea~invdat$Behandling-1)
summary(m1)

ggplot(invdat,(x=X), y=Medel.diameter) +
  #geom_bar(width = 0.75, stat = "identity", position ="dodge", alpha = 0.8) +
  #geom_errorbar(data=FDAMeans, aes(ymin= mean.antal - se.FDA, ymax=mean.antal + se.FDA),
                #width = 0.13, alpha = 1, position=position_dodge(0.75)) +
  theme_classic() + 
  #scale_fill_manual(values=c("#73ba2e", "#034c77")) +
  #scale_color_manual(values=c("#000000", "#000000")) +
  scale_y_continuous(limits = c(0,150), expand = c(0,0)) +
  labs(y="Medelantal per prov", x="", title = "Levande celler per tvärsnitt från varje lokal") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 20, angle = 60,
                                   hjust = 1, color = "grey1")) +
  theme(axis.ticks.length=unit(.25, "cm"))


