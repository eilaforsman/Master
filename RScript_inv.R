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

invdat_sub$Behandling = factor(invdat_sub$Behandling, levels=c("Kontroll", "Klippt", "Heatweed"))
#Change order so Kontroll is intercept

#Height
m = lm(invdat_sub$Höjd~invdat_sub$Behandling)
anova(m)

summary(m)#48.1% of variance described
hist(invdat_sub$Höjd)

m1 = lm(invdat_sub$Höjd~invdat_sub$Behandling-1)
summary(m1)

#Mean diameter
m2 = lm(invdat_sub$Medel.diameter~invdat_sub$Behandling)
anova(m2)
summary(m2)#53.9% of variance described

hist(invdat_sub$Medel.diameter)

m1 = lm(invdat_sub$Medel.diameter~invdat_sub$Behandling-1)
summary(m1)

#Mean number
library(lme4)
library(MASS)

hist(invdat_sub$Medelantal)
invdat_sub$logantal <- log(invdat_sub$Medelantal)
hist(invdat_sub$logantal)

m3 = lm(logantal~Behandling, data=invdat_sub) 
anova(m3)
summary(m3) #OBS log scale, #35.6% of variance is explained 

m1 = lm(invdat_sub$logantal~invdat_sub$Behandling-1) 
summary(m1) 

coefs = summary(m1)$coef
exp(coefs[1,1]) # [1] 17.22322 Control
exp(coefs[2,1]) # [1] 4.543839 Mowed
exp(coefs[3,1]) # [1] 3.05623 HW

#Area
hist(invdat_sub$Area)
invdat_sub$logarea <- log(invdat_sub$Area)
hist(invdat_sub$logarea)

invdat_sub$Behandling = factor(invdat_sub$Behandling, levels=c("Kontroll", "Klippt", "Heatweed"))
#Change order so kontroll is intercept
m4 = lm(invdat_sub$logarea~invdat_sub$Behandling) 
anova(m4)
summary(m4) #25.2% of variance explained

m1 = lm(invdat_sub$logarea~invdat_sub$Behandling-1)
summary(m1)

#Plotting####
library(ggplot2)
library(plyr)

#Mean diameter
stats <- ddply(invdat_sub,"Behandling", summarize, N=length(Medel.diameter),
                     mean.diameter=mean(na.omit(Medel.diameter)),
                     sd=sd(na.omit(Medel.diameter)),
                     se=sd/sqrt(N)) 

ggplot(data = stats, aes(x=Behandling, y=mean.diameter, 
                       fill=Behandling)) +
  geom_col() +
  scale_fill_grey() +
  geom_errorbar(data=stats, 
                aes(ymin = mean.diameter - se, 
                    ymax = mean.diameter + se),
                width = 0.13, alpha = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
  labs(y="Medeldiameter (mm)", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black"))+
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("medeldiameter.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Mean height

stats <- ddply(invdat_sub,"Behandling", summarize, N=length(Höjd),
               mean.height=mean(na.omit(Höjd)),
               sd=sd(na.omit(Höjd)),
               se=sd/sqrt(N)) 


ggplot(data = stats, aes(x=Behandling, y=mean.height, 
                         fill=Behandling)) +
  geom_col() +
  scale_fill_grey() +
  geom_errorbar(data=stats, 
                aes(ymin = mean.height - se, 
                    ymax = mean.height + se),
                width = 0.13, alpha = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,250), expand = c(0,0)) +
  labs(y="Medelhöjd (cm)", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("medelhöjd.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Mean number of shoots


stats <- ddply(invdat_sub,"Behandling", summarize, N=length(logantal),
               mean.number=mean(na.omit(logantal)),
               sd=sd(na.omit(logantal)),
               se=sd/sqrt(N)) 

ggplot(data = stats, aes(x=Behandling, y=mean.number, 
                         fill=Behandling)) +
  geom_col() +
  scale_fill_grey() +
  geom_errorbar(data=stats, 
                aes(ymin = mean.number - se, 
                    ymax = mean.number + se),
                width = 0.13, alpha = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
  labs(y="Täthet av skott", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("logantal.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Mean area

stats <- ddply(invdat_sub,"Behandling", summarize, N=length(logarea),
               mean.area=mean(na.omit(logarea)),
               sd=sd(na.omit(logarea)),
               se=sd/sqrt(N)) 

ggplot(data = stats, aes(x=Behandling, y=mean.area, 
                         fill=Behandling)) +
  geom_col() +
  scale_fill_grey() +
  geom_errorbar(data=stats, 
                aes(ymin = mean.area - se, 
                    ymax = mean.area + se),
                width = 0.13, alpha = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  labs(y="Log Medelarea (m^2)", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("medelarea.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#English plots####

invdat_sub$treatment = c("Heatweed","Heatweed","Heatweed","Heatweed",
                         "Control", "Mowed", "Mowed", "Heatweed",
                         "Control", "Mowed", "Control", "Mowed",
                         "Control","Heatweed","Heatweed","Control")

invdat_sub$treatment = factor(invdat_sub$treatment, levels=c("Control", "Mowed", "Heatweed"))

#Mean diameter
stats <- ddply(invdat_sub,"treatment", summarize, N=length(Medel.diameter),
               mean.diameter=mean(na.omit(Medel.diameter)),
               sd=sd(na.omit(Medel.diameter)),
               se=sd/sqrt(N)) 

ggplot(data = stats, aes(x=treatment, y=mean.diameter, 
                         fill=treatment)) +
  geom_col() +
  scale_fill_grey() +
  geom_errorbar(data=stats, 
                aes(ymin = mean.diameter - se, 
                    ymax = mean.diameter + se),
                width = 0.13, alpha = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
  labs(y="Mean diameter (mm)", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black"))+
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("meandiameter.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Mean height

stats <- ddply(invdat_sub,"treatment", summarize, N=length(Höjd),
               mean.height=mean(na.omit(Höjd)),
               sd=sd(na.omit(Höjd)),
               se=sd/sqrt(N)) 


ggplot(data = stats, aes(x=treatment, y=mean.height, 
                         fill=treatment)) +
  geom_col() +
  scale_fill_grey() +
  geom_errorbar(data=stats, 
                aes(ymin = mean.height - se, 
                    ymax = mean.height + se),
                width = 0.13, alpha = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,250), expand = c(0,0)) +
  labs(y="Mean height (cm)", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("meanheight.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Mean number of shoots log scale


stats <- ddply(invdat_sub,"treatment", summarize, N=length(logantal),
               mean.number=mean(na.omit(logantal)),
               sd=sd(na.omit(logantal)),
               se=sd/sqrt(N)) 

ggplot(data = stats, aes(x=treatment, y=mean.number, 
                         fill=treatment)) +
  geom_col() +
  scale_fill_grey() +
  geom_errorbar(data=stats, 
                aes(ymin = mean.number - se, 
                    ymax = mean.number + se),
                width = 0.13, alpha = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
  labs(y="Shoot density (log)", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("logdensity.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Mean number of shoots normal scale


stats <- ddply(invdat_sub,"treatment", summarize, N=length(Medelantal),
               mean.number=mean(na.omit(Medelantal)),
               sd=sd(na.omit(Medelantal)),
               se=sd/sqrt(N)) 

ggplot(data = stats, aes(x=treatment, y=mean.number, 
                         fill=treatment)) +
  geom_col() +
  scale_fill_grey() +
  geom_errorbar(data=stats, 
                aes(ymin = mean.number - se, 
                    ymax = mean.number + se),
                width = 0.13, alpha = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,40), expand = c(0,0)) +
  labs(y="Shoot density", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("shootdensity.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Mean area

stats <- ddply(invdat_sub,"treatment", summarize, N=length(logarea),
               mean.area=mean(na.omit(logarea)),
               sd=sd(na.omit(logarea)),
               se=sd/sqrt(N)) 

ggplot(data = stats, aes(x=treatment, y=mean.area, 
                         fill=treatment)) +
  geom_col() +
  scale_fill_grey() +
  geom_errorbar(data=stats, 
                aes(ymin = mean.area - se, 
                    ymax = mean.area + se),
                width = 0.13, alpha = 1) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  labs(y="Log mean area (m^2)", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("meanarea.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)
