#Andel_VH####

rm(list=ls())

getwd()
VH <- read.csv("Andel_VH.csv", sep=';')
str(VH)

head(VH)

#Data sorting####

VH$Procent.levande <- gsub(",",".", VH$Procent.levande)
VH$Procent.levande <- as.numeric(VH$Procent.levande)

VH$Procent.död = gsub(",",".", VH$Procent.död)
VH$Procent.död<-as.numeric(VH$Procent.död)

VH$Plats = gsub(":", "_", VH$Plats)
VH$Plats<- as.factor(VH$Plats)
str(VH)

VH = subset(VH, select = c(Plats, Grupp, levande, död, Procent.levande, 
                           Procent.död, Totala.antal.behandlingar))
VH <- VH[!apply(is.na(VH) | VH == "", 1, all),]

VH_sub = subset(VH, VH$Plats != "H_ Pålsjö efter HW")

VH_sub$Plats = gsub("Helsingborg kontroll","Helsingborg control", VH_sub$Plats)
VH_sub$Plats = gsub("Flommen kontroll","Vellinge control", VH_sub$Plats)
VH_sub$Plats = gsub("Motala kontroll","Motala control", VH_sub$Plats)

#Plotting live samples####
library(ggplot2)

ggplot(VH_sub, aes(x=reorder(Plats, Grupp), y=Procent.levande)) +
  geom_bar(width = 0.75, stat = "identity", position ="dodge", alpha = 0.8) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  labs(y="Proportion of shoots emerged (%)", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 20, angle = 60,
                                   hjust = 1, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("VH_live.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Difference in sprouting between sites####

VH_sub_HW = subset(VH_sub, VH_sub$Plats != "Helsingborg control" &
                VH_sub$Plats != "Vellinge control" &
                VH_sub$Plats != "Motala control")

#Mean then create deviance vector and do one sample t test

deviance = VH_sub_HW$levande - mean(VH_sub_HW$levande)

t.test(deviance, mu = 0, alternative = "two.sided")
#t = 3.2749e-16, df = 10, p-value = 1

#Difference between treatment####
library(glmmTMB)
library(lme4)
library(lmerTest)
library(MASS)

VH_sub$Treatment = c("Control","Control","Control","Heatweed","Heatweed",
                     "Heatweed","Heatweed","Heatweed","Heatweed","Heatweed",
                     "Heatweed","Heatweed","Heatweed","Heatweed")

hist(VH_sub$levande)
qqnorm(VH_sub$levande)

m = glm.nb(levande ~ Treatment, data=VH_sub)
summary(m)

coef = summary(m)$coef
exp(coef[1,1]) #8 mean Control
exp(coef[2,1]) #0.3 mean HW

exp(coef[1,2]) #1.3 SE Control
exp(coef[2,2]) #1.4 SE HW

#Plotting treatment####

ggplot(VH_sub, aes(x=Treatment, y=levande, fill=Treatment)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(0,8), expand = c(0,0)) +
  labs(y="Number of shoots", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("VH_treat.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)

#Municipalites####

#Kruskal wallis!!!!

VH_sub_HW$Kommun = c("Helsingborg","Helsingborg","Vellinge","Vellinge","Gothenburg",
                     "Motala","Motala","Motala","Motala","Motala","Motala")

VH_sub_HW$Kommun = factor(VH_sub_HW$Kommun, levels=c("Vellinge", "Gothenburg","Motala","Helsingborg"))

hist(VH_sub_HW$levande)

m2 = kruskal.test(levande ~ Kommun, data=VH_sub_HW)
m2 #Kruskal-Wallis chi-squared = 5.8605, df = 3, p-value = 0.1186

m3 = glm(levande ~ Kommun, data=VH_sub_HW, family="nbinom1")
summary(m3)


#exp(coef[2,1]) # [1] 1.35 mean Helsingborg
#exp(coef[3,1]) # [1] 1.241981e-09 mean Vellinge
#exp(coef[4,1]) # [1] 0.3 mean Gothenburg

#exp(coef[1,2]) # [1] 1.269107 SE Motala
#exp(coef[2,2]) # [1] 1.543291 SE Helsingborg
#exp(coef[3,2]) # [1] Inf SE Vellinge
#exp(coef[4,2]) # [1] 2.85092 SE Gothenburg

#Plotting municipalities

ggplot(VH_sub_HW, aes(x=Kommun, y=levande, fill=Kommun)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(0,8), expand = c(0,0)) +
  labs(y="Number of shoots", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("VH_Kommun.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)

#Fixing greenhouse live samples plot

VH_sub$Plats = gsub("Helsingborg control","Control Helsingborg", VH_sub$Plats)
VH_sub$Plats = gsub("Vellinge control","Control Vellinge", VH_sub$Plats)
VH_sub$Plats = gsub("Motala control","Control Motala", VH_sub$Plats)

ggplot(VH_sub, aes(x=reorder(Plats, Grupp), y=Procent.levande)) +
  geom_bar(width = 0.75, stat = "identity", position ="dodge", alpha = 0.8) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  labs(y="Proportion of shoots emerged (%)", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 20, angle = 60,
                                   hjust = 1, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("VH_live.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)
