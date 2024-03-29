rm(list = ls())
list.files ()

#Stapeldiagram medelantal per plats, medel på enskilda snitten först 
#Antal HW vs kontroll
#Antal skillnad i lokaler beror på region?
#Nested HW vs kontroll och skillnad mellan lokal
#Medelvärde antal per lokal (Pålsjö 1,2,3,4,5)
#5vs6ggr
#Pålsjö efter HW
#Skillnad i död/levande eller antal mellan strukturer

#Basic####
getwd()

setwd("~/Documents/Documents/R/Master/Master")

dataFDA <- read.csv("dataFDA.csv", sep=";")
head(dataFDA)
str(dataFDA)
dataFDA$Struktur <- as.factor(dataFDA$Struktur)
dataFDA$Lokal <- as.factor(dataFDA$Lokal)
dataFDA$Kommun <- as.factor(dataFDA$Kommun)
table(dataFDA$Lokal)
summary(dataFDA)


library(plyr)
dataFDA_sub <- subset(dataFDA, Lokal != "H_Pålsjö_efter_HW") #Remove "efter_HW" because its a different treatment
str(dataFDA_sub)
head(dataFDA_sub)

#install.packages("stringr")
library(stringr)

#Sort_data####

ID_split <- as.data.frame(str_split_fixed(dataFDA_sub$ID, "_", 2))
head(ID_split)


dataFDA_sub$ID_spec <- ID_split$V1
paste(dataFDA_sub$Lokal, dataFDA_sub$ID_spec, sep ="_")
dataFDA_sub$Prov<-paste(dataFDA_sub$Lokal, dataFDA_sub$ID_spec, sep ="_")
head(dataFDA_sub)

#Created new column of sample without sample ID

Prov<-ddply(dataFDA_sub, "Prov", summarise, mean_antal=mean(Antal))
head(Prov)

nchar(dataFDA_sub$Prov)
dataFDA_sub$newprov <- substr(dataFDA_sub$Prov, 1, nchar(dataFDA_sub$Prov) - 2) #New column of sites

Meantot<-ddply(dataFDA_sub, "newprov", summarise, mean_tot=mean(Antal)) #Calculate mean number per site
Meantot$order <- c(4,5,1,9,10,11,12,13,14,2,8,6,7,3) #Set order of sites for later plotting


#Data exploring####

hist(dataFDA_sub$Antal, xlab="Number of cells", las=1, main=NULL)

dataFDA_sub$Behandling = factor(dataFDA_sub$Behandling, levels=c("Kontroll", "Heatweed"))
m = lm(Antal~Behandling, data=dataFDA_sub)
anova(m)
summary(m)

#Not allowed to use ANOVA to draw conclusions because data is not normally distributed

#Site test#####

#Is there a difference in number of live cells between the Heatweed treated sites?
#Do we need to take site into account in the model?

kruskal.test(dataFDA_sub$Antal ~ dataFDA_sub$Prov) #Big difference between root samples

kruskal.test(dataFDA_sub$Antal ~ dataFDA_sub$newprov) #Big difference between sites if control included

HWsites <- subset(dataFDA_sub, dataFDA_sub$newprov != "Helsingborg_kontroll" & dataFDA_sub$newprov != "Motala_kontroll" & dataFDA_sub$newprov != "Vellinge_kontroll")
# Removes control samples to check if test sites has an effect

kruskal.test(HWsites$Antal ~ HWsites$newprov) #Big difference between HWsites, (pseudo rep)!

#Mixed Model#####

library(glmmTMB)
library(MASS)
library(car)
library(lme4)
library(lmerTest)

?glmmTMB
?lme4
?lmerTest

#LMM assumes normal distribution
#m3 = lmer(response ~ fixed factor + ev fixed 2+ fixed factor 1* fixed 2 + (1|random), data=data_sub)


#GLMM
#m3 = glmer(response ~ fixed factor + ev fixed 2+ fixed factor 1* fixed 2 + (1|random), family= "poisson", data=data_sub)

#Assumes equal variances between groups, (levene.test()), sig is problem 
#If anova doesn't work -> try Anova


m2 = glm.nb(Antal ~ Behandling, data=dataFDA_sub)
summary(m2) #Fitted a generalized linear model with negative binomial distribution
#to account for overdispersion. Big difference but doesn't take site into account

1-(m2$deviance/m2$null.deviance)
#[1] 0.04534703

m3 = glmer(dataFDA_sub$Antal ~ dataFDA_sub$Behandling + (1|dataFDA_sub$Lokal), family="poisson", data=dataFDA_sub)

Anova(m3) #assumes equal variance between groups

leveneTest(dataFDA_sub$Antal ~ dataFDA_sub$Behandling) #was significant -> not equal variances

library(nlme) #for weights


m4 = lme(Antal ~ Behandling, random = ~ 1|Lokal, weights = varIdent(form= ~1|Behandling), data=dataFDA_sub)

anova(m4) #Assumes normal distribution

m5 = glmmTMB(Antal ~ Behandling + (1|Lokal), dispformula = ~ Behandling, family="poisson", data=dataFDA_sub)
summary(m5)

Anova(m5) #Correct data distribution and with non equal variance, site as a random factor

boxplot (sqrt(dataFDA_sub$Antal) ~ dataFDA_sub$Behandling)

resid(m5)
boxplot(resid(m5) ~ dataFDA_sub$Behandling)

m6 = glmer.nb(Antal ~ Behandling + (1|Lokal), family="poisson", data=dataFDA_sub)
m7 = glmmTMB(Antal ~ Behandling + (1|Lokal/Prov), family="nbinom1", data=dataFDA_sub)
summary(m6) 
summary(m7)#Includes site and sample as random factor and corrects for overdispersion, 
#correct data distribution

boxplot(resid(m6) ~ dataFDA_sub$Behandling)
boxplot(resid(m7) ~ dataFDA_sub$Behandling)

plot(Antal ~ Behandling, data=dataFDA_sub) #Visualize data distribution between treatments

m10 = glmmTMB(Antal ~ Behandling, family="nbinom1", data=dataFDA_sub) #Trying wothout random effects to see how they affect the result
summary(m10)

#Results####

coef = summary(m7)
coef = coef[["coefficients"]][["cond"]]
exp(coef[1,1])
exp(coef[1,1] + coef[2,1])

#Create data frame with statistics for each treatment
stats = ddply(dataFDA_sub, "Behandling", summarize, 
              median_treat = median(Antal),
              mean_treat = mean(Antal),
              sd_treat = sd(Antal),
              se_treat = sd_treat/sqrt(length(Antal))) 


#Fancy plotting######

#install.packages("ggplot2")
library(ggplot2)

head(Meantot)
par(mfrow=c(1,1))
?ggplot()

#Mean number of live cells per micrograph####

FDAMeans<-ddply(dataFDA_sub,"Lokal", summarize, N=length(Antal),
                mean.antal=mean(na.omit(Antal)),
                sd.FDA=sd(na.omit(Antal)),
                se.FDA=sd.FDA/sqrt(N)) #Summarize data per site with sd and se for creating error bars

Meantot$newprov = gsub("Helsingborg_kontroll", "Helsingborg_control", Meantot$newprov)
Meantot$newprov = gsub("Motala_kontroll", "Motala_control", Meantot$newprov)
Meantot$newprov = gsub("Vellinge_kontroll", "Vellinge_control", Meantot$newprov)

ggplot(Meantot, aes(x=reorder(Meantot$newprov, Meantot$order), y=Meantot$mean_tot)) +
  geom_bar(width = 0.75, stat = "identity", position ="dodge", alpha = 0.8) +
  geom_errorbar(data=FDAMeans, aes(ymin= mean.antal - se.FDA, ymax=mean.antal + se.FDA),
                width = 0.13, alpha = 1, position=position_dodge(0.75)) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,150), expand = c(0,0)) +
  labs(y="Mean number of live cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 20, angle = 60,
                                   hjust = 1, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("Celler_FDA_plot.png", plot = last_plot(), device = "png",
       scale = 1, width = 12, height = 8,
       dpi = 600)

#Mean number of cells per sample####


SampleMeans <- ddply(dataFDA_sub,"Prov", summarize, N=length(Antal),
                     mean.antal=mean(na.omit(Antal)),
                     sd.FDA=sd(na.omit(Antal)),
                     se.FDA=sd.FDA/sqrt(N)) #Create summarize of number of cells per sample with means, sd and se.
str(SampleMeans)
SampleMeans$Prov = as.character(SampleMeans$Prov)


SampleMeans$newprov <- substr(SampleMeans$Prov, 1, nchar(SampleMeans$Prov) - 2)#New column of sites without ID

#New data frame with mean, sd and se and same length of data
Sample = aggregate(cbind(SampleMeans$mean.antal, SampleMeans$sd.FDA, 
                         SampleMeans$se.FDA), 
                   by=list(newprov=SampleMeans$newprov), FUN =sum) #V1 =mean, V2=sd, V3=se


#Add sample summary data to Meantot (Total means) data frame
Meantot$mean_sample = Sample$V1
Meantot$sd = Sample$V2
Meantot$se = Sample$V3

#Changing names of control samples
Meantot$newprov = gsub("Helsingborg_control","Control Helsingborg", Meantot$newprov)
Meantot$newprov = gsub("Motala_control","Control Motala", Meantot$newprov)
Meantot$newprov = gsub("Vellinge_control","Control Vellinge", Meantot$newprov)

ggplot(Meantot, aes(x=reorder(newprov, order), y=mean_sample)) +
  geom_bar(width = 0.75, stat = "identity", position ="dodge", alpha = 0.8) +
  geom_errorbar(data=Meantot, aes(ymin= mean_sample - se, 
                                  ymax=mean_sample + se),
                width = 0.13, alpha = 1, position=position_dodge(0.75)) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,500), expand = c(0,0)) +
  labs(y="Mean number of cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 20, angle = 60,
                                   hjust = 1, color = "grey1")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("Celler_prov_plot.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Variance####

#Split data according to treatment
vardat = split(dataFDA_sub, f=dataFDA_sub$Behandling) 
datC = vardat$Kontroll
datHW = vardat$Heatweed

hist(datC$Antal) #Visualize Control data

m8 = glm.nb(Antal ~ 1, data=datC) #Fit model to control data
summary(m8)

hist(datHW$Antal)#Visualize HW data

m9 = glm.nb(Antal ~ 1, data=datHW) #Fit model to HW data
summary(m9)


library(performance)
r2(m7)

#Variance explained by sites for control samples
m10 = glmmTMB(Antal ~ 1 + (1|Lokal), family="nbinom1", data=datC)
summary(m10)
r2(m10) #Conditional takes fixed and random effects into account. Marginal only fixed


#Variance explained by sites for HW samples
m11 = glmmTMB(Antal ~ 1 + (1|Lokal), family="nbinom1", data=datHW)
summary(m11)
r2(m11)
#High variance within data not the same as R^2 which shows how much of that variance that the model explains

#Variance explained by samples within sites for control samples
m12 = glmmTMB(Antal ~ 1 + (1|Lokal/Prov), family="nbinom1", data=datC)
summary(m12)
r2(m12)

#Variance explained by samples within sites for HW samples
m13 = glmmTMB(Antal ~ 1 + (1|Lokal/Prov), family="nbinom1", data=datHW)
summary(m13)
r2(m13)

#Plotting treatment difference####

plot(dataFDA_sub$Antal ~ dataFDA_sub$Behandling, las=1)

dataFDA_sub$Behandling = gsub("Kontroll", "Control", dataFDA_sub$Behandling)
dataFDA_sub$Behandling = gsub("Heatweed", "Hot water", dataFDA_sub$Behandling)

ggplot(dataFDA_sub, aes(x=Behandling, y=Antal, fill=Behandling)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(0,500), expand = c(0,0)) +
  labs(y="Number of live cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("Behandling_bild_boxplot.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)

#Same plot in English

ggplot(dataFDA_sub, aes(x=Behandling, y=Antal, fill=Behandling)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(0,500), expand = c(0,0)) +
  labs(y="Number of cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("Treatment_bild_boxplot.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)


#Difference between municipalities####

plot(datHW$Antal ~ datHW$Kommun, las=1) #visualize data

mm = glm.nb(Antal ~ Kommun, data=datHW) #Fit model for number of cells per micrograph
anova(mm)
summary(mm) #Vellinge smaller

mm1 = glm.nb(Antal ~ Kommun-1, data=datHW)
summary(mm1)

coef2 = summary(mm1)$coef
exp(coef2[1,1]) #24.69804 mean Gothenburg
exp(coef2[2,1]) #21.44664 mean Helsingborg
exp(coef2[3,1]) #23.43799 mean Motala
exp(coef2[4,1]) #4.520792 mean Vellinge

datHW$Kommun = factor(datHW$Kommun, levels=c("Vellinge", "Göteborg", "Helsingborg", "Motala"))

mm = glmmTMB(Antal ~ Kommun + (1|Lokal/Prov), family="nbinom1", datHW)
summary(mm) #Fit model where site and sample are random factors
#Vellinge fewest cells, differ from Motala.

mm2 = glmmTMB(Antal ~ Kommun + (1|Lokal/Prov), family="nbinom1", datC)
summary(mm2) #Model difference between municipalities in control site, no difference

boxplot(resid(mm) ~ datHW$Kommun)

coef4 = summary(mm)
coef4 = coef4[["coefficients"]][["cond"]] # Mean cells per micrograph, log scale
exp(coef4[1,1]) # [1] 1.120687 mean V
exp(coef4[1,1] + coef4[2,1]) # [1] 6.475923 mean G
exp(coef4[1,1] + coef4[3,1]) # [1] 5.773852 mean H
exp(coef4[1,1] + coef4[4,1]) # [1] 6.131035 mean M

exp(coef4[1,2]) # [1] 1.901537 SE V
exp(coef4[1,2] + coef4[2,2]) # [1] 5.541874 SE G
exp(coef4[1,2] + coef4[3,2]) # [1] 4.608651 SE H
exp(coef4[1,2] + coef4[4,2]) # [1] 3.949547 SE M


#Sample modelling
Meantot$Kommun = c("Helsingborg", "Helsingborg", "Helsingborg", "Motala",
                   "Motala","Motala","Motala","Motala","Motala","Motala",
                   "Göteborg", "Vellinge","Vellinge","Vellinge")
meantot_sub = subset(Meantot, Meantot$newprov != "Helsingborg_control" & Meantot$newprov != "Motala_control" & Meantot$newprov != "Vellinge_control")
meantot_sub$Kommun = as.factor(meantot_sub$Kommun)
str(meantot_sub)

meantot_sub$Kommun = factor(meantot_sub$Kommun, levels=c("Vellinge", "Göteborg", "Helsingborg", "Motala"))


hist(meantot_sub$mean_sample)
meantot_sub$logsample = log(meantot_sub$mean_sample)
hist(meantot_sub$logsample)

mm2 = lm(logsample ~ Kommun, data=meantot_sub) #Fit model to log data 
summary(mm2) 

mm2 = lm(logsample ~ Kommun-1, data=meantot_sub) #Get estimated means (log scale)
summary(mm2)

coef3 = summary(mm2)$coef
exp(coef3[1,1]) #[1] 4.32938 Mean Vellinge
exp(coef3[2,1]) #[1] 122.446 Mean Göteborg
exp(coef3[3,1]) # [1] 62.73146 Mean Helsingborg
exp(coef3[4,1]) # [1] 90.49646 Mean Motala

exp(coef3[1,2]) # [1] 2.960874 SE V
exp(coef3[2,2]) # [1] 4.641823 SE G
exp(coef3[3,2]) # [1] 2.960874 SE H
exp(coef3[4,2]) # [1] 1.871434 SE M


#Plotting municipalities sample####

plot(meantot_sub$mean_sample ~ meantot_sub$Kommun)

ggplot(meantot_sub, aes(x=Kommun, y=mean_sample, fill=Kommun)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(0,400), expand = c(0,0)) +
  labs(y="Number of cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("Kommun_prov_boxplot.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)

#Plotting municipalities micrograph

meantot_sub$Kommun = gsub("Göteborg","Gothenburg", meantot_sub$Kommun)
meantot_sub$Kommun = factor(meantot_sub$Kommun, levels=c("Vellinge", "Gothenburg", "Helsingborg", "Motala"))


ggplot(meantot_sub, aes(x=Kommun, y=mean_tot, fill=Kommun)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(-5,60), expand = c(0,0)) +
  labs(y="Number of cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("Kommun_bild_boxplot.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)


#Structure####

plot(Antal ~ Struktur, data=datHW)

HW_str <- subset(datHW, datHW$Struktur != "pith, cortex")

HW_str$Struktur <- as.factor(HW_str$Struktur)

str(HW_str)
summary(HW_str)

hist(HW_str$Antal)

mod1 = glm.nb(Antal ~ Struktur, data=HW_str)
anova(mod1)
summary(mod1) # sig dif betweeen structures, vt fewest cells, pith most, not great model fit

coef5 = summary(mod1)$coef

exp(coef5[1,1]) # [1] 28.75203 mean cortex
exp(coef5[2,1] + coef5[1,1]) # [1] 43.55789 mean pith
exp(coef5[3,1] + coef5[1,1]) # [1] 0.2593235 mean vt

exp(coef5[1,2]) # [1] 1.086046 SE cortex
exp(coef5[2,2] + coef5[1,2]) # [1] 1.287328 SE pith
exp(coef5[3,2] + coef5[1,2]) # [1] 1.239278 SE vt

mod2 = glmmTMB(Antal ~ Struktur + (1|Lokal/Prov), family="nbinom1", data=HW_str)
summary(mod2) #sig dif between str, takes site and sample into account

coef6 = summary(mod2)$coef
coef6 = coef6[["cond"]]

exp(coef6[1,1]) # [1] 5.541543 mean cortex
exp(coef6[1,1] + coef6[2,1]) # [1] 8.660902 mean pith
exp(coef6[1,1] + coef6[3,1]) # [1] 0.6752093 mean vt

exp(coef6[1,2]) #[1] 1.368077 SE cortex
exp(coef6[1,2] + coef6[2,2]) # [1] 1.462484 SE pith
exp(coef6[1,2] + coef6[3,2]) #[1] 1.55402 SE vt

#Means don't make sense should be 20 to 80 something per micrograph

kruskal.test(Antal ~ Struktur, data=HW_str)
#Kruskal-Wallis chi-squared = 638.68, df = 2, p-value < 2.2e-16
#sig dif between str

str_means = ddply(HW_str, "Struktur", summarize, 
                          mean = mean(Antal),
                          sd = sd(Antal),
                          se = sd /sqrt(length(Antal)))

#Plotting str####

plot(Antal ~ Struktur, data=HW_str)

ggplot(HW_str, aes(x=Struktur, y=Antal, fill=Struktur)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(0,500), expand = c(0,0)) +
  labs(y="Number of cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("STR_boxplot.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)

#Difference in control vs HW and structure#####

#data exploring 
data = dataFDA_sub 

data = subset(dataFDA_sub, dataFDA_sub$Struktur!="pith, cortex")

data$Struktur = as.factor(data$Struktur)
data$Behandling = as.factor(data$Behandling)
str(data)

plot(data$Struktur)
hist(data$Antal)

means = tapply(data$Antal, list(data$Struktur, data$Behandling), mean)

means = na.omit(means)

se = tapply(data$Antal, 
            list(data$Struktur, data$Behandling), 
            function(x) sd(x)/sqrt(sum(!is.na(x))))
se = na.omit(se)

colMeans(means)
#Control Heatweed 
#56.03020 24.18975 

rowMeans(means)
#  cortex            pith vascular tissue 
#48.7601593      71.3133389       0.2564223

colMeans(se)
# Control Heatweed 
#2.515320 1.730181 

rowMeans(se)
#   cortex            pith vascular tissue 
#2.17440696      4.10192672      0.09191741 

#mann whitney U compare seperatly

d = subset(data, data$Struktur =="cortex")
mu = wilcox.test(d$Antal ~ d$Behandling)
mu # W = 481916, p-value < 2.2e-16

d = subset(data, data$Struktur =="pith")
mu = wilcox.test(d$Antal ~ d$Behandling)
mu # W = 56851, p-value < 2.2e-16

d = subset(data, data$Struktur =="vascular tissue")
mu = wilcox.test(d$Antal ~ d$Behandling)
mu # W = 41317, p-value = 0.7502

mu_treat = wilcox.test(data$Antal ~ data$Behandling)
mu_treat # W = 1719611, p-value < 2.2e-16
summary(mu_treat)

#Plotting

library(dplyr)
#data %>% 
  #group_by(Struktur, Behandling) %>% 
  #summarise(data_groups = mean(Antal)) -> data2

data %>% 
  ggplot() +
  aes(x = Struktur, color = Behandling, group = Behandling, y = Antal) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  theme_classic() +
  scale_color_grey() +
  labs(y="Mean number of cells", x="", 
       title = "", color = "Treatment") +
  theme(text = element_text(size = 28, family="Times"),
        legend.position = c(0.9,0.9),
        legend.title = element_text("Times"),
        axis.text.x = element_text(size = 28, color = "black"))

ggsave("str_treat_interaction.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)




