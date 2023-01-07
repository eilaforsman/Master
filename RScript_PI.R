rm(list = ls())


#Stapeldiagram medelantal per plats, medel på enskilda snitten först 
#Antal HW vs kontroll
#Antal skillnad i lokaler beror på region?
#Nested HW vs kontroll och skillnad mellan lokal
#Medelvärde antal per lokal (Pålsjö 1,2,3,4,5)
#5vs6ggr
#Skillnad i död/levande eller antal mellan strukturer

#Basic data exploration####
getwd()

setwd("~/Documents/Documents/R/Master/Master")
list.files ()

dataPI <- read.csv("dataPI.csv", sep=";")
head(dataPI)
str(dataPI)

#Specify factors
dataPI$Struktur <- as.factor(dataPI$Struktur)
dataPI$Lokal <- as.factor(dataPI$Lokal)
dataPI$Kommun <- as.factor(dataPI$Kommun)
dataPI$Behandling <- as.factor(dataPI$Behandling)

table(dataPI$Lokal)
summary(dataPI)


library(plyr)
dataPI_sub <- subset(dataPI, Lokal != "H_Pålsjö_efter_HW") #Remove "efter_HW" because its a different treatment
str(dataPI_sub)
head(dataPI_sub)

#install.packages("stringr")
library(stringr)

#Sort_data####

ID_split <- as.data.frame(str_split_fixed(dataPI_sub$ID, "_", 2))
head(ID_split)

#Create coloumn for samples
dataPI_sub$ID_spec <- ID_split$V1
paste(dataPI_sub$Lokal, dataPI_sub$ID_spec, sep ="_")
dataPI_sub$Prov<-paste(dataPI_sub$Lokal, dataPI_sub$ID_spec, sep ="_")
head(dataPI_sub)

#Created new df of mean cells per sample 

Prov<-ddply(dataPI_sub, "Prov", summarise, mean_antal=mean(Antal)) #Mean cells per sample
head(Prov)

#Calculate mean number per site

nchar(dataPI_sub$Prov)
dataPI_sub$newprov <- substr(dataPI_sub$Prov, 1, nchar(dataPI_sub$Prov) - 2) #New column of sites

Meantot<-ddply(dataPI_sub, "newprov", summarise, mean_tot=mean(Antal)) 
Meantot$order <- c(4,5,1,8,10,11,12,13,14,9,2,6,7,3) #Set order of sites for later plotting


#Data exploring####

hist(dataPI_sub$Antal, xlab="Number of cells", las=1, main=NULL)

dataPI_sub$Behandling = factor(dataPI_sub$Behandling, levels=c("Control", "Heatweed"))
m = lm(Antal~Behandling, data=dataPI_sub)
anova(m)
summary(m)

#Not allowed to use ANOVA to draw conclusions because data is not normally distributed

#Site test#####

#Is there a difference in number of live cells between the Heatweed treated sites?
#Do we need to take site into account in the model?

kruskal.test(dataPI_sub$Antal ~ dataPI_sub$Prov) #Big difference between root samples

kruskal.test(dataPI_sub$Antal ~ dataPI_sub$newprov) #Big difference between sites if control included

HWsites <- subset(dataPI_sub, dataPI_sub$newprov != "Helsingborg_kontroll" & dataPI_sub$newprov != "Motala_kontroll" & dataPI_sub$newprov != "Vellinge_kontroll")
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


m2 = glm.nb(Antal ~ Behandling, data=dataPI_sub)
summary(m2) #Fitted a generalized linear model with negative binomial distribution
#to account for overdispersion. Big difference but doesn't take site into account

1-(m2$deviance/m2$null.deviance)
#[1] 0.02115497

m3 = glmer(dataPI_sub$Antal ~ dataPI_sub$Behandling + (1|dataPI_sub$Lokal), family="poisson", data=dataPI_sub)

Anova(m3) #assumes equal variance between groups

leveneTest(dataPI_sub$Antal ~ dataPI_sub$Behandling) #was significant -> not equal variances

library(nlme) #for weights


m4 = lme(Antal ~ Behandling, random = ~ 1|Lokal, weights = varIdent(form= ~1|Behandling), data=dataPI_sub)

anova(m4) #Assumes normal distribution

m5 = glmmTMB(Antal ~ Behandling + (1|Lokal), dispformula = ~ Behandling, family="poisson", data=dataPI_sub)
summary(m5)

Anova(m5) #Correct data distribution and with non equal variance, site as a random factor

boxplot (sqrt(dataPI_sub$Antal) ~ dataPI_sub$Behandling)

resid(m5)
boxplot(resid(m5) ~ dataPI_sub$Behandling)

m6 = glmer.nb(Antal ~ Behandling + (1|Lokal), family="poisson", data=dataPI_sub)
m7 = glmmTMB(Antal ~ Behandling + (1|Lokal/Prov), family="nbinom1", data=dataPI_sub)
summary(m6) 
summary(m7)#Includes site and sample as random factor and corrects for overdispersion, 
#correct data distribution

#No difference between control and Heatweed

boxplot(resid(m7) ~ dataPI_sub$Behandling)

plot(Antal ~ Behandling, data=dataPI_sub) #Visualize data distribution between treatments

#Results####

coef = summary(m7)
coef = coef[["coefficients"]][["cond"]]
exp(coef[1,1]) #[1] 0.01823506
exp(coef[1,1] + coef[2,1]) #[1] 0.0962816

#Create data frame with statistics for each treatment
stats = ddply(dataPI_sub, "Behandling", summarize, 
              median_treat = median(Antal),
              mean_treat = mean(Antal),
              sd_treat = sd(Antal),
              se_treat = sd_treat/sqrt(length(Antal))) 


#Fancy plotting######

library(ggplot2)

head(Meantot)
par(mfrow=c(1,1))

#Mean number of live cells per micrograph####

PIMeans<-ddply(dataPI_sub,"Lokal", summarize, N=length(Antal),
                mean.antal=mean(na.omit(Antal)),
                sd.PI=sd(na.omit(Antal)),
                se.PI=sd.PI/sqrt(N)) #Summarize data per site with sd and se for creating error bars

#Rename control sites to english names
Meantot$newprov = gsub("Helsingborg_kontroll", "Helsingborg_control", Meantot$newprov)
Meantot$newprov = gsub("Motala_kontroll", "Motala_control", Meantot$newprov)
Meantot$newprov = gsub("Vellinge_kontroll", "Vellinge_control", Meantot$newprov)

ggplot(Meantot, aes(x=reorder(Meantot$newprov, Meantot$order), y=Meantot$mean_tot)) +
  geom_bar(width = 0.75, stat = "identity", position ="dodge", alpha = 0.8) +
  geom_errorbar(data=PIMeans, aes(ymin= mean.antal - se.PI, ymax=mean.antal + se.PI),
                width = 0.13, alpha = 1, position=position_dodge(0.75)) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
  labs(y="Mean number of dead cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 20, angle = 60,
                                   hjust = 1, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("Mean_cells_PI_plot.png", plot = last_plot(), device = "png",
       scale = 1, width = 12, height = 8,
       dpi = 600)

#Mean number of cells per sample####


SampleMeans <- ddply(dataPI_sub,"Prov", summarize, N=length(Antal),
                     mean.antal=mean(na.omit(Antal)),
                     sd.PI=sd(na.omit(Antal)),
                     se.PI=sd.PI/sqrt(N)) #Create summarize of number of cells per sample with means, sd and se.
str(SampleMeans)
SampleMeans$Prov = as.character(SampleMeans$Prov)

SampleMeans$newprov <- substr(SampleMeans$Prov, 1, nchar(SampleMeans$Prov) - 2)#New column of sites without ID

#New data frame with mean, sd and se and same length of data
Sample = aggregate(cbind(SampleMeans$mean.antal, SampleMeans$sd.PI, 
                         SampleMeans$se.PI), 
                   by=list(newprov=SampleMeans$newprov), FUN =sum) #V1 =mean, V2=sd, V3=se


#Add sample summary data to Meantot (Total means) data frame
Meantot$mean_sample = Sample$V1
Meantot$sd = Sample$V2
Meantot$se = Sample$V3

ggplot(Meantot, aes(x=reorder(newprov, order), y=mean_sample)) +
  geom_bar(width = 0.75, stat = "identity", position ="dodge", alpha = 0.8) +
  geom_errorbar(data=Meantot, aes(ymin= mean_sample - se, 
                                  ymax=mean_sample + se),
                width = 0.13, alpha = 1, position=position_dodge(0.75)) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
  labs(y="Mean number of cells per sample", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 20, angle = 60,
                                   hjust = 1, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("CellsPI_sample_plot.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)

#Variance####

#Split data according to treatment
vardat = split(dataPI_sub, f=dataPI_sub$Behandling) 
datC = vardat$Control
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
#My estimates for variance are unreliable because the data is too close to 0

#Plotting treatment difference####

plot(dataPI_sub$Antal ~ dataPI_sub$Behandling, las=1)

ggplot(dataPI_sub, aes(x=Behandling, y=Antal, fill=Behandling)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(0,200), expand = c(0,0)) +
  labs(y="Number of dead cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("Treatment_PI_boxplot.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)

#Difference between municipalities####

#Set order of municipalities
datHW$Kommun = factor(datHW$Kommun, levels=c("Vellinge", "Göteborg", "Helsingborg", "Motala"))

plot(datHW$Antal ~ datHW$Kommun, las=1) #visualize data

mm = glm.nb(Antal ~ Kommun, data=datHW) #Fit model for number of cells per micrograph
anova(mm)
summary(mm) #Sig dif between GB, HB and Vellinge but not between Vellinge and Motala

mm1 = glm.nb(Antal ~ Kommun-1, data=datHW)
summary(mm1)

coef2 = summary(mm)$coef
exp(coef2[1,1]) #[1] 0.2840237 mean Vellinge
exp(coef2[2,1]) #[1] 7.055474 mean Gothenburg
exp(coef2[3,1]) #[1] 5.82589 mean Helsingborg
exp(coef2[4,1]) #[1] 0.9714246 mean Motala

#Take site and sample into account

mm = glmmTMB(Antal ~ Kommun + (1|Lokal/Prov), family="nbinom1", datHW)
summary(mm) #Fit model where site and sample are random factors
#No difference.

plot(Antal ~ Kommun, data=datC)
mm2 = glmmTMB(Antal ~ Kommun + (1|Lokal/Prov), family="nbinom1", datC)
summary(mm2) #Model difference between municipalities in control site
#Doesn't work, only 1 group with values other than 0

boxplot(resid(mm) ~ datHW$Kommun)

coef3 = summary(mm)
coef3 = coef3[["coefficients"]][["cond"]] # Mean cells per micrograph, log scale
exp(coef3[1,1]) # [1] 0.09597268 mean V
exp(coef3[1,1] + coef3[2,1]) # [1] 0.7730085 mean G
exp(coef3[1,1] + coef3[3,1]) # [1] 0.40313 mean H
exp(coef3[1,1] + coef3[4,1]) # [1] 0.05495027 mean M

exp(coef3[1,2]) # [1] 2.088992 SE V
exp(coef3[1,2] + coef3[2,2]) # [1] 6.328839 SE G
exp(coef3[1,2] + coef3[3,2]) # [1] 5.350443 SE H
exp(coef3[1,2] + coef3[4,2]) # [1] 4.777567 SE M

#Add municipalities to mean df
Meantot$Kommun = c("Helsingborg", "Helsingborg", "Helsingborg", "Göteborg",
                   "Motala","Motala","Motala","Motala","Motala","Motala",
                   "Motala", "Vellinge","Vellinge","Vellinge")

#Remove controls
meantot_sub = subset(Meantot, Meantot$newprov != "Helsingborg_control" & Meantot$newprov != "Motala_control" & Meantot$newprov != "Vellinge_control")
meantot_sub$Kommun = as.factor(meantot_sub$Kommun) #Set municipality as factor
str(meantot_sub)

meantot_sub$Kommun = factor(meantot_sub$Kommun, levels=c("Vellinge", "Göteborg", "Helsingborg", "Motala")) #Set plotting order

#Unnecessary? sample modelling
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
  scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
  labs(y="Number of cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("PI_Kommun_prov_boxplot.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)

#Plotting municipalities micrograph

#Change name of municipalities to english
meantot_sub$Kommun = gsub("Göteborg","Gothenburg", meantot_sub$Kommun)
meantot_sub$Kommun = factor(meantot_sub$Kommun, levels=c("Vellinge", "Gothenburg", "Helsingborg", "Motala"))


ggplot(meantot_sub, aes(x=Kommun, y=mean_tot, fill=Kommun)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  labs(y="Number of dead cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("Kommun_PI_boxplot.png", plot = last_plot(), device = "png",
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
summary(mod1) # sig dif betweeen structures, vt fewest cells, pith most


coef5 = summary(mod1)$coef

exp(coef5[1,1]) # [1] 0.8013621 mean cortex
exp(coef5[2,1] + coef5[1,1]) # [1] 2.276029 mean pith
exp(coef5[3,1] + coef5[1,1]) # [1] 0.1378378 mean vt

exp(coef5[1,2]) # [1] 1.265811 SE cortex
exp(coef5[2,2] + coef5[1,2]) # [1] 1.916672 SE pith
exp(coef5[3,2] + coef5[1,2]) # [1] 1.716296 SE vt

#take site and sample into account
mod2 = glmmTMB(Antal ~ Struktur + (1|Lokal/Prov), family="nbinom1", data=HW_str)
summary(mod2) #sig dif between str, 

coef6 = summary(mod2)$coef
coef6 = coef6[["cond"]]

exp(coef6[1,1]) # [1] 0.1436847 mean cortex
exp(coef6[1,1] + coef6[2,1]) # [1] 0.2765676
exp(coef6[1,1] + coef6[3,1]) # [1] 0.02647525

exp(coef6[1,2]) #[1] 1.723259 SE cortex
exp(coef6[1,2] + coef6[2,2]) # [1] 2.052546 SE pith
exp(coef6[1,2] + coef6[3,2]) #[1] 2.199552 SE vt

#Plotting str####

plot(Antal ~ Struktur, data=HW_str)

ggplot(HW_str, aes(x=Struktur, y=Antal, fill=Struktur)) +
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() + 
  scale_y_continuous(limits = c(0,200), expand = c(0,0)) +
  labs(y="Number of cells", x="", 
       title = "") +
  theme(legend.position = c(0.9,0.9), 
        legend.title = element_blank(),
        plot.title = element_text (hjust = 0.5),
        text = element_text(size=28, family= "Times"),
        axis.text.x = element_text(size = 28, angle = 0,
                                   hjust = 0.5, color = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))

ggsave("STR__PI_boxplot.png", plot = last_plot(), device = "png",
       scale = 1, width = 13, height = 8,
       dpi = 600)

#Difference in control vs HW and structure#####

#data exploring 
data = dataPI_sub 

data = subset(dataPI_sub, dataPI_sub$Struktur!="pith, cortex")
#Why is "pith, cortex" still included in graphs and tapply?

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
#Control  Heatweed 
#0.2232044 1.0717430 

rowMeans(means)
#      cortex            pith vascular tissue 
#0.57902499      1.29447711      0.06891892 

colMeans(se)
#  Control  Heatweed 
#0.1024732 0.2340615  

rowMeans(se)
#       cortex            pith vascular tissue 
#0.16939395      0.31242267      0.02298541 

#model fitting
mod3 = glm.nb(Antal ~ Behandling*Struktur, data=data)
anova(mod3)
summary(mod3)

mod4 = glm.nb(Antal ~ Behandling*Struktur -1, data=data)
anova(mod4)
summary(mod4) #?? HELP what happened to cortex?

data$log = log(data$Antal +1) #Try to get normal data for 2-way anova
hist(data$log) #Still not normal
model = aov(log ~ Behandling*Struktur, data=data)
summary(model) #Gives sig. but not really allowed to use this test

model2 = glmmTMB(log ~ Behandling*Struktur + (1|Lokal), data=data)
summary(model2)

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

ggsave("str_treat_interaction_PI.png", plot = last_plot(), device = "png",
       scale = 1, width = 10, height = 8,
       dpi = 600)





