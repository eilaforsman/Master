rm(list = ls())
list.files ()

#Stapeldiagram medelantal per plats, medel på enskilda snitten först 
#Antal HW vs kontroll
#Antal skillnad i lokaler beror på region?
#Two-way ANOVA HW vs kontroll och antal per lokal
#Medelvärde antal per lokal (Pålsjö 1,2,3,4,5)
#5vs6ggr
#Skillnad i död/levande eller antal mellan strukturer

#Basic####
getwd()

setwd("~/Documents/Documents/R/Master")

dataFDA <- read.csv("dataFDA.csv", sep = ';')
head(dataFDA)
str(dataFDA)
dataFDA$Struktur <- as.factor(dataFDA$Struktur)
dataFDA$Lokal <- as.factor(dataFDA$Lokal)
table(dataFDA$Lokal)
summary(dataFDA)
boxplot(Antal~Lokal,data=dataFDA) #Testar boxplot


mean_LokalFDA<-ddply(dataFDA,"Lokal",summarise,
                  mean_antal=mean(Antal))

mean_ID_FDA<-ddply(dataFDA,"ID",summarise,
                     mean_antal=mean(Antal))
head(mean_ID_FDA)


head(dataFDA)

#install.packages("stringr")
library(stringr)

#Sort_data####

ID_split <- as.data.frame(str_split_fixed(dataFDA$ID, "_", 2))
head(ID_split)

nchar(dataFDA$Prov)
dataFDA$newprov <- substr(dataFDA$Prov, 1, nchar(dataFDA$Prov) - 2)

dataFDA$ID_spec <- ID_split$V1
paste(dataFDA$Lokal, dataFDA$ID_spec, sep ="_")
dataFDA$Prov<-paste(dataFDA$Lokal, dataFDA$ID_spec, sep ="_")
head(dataFDA)

Prov<-ddply(dataFDA, "Prov", summarise, mean_antal=mean(Antal))
head(Prov)

Meantot<-ddply(Prov, "newprov", summarise, mean_tot=mean(mean_antal))

#Plotting####

library(plyr)
barplot(mean_antal~Lokal,data=mean_LokalFDA,srt=35)
str(mean_Antal)
axis(side=1,labels = FALSE)

barplot(mean_antal~Lokal,data=mean_LokalFDA, xaxt="n")
text(x = 1:length(levels(mean_LokalFDA$Lokal)),
     y = par("usr")[3] - 0.45,
     labels = levels(mean_LokalFDA$Lokal),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 50,
     cex = 1.0)

#install.packages("ggplot2")
library(ggplot2)

head(MeanProv)

?ggplot()

FDAMeans<-ddply(dataFDA,"Lokal",summarize,N=length(Antal),
                         mean.antal=mean(Antal)
                         ,sd.FDA=sd(mean.antal),
                         se.FDA=sd.FDA/sqrt(N))

ggplot(Meantot, aes(x= reorder(newprov, -mean_tot), y=mean_tot)) +
        geom_bar(width = 0.75, stat = "identity", position ="dodge", alpha = 0.8) +
        #geom_errorbar(aes(ymin= Meantot - se_Meantot, ymax=Meantot + se_Meantot),
         #             width = 0.13, alpha = 1, position=position_dodge(0.75)) +
        theme_classic() + 
        #scale_fill_manual(values=c("#73ba2e", "#034c77")) +
        #scale_color_manual(values=c("#000000", "#000000")) +
        scale_y_continuous(limits = c(0,150), expand = c(0,0)) +
        labs(y="Antal per prov", x="", title = "") +
        theme(legend.position = c(0.9,0.9), legend.title = element_blank(),
              plot.title = element_text(hjust = -0.15),
              text = element_text(size=20, family= "Times"),
              axis.text.x = element_text(size = 12, angle = 45,
                                         hjust = 1, color = "grey1")) +
        theme(axis.ticks.length=unit(.25, "cm"))









#Växthus####

dataVH <- read.csv("Växthus_data.csv", sep = ';')

head(dataVH)
str(dataVH)
dataVH$Lokal<- as.factor(dataVH$Lokal)
table(dataVH$Lokal)
mean_LokalVH<-ddply(dataVH,"Lokal",summarise,mean_kommitupp=mean(Kommit.upp))
barplot(mean_kommitupp~Lokal, data=mean_LokalVH)


#VH_andel_död####


dataAndel<- read.csv("Andel_VH.csv" , sep= ';')

head(dataAndel)
str(dataAndel)
dataAndel$Plats<- as.factor(dataAndel$Plats)
table(dataAndel$Plats)
dataAndel$Procent.död<-as.numeric(data$Procent.död)
mean_LokalAndel<-ddply(dataAndel,"Plats",summarise,mean_död=mean(Procent.död))
barplot(mean_kommitupp~Lokal, data=mean_LokalVH)




