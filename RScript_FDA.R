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

setwd("~/Documents/Documents/R/Master/Master")

dataFDA <- read.csv("dataFDA.csv", sep = ';')
head(dataFDA)
str(dataFDA)
dataFDA$Struktur <- as.factor(dataFDA$Struktur)
dataFDA$Lokal <- as.factor(dataFDA$Lokal)
table(dataFDA$Lokal)
summary(dataFDA)
boxplot(Antal~Lokal,data=dataFDA) #Testar boxplot

library(plyr)
dataFDA_sub <- subset(dataFDA, Lokal != "H_Pålsjö_efter_HW")
mean_LokalFDA<-ddply(dataFDA_sub,"Lokal",summarise,
                  mean_antal=mean(Antal))

mean_ID_FDA<-ddply(dataFDA_sub,"ID",summarise,
                     mean_antal=mean(Antal))
head(mean_ID_FDA)


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


Prov<-ddply(dataFDA_sub, "Prov", summarise, mean_antal=mean(Antal))
head(Prov)

nchar(dataFDA_sub$Prov)
dataFDA_sub$newprov <- substr(dataFDA_sub$Prov, 1, nchar(dataFDA_sub$Prov) - 2)
Meantot<-ddply(dataFDA_sub, "newprov", summarise, mean_tot=mean(Antal))
Meantot$order <- c(4,5,1,6,7,8,9,10,11,2,12,13,14,3)

#Basic Plotting####

library(plyr)
barplot(mean_antal~Lokal,data=mean_LokalFDA,srt=35)
str(mean_antal)
axis(side=1,labels = FALSE)

barplot(mean_antal~Lokal,data=mean_LokalFDA, xaxt="n")
text(x = 1:length(levels(mean_LokalFDA$Lokal)),
     y = par("usr")[3] - 0.45,
     labels = levels(mean_LokalFDA$Lokal),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 50,
     cex = 1.0)

#Fancy plotting######

#install.packages("ggplot2")
library(ggplot2)

head(Meantot)

?ggplot()



FDAMeans<-ddply(dataFDA_sub,"Lokal", summarize, N=length(Antal),
                         mean.antal=mean(na.omit(Antal)),
                        sd.FDA=sd(na.omit(Antal)),
                         se.FDA=sd.FDA/sqrt(N))

ggplot(Meantot, aes(x=reorder(Meantot$newprov, Meantot$order), y=Meantot$mean_tot)) +
        geom_bar(width = 0.75, stat = "identity", position ="dodge", alpha = 0.8) +
        geom_errorbar(data=FDAMeans, aes(ymin= mean.antal - se.FDA, ymax=mean.antal + se.FDA),
                   width = 0.13, alpha = 1, position=position_dodge(0.75)) +
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

ggsave("Celler_FDA_plot", plot = last_plot(), device = "png",
       scale = 1, width = 12, height = 8,
       dpi = 600)



#ANOVA#####

dataFDA_sub$Behandling = factor(dataFDA_sub$Behandling, levels=c("Kontroll", "Heatweed"))
m = lm(dataFDA_sub$Antal~dataFDA_sub$Behandling)
anova(m)
summary(m)

hist(dataFDA_sub$Antal)
