#Andel_VH####

rm(list=ls())

getwd()
VH <- read.csv("Andel_VH.csv", sep=';')
str(VH)

head(dataAndel)
str(dataAndel)
dataAndel$Plats<- as.factor(dataAndel$Plats)
table(dataAndel$Plats)
dataAndel$Procent.död<-as.numeric(data$Procent.död)
mean_LokalAndel<-ddply(dataAndel,"Plats",summarise,mean_död=mean(Procent.död))
barplot(mean_kommitupp~Lokal, data=mean_LokalVH)