#Andel_VH####

rm(list=ls())

getwd()
VH <- read.csv("Andel_VH.csv", sep=';')
str(VH)

head(VH)
str(VH)
VH$Plats<- as.factor(VH$Plats)
table(VH$Plats)
VH$Procent.död = gsub(",",".", VH$Procent.död)
VH$Procent.död<-as.numeric(VH$Procent.död)
mean_LokalAndel<-ddply(VH,"Plats",summarise,mean_död=mean(Procent.död))
barplot(mean_död~Plats, data=mean_LokalAndel)

qqnorm(VH$död)
hist(VH$död)

