#### Statistics ####
# 11.07.2023
# Helmut Hillebrand & Charlotte Kunze

library(tidyverse)

#### data ####
netdiv <- read.csv('~/Desktop/Exp22/MicrocosmExp22/Data/netEffect.csv')%>%
  mutate(absNBE = abs(NBE))


names(netdiv)
netdiv$combination<-as.factor(netdiv$combination)
summary(netdiv)

aov1<-aov(NBE~temp*combination, netdiv)
summary(aov1)
TukeyHSD(aov1)

ggplot(netdiv, aes(y = NBE, x = as.factor(temp), fill=as.factor(N)))+
  geom_boxplot()

# test for normality
fligner.test(NBE~interaction(temp,combination), netdiv)


#test against zero
t.test(netdiv$NBE, mu = 0, alternative = "two.sided")


#### explore the interaction of temperature and species combinations ####
# analyse if presence absence of species affects NBE
combiEffect <- netdiv %>%
  mutate(composition = combination) %>%
  mutate(A = ifelse(str_detect(combination, 'A'), 1, 0), # for presence 1 for absence 0 
         D = ifelse(str_detect(combination, 'D'), 1, 0),
         G = ifelse(str_detect(combination, 'G'), 1, 0),
         R = ifelse(str_detect(combination, 'R'), 1, 0), 
         T = ifelse(str_detect(combination, 'T'), 1, 0)) 
  
fluct <- combiEffect %>% filter(temp == 'Fluctuation')
aov2<-aov(NBE~A+D+G+R+T, fluct)
summary(aov2)
TukeyHSD(aov2)

incfluct <- combiEffect %>% filter(temp == 'IncreaseFluctuation')
aov3<-aov(NBE~A+D+G+R+T, incfluct)
summary(aov3)

inc <- combiEffect %>% filter(temp == 'Increase')
aov4<-aov(NBE~A+D+G+R+T, inc)
summary(aov4)

#### contrasts ####
# contrast 1: 2 versus 4 species

netdiv$con1<-NA
netdiv$con1[netdiv$N==2]<- "A"
netdiv$con1[netdiv$N==4]<- "B"

lm1<-aov(NBE~temp*con1, netdiv)
summary(lm1)
TukeyHSD(lm1)

# contrast 2: 2 versus 4 species

netdiv$con2<-NA
netdiv$con2[netdiv$N==2]<- "A"
netdiv$con2[netdiv$N==5]<- "B"

lm2<-aov(NBE~temp*con2, netdiv)
summary(lm2)
TukeyHSD(lm2)

# contrast 3: 4 versus 5 species

netdiv$con3<-NA
netdiv$con3[netdiv$N==4]<-"A"
netdiv$con3[netdiv$N==5]<-"B"

lm3<-aov(NBE~temp*con3, netdiv)
summary(lm3)


#### NBE data ####

HectorRaw <- read.csv('~/Desktop/Exp22/MicrocosmExp22/Data/NBEonFunctioning.csv')
str(HectorRaw)

HectorRaw$combination<-as.factor(HectorRaw$combination)
summary(HectorRaw)


aov1<-aov(NetEffect~temp*combination, HectorRaw)
summary(aov1)
TukeyHSD(aov1)

ggplot(HectorRaw, aes(y = NetEffect, x = as.factor(temp), fill=as.factor(N)))+
  geom_boxplot()

# test for normality
maxNE <- max(abs(HectorRaw$NetEffect))
HectorRaw$trans = log(HectorRaw$NetEffect + maxNE)
fligner.test(trans~interaction(temp,combination), HectorRaw)

plot(resid(aov1))
hist(resid(aov1))
hist(((HectorRaw$NetEffect)^2))


hist(cuberoot(HectorRaw$NetEffect))
dummy <- (HectorRaw$NetEffect^2)
summary(dummy)
qqnorm(log(HectorRaw$NetEffect^2))

#test against zero
t.test(HectorRaw$NetEffect, mu = 0, alternative = "two.sided")



