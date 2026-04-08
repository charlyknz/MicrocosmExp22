#### Start Rscript ####
# 11.07.2023
# Helmut Hillebrand & Charlotte Kunze

library(tidyverse)
library(ggpubr)
library(here)
library(mgcv)
library(rr2)


#### data NBES ####

# run R scripts 
source(here("code/03_NBES_calculation.R"))
str(d3)
netdiv <-d3%>%
  mutate(absNBE = abs(NBE))

# or import as csv
#netdiv <-  read.csv('Data/NBES.csv')

## check data 
summary(netdiv)
names(netdiv)

# combination as factor for analysis
netdiv$combination<-as.factor(netdiv$combination)
netdiv$Nfac<-as.factor(netdiv$N)

##### Mixed Model: NBES ~ SR * temperature #####
#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
M1 = lme(NBE ~ temp*Nfac, random = ~1|combination, method = 'REML', data = netdiv)
M2 = lme(NBE ~ temp*Nfac, random = ~0+Nfac|combination, method = 'REML', data = netdiv)
M3 = lme(NBE ~ temp*Nfac, random = ~Nfac|combination, method = 'REML',
         control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = netdiv)

#compare models
anova(M1,M2, M3)

#save fitted values
netdiv$fit_InterceptOnly2 <- predict(M2)


# mixed model vs. model without random component: gls 
M0=gls(NBE~temp*Nfac, method="REML",data =netdiv, na.action=na.omit)
anova(M2, M0) #gls is better

#Autocorrelation test 
plot(ACF(M2), alpha=0.05)

#Residuals
par(mfrow=c(1,1),cex.axis=1.2, cex.lab=1.5)
plot(resid(M2, type = "normalized"), ylab="residuales")
hist(resid(M2, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(M2),resid(M2, type = "normalized"),ylab="residuales")
qqnorm(resid(M2, type = "normalized"), main=""); qqline(resid(M2, type = "normalized"))

anova(M2)
summary(M2)

# posthoc test
emmeans::emmeans(M2, ~Nfac)
emmeans::emmeans(M2, ~temp)


##### NBES: T-test #####
#test against zero
test1 <- t.test(netdiv$NBE, mu = 0, alternative = "two.sided")
test1


##### NBES: Contrasts #####
# contrast 1 (con): 2 versus 4 species

netdiv$con1<-NA
netdiv$con1[netdiv$N==2]<- "A"
netdiv$con1[netdiv$N==4]<- "B"

lm1<-aov(NBE~temp*con1, netdiv)
summary(lm1)
TukeyHSD(lm1)
emmeans::emmeans(lm1, ~pairwise ~ con1 | temp)

# contrast 2: 2 versus 5 species

netdiv$con2<-NA
netdiv$con2[netdiv$N==2]<- "A"
netdiv$con2[netdiv$N==5]<- "B"

lm2<-aov(NBE~temp*con2, netdiv)
summary(lm2)
TukeyHSD(lm2)
emmeans::emmeans(lm2, ~pairwise ~ con2 | temp)

# contrast 3: 4 versus 5 species

netdiv$con3<-NA
netdiv$con3[netdiv$N==4]<-"A"
netdiv$con3[netdiv$N==5]<-"B"

lm3<-aov(NBE~temp*con3, netdiv)
summary(lm3)
emmeans::emmeans(lm3, ~pairwise ~ con3 | temp)

#### data NBE Functioning ####

source(here("code/04_NBE_HectorLoreau_NetBiodivEffect.R"))
str(allNetBiodiv)
HectorRaw <- allNetBiodiv 

#or
#import as csv
#HectorRaw <- read.csv('Data/NBEonFunctioning.csv') %>%
#  select(-X)
str(HectorRaw)

HectorRaw$combination<-as.factor(HectorRaw$combination)
summary(HectorRaw)


##### Mixed Model: NBEF ~ SR * temperature #####
#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
H1 = lme(NetEffect ~ temp*N, random = ~1|combination, method = 'REML', data = HectorRaw)
H2 = lme(NetEffect ~ temp*N, random = ~0+N|combination, method = 'REML', data = HectorRaw)
H3 = lme(NetEffect ~ temp*N, random = ~N|combination, method = 'REML',
         control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = HectorRaw)

#compare models
anova(H1,H2, H3)

#save fitted values
HectorRaw$fit_InterceptOnly2 <- predict(H1)


# mixed model vs. model without random component: gls 
M0=gls(NetEffect~temp*N, method="REML",data =HectorRaw, na.action=na.omit)
anova(H1, M0) #gls is better

#Autocorrelation test 
plot(ACF(H1), alpha=0.05)

#Residuals
par(mfrow=c(1,1),cex.axis=1.2, cex.lab=1.5)
plot(resid(H1, type = "normalized"), ylab="residuales")
hist(resid(H1, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(H1),resid(H2, type = "normalized"),ylab="residuales")
qqnorm(resid(H1, type = "normalized"), main=""); qqline(resid(H1, type = "normalized"))

anova(H1)
summary(H1)

##### NBE on F: T-Test #####
#test against zero
test2<-t.test(HectorRaw$NetEffect, mu = 0, alternative = "two.sided")
test2


##### NBES - NBE on F: Correlation #####
Corr.data <- netdiv %>%
  select(combination, temp, rep, NBE, N) %>%
  mutate(N = as.factor(N)) %>%
  left_join(.,HectorRaw, by = c('combination', 'temp', 'rep', 'N'))


ggplot(Corr.data, aes(x= NetEffect, y = NBE, color = N, shape = temp))+
 geom_abline (slope=1, linetype = "dashed", color="black")+
 geom_point()+
  labs(x ='NBEF', y= 'NBES', shape ='Treatment')+
  facet_wrap(~temp)+
geom_smooth(method = 'lm', se =F)+
  theme_bw()
#ggsave(plot = last_plot(), file = here('NBES_NBEF_correlation.png'), width = 10, height = 4)
ggscatter(Corr.data, x = 'NetEffect', y='NBE', add = 'reg.line', col = 'temp',cor.coef = T, xlab = 'NBE on Functioning', ylab = 'NBES')+
  stat_cor(aes(color = temp), label.x = 3)

ggsave(plot = last_plot(), file = here('output/Correlation_NBE_NBES.png'), width = 6, height = 5)


#### Additional statistical tests ####

##### NBES: ANOVA temp combination #####
aov1<-aov(NBE~temp*combination, netdiv)
summary(aov1)

ggplot(netdiv, aes(y = NBE, x = as.factor(temp), fill=as.factor(N)))+
  geom_boxplot()

# test for normality
hist(netdiv$NBE)
qqnorm(netdiv$NBE)
qqline(netdiv$NBE)

##### NBE on F: ANOVA temp, combination #####
aov1<-aov(NetEffect~temp*combination, HectorRaw)
summary(aov1)
TukeyHSD(aov1)


ggplot(HectorRaw, aes(y = NetEffect, x = as.factor(temp), fill=as.factor(N)))+
  geom_boxplot()

# test for normality
fligner.test(NetEffect~interaction(temp,combination), HectorRaw)

plot(resid(aov1))
hist(resid(aov1))
hist(((HectorRaw$NetEffect)))
qqnorm(HectorRaw$NetEffect)
qqline(HectorRaw$NetEffect)

##### NBES: Presence/absence of species ~NBES #####
# explore the interaction of temperature and species combinations
# analyse if presence absence of species affects NBE
# separately for each temperature treatment

combiEffect <- netdiv %>%
  mutate(composition = combination) %>%
  mutate(A = ifelse(str_detect(combination, 'A'), 1, 0), # for presence 1 for absence 0 
         D = ifelse(str_detect(combination, 'D'), 1, 0),
         G = ifelse(str_detect(combination, 'G'), 1, 0),
         R = ifelse(str_detect(combination, 'R'), 1, 0), 
         Th = ifelse(str_detect(combination, 'T'), 1, 0)) 

fluct <- combiEffect %>% filter(temp == 'Fluctuation')
aov2<-aov(NBE~A+D+G+R+Th, fluct) 
summary(aov2)

incfluct <- combiEffect %>% filter(temp == 'Increase + Fluctuation')
aov3<-aov(NBE~A+D+G+R+Th, incfluct)
summary(aov3)

inc <- combiEffect %>% filter(temp == 'Increase')
aov4<-aov(NBE~A+D+G+R+Th, inc)
summary(aov4)


### richness levels ###
library(rstatix )
fluct2 <- combiEffect %>% filter(temp == 'Fluctuation' & N == 2)
aov5<-aov(NBE~A+D+G+R+Th, fluct2) 
anova(aov5)
rstatix::eta_squared(aov5)

fluct4 <- combiEffect %>% filter(temp == 'Fluctuation' & N == 4)
aov6<-aov(NBE~A+D+G+R+Th, fluct4) 
summary(aov6)

incfluct2 <- combiEffect %>% filter(temp == 'Increase + Fluctuation'& N == 2)
aov7<-aov(NBE~A+D+G+R+Th, incfluct2)
summary(aov7)

incfluct4 <- combiEffect %>% filter(temp == 'Increase + Fluctuation'& N == 4)
aov8<-aov(NBE~A+D+G+R+Th, incfluct4)
summary(aov8)

inc2 <- combiEffect %>% filter(temp == 'Increase'& N == 2)
aov9<-aov(NBE~A+D+G+R+Th, inc2)
summary(aov9)

inc4 <- combiEffect %>% filter(temp == 'Increase'& N == 4)
aov10<-aov(NBE~A+D+G+R+Th, inc4)
summary(aov10)

summary.table<-tibble()
summary.table  <- anova(aov10)["Sum Sq"]

cor(inc4[,12:16])
