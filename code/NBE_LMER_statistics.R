#### Statistics ####

# 03.04.2023
# Charlotte Kunze

# packages
library(lme4)
library(tidyverse)
library(mgcv)
library(here)
library(rr2)
library(performance)

# data
data <- read.csv('~/Desktop/Exp22/MicrocosmExp22/output/netEffect.csv')
str(data)
tempPalette1 <- c("#E41A1C" ,"#377EB8" ,"#4DAF4A" )

#first look at data
dev.off()
data %>%
  group_by(N, temp)%>%
  mutate(mean = mean(AUC.ges.RR, na.rm = T)) %>%
  ggplot(.,)+
  geom_point(aes(N, AUC.ges.RR,colour = temp),size=1.5, alpha=0.5, position = position_dodge(width = 0.2))+
  geom_point(aes(N, mean,shape = temp),alpha = 0.5,size = 2,color = 'black', position = position_dodge(width = 0.2))+
 # geom_errorbar(aes(N, ymin = mean-se, ymax = mean+se), position = position_dodge(width = 0.2),width = .2)+
  labs(y = 'Net Effect (Mean in black)', colour = 'Treatment', shape = 'Treatment')+
  scale_color_manual(values=tempPalette1)+
  scale_shape_manual(values = c(1,2,7))+
  theme_bw()+
  theme(legend.position= 'bottom')
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/rawData.png'), width = 4, height= 3 )
str(master) # structure check


#### Statistics ####

### Nested model: combination as random effect nested in N 
### 3 singular models instead of 1 as the two fixed effects are correlated --> to decipher underlying mechanisms


### t-test ###
t.test(data$NBE)


##### LMM: Interaction// Fixed and random effects ######
# Response variable: ln of zooplankton Carbon mg (measured over time)
# Explanatory variables: treatments (fluctuation hours) + time
# Random: MC number

# If you want to estimate a trend then you need to include TEMP as a fixed effect as mentioned above - 
## and you could allow this to differ between richness levels by including it as a random slope, 

#read ZUUR et al. 2008

# example: combination is nested within Nfax 

## format variable
names(data)
data$temp <- as.factor(data$temp)
data$N <- as.numeric(data$N)
data$Nfax <- as.factor(data$N)

#subset of data 
fluc <- subset(data, temp == 'Fluctuation')
inc <- subset(data, temp == 'Increase')
incfluc <- subset(data, temp == 'IncreaseFluctuation')


### fluctuation ###

#random intercept
flucM1 <- lme(AUC.ges.RR ~ Nfax, random = ~1|combination, method = 'REML', data = fluc)

#random intercept and slope
ctrl <- lmeControl(opt='optim'); 
flucM2 <- lme(AUC.ges.RR ~ Nfax, random = ~1+Nfax|combination, method = 'REML', control = ctrl, data = fluc)
flucM3  = lme(AUC.ges.RR ~ Nfax, random = ~Nfax|combination, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = fluc)

#compare models
anova(flucM1, flucM2, flucM3)

#check output
anova(flucM1)
summary(flucM1)
check_model(flucM1)
ggsave(plot = last_plot(), file = 'flucPlot.png')

# to check explanatory power of fixed effects, use simple lm 
flucM0 <- gls(AUC.ges.RR ~ Nfax, method = 'REML', data = fluc) # 0 Model

rr2::R2_lik(flucM0)

### inc ###

#random intercept
incM1 <- lme(AUC.ges.RR ~ Nfax, random = ~1|combination, method = 'REML', data = inc)

#random intercept and slope
ctrl <- lmeControl(opt='optim'); 
incM2 <- lme(AUC.ges.RR ~ Nfax, random = ~1+Nfax|combination, method = 'REML', control = ctrl, data = inc)
incM3  = lme(AUC.ges.RR ~ Nfax, random = ~Nfax|combination, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = inc)

#compare models
anova(incM1, incM2, incM3)

#check output
anova(incM1)
summary(incM1)
check_model(incM1)
ggsave(plot = last_plot(), file = 'incPlot.png')

# to check explanatory power of fixed effects, use simple lm 
incM0 <- gls(AUC.ges.RR ~ Nfax, method = 'REML', data = inc) # 0 Model

rr2::R2_lik(incM0)

### incfluc ###

#random intercept
incflucM1 <- lme(AUC.ges.RR ~ Nfax, random = ~1|combination, method = 'REML', data = incfluc)

#random intercept and slope
ctrl <- lmeControl(opt='optim'); 
incflucM2 <- lme(AUC.ges.RR ~ Nfax, random = ~1+Nfax|combination, method = 'REML', control = ctrl, data = incfluc)
incflucM3  = lme(AUC.ges.RR ~ Nfax, random = ~Nfax|combination, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = incfluc)

#compare models
anova(incflucM1, incflucM2, incflucM3)

#check output
anova(incflucM1)
summary(incflucM1)
check_model(incflucM1)
ggsave(plot = last_plot(), file = 'incFlucPlot.png')


# to check explanatory power of fixed effects, use simple lm 
incflucM0 <- gls(AUC.ges.RR ~ Nfax, method = 'REML', data = incfluc) # 0 Model

rr2::R2_lik(incflucM0)



###############################################################################
#random intercept
NMM1 <- lme(AUC.ges.RR ~ Nfax, random = ~1|combination, method = 'REML', data = data)

#random intercept and slope
ctrl <- lmeControl(opt='optim'); 
NMM2 <- lme(AUC.ges.RR ~ Nfax, random = ~1+Nfax|temp, method = 'REML', control = ctrl, data = data)
NMM3 = lme(AUC.ges.RR ~ Nfax, random = ~Nfax|temp, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = data)

#compare models
anova(NMM1, NMM3, NMM2)

#check output
anova(NMM2)
check_model(NMM1)

#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
M1 = lme(AUC.ges.RR ~ Nfax*temp,random = ~1|combination, method = 'REML', data = data) #random intercept, fixed slope
M2 = lme(AUC.ges.RR ~ Nfax*temp, random = ~0+N|combination, method = 'REML', data = data)#random intercept and slope 
M3 = lme(AUC.ges.RR ~ Nfax*temp, random = ~N|combination, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = data) #random slope? fixed intercept


#compare models
anova(M1,M2, M3)


#save fitted values
data$fit_InterceptOnly2 <- predict(M2)


# mixed model vs. model without random component: gls 
M0=gls(AUC.ges.RR ~ Nfax*temp, method="REML",data =data, na.action=na.omit)
anova(M2, M0) #gls is better
summary(M2)
#Autocorrelation test 
plot(ACF(M2), alpha=0.05)

#Residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(M2, type = "normalized"), ylab="residuales")
hist(resid(M2, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(M2),resid(M1, type = "normalized"),ylab="residuales")
qqnorm(resid(M2, type = "normalized"), main=""); qqline(resid(M2, type = "normalized"))
check_model(M2)

anova(M2)
summary(M2)
# to check explanatory power of fixed effects, use simple lm 
rr2::R2_lik(M0)


data$Nfax <- as.factor(data$N)
### aditive ###
M4 = lme(AUC.ges.RR ~ Nfax+temp, random = ~1|combination, method = 'REML', data = data)
M5 = lme(AUC.ges.RR ~ Nfax+temp, random = ~0+N|combination, method = 'REML', data = data)
M6 = lme(AUC.ges.RR ~ Nfax+temp, random = ~N|combination, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = data)
anova(M4, M5, M6)
check_model(M5)
summary(M5)

M0.1=gls(AUC.ges.RR ~ Nfax+temp, method="REML",data =data, na.action=na.omit)

# to check explanatory power of fixed effects, use simple lm 
rr2::R2_lik(M0.1)


## GAM ##
#Additive mixed model 

M_total <- gamm(AUC.ges.RR ~ s(N, k = 3) + temp,random = list(combination =~ 1),  na.action=na.omit, method = "REML", data = data)

summary(M_total $gam) #This gives detailed output on the smoother and parametric
gam.check(M_total$gam)
#terms in the models.
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)

anova(M_total $gam) #This command gives a more compact presentation of the
#results. The anova table is not doing sequential testing!
plot(M_total $gam) #This command plots the smoothers.
plot(M_total $lme) #This command plots the normalised residuals versus fitted

#values and can be used to assess homogeneity.
summary(M_total $lme) #Detailed output on the estimated variances. 
summary(M_total$gam)


