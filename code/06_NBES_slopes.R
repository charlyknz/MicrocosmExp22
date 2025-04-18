#### look at NBES slopes ####
# 09.12.2024 by Charlotte Kunze

library(tidyverse)
library(here)
library(readxl)

### import data ###
data <- read.csv(('Data/NBES_revisited.csv')) %>%
  select(-X) %>%
  mutate(communityID = paste(rep, temp,sep ='_'))

#### Resilience slope ####
USI <- unique(data$communityID)

#create an empty data frame
slope<-tibble()

for(i in 1:length(USI)){
  temp<-data[data$communityID==USI[i], ]#creates a temporary data frame for each case
  lm1<-lm(NBE~log(N), temp)#makes a linear regression
  intcp.lm <- coef(summary(lm1))[1, 1]#selects the intercept
  se.intcp.lm<- coef(summary(lm1))[1, 2]#selects its standard error
  resil.lm <- coef(summary(lm1))[2, 1]#selects the slope
  se.slp.lm<- coef(summary(lm1))[2, 2]#selects its standard error
  sd.res.lm<- sd(resid(lm1)) #selects the standard deviation of the residuals
  slope<-rbind(slope,data.frame(temp[1,c(2,3,8)],resil.lm,intcp.lm,se.intcp.lm,sd.res.lm,se.slp.lm))
  rm(temp)
}
names(slope)


#### Start ####
ggplot(slope, aes(x = temp, y = resil.lm))+
  geom_point()

data1 <- data%>%
 left_join(., slope)

data1%>%
  #filter(resil.lm <0) %>%
  group_by(N, temp, communityID) %>%
  reframe(mean = mean(NBE))%>%
  ggplot(., aes (x = N, y = mean, color = communityID, group = communityID))+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_line()+
  facet_grid(~temp)
 # theme(legend.position = 'none')


##### TPC- Grand Mean ####
TPC<- read_excel("~/Desktop/phD/Exp22/Experiments/CharlyTPC2021/createTPC/Species_TPC_maxBiom.xlsx")
str(TPC)

all_tpc_output <- left_join(TPC, anovaOutput) %>%
  mutate(q= ifelse(p.value<0.05, 1,0)) %>%
  drop_na(p.value)
ggplot(all_tpc_output, aes(x = topt, y = sum_sq, color = treatment, group = treatment))+
  geom_point(size =2)+
  facet_grid(~richness)+
  theme_bw()


grandMeanA <- data1 %>%
  filter(str_detect(combination, 'A'))%>%
  group_by(N, temp)%>%
  reframe(mean_NBES = mean(NBE),
          species = paste('A'))
  
grandMeanD <- data1 %>%
  filter(str_detect(combination, 'D'))%>%
  group_by(N, temp)%>%
  reframe(mean_NBES = mean(NBE),
          species = paste('D'))

grandMeanG <- data1 %>%
  filter(str_detect(combination, 'G'))%>%
  group_by(N, temp)%>%
  reframe(mean_NBES = mean(NBE),
          species = paste('G'))

grandMeanR <- data1 %>%
  filter(str_detect(combination, 'R'))%>%
  group_by(N, temp)%>%
  reframe(mean_NBES = mean(NBE),
          species = paste('R'))

grandMeanT <- data1 %>%
  filter(str_detect(combination, 'T'))%>%
  group_by(N, temp)%>%
  reframe(mean_NBES = mean(NBE),
          species = paste('T'))

all_grandMean <- grandMeanA %>%
  bind_rows(., grandMeanD)%>%
  bind_rows(., grandMeanG)%>%
  bind_rows(., grandMeanR)%>%
  bind_rows(., grandMeanT) 

GrandMean <- data1 %>%
  group_by(N, temp)%>%
  reframe(gMean = mean(NBE)) %>%
  right_join(., all_grandMean) %>%
  mutate(devFromGrandMean = mean_NBES-gMean) %>%
  filter(N != 5) %>%
  left_join(., TPC)

GrandMean$temp<-factor(GrandMean$temp, levels = c("Increase" ,'Fluctuation', "Increase + Fluctuation"))
GrandMean$interaction <- NA
GrandMean$interaction <- paste('Experimental communities')

####plot grand mean ####
# palette
colours<-c("darkblue", "skyblue", '#EFC000FF', 'darkorange')
pa <- GrandMean%>%
  filter(temp == 'Increase')%>%
  ggplot(., aes( y = topt, x = devFromGrandMean, color = N))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2)+
  labs(x = 'Influence on NBES',  title = 'Increase', color = 'Richness')+
  ylab(bquote(b[opt]))+
  facet_grid(~interaction)+
  scale_color_gradientn(colours = colours)+
  theme_bw()+
  theme(legend.position = 'none')+
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.text.x  = element_text(size = 12))


pb <-GrandMean%>%
  filter(temp == 'Fluctuation')%>%
  ggplot(., aes( y = topt, x = devFromGrandMean, color = N))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2)+
  labs(x = 'Influence on NBES',  title = 'Fluctuation', color = 'Richness')+
  scale_color_gradientn(colours = colours)+
  ylab(bquote(b[opt]))+
  facet_grid(~interaction)+
  theme_bw()+
  theme(legend.position = 'none')+
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.text.x  = element_text(size = 12))
  

pc <-GrandMean%>%
  filter(temp == 'Increase + Fluctuation')%>%
  ggplot(., aes( y = topt, x = devFromGrandMean, color = N))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2)+
  labs(x = 'Influence on NBES',  title = 'Increase + Fluctuation', color = 'Richness')+
  ylab(bquote(b[opt]))+
  facet_grid(~interaction)+
  scale_color_gradientn(colours = colours)+
  theme_bw()+
  theme(legend.position = 'none')+
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.text.x  = element_text(size = 12))


cowplot::plot_grid(pa, pb, pc, ncol = 1,labels = c('(d)', '(e)', '(f)'))
ggsave(plot=last_plot(), file = here('output/Figure3_topt_NBESeffect-temp.png'), width = 3, height = 8)
